#!/usr/bin/env python3
"""
Build and freeze the expanded ClinicalTrials.gov bacterial/AMR diagnostic
analytic release v3.3.0.

This builder deliberately DOES NOT run H1-H4.

Inputs
------
1. Frozen historic 360-study analytic cohort v3.2.7.
2. Validated/frozen neutral descriptor packet for the 213 newly eligible
   v3.2.9 rescue/expansion studies.
3. Frozen P35 full-cohort imaging-sensitivity exclusion NCT list.

Core guarantees
---------------
- Historic 360 rows are preserved exactly on every pre-existing historic column.
- New 213 eligibility, BROAD/CORE stratum, AMR depth, and depth-derived output
  type are not altered.
- New neutral descriptor signoff is standardized into the historic release
  descriptor-adjudication column names.
- P35 is represented as a separate analysis-stage flag and does not alter
  primary-cohort eligibility.
- The primary release must contain exactly 573 unique eligible NCTs.
- Hard-coded distribution gates protect the frozen 573-study counts.
- Existing historic special_sensitivity_flags are preserved.
- New H2 rare-pathogen sensitivity flags are added only if a new
  Enterobacterales study explicitly contains a prespecified typhoid/plague
  term in the frozen registry evidence.
- A heuristic H3 quantitative-mechanism near-miss candidate queue is generated
  for later descriptive QC only; it does NOT change any analytic flag or depth.
- Outputs are checksummed and versioned. Existing non-empty output directories
  are never overwritten.

The v3.3.1 statistical analysis should be run only after this release passes
all checks and is made read-only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


VERSION = "v3.3.0"

EXPECTED_HISTORIC_SHA256 = (
    "a59a9ec30d188533c2e4508bb8044150fb86e563208273cf10c4563b1543bda6"
)

EXPECTED = {
    "historic_n": 360,
    "new_n": 213,
    "total_n": 573,
    "depth": {"0": 435, "1": 83, "2": 52, "3": 3, "4": 0},
    "h2": {
        "GRAM_POSITIVE": 62,
        "ENTEROBACTERALES": 25,
        "OTHER_EXCLUDED": 486,
    },
    "utility": {"NO": 312, "YES": 261},
    "status": {"FINAL": 521, "FINAL_WITH_UNCERTAINTY": 52},
    "p35_excluded_n": 8,
    "p35_remaining_n": 565,
}

EXPECTED_NEW = {
    "depth": {"0": 177, "1": 23, "2": 13, "3": 0, "4": 0},
    "stratum": {
        "BROAD_BACTERIAL_DIAGNOSTIC": 177,
        "CORE_AMR_DIAGNOSTIC": 36,
    },
    "h2": {
        "GRAM_POSITIVE": 19,
        "ENTEROBACTERALES": 8,
        "OTHER_EXCLUDED": 186,
    },
    "utility": {"NO": 97, "YES": 116},
    "status": {"FINAL": 184, "FINAL_WITH_UNCERTAINTY": 29},
}

EXPECTED_P35_IDS = {
    "NCT01378728",
    "NCT02450942",
    "NCT02491164",
    "NCT02558062",
    "NCT03091361",
    "NCT03290690",
    "NCT05285072",
    "NCT06986512",
}

DEPTH_OUTPUT = {
    "0": "ORGANISM_ONLY",
    "1": "BINARY_OR_CATEGORICAL_RESISTANCE",
    "2": "PHENOTYPIC_AST_MIC_ZONE",
    "3": "INTEGRATED_MULTIMECHANISM",
    "4": "QUANTITATIVE_AMR_MECHANISM",
}

FINAL_DESCRIPTOR_FIELDS = [
    "final_primary_diagnostic_modality",
    "final_all_diagnostic_modalities",
    "final_organism_group",
    "final_gram_group",
    "final_h2_comparison_group",
    "final_analytical_endpoint_categories",
    "final_clinical_utility_endpoint_categories",
    "final_clinical_utility_any",
    "final_preanalytical_flag",
    "final_amr_reporting_intervention_flag",
    "final_mixed_viral_bacterial_panel_flag",
    "final_direct_patient_specimen_flag",
    "final_index_test_output_type",
]

REQUIRED_ANALYTIC_FIELDS = [
    "nct_id",
    "final_primary_eligible",
    "final_stratum",
    "final_amr_depth",
    "start_year",
    "study_type",
    *FINAL_DESCRIPTOR_FIELDS,
    "descriptor_adjudication_notes",
    "descriptor_adjudicator_initials",
    "descriptor_adjudication_status",
    "special_sensitivity_flags",
]

NEW_EVIDENCE_FIELDS_FOR_POLICY_AUDIT = [
    "brief_title",
    "official_title",
    "conditions",
    "keywords",
    "intervention_names",
    "primary_outcomes",
    "secondary_outcomes",
    "other_outcomes",
    "summary",
]

PROVENANCE_FIELDS = [
    "analytic_release_version",
    "cohort_origin",
    "source_descriptor_release",
    "p35_imaging_sensitivity_exclude",
    "p35_sensitivity_frozen_version",
]

# Prespecified H2 rare-pathogen sensitivity concept: typhoid/plague.
# This is deliberately constrained to rows already adjudicated as
# ENTEROBACTERALES, and it uses only the supplied frozen evidence.
RARE_H2_PATTERNS = [
    (
        "TYPHOID",
        re.compile(
            r"\btyphoid\b"
            r"|\bsalmonella\s+(?:enterica\s+)?(?:serovar\s+)?typhi\b"
            r"|\bs\.?\s*typhi\b",
            re.IGNORECASE,
        ),
    ),
    (
        "PLAGUE",
        re.compile(
            r"\bplague\b|\byersinia\s+pestis\b",
            re.IGNORECASE,
        ),
    ),
]

# Nonbinding H3 near-miss candidate scan. This is intentionally broad because
# it only makes a QC queue and does not modify depth or special flags.
H3_MECHANISM_TERMS = re.compile(
    r"\b("
    r"copy\s*number|gene\s*dosage|transcript(?:ion|omic)?|expression|"
    r"porin|permeability|efflux|enzyme\s*activity|beta[- ]?lactamase|"
    r"carbapenemase|resistance\s+(?:gene|mechanism|determinant)|"
    r"quantitative\s+(?:pcr|pcr|assay|measurement)|qPCR"
    r")\b",
    re.IGNORECASE,
)
H3_QUANT_TERMS = re.compile(
    r"\b(quantitative|quantify|quantification|copy\s*number|level|abundance|"
    r"concentration|expression|activity)\b",
    re.IGNORECASE,
)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def text(v: Any) -> str:
    if pd.isna(v):
        return ""
    return str(v).strip()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        encoding="utf-8-sig",
    )


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, sep="\t", index=False, lineterminator="\n")


def parse_flags(v: Any) -> list[str]:
    vals = [x.strip() for x in text(v).split("|") if x.strip()]
    if not vals:
        return []
    if "NONE" in vals and len(vals) > 1:
        vals = [x for x in vals if x != "NONE"]
    return vals


def add_flag(v: Any, flag: str) -> str:
    vals = parse_flags(v)
    if flag not in vals:
        vals.append(flag)
    return "|".join(sorted(set(vals))) if vals else "NONE"


def normalized_counter(series: pd.Series) -> dict[str, int]:
    c = Counter(text(x) for x in series)
    return dict(sorted(c.items()))


def depth_counter(series: pd.Series) -> dict[str, int]:
    c = Counter(text(x) for x in series)
    return {d: c.get(d, 0) for d in ["0", "1", "2", "3", "4"]}


def require_exact_counter(
    observed: dict[str, int],
    expected: dict[str, int],
    label: str,
    errors: list[str],
) -> None:
    if observed != expected:
        errors.append(f"{label} mismatch: observed={observed}; expected={expected}")


def read_p35_ids(path: Path) -> set[str]:
    ids = set()
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        # Permit a simple NCT-only file or a one-column TSV with a header.
        first = s.split("\t", 1)[0].strip()
        if re.fullmatch(r"NCT\d{8}", first, re.IGNORECASE):
            ids.add(first.upper())
    return ids


def ensure_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for c in columns:
        if c not in out.columns:
            out[c] = ""
    return out


def frozen_evidence_blob(row: pd.Series) -> str:
    vals = []
    for f in NEW_EVIDENCE_FIELDS_FOR_POLICY_AUDIT:
        if f in row.index:
            v = text(row[f])
            if v:
                vals.append(v)
    return "\n".join(vals)


def validate_depth_output(df: pd.DataFrame, label: str, errors: list[str]) -> None:
    for _, r in df.iterrows():
        nct = text(r["nct_id"])
        d = text(r["final_amr_depth"])
        out = text(r["final_index_test_output_type"])
        expected = DEPTH_OUTPUT.get(d)
        if expected is None:
            errors.append(f"{label} {nct}: invalid depth {d!r}")
        elif out != expected:
            errors.append(
                f"{label} {nct}: depth/output mismatch: depth={d}, "
                f"output={out!r}, expected={expected!r}"
            )


def make_manifest(paths: list[Path], out_path: Path) -> None:
    rows = []
    for p in sorted(paths, key=lambda x: x.name):
        parsed_rows = ""
        if p.suffix == ".tsv":
            try:
                parsed_rows = str(len(read_tsv(p)))
            except Exception:
                parsed_rows = "PARSE_ERROR"
        rows.append(
            {
                "filename": p.name,
                "sha256": sha256(p),
                "bytes": p.stat().st_size,
                "physical_lines": sum(
                    1 for _ in p.open("r", encoding="utf-8-sig", errors="replace")
                ),
                "parsed_rows_if_tsv": parsed_rows,
            }
        )
    write_tsv(pd.DataFrame(rows), out_path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--historic-eligible", required=True, type=Path)
    ap.add_argument("--new-neutral-frozen", required=True, type=Path)
    ap.add_argument("--p35-exclusion-ids", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    for p in [
        args.historic_eligible,
        args.new_neutral_frozen,
        args.p35_exclusion_ids,
    ]:
        if not p.exists():
            raise SystemExit(f"Required input missing: {p}")

    out = args.output_dir.expanduser().resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(
            f"Output directory is not empty: {out}\n"
            "Refusing to overwrite a prior release."
        )
    out.mkdir(parents=True, exist_ok=True)

    errors: list[str] = []

    # ------------------------------------------------------------------
    # 1. Validate authoritative inputs.
    # ------------------------------------------------------------------
    hist_hash = sha256(args.historic_eligible)
    if hist_hash != EXPECTED_HISTORIC_SHA256:
        errors.append(
            "Historic eligible SHA-256 mismatch: "
            f"observed={hist_hash}; expected={EXPECTED_HISTORIC_SHA256}"
        )

    hist = read_tsv(args.historic_eligible)
    new = read_tsv(args.new_neutral_frozen)

    if len(hist) != EXPECTED["historic_n"]:
        errors.append(f"Historic rows: {len(hist)} != {EXPECTED['historic_n']}")
    if hist["nct_id"].nunique() != EXPECTED["historic_n"]:
        errors.append("Historic NCT IDs are not unique")

    if len(new) != EXPECTED["new_n"]:
        errors.append(f"New rows: {len(new)} != {EXPECTED['new_n']}")
    if "nct_id" not in new.columns or new["nct_id"].nunique() != EXPECTED["new_n"]:
        errors.append("New NCT IDs are missing or not unique")

    hist_ids = set(hist["nct_id"].map(lambda x: text(x).upper()))
    new_ids = set(new["nct_id"].map(lambda x: text(x).upper()))
    overlap = sorted(hist_ids & new_ids)
    if overlap:
        errors.append(f"Historic/new NCT overlap detected: {overlap[:20]}")

    if set(hist["final_primary_eligible"].map(text)) != {"YES"}:
        errors.append("Historic eligible file contains non-YES eligibility")

    # Validate expected new-study frozen fields before transformation.
    new_required = {
        "nct_id",
        "final_stratum",
        "final_amr_depth",
        *FINAL_DESCRIPTOR_FIELDS,
        "neutral_adjudication_notes",
        "neutral_adjudicator_initials",
        "neutral_adjudication_status",
    }
    missing_new = sorted(new_required - set(new.columns))
    if missing_new:
        errors.append(f"New frozen packet missing columns: {missing_new}")

    if not missing_new:
        require_exact_counter(
            depth_counter(new["final_amr_depth"]),
            EXPECTED_NEW["depth"],
            "New depth",
            errors,
        )
        require_exact_counter(
            normalized_counter(new["final_stratum"]),
            EXPECTED_NEW["stratum"],
            "New stratum",
            errors,
        )
        require_exact_counter(
            normalized_counter(new["final_h2_comparison_group"]),
            EXPECTED_NEW["h2"],
            "New H2",
            errors,
        )
        require_exact_counter(
            normalized_counter(new["final_clinical_utility_any"]),
            EXPECTED_NEW["utility"],
            "New utility",
            errors,
        )
        require_exact_counter(
            normalized_counter(new["neutral_adjudication_status"]),
            EXPECTED_NEW["status"],
            "New adjudication status",
            errors,
        )

        if set(new["neutral_adjudicator_initials"].map(text)) != {"GB"}:
            errors.append("New neutral adjudicator initials are not uniformly GB")

        blank_final = []
        for _, r in new.iterrows():
            if any(not text(r.get(f, "")) for f in FINAL_DESCRIPTOR_FIELDS):
                blank_final.append(text(r["nct_id"]))
        if blank_final:
            errors.append(
                f"{len(blank_final)} new rows have missing final descriptors: "
                f"{blank_final[:20]}"
            )

        validate_depth_output(new, "NEW", errors)

    p35_ids = read_p35_ids(args.p35_exclusion_ids)
    if p35_ids != EXPECTED_P35_IDS:
        errors.append(
            "Frozen P35 exclusion ID set mismatch: "
            f"observed={sorted(p35_ids)}; expected={sorted(EXPECTED_P35_IDS)}"
        )

    if errors:
        (out / "ANALYTIC_RELEASE_BUILD_ERRORS_v3_3_0.txt").write_text(
            "\n".join(errors) + "\n", encoding="utf-8"
        )
        raise SystemExit(
            "ANALYTIC RELEASE v3.3.0 BUILD: FAIL\n" + "\n".join(errors)
        )

    # ------------------------------------------------------------------
    # 2. Standardize the new 213 rows into historic analytic field names.
    # ------------------------------------------------------------------
    new2 = new.copy()
    new2["nct_id"] = new2["nct_id"].map(lambda x: text(x).upper())
    new2["final_primary_eligible"] = "YES"
    new2["descriptor_adjudication_notes"] = new2[
        "neutral_adjudication_notes"
    ].map(text)
    new2["descriptor_adjudicator_initials"] = new2[
        "neutral_adjudicator_initials"
    ].map(text)
    new2["descriptor_adjudication_status"] = new2[
        "neutral_adjudication_status"
    ].map(text)
    new2["special_sensitivity_flags"] = "NONE"

    # ------------------------------------------------------------------
    # 3. Apply the pre-existing H2 typhoid/plague sensitivity policy
    #    symmetrically to the new Enterobacterales rows.
    # ------------------------------------------------------------------
    rare_audit_rows = []
    for idx, r in new2.iterrows():
        nct = text(r["nct_id"])
        is_entero = text(r["final_h2_comparison_group"]) == "ENTEROBACTERALES"
        blob = frozen_evidence_blob(r)
        matches = []
        if is_entero:
            for label, pat in RARE_H2_PATTERNS:
                if pat.search(blob):
                    matches.append(label)

        apply_flag = bool(matches)
        if apply_flag:
            new2.at[idx, "special_sensitivity_flags"] = add_flag(
                new2.at[idx, "special_sensitivity_flags"],
                "H2_EXCLUDE_TYPHOID_PLAGUE",
            )

        if is_entero:
            rare_audit_rows.append(
                {
                    "nct_id": nct,
                    "final_h2_comparison_group": text(
                        r["final_h2_comparison_group"]
                    ),
                    "matched_prespecified_rare_pathogen_terms": "|".join(matches),
                    "H2_EXCLUDE_TYPHOID_PLAGUE_applied": (
                        "YES" if apply_flag else "NO"
                    ),
                    "audit_basis": (
                        "Frozen registration evidence only; flag applied only "
                        "to new rows already adjudicated ENTEROBACTERALES."
                    ),
                }
            )

    rare_audit = pd.DataFrame(
        rare_audit_rows,
        columns=[
            "nct_id",
            "final_h2_comparison_group",
            "matched_prespecified_rare_pathogen_terms",
            "H2_EXCLUDE_TYPHOID_PLAGUE_applied",
            "audit_basis",
        ],
    )

    # ------------------------------------------------------------------
    # 4. Create a nonbinding H3 near-miss QC queue for the new 213.
    #    No analytic field is changed by this scan.
    # ------------------------------------------------------------------
    h3_rows = []
    for _, r in new2.iterrows():
        blob = frozen_evidence_blob(r)
        mech = sorted(set(m.group(0) for m in H3_MECHANISM_TERMS.finditer(blob)))
        quant = sorted(set(m.group(0) for m in H3_QUANT_TERMS.finditer(blob)))
        if mech and quant:
            h3_rows.append(
                {
                    "nct_id": text(r["nct_id"]),
                    "final_amr_depth": text(r["final_amr_depth"]),
                    "matched_mechanism_terms": "|".join(mech),
                    "matched_quantitative_terms": "|".join(quant),
                    "qc_status": "HEURISTIC_CANDIDATE_ONLY",
                    "analytic_effect": "NONE",
                }
            )

    h3_qc = pd.DataFrame(
        h3_rows,
        columns=[
            "nct_id",
            "final_amr_depth",
            "matched_mechanism_terms",
            "matched_quantitative_terms",
            "qc_status",
            "analytic_effect",
        ],
    )

    # ------------------------------------------------------------------
    # 5. Harmonize columns while preserving every historic value exactly.
    # ------------------------------------------------------------------
    hist_original_columns = list(hist.columns)

    # Add all analytically useful columns from the new frozen packet without
    # importing SB/ZB review-context columns into the final analytic release.
    new_release_allow = [
        "brief_title",
        "official_title",
        "conditions",
        "keywords",
        "intervention_names",
        "intervention_types",
        "primary_outcomes",
        "secondary_outcomes",
        "other_outcomes",
        "summary",
        "clinicaltrials_url",
        "overall_status",
        "start_year",
        "study_type",
        "has_results",
        "enrollment_count",
        "enrollment_type",
        "lead_sponsor_name",
        "lead_sponsor_class",
        "countries",
        "source_registry_shard",
        "screening_final_decision_source",
        "screening_final_decision_basis",
        "final_primary_eligible",
        "final_stratum",
        "final_amr_depth",
        *FINAL_DESCRIPTOR_FIELDS,
        "descriptor_adjudication_notes",
        "descriptor_adjudicator_initials",
        "descriptor_adjudication_status",
        "special_sensitivity_flags",
    ]

    final_columns = list(hist_original_columns)
    for c in new_release_allow:
        if c not in final_columns:
            final_columns.append(c)
    for c in PROVENANCE_FIELDS:
        if c not in final_columns:
            final_columns.append(c)

    hist2 = ensure_columns(hist, final_columns)
    new3 = ensure_columns(new2, final_columns)

    # Provenance is added, not substituted for any historic source column.
    hist2["analytic_release_version"] = VERSION
    hist2["cohort_origin"] = "HISTORIC_V3_2_7"
    hist2["source_descriptor_release"] = "v3.2.7"

    new3["analytic_release_version"] = VERSION
    new3["cohort_origin"] = "NEW_V3_2_9_RESCUE"
    new3["source_descriptor_release"] = "v3.3.0_NEUTRAL_GB"

    for df in [hist2, new3]:
        df["p35_imaging_sensitivity_exclude"] = df["nct_id"].map(
            lambda x: "YES" if text(x).upper() in p35_ids else "NO"
        )
        df["p35_sensitivity_frozen_version"] = "v3.3.0"

    # Force same column order.
    hist2 = hist2[final_columns].copy()
    new3 = new3[final_columns].copy()

    merged = pd.concat([hist2, new3], ignore_index=True)
    merged["nct_id"] = merged["nct_id"].map(lambda x: text(x).upper())
    merged = merged.sort_values("nct_id", kind="stable").reset_index(drop=True)

    # ------------------------------------------------------------------
    # 6. Historic preservation audit.
    # ------------------------------------------------------------------
    hist_source = hist.copy()
    hist_source["nct_id"] = hist_source["nct_id"].map(lambda x: text(x).upper())
    merged_hist = merged[
        merged["cohort_origin"] == "HISTORIC_V3_2_7"
    ].copy()

    hx = hist_source.set_index("nct_id", drop=False)
    mx = merged_hist.set_index("nct_id", drop=False)

    preservation_rows = []
    historic_changed_n = 0
    for nct in sorted(hx.index):
        changed = []
        for c in hist_original_columns:
            if text(hx.at[nct, c]) != text(mx.at[nct, c]):
                changed.append(c)
        if changed:
            historic_changed_n += 1
        preservation_rows.append(
            {
                "nct_id": nct,
                "historic_columns_checked": len(hist_original_columns),
                "changed_historic_columns": "|".join(changed),
                "preserved": "YES" if not changed else "NO",
            }
        )
    preservation = pd.DataFrame(preservation_rows)

    # ------------------------------------------------------------------
    # 7. Full release validation.
    # ------------------------------------------------------------------
    errors = []

    missing_release = sorted(set(REQUIRED_ANALYTIC_FIELDS) - set(merged.columns))
    if missing_release:
        errors.append(f"Merged release missing analytic columns: {missing_release}")

    if len(merged) != EXPECTED["total_n"]:
        errors.append(f"Total rows {len(merged)} != {EXPECTED['total_n']}")
    if merged["nct_id"].nunique() != EXPECTED["total_n"]:
        errors.append("Merged NCT IDs are not unique")
    if set(merged["final_primary_eligible"].map(text)) != {"YES"}:
        errors.append("Merged primary cohort contains non-YES eligibility")

    if historic_changed_n != 0:
        errors.append(
            f"Historic preservation failed for {historic_changed_n} row(s)"
        )

    require_exact_counter(
        depth_counter(merged["final_amr_depth"]),
        EXPECTED["depth"],
        "573 depth",
        errors,
    )
    require_exact_counter(
        normalized_counter(merged["final_h2_comparison_group"]),
        EXPECTED["h2"],
        "573 H2",
        errors,
    )
    require_exact_counter(
        normalized_counter(merged["final_clinical_utility_any"]),
        EXPECTED["utility"],
        "573 utility",
        errors,
    )
    require_exact_counter(
        normalized_counter(merged["descriptor_adjudication_status"]),
        EXPECTED["status"],
        "573 descriptor status",
        errors,
    )

    blank_descriptor_rows = []
    for _, r in merged.iterrows():
        if any(not text(r.get(f, "")) for f in FINAL_DESCRIPTOR_FIELDS):
            blank_descriptor_rows.append(text(r["nct_id"]))
    if blank_descriptor_rows:
        errors.append(
            f"{len(blank_descriptor_rows)} merged rows have missing final descriptors: "
            f"{blank_descriptor_rows[:20]}"
        )

    blank_flags = [
        text(r["nct_id"])
        for _, r in merged.iterrows()
        if not text(r.get("special_sensitivity_flags", ""))
    ]
    if blank_flags:
        errors.append(
            f"{len(blank_flags)} merged rows have blank special_sensitivity_flags"
        )

    validate_depth_output(merged, "MERGED", errors)

    observed_p35 = set(
        merged.loc[
            merged["p35_imaging_sensitivity_exclude"] == "YES", "nct_id"
        ].map(text)
    )
    if observed_p35 != EXPECTED_P35_IDS:
        errors.append(
            f"P35 merged flag set mismatch: {sorted(observed_p35)}"
        )

    p35_n = int((merged["p35_imaging_sensitivity_exclude"] == "YES").sum())
    if p35_n != EXPECTED["p35_excluded_n"]:
        errors.append(
            f"P35 excluded n {p35_n} != {EXPECTED['p35_excluded_n']}"
        )

    p35_sens = merged[
        merged["p35_imaging_sensitivity_exclude"] != "YES"
    ].copy()
    if len(p35_sens) != EXPECTED["p35_remaining_n"]:
        errors.append(
            f"P35 sensitivity cohort n {len(p35_sens)} "
            f"!= {EXPECTED['p35_remaining_n']}"
        )

    # Legacy historic flags must survive exactly.
    legacy_flag_mismatches = []
    hs = hist_source.set_index("nct_id")
    mm = merged_hist.set_index("nct_id")
    if "special_sensitivity_flags" not in hist_source.columns:
        errors.append("Historic release lacks special_sensitivity_flags")
    else:
        for nct in hs.index:
            if text(hs.at[nct, "special_sensitivity_flags"]) != text(
                mm.at[nct, "special_sensitivity_flags"]
            ):
                legacy_flag_mismatches.append(nct)
        if legacy_flag_mismatches:
            errors.append(
                "Historic special_sensitivity_flags changed for: "
                + ", ".join(legacy_flag_mismatches[:20])
            )

    if errors:
        (out / "ANALYTIC_RELEASE_BUILD_ERRORS_v3_3_0.txt").write_text(
            "\n".join(errors) + "\n", encoding="utf-8"
        )
        raise SystemExit(
            "ANALYTIC RELEASE v3.3.0 BUILD: FAIL\n" + "\n".join(errors)
        )

    # ------------------------------------------------------------------
    # 8. Write release artifacts.
    # ------------------------------------------------------------------
    primary_path = out / "eligible_primary_cohort_final_v3_3_0.tsv"
    p35_path = out / "eligible_primary_cohort_P35_sensitivity_v3_3_0.tsv"
    p35_excluded_path = out / "P35_imaging_sensitivity_excluded_records_v3_3_0.tsv"
    preserve_path = out / "historic_preservation_audit_v3_3_0.tsv"
    rare_path = out / "new_H2_rare_pathogen_flag_audit_v3_3_0.tsv"
    h3_path = out / "new_H3_quantitative_near_miss_candidate_QC_v3_3_0.tsv"

    write_tsv(merged, primary_path)
    write_tsv(p35_sens, p35_path)
    write_tsv(
        merged[merged["p35_imaging_sensitivity_exclude"] == "YES"].copy(),
        p35_excluded_path,
    )
    write_tsv(preservation, preserve_path)
    write_tsv(rare_audit, rare_path)
    write_tsv(h3_qc, h3_path)

    new_rare_flagged = sorted(
        new3.loc[
            new3["special_sensitivity_flags"].map(
                lambda x: "H2_EXCLUDE_TYPHOID_PLAGUE" in parse_flags(x)
            ),
            "nct_id",
        ].map(text)
    )

    summary = {
        "created_at": now_utc(),
        "release_version": VERSION,
        "build_pass": True,
        "primary_eligible_n": len(merged),
        "historic_n": int((merged["cohort_origin"] == "HISTORIC_V3_2_7").sum()),
        "new_n": int((merged["cohort_origin"] == "NEW_V3_2_9_RESCUE").sum()),
        "historic_input_sha256": hist_hash,
        "new_neutral_input_sha256": sha256(args.new_neutral_frozen),
        "p35_exclusion_list_sha256": sha256(args.p35_exclusion_ids),
        "historic_rows_with_changed_preexisting_values": historic_changed_n,
        "depth_distribution": depth_counter(merged["final_amr_depth"]),
        "stratum_distribution": normalized_counter(merged["final_stratum"]),
        "h2_group_distribution": normalized_counter(
            merged["final_h2_comparison_group"]
        ),
        "clinical_utility_distribution": normalized_counter(
            merged["final_clinical_utility_any"]
        ),
        "descriptor_status_distribution": normalized_counter(
            merged["descriptor_adjudication_status"]
        ),
        "p35_imaging_sensitivity": {
            "excluded_n": p35_n,
            "remaining_n": len(p35_sens),
            "excluded_nct_ids": sorted(observed_p35),
        },
        "new_h2_rare_pathogen_flagged_n": len(new_rare_flagged),
        "new_h2_rare_pathogen_flagged_nct_ids": new_rare_flagged,
        "new_h3_near_miss_heuristic_candidate_n": len(h3_qc),
        "new_h3_near_miss_note": (
            "Nonbinding QC queue only. No new H3 near-miss analytic flag or depth "
            "is assigned by this release builder."
        ),
        "analysis_gate": (
            "PASS. Freeze this directory read-only before generating or running "
            "the v3.3.1 H1-H4 analysis script."
        ),
    }

    summary_path = out / "analytic_release_summary_v3_3_0.json"
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    readme_path = out / "README_ANALYTIC_RELEASE_v3_3_0.txt"
    readme_path.write_text(
        f"""ClinicalTrials.gov bacterial/AMR diagnostic landscape
Analytic release {VERSION}

Primary cohort
--------------
573 eligible studies:
- 360 frozen historic studies carried forward from analytic release v3.2.7.
- 213 newly eligible studies from the v3.2.9 screening rescue/expansion with
  validated neutral descriptors frozen at v3.3.0.

Historic preservation
---------------------
Every column that existed in the historic v3.2.7 eligible file is required to
remain identical for all 360 historic rows. Added v3.3.0 provenance/P35 columns
do not alter historic scientific values.

P35 sensitivity
---------------
P35 is an analysis-stage exclusion flag only. It does not change eligibility.
Exactly eight frozen NCT IDs are marked YES in
p35_imaging_sensitivity_exclude, yielding a 565-study P35 sensitivity cohort.

Special sensitivity flags
-------------------------
Historic special_sensitivity_flags are preserved exactly. For the new 213,
H2_EXCLUDE_TYPHOID_PLAGUE is added only if a row already adjudicated as
ENTEROBACTERALES explicitly contains a prespecified typhoid/plague term in
the supplied frozen registry evidence. The audit is written separately.

H3 near-miss QC
---------------
new_H3_quantitative_near_miss_candidate_QC_v3_3_0.tsv is heuristic,
descriptive, and nonbinding. It does not change diagnostic depth,
special_sensitivity_flags, or the primary H3 zero-depth-4 analysis.

Analysis boundary
-----------------
Do not edit this release after statistical inspection. After the release
manifest and checksums pass, make the directory read-only and create the
v3.3.1 analysis as a new version using the already-established v3.2.8
inferential hierarchy.
""",
        encoding="utf-8",
    )

    output_files = [
        primary_path,
        p35_path,
        p35_excluded_path,
        preserve_path,
        rare_path,
        h3_path,
        summary_path,
        readme_path,
    ]

    manifest_path = out / "release_manifest_v3_3_0.tsv"
    make_manifest(output_files, manifest_path)
    output_files.append(manifest_path)

    sums_path = out / "SHA256SUMS.txt"
    with sums_path.open("w", encoding="utf-8") as h:
        for p in sorted(output_files, key=lambda x: x.name):
            h.write(f"{sha256(p)}  {p.name}\n")

    print("ANALYTIC RELEASE v3.3.0 BUILD: PASS")
    print(json.dumps(summary, indent=2))
    print(f"Output directory: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
