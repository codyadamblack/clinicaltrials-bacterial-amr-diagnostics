#!/usr/bin/env python3
"""
Validate and freeze the completed v3.3.0 neutral canonical descriptor
adjudication for the 213 newly eligible ClinicalTrials.gov studies.

This validator:
- requires exactly 213 unique NCT IDs;
- requires identical headers and row order versus the original prefilled packet;
- protects all frozen/evidence/reviewer-context fields;
- permits edits only to final descriptor fields (except the depth-derived output
  type, which remains protected) and GB signoff fields;
- validates codebook vocabularies and canonical multiselect ordering;
- validates depth <-> output type, primary modality subset, H2 derivation,
  clinical-utility consistency, and uncertainty/status consistency;
- audits the four GB corrections to previously agreement-prefilled cells;
- writes an immutable full completed packet, a compact harmonized descriptor
  table for the 213 new studies, an audit of corrected prefills, a JSON summary,
  and SHA-256 checksums.

It does NOT reopen eligibility, stratum, diagnostic depth, or P35.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

EXPECTED_N = 213
EXPECTED_DEPTHS = {"0": 177, "1": 23, "2": 13, "3": 0, "4": 0}
EXPECTED_CODEBOOK_SHA256 = (
    "a984f470979f0b914de67a4d989573dd84725f6f74fe9a3801530487a1387d6c"
)

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

GB_EDITABLE_FIELDS = {
    *FINAL_DESCRIPTOR_FIELDS,
    "neutral_adjudication_notes",
    "neutral_adjudicator_initials",
    "neutral_adjudication_status",
}

# The output type is depth-derived and explicitly protected by the protocol.
GB_EDITABLE_FIELDS.remove("final_index_test_output_type")

MULTI_FIELDS = {
    # The v3.2.5 codebook represents all_diagnostic_modalities as an instruction
    # ("Pipe-separated values from primary_diagnostic_modality allowed set"),
    # so validate its tokens/order against the primary-modality vocabulary.
    "final_all_diagnostic_modalities": "primary_diagnostic_modality",
    "final_analytical_endpoint_categories": "analytical_endpoint_categories",
    "final_clinical_utility_endpoint_categories": "clinical_utility_endpoint_categories",
}

SINGLE_CODEBOOK_FIELDS = {
    "final_primary_diagnostic_modality": "primary_diagnostic_modality",
    "final_organism_group": "organism_group",
    "final_gram_group": "gram_group",
    "final_h2_comparison_group": "h2_comparison_group",
    "final_clinical_utility_any": "clinical_utility_any",
    "final_preanalytical_flag": "preanalytical_flag",
    "final_amr_reporting_intervention_flag": "amr_reporting_intervention_flag",
    "final_mixed_viral_bacterial_panel_flag": "mixed_viral_bacterial_panel_flag",
    "final_direct_patient_specimen_flag": "direct_patient_specimen_flag",
    "final_index_test_output_type": "index_test_output_type",
}

DEPTH_OUTPUT = {
    "0": "ORGANISM_ONLY",
    "1": "BINARY_OR_CATEGORICAL_RESISTANCE",
    "2": "PHENOTYPIC_AST_MIC_ZONE",
    "3": "INTEGRATED_MULTIMECHANISM",
    "4": "QUANTITATIVE_AMR_MECHANISM",
}

EXPECTED_PREFILL_CORRECTIONS = {
    (
        "NCT02142933",
        "final_clinical_utility_endpoint_categories",
        "NONE_REGISTERED",
        "COST_RESOURCE_USE",
    ),
    (
        "NCT02142933",
        "final_clinical_utility_any",
        "NO",
        "YES",
    ),
    (
        "NCT03759470",
        "final_all_diagnostic_modalities",
        "CULTURE_OR_MICROSCOPY",
        "CULTURE_OR_MICROSCOPY|ANTIGEN_OR_IMMUNOASSAY",
    ),
    (
        "NCT07675018",
        "final_analytical_endpoint_categories",
        "NONE_REGISTERED",
        "ACCURACY|OTHER",
    ),
}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def text(v: Any) -> str:
    return str(v or "").strip()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def set_field_size_limit() -> None:
    value = sys.maxsize
    while True:
        try:
            csv.field_size_limit(value)
            return
        except OverflowError:
            value //= 10


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    set_field_size_limit()
    with path.open("r", encoding="utf-8-sig", newline="") as h:
        reader = csv.DictReader(h, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_tsv(path: Path, rows, fields) -> None:
    with path.open("w", encoding="utf-8", newline="") as h:
        writer = csv.DictWriter(
            h,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def token_list(v: Any) -> list[str]:
    return [x.strip() for x in text(v).split("|") if x.strip()]


def token_set(v: Any) -> set[str]:
    return set(token_list(v))


def has_uncertainty(row: dict[str, str]) -> bool:
    for field in FINAL_DESCRIPTOR_FIELDS:
        if "UNCERTAIN" in token_set(row.get(field, "")):
            return True
    return False


def load_codebook(path: Path):
    observed_hash = sha256(path)
    if observed_hash != EXPECTED_CODEBOOK_SHA256:
        raise SystemExit(
            "Descriptor codebook SHA-256 mismatch.\n"
            f"Expected: {EXPECTED_CODEBOOK_SHA256}\n"
            f"Observed: {observed_hash}"
        )

    rows, fields = read_tsv(path)
    required = {"field", "position", "allowed_value_or_instruction", "multi_select"}
    if not required <= set(fields):
        raise SystemExit(
            f"Descriptor codebook missing columns: {sorted(required - set(fields))}"
        )

    grouped: dict[str, list[tuple[int, str]]] = {}
    multi: dict[str, str] = {}
    for r in rows:
        f = text(r["field"])
        tok = text(r["allowed_value_or_instruction"])
        try:
            pos = int(text(r["position"]))
        except Exception:
            raise SystemExit(f"Invalid codebook position for {f}: {r['position']!r}")
        grouped.setdefault(f, []).append((pos, tok))
        multi[f] = text(r["multi_select"])

    order = {
        f: [tok for _, tok in sorted(vals)]
        for f, vals in grouped.items()
    }
    return order, multi, observed_hash


def canonical_multiselect(value: str, order: list[str]) -> str:
    tokens = token_list(value)
    rank = {tok: i for i, tok in enumerate(order)}
    unknown = [t for t in tokens if t not in rank]
    if unknown:
        raise ValueError(f"unknown token(s) {unknown}")
    if len(tokens) != len(set(tokens)):
        raise ValueError("duplicate multiselect token")
    return "|".join(sorted(tokens, key=lambda t: rank[t]))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--original-prefilled", required=True, type=Path)
    ap.add_argument("--completed", required=True, type=Path)
    ap.add_argument("--descriptor-codebook", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    for p in [args.original_prefilled, args.completed, args.descriptor_codebook]:
        if not p.exists():
            raise SystemExit(f"File not found: {p}")

    outdir = args.output_dir.expanduser().resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"Output directory is not empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    codebook_order, codebook_multi, codebook_hash = load_codebook(
        args.descriptor_codebook
    )

    before, before_fields = read_tsv(args.original_prefilled)
    after, after_fields = read_tsv(args.completed)

    errors: list[str] = []

    if before_fields != after_fields:
        errors.append("Headers or column order changed")

    if len(before) != EXPECTED_N:
        errors.append(f"Original packet expected {EXPECTED_N} rows; found {len(before)}")
    if len(after) != EXPECTED_N:
        errors.append(f"Completed packet expected {EXPECTED_N} rows; found {len(after)}")

    before_ids = [text(r.get("nct_id")) for r in before]
    after_ids = [text(r.get("nct_id")) for r in after]

    if before_ids != after_ids:
        errors.append("NCT IDs or row order changed")
    if len(set(after_ids)) != EXPECTED_N:
        errors.append("Completed packet NCT IDs are not unique")
    if any(not x for x in after_ids):
        errors.append("Blank NCT ID present")

    # Only final descriptors (other than output type) and GB signoff fields
    # are permitted to change.
    protected_fields = [
        f for f in before_fields
        if f not in GB_EDITABLE_FIELDS
    ]

    prefill_corrections = []
    uncertain_ncts = []

    for idx, (b, a) in enumerate(zip(before, after), start=2):
        nct = text(a.get("nct_id"))
        if text(b.get("nct_id")) != nct:
            continue

        for f in protected_fields:
            if str(b.get(f, "")) != str(a.get(f, "")):
                errors.append(f"{nct}: protected field changed: {f}")

        for f in FINAL_DESCRIPTOR_FIELDS:
            if not text(a.get(f)):
                errors.append(f"{nct}: blank final descriptor: {f}")

        if not text(a.get("neutral_adjudication_notes")):
            errors.append(f"{nct}: neutral_adjudication_notes is blank")
        if text(a.get("neutral_adjudicator_initials")) != "GB":
            errors.append(
                f"{nct}: neutral_adjudicator_initials must be GB, found "
                f"{a.get('neutral_adjudicator_initials')!r}"
            )
        status = text(a.get("neutral_adjudication_status"))
        if status not in {"FINAL", "FINAL_WITH_UNCERTAINTY"}:
            errors.append(f"{nct}: invalid neutral status {status!r}")

        # Codebook validation: single-select fields.
        for final_field, cb_field in SINGLE_CODEBOOK_FIELDS.items():
            value = text(a.get(final_field))
            allowed = codebook_order.get(cb_field)
            if allowed is None:
                errors.append(f"Codebook lacks field {cb_field}")
            elif value not in allowed:
                errors.append(
                    f"{nct}: {final_field} has non-codebook value {value!r}"
                )

        # Codebook validation: multiselect fields and canonical ordering.
        for final_field, cb_field in MULTI_FIELDS.items():
            value = text(a.get(final_field))
            allowed_order = codebook_order.get(cb_field)
            if allowed_order is None:
                errors.append(f"Codebook lacks field {cb_field}")
                continue
            try:
                canon = canonical_multiselect(value, allowed_order)
            except ValueError as exc:
                errors.append(f"{nct}: {final_field}: {exc}")
                continue
            if value != canon:
                errors.append(
                    f"{nct}: {final_field} not in canonical codebook order: "
                    f"{value!r} != {canon!r}"
                )

            toks = token_set(value)
            if "NONE_REGISTERED" in toks and len(toks) != 1:
                errors.append(
                    f"{nct}: NONE_REGISTERED must stand alone in {final_field}"
                )

        # Primary modality must appear in the study's all-modality set.
        primary = text(a.get("final_primary_diagnostic_modality"))
        all_mods = token_set(a.get("final_all_diagnostic_modalities"))
        if primary not in all_mods:
            errors.append(
                f"{nct}: primary modality {primary!r} absent from all modalities"
            )

        # Depth-derived output is frozen and must match the frozen depth.
        depth = text(a.get("final_amr_depth"))
        expected_output = DEPTH_OUTPUT.get(depth)
        observed_output = text(a.get("final_index_test_output_type"))
        if expected_output is None:
            errors.append(f"{nct}: invalid frozen depth {depth!r}")
        elif observed_output != expected_output:
            errors.append(
                f"{nct}: depth/output mismatch {depth} -> {observed_output!r}; "
                f"expected {expected_output!r}"
            )

        # H2 is derived conservatively from organism grouping.
        organism = text(a.get("final_organism_group"))
        h2 = text(a.get("final_h2_comparison_group"))
        if organism == "GRAM_POSITIVE" and h2 != "GRAM_POSITIVE":
            errors.append(f"{nct}: GRAM_POSITIVE organism must map to H2 GRAM_POSITIVE")
        elif organism == "ENTEROBACTERALES" and h2 != "ENTEROBACTERALES":
            errors.append(
                f"{nct}: ENTEROBACTERALES organism must map to H2 ENTEROBACTERALES"
            )
        elif organism == "UNCERTAIN":
            if h2 not in {"UNCERTAIN", "OTHER_EXCLUDED"}:
                errors.append(f"{nct}: unexpected H2 mapping for uncertain organism")
        elif organism not in {"GRAM_POSITIVE", "ENTEROBACTERALES"}:
            if h2 != "OTHER_EXCLUDED":
                errors.append(
                    f"{nct}: non-H2 organism {organism!r} must map to OTHER_EXCLUDED"
                )

        # Basic organism -> Gram consistency.
        gram = text(a.get("final_gram_group"))
        required_gram = {
            "GRAM_POSITIVE": "GRAM_POSITIVE",
            "ENTEROBACTERALES": "GRAM_NEGATIVE",
            "OTHER_GRAM_NEGATIVE": "GRAM_NEGATIVE",
            "MIXED_OR_PAN_BACTERIAL": "MIXED",
            "BACTERIAL_STI": "NOT_APPLICABLE",
            "NOT_SPECIFIED": "NOT_SPECIFIED",
        }.get(organism)
        if required_gram is not None and gram != required_gram:
            errors.append(
                f"{nct}: organism {organism!r} requires gram {required_gram!r}, "
                f"found {gram!r}"
            )

        # clinical_utility_any must agree with the registered utility categories.
        util_tokens = token_set(a.get("final_clinical_utility_endpoint_categories"))
        util_any = text(a.get("final_clinical_utility_any"))
        if util_tokens == {"NONE_REGISTERED"}:
            if util_any != "NO":
                errors.append(
                    f"{nct}: NONE_REGISTERED utility categories require utility_any=NO"
                )
        elif util_tokens == {"UNCERTAIN"}:
            if util_any != "UNCERTAIN":
                errors.append(
                    f"{nct}: UNCERTAIN-only utility category requires utility_any=UNCERTAIN"
                )
        else:
            known_utility = util_tokens - {"UNCERTAIN"}
            if known_utility and util_any != "YES":
                errors.append(
                    f"{nct}: registered utility category requires utility_any=YES"
                )

        # Status must correspond exactly to whether any final descriptor retains
        # a codebook UNCERTAIN token.
        uncertain = has_uncertainty(a)
        if uncertain:
            uncertain_ncts.append(nct)
            if status != "FINAL_WITH_UNCERTAINTY":
                errors.append(
                    f"{nct}: contains UNCERTAIN but status is {status!r}"
                )
        else:
            if status != "FINAL":
                errors.append(
                    f"{nct}: contains no UNCERTAIN but status is {status!r}"
                )

        # Audit modifications to fields that were already prefilled from exact
        # SB/ZB agreement. Blank original final fields were the 612 intended
        # adjudication cells, not prefill corrections.
        for f in FINAL_DESCRIPTOR_FIELDS:
            old = text(b.get(f))
            new = text(a.get(f))
            if old and old != new:
                prefill_corrections.append(
                    {
                        "nct_id": nct,
                        "field": f,
                        "prefilled_value": old,
                        "GB_final_value": new,
                        "neutral_adjudication_notes": text(
                            a.get("neutral_adjudication_notes")
                        ),
                    }
                )

    # Expected frozen depth profile for the 213 new studies.
    depth_counts = Counter(text(r.get("final_amr_depth")) for r in after)
    for d, expected_n in EXPECTED_DEPTHS.items():
        if depth_counts.get(d, 0) != expected_n:
            errors.append(
                f"Frozen depth {d}: expected {expected_n}, observed "
                f"{depth_counts.get(d, 0)}"
            )

    # The completed return was independently inspected before this validator was
    # issued. Protect the documented four audit corrections from accidental loss.
    observed_corrections = {
        (
            r["nct_id"],
            r["field"],
            r["prefilled_value"],
            r["GB_final_value"],
        )
        for r in prefill_corrections
    }
    if observed_corrections != EXPECTED_PREFILL_CORRECTIONS:
        missing = sorted(EXPECTED_PREFILL_CORRECTIONS - observed_corrections)
        extra = sorted(observed_corrections - EXPECTED_PREFILL_CORRECTIONS)
        errors.append(
            "Agreement-prefill correction audit mismatch. "
            f"Missing={missing}; Extra={extra}"
        )

    if errors:
        err_path = outdir / "NEUTRAL_DESCRIPTOR_VALIDATION_ERRORS_v3_3_0.txt"
        err_path.write_text("\n".join(errors) + "\n", encoding="utf-8")
        raise SystemExit(
            "V3.3.0 NEUTRAL DESCRIPTOR VALIDATION: FAIL\n"
            + "\n".join(errors[:40])
            + (f"\n... plus {len(errors)-40} more errors" if len(errors) > 40 else "")
        )

    # Freeze the authoritative full adjudication packet.
    frozen_path = (
        outdir / "neutral_descriptor_adjudication_v3_3_0_completed_frozen.tsv"
    )
    write_tsv(frozen_path, after, after_fields)

    # Compact harmonized 213-study descriptor artifact. The signoff column names
    # intentionally match the historic descriptor-release convention.
    compact_fields = [
        "nct_id",
        "brief_title",
        "final_stratum",
        "final_amr_depth",
        *FINAL_DESCRIPTOR_FIELDS,
        "descriptor_adjudication_notes",
        "descriptor_adjudicator_initials",
        "descriptor_adjudication_status",
    ]
    compact_rows = []
    for r in after:
        x = {
            "nct_id": text(r["nct_id"]),
            "brief_title": text(r.get("brief_title")),
            "final_stratum": text(r.get("final_stratum")),
            "final_amr_depth": text(r.get("final_amr_depth")),
        }
        for f in FINAL_DESCRIPTOR_FIELDS:
            x[f] = text(r.get(f))
        x["descriptor_adjudication_notes"] = text(
            r.get("neutral_adjudication_notes")
        )
        x["descriptor_adjudicator_initials"] = text(
            r.get("neutral_adjudicator_initials")
        )
        x["descriptor_adjudication_status"] = text(
            r.get("neutral_adjudication_status")
        )
        compact_rows.append(x)

    compact_path = outdir / "final_new_canonical_descriptors_v3_3_0.tsv"
    write_tsv(compact_path, compact_rows, compact_fields)

    correction_path = outdir / "prefill_corrections_audit_v3_3_0.tsv"
    write_tsv(
        correction_path,
        prefill_corrections,
        [
            "nct_id",
            "field",
            "prefilled_value",
            "GB_final_value",
            "neutral_adjudication_notes",
        ],
    )

    # Summary distributions.
    distributions = {
        "depth": dict(sorted(Counter(
            text(r["final_amr_depth"]) for r in after
        ).items())),
        "stratum": dict(sorted(Counter(
            text(r["final_stratum"]) for r in after
        ).items())),
        "primary_modality": dict(sorted(Counter(
            text(r["final_primary_diagnostic_modality"]) for r in after
        ).items())),
        "organism_group": dict(sorted(Counter(
            text(r["final_organism_group"]) for r in after
        ).items())),
        "h2_group": dict(sorted(Counter(
            text(r["final_h2_comparison_group"]) for r in after
        ).items())),
        "clinical_utility_any": dict(sorted(Counter(
            text(r["final_clinical_utility_any"]) for r in after
        ).items())),
        "adjudication_status": dict(sorted(Counter(
            text(r["neutral_adjudication_status"]) for r in after
        ).items())),
    }

    def multiselect_other_count(field: str) -> int:
        return sum(
            "OTHER" in token_set(r.get(field, ""))
            for r in after
        )

    summary = {
        "created_at": now_utc(),
        "version": "v3.3.0-neutral-descriptor-validated",
        "validation_pass": True,
        "rows": len(after),
        "unique_nct_ids": len(set(after_ids)),
        "headers_identical": True,
        "row_order_identical": True,
        "protected_fields_identical": True,
        "depth_output_consistent": True,
        "codebook_legal": True,
        "cross_field_logic_errors": 0,
        "prefill_correction_count": len(prefill_corrections),
        "prefill_correction_ncts": sorted(
            {r["nct_id"] for r in prefill_corrections}
        ),
        "uncertainty_row_count": len(uncertain_ncts),
        "uncertainty_nct_ids": uncertain_ncts,
        "other_token_usage": {
            "all_diagnostic_modalities": multiselect_other_count(
                "final_all_diagnostic_modalities"
            ),
            "analytical_endpoint_categories": multiselect_other_count(
                "final_analytical_endpoint_categories"
            ),
            "clinical_utility_endpoint_categories": multiselect_other_count(
                "final_clinical_utility_endpoint_categories"
            ),
        },
        "distributions": distributions,
        "input_sha256": {
            "original_prefilled": sha256(args.original_prefilled),
            "completed": sha256(args.completed),
            "descriptor_codebook": codebook_hash,
        },
        "outputs": {
            "frozen_completed_packet": str(frozen_path),
            "compact_new_descriptors": str(compact_path),
            "prefill_corrections_audit": str(correction_path),
        },
        "next_gate": (
            "Merge these frozen 213 descriptors with the frozen historic 360 "
            "eligible studies and the already-frozen P35 sensitivity set to "
            "construct analytic release v3.3.0. Do not run H1-H4 until that "
            "release passes integrity checks."
        ),
    }

    summary_path = outdir / "neutral_descriptor_validation_summary_v3_3_0.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    readme_path = outdir / "README_FROZEN_NEW_DESCRIPTORS_v3_3_0.txt"
    readme_path.write_text(
        """ClinicalTrials.gov bacterial/AMR diagnostic landscape
Frozen neutral canonical descriptors for 213 newly eligible studies
Version v3.3.0

This directory is the validated descriptor freeze for the 213 studies added
after the v3.2.9 screening-coverage rescue.

Eligibility, BROAD/CORE stratum, AMR depth, and depth-derived index-test output
type remain frozen. Four agreement-prefilled descriptor cells were corrected
by the neutral adjudicator from restored frozen evidence and are recorded in
prefill_corrections_audit_v3_3_0.tsv.

Codebook OTHER values are retained as valid versioned codes. No vocabulary
expansion or recoding is performed at this stage.

Do not overwrite this freeze. Any later genuine correction requires a new
version and an explicit provenance record.
""",
        encoding="utf-8",
    )

    files = [
        frozen_path,
        compact_path,
        correction_path,
        summary_path,
        readme_path,
    ]
    checksum_path = outdir / "SHA256SUMS.txt"
    with checksum_path.open("w", encoding="utf-8") as h:
        for p in sorted(files, key=lambda x: x.name):
            h.write(f"{sha256(p)}  {p.name}\n")

    print("V3.3.0 NEUTRAL DESCRIPTOR VALIDATION: PASS")
    print(json.dumps(summary, indent=2))
    print(f"Output directory: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
