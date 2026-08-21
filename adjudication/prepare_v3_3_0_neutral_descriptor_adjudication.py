#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

EXPECTED_N = 213

CODEBOOK_ORDER = {
    "all_diagnostic_modalities": [
        "CULTURE_OR_MICROSCOPY",
        "ANTIGEN_OR_IMMUNOASSAY",
        "SINGLEPLEX_NAAT_PCR",
        "MULTIPLEX_SYNDROMIC_PANEL",
        "SEQUENCING_OR_METAGENOMICS",
        "MASS_SPECTROMETRY_OR_PROTEOMICS",
        "PHENOTYPIC_AST_OR_MIC",
        "BREATH_VOC_OR_METABOLIC",
        "PREANALYTICAL_SPECIMEN_METHOD",
        "INFORMATICS_OR_ALGORITHM",
        "OTHER",
        "UNCERTAIN",
    ],
    "analytical_endpoint_categories": [
        "ACCURACY",
        "SENSITIVITY_SPECIFICITY",
        "CONCORDANCE",
        "YIELD_INCREMENTAL_DETECTION",
        "TURNAROUND_TIME",
        "LIMIT_OF_DETECTION",
        "REPRODUCIBILITY_PRECISION",
        "FEASIBILITY_USABILITY",
        "NONE_REGISTERED",
        "OTHER",
        "UNCERTAIN",
    ],
    "clinical_utility_endpoint_categories": [
        "TIME_TO_APPROPRIATE_THERAPY",
        "ANTIBIOTIC_ESCALATION_DEESCALATION",
        "ANTIBIOTIC_DURATION_OR_STOP",
        "ANTIBIOTIC_SPECTRUM_OR_TARGETING",
        "CLINICAL_OUTCOME_OTHER",
        "LENGTH_OF_STAY",
        "MORTALITY",
        "COST_RESOURCE_USE",
        "INFECTION_CONTROL_OR_ISOLATION",
        "NONE_REGISTERED",
        "OTHER",
        "UNCERTAIN",
    ],
}

FINAL_DESCRIPTOR_FIELDS = [
    "primary_diagnostic_modality",
    "all_diagnostic_modalities",
    "organism_group",
    "gram_group",
    "h2_comparison_group",
    "analytical_endpoint_categories",
    "clinical_utility_endpoint_categories",
    "clinical_utility_any",
    "preanalytical_flag",
    "amr_reporting_intervention_flag",
    "mixed_viral_bacterial_panel_flag",
    "direct_patient_specimen_flag",
    "index_test_output_type",
]

REVIEWER_CONTEXT_FIELDS = [
    "registry_lookup_performed",
    *FINAL_DESCRIPTOR_FIELDS,
    "descriptor_notes",
    "review_status",
]

PROTECTED_REVIEW_PACKET_FIELDS = [
    "brief_title",
    "official_title",
    "conditions",
    "keywords",
    "intervention_names",
    "intervention_types",
    "primary_outcomes",
    "secondary_outcomes",
    "summary",
    "clinicaltrials_url",
    "final_stratum",
    "final_amr_depth",
    "study_type",
    "overall_status",
    "start_year",
    "has_results",
    "enrollment_count",
    "lead_sponsor_class",
    "countries",
]

DEPTH_OUTPUT = {
    "0": "ORGANISM_ONLY",
    "1": "BINARY_OR_CATEGORICAL_RESISTANCE",
    "2": "PHENOTYPIC_AST_MIC_ZONE",
    "3": "INTEGRATED_MULTIMECHANISM",
    "4": "QUANTITATIVE_AMR_MECHANISM",
}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


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


def canonicalize(field: str, value: str) -> str:
    value = str(value or "").strip()
    if field not in CODEBOOK_ORDER:
        return value
    tokens = [x.strip() for x in value.split("|") if x.strip()]
    order = {tok: i for i, tok in enumerate(CODEBOOK_ORDER[field])}
    unknown = [x for x in tokens if x not in order]
    if unknown:
        raise SystemExit(f"Unknown token(s) for {field}: {unknown}")
    return "|".join(sorted(set(tokens), key=lambda x: order[x]))


def require_unique(df: pd.DataFrame, label: str) -> None:
    if len(df) != EXPECTED_N:
        raise SystemExit(f"{label}: expected {EXPECTED_N} rows, observed {len(df)}")
    if "nct_id" not in df.columns:
        raise SystemExit(f"{label}: nct_id missing")
    if df["nct_id"].nunique() != EXPECTED_N:
        raise SystemExit(f"{label}: NCT IDs are not unique")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sb", required=True, type=Path)
    ap.add_argument("--zb", required=True, type=Path)
    ap.add_argument("--descriptor-source", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    sb = read_tsv(args.sb)
    zb = read_tsv(args.zb)
    src = read_tsv(args.descriptor_source)

    require_unique(sb, "SB")
    require_unique(zb, "ZB")
    require_unique(src, "descriptor source")

    if set(sb["nct_id"]) != set(zb["nct_id"]) or set(sb["nct_id"]) != set(src["nct_id"]):
        raise SystemExit("SB, ZB, and descriptor-source NCT sets differ")

    if "other_outcomes" not in src.columns:
        raise SystemExit(
            "descriptor source lacks other_outcomes; use the frozen enriched source from "
            "v3_3_0_descriptor_stage_preparation/descriptor_source"
        )

    sbx = sb.set_index("nct_id", drop=False)
    zbx = zb.set_index("nct_id", drop=False)
    srcx = src.set_index("nct_id", drop=False)

    # Reviewer packets must be identical on all protected fields they share.
    for field in PROTECTED_REVIEW_PACKET_FIELDS:
        if field not in sb.columns or field not in zb.columns:
            raise SystemExit(f"protected field missing from reviewer packet: {field}")
        bad = [nct for nct in sbx.index if sbx.at[nct, field] != zbx.at[nct, field]]
        if bad:
            raise SystemExit(f"protected-field mismatch in {field}: {bad[:5]}")

    # Depth/output consistency is frozen.
    for nct in sbx.index:
        d = sbx.at[nct, "final_amr_depth"]
        expected = DEPTH_OUTPUT.get(d)
        if expected is None:
            raise SystemExit(f"{nct}: invalid frozen depth {d!r}")
        if sbx.at[nct, "index_test_output_type"] != expected:
            raise SystemExit(f"{nct}: SB output type inconsistent with depth")
        if zbx.at[nct, "index_test_output_type"] != expected:
            raise SystemExit(f"{nct}: ZB output type inconsistent with depth")

    output_rows = []
    agreement_rows = []

    # Use frozen source NCT order for neutral adjudication, not either randomized reviewer order.
    for nct in src["nct_id"].tolist():
        a = sbx.loc[nct]
        b = zbx.loc[nct]
        source = srcx.loc[nct]

        row = {c: source.get(c, "") for c in src.columns}
        row["SB_descriptor_review_id"] = a["descriptor_review_id"]
        row["ZB_descriptor_review_id"] = b["descriptor_review_id"]

        for field in REVIEWER_CONTEXT_FIELDS:
            row[f"SB_{field}"] = a.get(field, "")
            row[f"ZB_{field}"] = b.get(field, "")

        remaining = []
        for field in FINAL_DESCRIPTOR_FIELDS:
            av = canonicalize(field, a.get(field, ""))
            bv = canonicalize(field, b.get(field, ""))
            final_field = f"final_{field}"
            if av == bv:
                row[final_field] = av
                agree = "YES"
            else:
                row[final_field] = ""
                remaining.append(final_field)
                agree = "NO"

            agreement_rows.append(
                {
                    "nct_id": nct,
                    "field": field,
                    "SB_value": av,
                    "ZB_value": bv,
                    "agreement": agree,
                }
            )

        row["remaining_final_fields"] = "|".join(remaining)
        row["remaining_final_field_count"] = str(len(remaining))
        row["adjudication_priority"] = (
            "EXACT_AGREEMENT_AUDIT" if not remaining else "REQUIRED_ADJUDICATION"
        )
        row["reviewer_any_registry_lookup"] = (
            "YES"
            if a.get("registry_lookup_performed", "") == "YES"
            or b.get("registry_lookup_performed", "") == "YES"
            else "NO"
        )
        row["reviewer_any_needs_discussion"] = (
            "YES"
            if a.get("review_status", "") == "NEEDS_DISCUSSION"
            or b.get("review_status", "") == "NEEDS_DISCUSSION"
            else "NO"
        )
        row["neutral_adjudication_notes"] = ""
        row["neutral_adjudicator_initials"] = ""
        row["neutral_adjudication_status"] = ""
        output_rows.append(row)

    packet = pd.DataFrame(output_rows)
    agreement = pd.DataFrame(agreement_rows)

    # Summary metrics.
    field_summary = []
    for field in FINAL_DESCRIPTOR_FIELDS:
        x = agreement[agreement["field"] == field]
        n_agree = int((x["agreement"] == "YES").sum())
        field_summary.append(
            {
                "field": field,
                "records": EXPECTED_N,
                "agreement_n": n_agree,
                "disagreement_n": EXPECTED_N - n_agree,
                "agreement_percent": round(100 * n_agree / EXPECTED_N, 2),
            }
        )
    field_summary_df = pd.DataFrame(field_summary)

    unresolved_cells = int(packet["remaining_final_field_count"].astype(int).sum())
    required_rows = int((packet["adjudication_priority"] == "REQUIRED_ADJUDICATION").sum())
    exact_rows = EXPECTED_N - required_rows

    summary = {
        "version": "v3.3.0-neutral-descriptor-preparation",
        "created_at": now_utc(),
        "rows": EXPECTED_N,
        "unique_nct_ids": EXPECTED_N,
        "protected_reviewer_fields_identical": True,
        "other_outcomes_restored_from_frozen_descriptor_source": True,
        "descriptor_cells_total": EXPECTED_N * len(FINAL_DESCRIPTOR_FIELDS),
        "prefilled_agreement_cells": EXPECTED_N * len(FINAL_DESCRIPTOR_FIELDS) - unresolved_cells,
        "unresolved_descriptor_cells": unresolved_cells,
        "required_adjudication_rows": required_rows,
        "exact_agreement_audit_rows": exact_rows,
        "SB_review_status": dict(Counter(sb["review_status"])),
        "ZB_review_status": dict(Counter(zb["review_status"])),
        "SB_registry_lookup": dict(Counter(sb["registry_lookup_performed"])),
        "ZB_registry_lookup": dict(Counter(zb["registry_lookup_performed"])),
        "input_sha256": {
            "SB": sha256(args.sb),
            "ZB": sha256(args.zb),
            "descriptor_source": sha256(args.descriptor_source),
        },
        "next_gate": (
            "Neutral adjudicator signs off all 213 rows, resolving blank final_* cells from frozen "
            "evidence (including other_outcomes), SB/ZB interpretations, and binding descriptor policies."
        ),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    packet_path = args.output_dir / "neutral_descriptor_adjudication_packet_v3_3_0_prefilled.tsv"
    agreement_path = args.output_dir / "descriptor_cell_comparison_v3_3_0.tsv"
    field_summary_path = args.output_dir / "descriptor_field_agreement_v3_3_0.tsv"
    summary_path = args.output_dir / "neutral_descriptor_preparation_summary_v3_3_0.json"
    instructions_path = args.output_dir / "README_NEUTRAL_DESCRIPTOR_ADJUDICATION_v3_3_0.txt"

    packet.to_csv(packet_path, sep="\t", index=False, lineterminator="\n")
    agreement.to_csv(agreement_path, sep="\t", index=False, lineterminator="\n")
    field_summary_df.to_csv(field_summary_path, sep="\t", index=False, lineterminator="\n")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    instructions_text = """Neutral canonical descriptor adjudication v3.3.0

Purpose: finalize descriptors for the 213 newly eligible studies without reopening eligibility, BROAD/CORE stratum, or frozen AMR depth.

All 213 rows require signoff. Exact reviewer agreements are prefilled; disagreements are blank. Review frozen evidence first, including other_outcomes restored from the source file. SB and ZB values/notes are competing interpretations, not authoritative labels.

Binding clarifications:
1. AMR-reporting intervention = YES only when the evaluated intervention changes notification, selective/cascade reporting, stewardship interpretation, or communication of AST/AMR results; an assay merely producing a resistance result is not sufficient.
2. Do not infer assay chemistry, organism coverage, specimen matrix, or AMR capability from product knowledge.
3. Taxonomic Enterobacterales includes Salmonella and Yersinia, subject to the prespecified rare-pathogen sensitivity policy.
4. direct_patient_specimen = YES if any registered index-test stream is directly on a patient specimen; positive blood-culture broth is not direct.
5. clinical_utility_any = YES only when a registered endpoint evaluates patient management, treatment, clinical outcome, infection control, or resource use.
6. Never edit final_stratum, final_amr_depth, or the depth-derived output type.

Complete neutral_adjudication_notes, neutral_adjudicator_initials, and neutral_adjudication_status for every row.
"""
    instructions_path.write_text(instructions_text, encoding="utf-8")

    manifest_files = [packet_path, agreement_path, field_summary_path, summary_path, instructions_path]
    manifest = args.output_dir / "SHA256SUMS.txt"
    with manifest.open("w", encoding="utf-8") as f:
        for p in sorted(manifest_files, key=lambda x: x.name):
            f.write(f"{sha256(p)}  {p.name}\n")

    print("V3.3.0 NEUTRAL DESCRIPTOR PREPARATION: PASS")
    print(json.dumps(summary, indent=2))
    print(f"Packet: {packet_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
