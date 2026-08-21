#!/usr/bin/env python3
"""
Validate the v3.3.0 descriptor-review dispatch after preparation.

This validator is intentionally CSV/TSV-parser based. Physical newline counts
are not row counts because ClinicalTrials.gov text fields may contain embedded
literal newlines inside quoted TSV cells.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import tarfile
import tempfile
from pathlib import Path
from collections import Counter

import pandas as pd


EXPECTED_N = 213

REVIEW_FIELDS = [
    "registry_lookup_performed",
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
    "descriptor_notes",
    "reviewer_initials",
    "review_status",
]

OUTPUT_BY_DEPTH = {
    "0": "ORGANISM_ONLY",
    "1": "BINARY_OR_CATEGORICAL_RESISTANCE",
    "2": "PHENOTYPIC_AST_MIC_ZONE",
    "3": "INTEGRATED_MULTIMECHANISM",
    "4": "QUANTITATIVE_AMR_MECHANISM",
}


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


def verify_internal_archive(archive: Path) -> dict:
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        with tarfile.open(archive, "r:gz") as tf:
            names = tf.getnames()
            tf.extractall(root)

        prohibited = [
            n for n in names
            if re.search(
                r"private|carry_forward|final_decision_ledger|"
                r"screening_v3_2_9_full_review_master",
                n,
                re.I,
            )
        ]
        if prohibited:
            raise SystemExit(
                f"{archive.name}: prohibited content: {prohibited}"
            )

        manifests = list(root.rglob("SHA256SUMS.txt"))
        if len(manifests) != 1:
            raise SystemExit(
                f"{archive.name}: expected one internal SHA256SUMS.txt, "
                f"found {len(manifests)}"
            )

        manifest = manifests[0]
        package_dir = manifest.parent
        checked = 0
        for line in manifest.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            digest, rel = line.split(None, 1)
            rel = rel.strip()
            p = package_dir / rel
            if not p.is_file():
                raise SystemExit(
                    f"{archive.name}: manifest file missing: {rel}"
                )
            observed = sha256(p)
            if observed != digest:
                raise SystemExit(
                    f"{archive.name}: internal hash mismatch for {rel}\n"
                    f"expected={digest}\nobserved={observed}"
                )
            checked += 1

        return {
            "archive": archive.name,
            "archive_sha256": sha256(archive),
            "members": len([n for n in names if not n.endswith("/")]),
            "internal_manifest_entries_verified": checked,
            "prohibited_members": prohibited,
        }


def validate_packet(
    reviewer: str,
    packet: pd.DataFrame,
    source: pd.DataFrame,
) -> dict:
    if len(packet) != EXPECTED_N:
        raise SystemExit(
            f"{reviewer}: expected {EXPECTED_N} parsed rows, "
            f"found {len(packet)}"
        )
    if packet["nct_id"].nunique() != EXPECTED_N:
        raise SystemExit(f"{reviewer}: duplicate NCT IDs.")
    if packet["descriptor_review_id"].nunique() != EXPECTED_N:
        raise SystemExit(f"{reviewer}: duplicate descriptor_review_id.")
    if not packet["descriptor_review_id"].str.startswith(
        f"{reviewer}NEWDR-"
    ).all():
        raise SystemExit(f"{reviewer}: unexpected review ID prefix.")

    source_by = source.set_index("nct_id")
    packet_by = packet.set_index("nct_id")

    if set(packet_by.index) != set(source_by.index):
        missing = sorted(set(source_by.index) - set(packet_by.index))
        extra = sorted(set(packet_by.index) - set(source_by.index))
        raise SystemExit(
            f"{reviewer}: source/packet NCT mismatch. "
            f"missing={missing} extra={extra}"
        )

    for frozen in ["final_stratum", "final_amr_depth"]:
        left = packet_by.loc[source_by.index, frozen]
        right = source_by[frozen]
        if not left.equals(right):
            bad = source_by.index[left != right].tolist()
            raise SystemExit(
                f"{reviewer}: frozen field {frozen} differs: {bad}"
            )

    # All review fields blank except the depth-derived output type.
    for field in REVIEW_FIELDS:
        if field == "index_test_output_type":
            continue
        if packet[field].ne("").any():
            bad = packet.loc[packet[field].ne(""), "nct_id"].tolist()
            raise SystemExit(
                f"{reviewer}: reviewer field {field} unexpectedly prefilled: "
                f"{bad[:10]}"
            )

    expected_output = packet["final_amr_depth"].map(OUTPUT_BY_DEPTH)
    if not packet["index_test_output_type"].equals(expected_output):
        bad = packet.loc[
            packet["index_test_output_type"].ne(expected_output),
            "nct_id",
        ].tolist()
        raise SystemExit(
            f"{reviewer}: output-type prefill inconsistent with depth: {bad}"
        )

    return {
        "reviewer": reviewer,
        "parsed_rows": len(packet),
        "unique_nct_ids": packet["nct_id"].nunique(),
        "unique_descriptor_review_ids": packet[
            "descriptor_review_id"
        ].nunique(),
        "depth_counts": dict(
            Counter(packet["final_amr_depth"])
        ),
        "all_non_output_review_fields_blank": True,
        "output_type_prefill_valid": True,
        "frozen_fields_match_source": True,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prep-dir", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    prep = args.prep_dir
    out = args.output_dir

    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"Output directory not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)

    summary_path = prep / "PREP_SUMMARY_v3_3_0.json"
    source_path = (
        prep
        / "descriptor_source"
        / "newly_eligible_descriptor_source_enriched_v3_3_0.tsv"
    )
    sb_path = (
        prep
        / "reviewer_dispatch"
        / "SB_new_eligible_descriptor_review_v3_3_0"
        / "SB_new_eligible_canonical_descriptor_review_blinded_v3_3_0.tsv"
    )
    zb_path = (
        prep
        / "reviewer_dispatch"
        / "ZB_new_eligible_descriptor_review_v3_3_0"
        / "ZB_new_eligible_canonical_descriptor_review_blinded_v3_3_0.tsv"
    )
    sb_archive = (
        prep / "SB_New_Eligible_Descriptor_Reviewer_Package_v3_3_0.tar.gz"
    )
    zb_archive = (
        prep / "ZB_New_Eligible_Descriptor_Reviewer_Package_v3_3_0.tar.gz"
    )

    for p in [
        summary_path,
        source_path,
        sb_path,
        zb_path,
        sb_archive,
        zb_archive,
    ]:
        if not p.exists():
            raise SystemExit(f"Missing required file: {p}")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    source = read_tsv(source_path)
    sb = read_tsv(sb_path)
    zb = read_tsv(zb_path)

    if len(source) != EXPECTED_N:
        raise SystemExit(
            f"Descriptor source should parse to 213 rows; got {len(source)}"
        )
    if source["nct_id"].nunique() != EXPECTED_N:
        raise SystemExit("Descriptor source has duplicate NCT IDs.")

    sb_result = validate_packet("SB", sb, source)
    zb_result = validate_packet("ZB", zb, source)

    if sb["nct_id"].tolist() == zb["nct_id"].tolist():
        raise SystemExit(
            "SB and ZB packets have identical order; independent shuffles "
            "were expected."
        )

    expected_sb_sha = summary["descriptor_review"]["packages"]["SB"][
        "archive_sha256"
    ]
    expected_zb_sha = summary["descriptor_review"]["packages"]["ZB"][
        "archive_sha256"
    ]
    if sha256(sb_archive) != expected_sb_sha:
        raise SystemExit("SB archive hash differs from PREP_SUMMARY.")
    if sha256(zb_archive) != expected_zb_sha:
        raise SystemExit("ZB archive hash differs from PREP_SUMMARY.")

    sb_archive_result = verify_internal_archive(sb_archive)
    zb_archive_result = verify_internal_archive(zb_archive)

    # Explain the wc -l discrepancy using the parsed source.
    multiline_fields = {}
    for col in source.columns:
        n = int(
            source[col]
            .astype(str)
            .str.contains(r"[\r\n]", regex=True)
            .sum()
        )
        if n:
            multiline_fields[col] = n

    result = {
        "version": "v3.3.0-descriptor-dispatch-validation",
        "status": "PASS",
        "descriptor_source": {
            "parsed_rows": len(source),
            "unique_nct_ids": source["nct_id"].nunique(),
            "depth_counts": dict(
                Counter(source["final_amr_depth"])
            ),
            "fields_containing_embedded_newlines": multiline_fields,
        },
        "SB": sb_result,
        "ZB": zb_result,
        "orders_are_independently_shuffled": True,
        "archives": {
            "SB": sb_archive_result,
            "ZB": zb_archive_result,
        },
        "wc_line_count_interpretation": (
            "Physical newline counts are not TSV record counts because "
            "ClinicalTrials.gov evidence fields contain embedded literal "
            "newlines inside quoted cells. Parsed record count is authoritative."
        ),
        "dispatch_authorized": True,
        "next_gate": (
            "Send SB and ZB their respective tar.gz packages. Do not send "
            "private keys or the other reviewer's packet. In parallel, freeze "
            "the full-cohort P35 sensitivity set before v3.3.1."
        ),
    }

    summary_out = out / "DESCRIPTOR_DISPATCH_VALIDATION_SUMMARY_v3_3_0.json"
    summary_out.write_text(
        json.dumps(result, indent=2) + "\n",
        encoding="utf-8",
    )

    note = out / "MULTILINE_TSV_NOTE_v3_3_0.txt"
    note.write_text(
        """Descriptor TSV physical-line-count clarification

The descriptor source and each blinded reviewer TSV contain 213 parsed
records. `wc -l` reports a much larger number because ClinicalTrials.gov
evidence fields contain embedded literal newline characters inside quoted
TSV cells.

This is valid RFC-style delimited-file behavior and does not represent extra
study rows. Record counts must be checked with a CSV/TSV parser such as
Python's csv module or pandas.read_csv, not with `wc -l`.

No reviewer-facing evidence was modified to remove or normalize these
embedded newlines.
""",
        encoding="utf-8",
    )

    manifest = out / "SHA256SUMS.txt"
    with manifest.open("w", encoding="utf-8") as f:
        for p in sorted([summary_out, note], key=lambda x: x.name):
            f.write(f"{sha256(p)}  {p.name}\n")

    print("V3.3.0 DESCRIPTOR DISPATCH VALIDATION: PASS")
    print(json.dumps(result, indent=2))
    print(f"Output: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
