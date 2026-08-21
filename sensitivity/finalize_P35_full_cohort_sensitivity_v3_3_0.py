#!/usr/bin/env python3
"""
Finalize the P35 full-cohort imaging sensitivity set after central completion.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from collections import Counter

import pandas as pd


INCLUDE = "INCLUDE_P35_IMAGING_SENSITIVITY"
EXCLUDE = "EXCLUDE_NOT_P35_IMAGING_SENSITIVITY"
ALLOWED = {INCLUDE, EXCLUDE}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path, sep="\t", dtype=str, keep_default_na=False,
        encoding="utf-8-sig"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--working-file", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    df = read_tsv(args.working_file)

    if len(df) != 104:
        raise SystemExit(f"Expected 104 candidates; found {len(df)}")
    if df["nct_id"].nunique() != 104:
        raise SystemExit("Duplicate NCT IDs.")

    bad = df.loc[
        ~df["central_p35_sensitivity_final"].isin(ALLOWED),
        ["nct_id", "central_p35_sensitivity_final"],
    ]
    if len(bad):
        raise SystemExit(
            "P35 sensitivity file still contains blank/invalid final "
            "classifications:\n" + bad.to_string(index=False)
        )

    bad = df.loc[
        df["central_p35_sensitivity_basis"].eq("")
        | df["central_adjudicator_initials"].eq("")
        | ~df["p35_sensitivity_status"].isin(
            {"FINAL", "FINAL_BY_FROZEN_P35"}
        ),
        [
            "nct_id",
            "central_p35_sensitivity_basis",
            "central_adjudicator_initials",
            "p35_sensitivity_status",
        ],
    ]
    if len(bad):
        raise SystemExit(
            "P35 sensitivity rows missing final signoff:\n"
            + bad.to_string(index=False)
        )

    if set(df["central_adjudicator_initials"]) != {"CB"}:
        raise SystemExit(
            "Every central P35 sensitivity row must have initials CB."
        )

    final_ids = sorted(
        df.loc[
            df["central_p35_sensitivity_final"].eq(INCLUDE),
            "nct_id",
        ].tolist()
    )

    # These known cases encode the pre-result sensitivity definition.
    for required in ["NCT05285072", "NCT06986512"]:
        if required not in final_ids:
            raise SystemExit(
                f"Known P35 pathogen-directed imaging case missing from "
                f"sensitivity set: {required}"
            )

    for excluded in ["NCT05667207", "NCT05872152"]:
        if excluded in final_ids:
            raise SystemExit(
                f"Non-target modality incorrectly included in sensitivity: "
                f"{excluded}"
            )

    out = args.output_dir
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"Output directory not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)

    final_tsv = out / "P35_full_cohort_imaging_sensitivity_FINAL_v3_3_0.tsv"
    df.to_csv(final_tsv, sep="\t", index=False)

    id_file = out / "P35_imaging_sensitivity_exclusion_nct_ids_v3_3_0.txt"
    id_file.write_text(
        "\n".join(final_ids) + ("\n" if final_ids else ""),
        encoding="utf-8",
    )

    summary = {
        "version": "v3.3.0-P35-imaging-sensitivity-final",
        "candidate_rows": len(df),
        "finalized_rows": len(df),
        "include_n": len(final_ids),
        "exclude_n": int(
            df["central_p35_sensitivity_final"].eq(EXCLUDE).sum()
        ),
        "include_nct_ids": final_ids,
        "classification_counts": dict(
            Counter(df["central_p35_sensitivity_final"])
        ),
        "status_counts": dict(
            Counter(df["p35_sensitivity_status"])
        ),
        "definition": (
            "Sensitivity excludes only already-eligible studies whose "
            "primary qualifying diagnostic is imaging and whose eligibility "
            "depends on the frozen registration explicitly stating that the "
            "imaging signal is bacteria/pathogen-derived. Direct conventional "
            "microscopy and pathogen-derived spectroscopy are not excluded."
        ),
        "primary_cohort_changed": False,
        "frozen_before_v3_3_1": True,
        "input_sha256": {
            "working_file": sha256(args.working_file)
        },
    }

    sj = out / "P35_IMAGING_SENSITIVITY_FREEZE_SUMMARY_v3_3_0.json"
    sj.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    freeze = out / "P35_IMAGING_SENSITIVITY_FREEZE_STATEMENT_v3_3_0.txt"
    freeze.write_text(
        """P35 imaging sensitivity definition is frozen before inspection of
v3.3.1 H1-H4 results.

This sensitivity does not alter the primary eligible cohort. It excludes
only the NCT IDs listed in
P35_imaging_sensitivity_exclusion_nct_ids_v3_3_0.txt for the prespecified
secondary sensitivity analysis.

Reopening requires a documented data-integrity or protocol issue.
""",
        encoding="utf-8",
    )

    manifest = out / "SHA256SUMS.txt"
    with manifest.open("w", encoding="utf-8") as f:
        for p in sorted([final_tsv, id_file, sj, freeze], key=lambda x: x.name):
            f.write(f"{sha256(p)}  {p.name}\n")

    print("P35 FULL-COHORT IMAGING SENSITIVITY FREEZE: PASS")
    print(json.dumps(summary, indent=2))
    print(f"Output: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
