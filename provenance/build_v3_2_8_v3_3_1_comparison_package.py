#!/usr/bin/env python3
"""
Build the frozen v3.2.8 -> v3.3.1 comparison and manuscript-transition package
for the ClinicalTrials.gov bacterial/AMR diagnostic landscape.

THIS SCRIPT DOES NOT RUN OR REFIT ANY STATISTICAL MODEL.

It:
1. verifies the frozen v3.2.8 and v3.3.1 analysis manifests;
2. verifies the frozen 360-study v3.2.7 and 573-study v3.3.0 releases;
3. extracts already-computed H1-H4 results from both analysis versions;
4. quantifies the 213-study cohort expansion and its depth/era composition;
5. compares the v3.3.1 primary 573-study analysis with the frozen P35
   565-study sensitivity analysis;
6. writes a versioned comparison table and an explicit manuscript-claim
   transition memo;
7. creates a manifest and SHA-256 checksums.

No eligibility, descriptor, diagnostic-depth, P35, sensitivity flag,
statistical model, estimate, confidence interval, or p-value is recomputed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


COMPARISON_VERSION = "v3.2.8-vs-v3.3.1-comparison"
COMPARISON_TAG = "v3_2_8_vs_v3_3_1"

EXPECTED_OLD_RELEASE_SHA256 = (
    "a59a9ec30d188533c2e4508bb8044150fb86e563208273cf10c4563b1543bda6"
)
EXPECTED_NEW_RELEASE_SHA256 = (
    "686a99d7e33b78822c7e402478589cfec086d8a58b7156468b43e7aaa609c4b2"
)

EXPECTED_OLD_N = 360
EXPECTED_NEW_N = 573
EXPECTED_ADDED_N = 213
EXPECTED_P35_N = 565

EXPECTED_OLD_DEPTH = {"0": 258, "1": 60, "2": 39, "3": 3}
EXPECTED_NEW_DEPTH = {"0": 435, "1": 83, "2": 52, "3": 3}
EXPECTED_ADDED_DEPTH = {"0": 177, "1": 23, "2": 13}

EXPECTED_NEW_P35_IDS = {
    "NCT01378728",
    "NCT02450942",
    "NCT02491164",
    "NCT02558062",
    "NCT03091361",
    "NCT03290690",
    "NCT05285072",
    "NCT06986512",
}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def text(v: Any) -> str:
    if pd.isna(v):
        return ""
    return str(v).strip()


def num(v: Any) -> float:
    s = text(v)
    return float(s) if s else math.nan


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


def require(cond: bool, message: str) -> None:
    if not cond:
        raise SystemExit(message)


def verify_sha_manifest(root: Path) -> dict[str, Any]:
    manifest = root / "SHA256SUMS.txt"
    require(manifest.exists(), f"Missing checksum manifest: {manifest}")

    checked = 0
    failures = []
    for raw in manifest.read_text(encoding="utf-8-sig").splitlines():
        raw = raw.strip()
        if not raw:
            continue
        try:
            expected, rel = raw.split(maxsplit=1)
        except ValueError:
            failures.append(f"Malformed checksum line: {raw!r}")
            continue
        rel = rel.lstrip("* ")
        p = root / rel
        if not p.exists():
            failures.append(f"Missing: {rel}")
            continue
        observed = sha256(p)
        if observed != expected:
            failures.append(
                f"SHA mismatch {rel}: observed={observed}, expected={expected}"
            )
        checked += 1

    require(
        not failures,
        "Checksum verification failed:\n" + "\n".join(failures[:30]),
    )
    return {
        "root": str(root),
        "manifest_sha256": sha256(manifest),
        "files_checked": checked,
        "failures": 0,
    }


def select_row(df: pd.DataFrame, column: str, value: str) -> pd.Series:
    m = df[column].map(text) == value
    require(
        int(m.sum()) == 1,
        f"Expected exactly one row where {column}={value!r}; found {int(m.sum())}",
    )
    return df.loc[m].iloc[0]


def select_contains(df: pd.DataFrame, column: str, phrase: str) -> pd.Series:
    m = df[column].map(text).str.contains(phrase, regex=False)
    require(
        int(m.sum()) == 1,
        f"Expected exactly one row containing {phrase!r} in {column}; "
        f"found {int(m.sum())}",
    )
    return df.loc[m].iloc[0]


def counter_nozero(series: pd.Series) -> dict[str, int]:
    c = Counter(series.map(text))
    return dict(sorted((k, v) for k, v in c.items() if k))


def depth_counts(df: pd.DataFrame) -> dict[str, int]:
    c = Counter(df["final_amr_depth"].map(text))
    return {str(d): c.get(str(d), 0) for d in range(5) if c.get(str(d), 0)}


def era(y: Any) -> str:
    try:
        yy = int(float(text(y)))
    except Exception:
        return "MISSING"
    if yy <= 2009:
        return "<=2009"
    if yy <= 2014:
        return "2010-2014"
    if yy <= 2019:
        return "2015-2019"
    return ">=2020"


def fmt(x: Any, digits: int = 3) -> str:
    try:
        y = float(x)
    except Exception:
        return text(x)
    if math.isnan(y):
        return "NA"
    if abs(y) < 0.001 and y != 0:
        return f"{y:.2e}"
    return f"{y:.{digits}f}"


def significance(p: float) -> str:
    if math.isnan(p):
        return "NA"
    return "p<0.05" if p < 0.05 else "p>=0.05"


def get_analysis_paths(old: Path, new: Path) -> dict[str, Path]:
    n = new / "primary_573"
    p35 = new / "P35_sensitivity_565"
    return {
        # H1
        "old_h1": old / "H1_primary_and_threshold_models_v3_2_8.tsv",
        "new_h1": n / "H1_primary_and_threshold_models_v3_3_1.tsv",
        "old_h1_diag": old / "H1_proportional_odds_diagnostic_v3_2_8.tsv",
        "new_h1_diag": n / "H1_proportional_odds_diagnostic_v3_3_1.tsv",
        "old_h1_sens": old / "H1_sensitivity_models_v3_2_8.tsv",
        "new_h1_sens": n / "H1_sensitivity_models_v3_3_1.tsv",
        "old_h1_era": old / "H1_depth_by_era_v3_2_8.tsv",
        "new_h1_era": n / "H1_depth_by_era_v3_3_1.tsv",
        # H2
        "old_h2_perm": old / "H2_primary_ordinal_permutation_v3_2_8.tsv",
        "new_h2_perm": n / "H2_primary_ordinal_permutation_v3_3_1.tsv",
        "old_h2_exact": old / "H2_primary_threshold_fisher_v3_2_8.tsv",
        "new_h2_exact": n / "H2_primary_threshold_fisher_v3_3_1.tsv",
        "old_h2_adj": old / "H2_year_adjusted_threshold_models_v3_2_8.tsv",
        "new_h2_adj": n / "H2_year_adjusted_threshold_models_v3_3_1.tsv",
        "old_h2_diag": old / "H2_proportional_odds_diagnostic_v3_2_8.tsv",
        "new_h2_diag": n / "H2_proportional_odds_diagnostic_v3_3_1.tsv",
        "old_h2_rare_perm": old / "H2_rare_pathogen_ordinal_permutation_v3_2_8.tsv",
        "new_h2_rare_perm": n / "H2_rare_pathogen_ordinal_permutation_v3_3_1.tsv",
        "old_h2_rare_exact": old / "H2_rare_pathogen_threshold_fisher_v3_2_8.tsv",
        "new_h2_rare_exact": n / "H2_rare_pathogen_threshold_fisher_v3_3_1.tsv",
        "old_h2_modperm": old / "H2_primary_modality_permutation_v3_2_8.tsv",
        "new_h2_modperm": n / "H2_primary_modality_permutation_v3_3_1.tsv",
        "old_h2_modtok": old / "H2_modality_token_fisher_FDR_v3_2_8.tsv",
        "new_h2_modtok": n / "H2_modality_token_fisher_FDR_v3_3_1.tsv",
        # H3
        "old_h3": old / "H3_depth4_exact_bound_v3_2_8.tsv",
        "new_h3": n / "H3_depth4_exact_bound_v3_3_1.tsv",
        "old_h3_records": old / "H3_depth3_and_quantitative_near_miss_records_v3_2_8.tsv",
        "new_h3_records": n / "H3_depth3_and_frozen_quantitative_near_miss_records_v3_3_1.tsv",
        # H4
        "old_h4_depth": old / "H4_utility_by_depth_v3_2_8.tsv",
        "new_h4_depth": n / "H4_utility_by_depth_v3_3_1.tsv",
        "old_h4_perm": old / "H4_primary_categorical_permutation_v3_2_8.tsv",
        "new_h4_perm": n / "H4_primary_categorical_permutation_v3_3_1.tsv",
        "old_h4_pair": old / "H4_pairwise_depth_fisher_v3_2_8.tsv",
        "new_h4_pair": n / "H4_pairwise_depth_fisher_v3_3_1.tsv",
        "old_h4_stratum": old / "H4_stratum_models_v3_2_8.tsv",
        "new_h4_stratum": n / "H4_stratum_models_v3_3_1.tsv",
        "old_h4_core": old / "H4_core_AMR_depth_fisher_v3_2_8.tsv",
        "new_h4_core": n / "H4_core_AMR_depth_fisher_v3_3_1.tsv",
        "old_h4_shape": old / "H4_depth_shape_diagnostic_v3_2_8.tsv",
        "new_h4_shape": n / "H4_depth_shape_diagnostic_v3_3_1.tsv",
        # P35
        "p35_main": p35 / "TABLE_MAIN_HYPOTHESIS_NUMERIC_RESULTS_v3_3_1.tsv",
    }


def required_files_exist(paths: dict[str, Path]) -> None:
    missing = [str(p) for p in paths.values() if not p.exists()]
    require(not missing, "Missing required comparison inputs:\n" + "\n".join(missing))


def old_new_result_row(
    hypothesis: str,
    analysis: str,
    hierarchy: str,
    metric: str,
    old_effect: float,
    old_ci_low: float,
    old_ci_high: float,
    old_p: float,
    new_effect: float,
    new_ci_low: float,
    new_ci_high: float,
    new_p: float,
    assessment: str,
    manuscript_action: str,
) -> dict[str, Any]:
    return {
        "hypothesis": hypothesis,
        "analysis": analysis,
        "hierarchy": hierarchy,
        "effect_metric": metric,
        "v3_2_8_effect": old_effect,
        "v3_2_8_ci_low": old_ci_low,
        "v3_2_8_ci_high": old_ci_high,
        "v3_2_8_p_value": old_p,
        "v3_2_8_significance": significance(old_p),
        "v3_3_1_effect": new_effect,
        "v3_3_1_ci_low": new_ci_low,
        "v3_3_1_ci_high": new_ci_high,
        "v3_3_1_p_value": new_p,
        "v3_3_1_significance": significance(new_p),
        "absolute_effect_change": (
            new_effect - old_effect
            if not math.isnan(old_effect) and not math.isnan(new_effect)
            else math.nan
        ),
        "assessment": assessment,
        "manuscript_action": manuscript_action,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--old-analysis", required=True, type=Path)
    ap.add_argument("--new-analysis", required=True, type=Path)
    ap.add_argument("--old-release", required=True, type=Path)
    ap.add_argument("--new-release", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    for p in [
        args.old_analysis,
        args.new_analysis,
        args.old_release,
        args.new_release,
    ]:
        require(p.exists(), f"Input does not exist: {p}")

    out = args.output_dir.expanduser().resolve()
    require(
        not out.exists() or not any(out.iterdir()),
        f"Output directory is not empty: {out}\nRefusing to overwrite.",
    )
    out.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Frozen-input validation.
    # ------------------------------------------------------------------
    old_manifest = verify_sha_manifest(args.old_analysis)
    new_manifest = verify_sha_manifest(args.new_analysis)

    require(
        sha256(args.old_release) == EXPECTED_OLD_RELEASE_SHA256,
        "Historic v3.2.7 release hash mismatch.",
    )
    require(
        sha256(args.new_release) == EXPECTED_NEW_RELEASE_SHA256,
        "Final v3.3.0 release hash mismatch.",
    )

    old_summary_path = args.old_analysis / "FINAL_ANALYSIS_RUN_SUMMARY_v3_2_8.json"
    new_summary_path = args.new_analysis / "FINAL_ANALYSIS_RUN_SUMMARY_v3_3_1.json"
    require(old_summary_path.exists(), f"Missing: {old_summary_path}")
    require(new_summary_path.exists(), f"Missing: {new_summary_path}")

    old_summary = json.loads(old_summary_path.read_text(encoding="utf-8"))
    new_summary = json.loads(new_summary_path.read_text(encoding="utf-8"))

    require(old_summary.get("analysis_version") == "v3.2.8", "Wrong old analysis version.")
    require(new_summary.get("analysis_version") == "v3.3.1", "Wrong new analysis version.")
    require(bool(old_summary.get("analysis_pass")), "v3.2.8 analysis did not PASS.")
    require(bool(new_summary.get("analysis_pass")), "v3.3.1 analysis did not PASS.")

    old_release = read_tsv(args.old_release)
    new_release = read_tsv(args.new_release)
    require(len(old_release) == EXPECTED_OLD_N, "Old release must have 360 rows.")
    require(len(new_release) == EXPECTED_NEW_N, "New release must have 573 rows.")
    require(old_release["nct_id"].nunique() == EXPECTED_OLD_N, "Old NCT IDs not unique.")
    require(new_release["nct_id"].nunique() == EXPECTED_NEW_N, "New NCT IDs not unique.")

    old_ids = set(old_release["nct_id"].map(text))
    new_ids = set(new_release["nct_id"].map(text))
    added_ids = new_ids - old_ids
    lost_ids = old_ids - new_ids
    require(len(added_ids) == EXPECTED_ADDED_N, f"Expected 213 additions; found {len(added_ids)}.")
    require(not lost_ids, f"Historic eligible records lost from v3.3.0: {sorted(lost_ids)}")
    require(depth_counts(old_release) == EXPECTED_OLD_DEPTH, "Old depth counts mismatch.")
    require(depth_counts(new_release) == EXPECTED_NEW_DEPTH, "New depth counts mismatch.")

    added = new_release[new_release["nct_id"].isin(added_ids)].copy()
    require(depth_counts(added) == EXPECTED_ADDED_DEPTH, "Added-study depth counts mismatch.")

    p35_ids = set(
        new_release.loc[
            new_release["p35_imaging_sensitivity_exclude"].map(text) == "YES",
            "nct_id",
        ].map(text)
    )
    require(p35_ids == EXPECTED_NEW_P35_IDS, "P35 exclusion set mismatch.")

    paths = get_analysis_paths(args.old_analysis, args.new_analysis)
    required_files_exist(paths)

    # ------------------------------------------------------------------
    # Cohort-expansion audit.
    # ------------------------------------------------------------------
    for frame in [old_release, new_release, added]:
        frame["comparison_era"] = frame["start_year"].map(era)

    cohort_rows = []
    cohort_rows.append({
        "domain": "COHORT",
        "category": "Eligible studies",
        "v3_2_8_n": len(old_release),
        "v3_3_1_n": len(new_release),
        "added_n": len(added),
        "absolute_change": len(new_release) - len(old_release),
    })

    for depth in ["0", "1", "2", "3", "4"]:
        o = int((old_release["final_amr_depth"] == depth).sum())
        n = int((new_release["final_amr_depth"] == depth).sum())
        a = int((added["final_amr_depth"] == depth).sum())
        cohort_rows.append({
            "domain": "DEPTH",
            "category": f"Depth {depth}",
            "v3_2_8_n": o,
            "v3_3_1_n": n,
            "added_n": a,
            "absolute_change": n - o,
        })

    for stratum in sorted(set(new_release["final_stratum"].map(text))):
        o = int((old_release["final_stratum"].map(text) == stratum).sum())
        n = int((new_release["final_stratum"].map(text) == stratum).sum())
        a = int((added["final_stratum"].map(text) == stratum).sum())
        cohort_rows.append({
            "domain": "STRATUM",
            "category": stratum,
            "v3_2_8_n": o,
            "v3_3_1_n": n,
            "added_n": a,
            "absolute_change": n - o,
        })

    for utility in ["NO", "YES"]:
        o = int((old_release["final_clinical_utility_any"].map(text) == utility).sum())
        n = int((new_release["final_clinical_utility_any"].map(text) == utility).sum())
        a = int((added["final_clinical_utility_any"].map(text) == utility).sum())
        cohort_rows.append({
            "domain": "CLINICAL_UTILITY",
            "category": utility,
            "v3_2_8_n": o,
            "v3_3_1_n": n,
            "added_n": a,
            "absolute_change": n - o,
        })

    for group in ["GRAM_POSITIVE", "ENTEROBACTERALES", "OTHER_EXCLUDED"]:
        o = int((old_release["final_h2_comparison_group"].map(text) == group).sum())
        n = int((new_release["final_h2_comparison_group"].map(text) == group).sum())
        a = int((added["final_h2_comparison_group"].map(text) == group).sum())
        cohort_rows.append({
            "domain": "H2_GROUP",
            "category": group,
            "v3_2_8_n": o,
            "v3_3_1_n": n,
            "added_n": a,
            "absolute_change": n - o,
        })

    cohort_compare = pd.DataFrame(cohort_rows)
    write_tsv(
        cohort_compare,
        out / f"cohort_expansion_summary_{COMPARISON_TAG}.tsv",
    )

    added_era_depth = (
        added.groupby(["comparison_era", "final_amr_depth"])
        .size()
        .reset_index(name="n")
    )
    era_totals = added.groupby("comparison_era").size().to_dict()
    added_era_depth["era_n"] = added_era_depth["comparison_era"].map(era_totals)
    added_era_depth["percent_within_added_era"] = (
        100 * added_era_depth["n"] / added_era_depth["era_n"]
    )
    era_order = ["<=2009", "2010-2014", "2015-2019", ">=2020", "MISSING"]
    added_era_depth["comparison_era"] = pd.Categorical(
        added_era_depth["comparison_era"],
        categories=era_order,
        ordered=True,
    )
    added_era_depth = added_era_depth.sort_values(
        ["comparison_era", "final_amr_depth"]
    )
    write_tsv(
        added_era_depth,
        out / f"added_213_by_era_and_depth_{COMPARISON_TAG}.tsv",
    )

    recent_added = added[added["comparison_era"] == ">=2020"]
    recent_added_n = len(recent_added)
    recent_added_d0 = int((recent_added["final_amr_depth"] == "0").sum())
    recent_added_d0_pct = (
        100 * recent_added_d0 / recent_added_n if recent_added_n else math.nan
    )

    # ------------------------------------------------------------------
    # H1 extraction.
    # ------------------------------------------------------------------
    old_h1 = read_tsv(paths["old_h1"])
    new_h1 = read_tsv(paths["new_h1"])
    old_h1_diag = read_tsv(paths["old_h1_diag"])
    new_h1_diag = read_tsv(paths["new_h1_diag"])

    oh1 = select_contains(old_h1, "analysis", "primary proportional odds")
    nh1 = select_contains(new_h1, "analysis", "primary proportional odds")
    oh1_ge1 = select_contains(old_h1, "analysis", "depth>=1")
    nh1_ge1 = select_contains(new_h1, "analysis", "depth>=1")
    oh1_ge2 = select_contains(old_h1, "analysis", "depth>=2")
    nh1_ge2 = select_contains(new_h1, "analysis", "depth>=2")
    oh1d = old_h1_diag.iloc[0]
    nh1d = new_h1_diag.iloc[0]

    # ------------------------------------------------------------------
    # H2 extraction.
    # ------------------------------------------------------------------
    old_h2p = read_tsv(paths["old_h2_perm"]).iloc[0]
    new_h2p = read_tsv(paths["new_h2_perm"]).iloc[0]
    old_h2e = read_tsv(paths["old_h2_exact"])
    new_h2e = read_tsv(paths["new_h2_exact"])
    old_h2a = read_tsv(paths["old_h2_adj"])
    new_h2a = read_tsv(paths["new_h2_adj"])
    old_h2d = read_tsv(paths["old_h2_diag"]).iloc[0]
    new_h2d = read_tsv(paths["new_h2_diag"]).iloc[0]
    old_h2rp = read_tsv(paths["old_h2_rare_perm"]).iloc[0]
    new_h2rp = read_tsv(paths["new_h2_rare_perm"]).iloc[0]
    old_h2re = read_tsv(paths["old_h2_rare_exact"])
    new_h2re = read_tsv(paths["new_h2_rare_exact"])
    old_h2mp = read_tsv(paths["old_h2_modperm"]).iloc[0]
    new_h2mp = read_tsv(paths["new_h2_modperm"]).iloc[0]

    oh2_ge1 = select_row(old_h2e, "outcome", "depth_ge1")
    nh2_ge1 = select_row(new_h2e, "outcome", "depth_ge1")
    oh2_ge2 = select_row(old_h2e, "outcome", "depth_ge2")
    nh2_ge2 = select_row(new_h2e, "outcome", "depth_ge2")
    oh2_util = select_row(old_h2e, "outcome", "utility_yes")
    nh2_util = select_row(new_h2e, "outcome", "utility_yes")
    oh2_adj2 = select_contains(old_h2a, "analysis", "depth>=2")
    nh2_adj2 = select_contains(new_h2a, "analysis", "depth>=2")

    # ------------------------------------------------------------------
    # H3 extraction.
    # ------------------------------------------------------------------
    old_h3 = read_tsv(paths["old_h3"]).iloc[0]
    new_h3 = read_tsv(paths["new_h3"]).iloc[0]

    # ------------------------------------------------------------------
    # H4 extraction.
    # ------------------------------------------------------------------
    old_h4d = read_tsv(paths["old_h4_depth"])
    new_h4d = read_tsv(paths["new_h4_depth"])
    old_h4p = read_tsv(paths["old_h4_perm"]).iloc[0]
    new_h4p = read_tsv(paths["new_h4_perm"]).iloc[0]
    old_h4pair = read_tsv(paths["old_h4_pair"])
    new_h4pair = read_tsv(paths["new_h4_pair"])
    old_h4str = read_tsv(paths["old_h4_stratum"])
    new_h4str = read_tsv(paths["new_h4_stratum"])
    old_h4core = read_tsv(paths["old_h4_core"])
    new_h4core = read_tsv(paths["new_h4_core"])
    old_h4shape = read_tsv(paths["old_h4_shape"]).iloc[0]
    new_h4shape = read_tsv(paths["new_h4_shape"]).iloc[0]

    oh4_21 = select_contains(old_h4pair, "analysis", "depth 2 vs depth 1")
    nh4_21 = select_contains(new_h4pair, "analysis", "depth 2 vs depth 1")
    oh4_core_broad = select_contains(
        old_h4str, "analysis", "structural exact: CORE_AMR"
    )
    nh4_core_broad = select_contains(
        new_h4str, "analysis", "structural exact: CORE_AMR"
    )
    oh4_core21 = select_contains(old_h4core, "analysis", "depth 2 vs depth 1")
    nh4_core21 = select_contains(new_h4core, "analysis", "depth 2 vs depth 1")

    # ------------------------------------------------------------------
    # Main old-vs-new comparison table.
    # ------------------------------------------------------------------
    comparisons = [
        old_new_result_row(
            "H1",
            "Diagnostic depth per 5-year increase in start year",
            "PRIMARY",
            "OR",
            num(oh1["odds_ratio"]),
            num(oh1["ci_low"]),
            num(oh1["ci_high"]),
            num(oh1["p_value"]),
            num(nh1["odds_ratio"]),
            num(nh1["ci_low"]),
            num(nh1["ci_high"]),
            num(nh1["p_value"]),
            "MATERIAL_CHANGE",
            "REVISE: do not claim statistically supported monotonic temporal increase.",
        ),
        old_new_result_row(
            "H1",
            "Cumulative-threshold slope interaction (>=1 vs >=2)",
            "DIAGNOSTIC",
            "interaction_OR_ratio",
            num(oh1d["interaction_or_ratio"]),
            math.nan,
            math.nan,
            num(oh1d["p_value"]),
            num(nh1d["interaction_or_ratio"]),
            math.nan,
            math.nan,
            num(nh1d["p_value"]),
            "MATERIAL_CHANGE",
            "REVISE: threshold-specific temporal slopes now differ; emphasize >=1 vs >=2 results.",
        ),
        old_new_result_row(
            "H1",
            "Depth >=1 per 5-year increase",
            "SUPPORTING",
            "OR",
            num(oh1_ge1["odds_ratio"]),
            num(oh1_ge1["ci_low"]),
            num(oh1_ge1["ci_high"]),
            num(oh1_ge1["p_value"]),
            num(nh1_ge1["odds_ratio"]),
            num(nh1_ge1["ci_low"]),
            num(nh1_ge1["ci_high"]),
            num(nh1_ge1["p_value"]),
            "ATTENUATED",
            "Report as modest/nonsignificant supporting signal.",
        ),
        old_new_result_row(
            "H1",
            "Depth >=2 per 5-year increase",
            "SUPPORTING",
            "OR",
            num(oh1_ge2["odds_ratio"]),
            num(oh1_ge2["ci_low"]),
            num(oh1_ge2["ci_high"]),
            num(oh1_ge2["p_value"]),
            num(nh1_ge2["odds_ratio"]),
            num(nh1_ge2["ci_low"]),
            num(nh1_ge2["ci_high"]),
            num(nh1_ge2["p_value"]),
            "NO_TEMPORAL_DEEPENING",
            "Report no evidence of increasing phenotypic-or-deeper resolution.",
        ),
        old_new_result_row(
            "H2",
            "Enterobacterales vs Gram-positive ordinal depth",
            "PRIMARY",
            "mean_depth_difference",
            num(old_h2p["mean_depth_difference_a_minus_b"]),
            math.nan,
            math.nan,
            num(old_h2p["permutation_p_value_two_sided"]),
            num(new_h2p["mean_depth_difference_a_minus_b"]),
            math.nan,
            math.nan,
            num(new_h2p["permutation_p_value_two_sided"]),
            "ROBUST",
            "RETAIN: organism-group depth difference persists.",
        ),
        old_new_result_row(
            "H2",
            "Depth >=2 Enterobacterales vs Gram-positive",
            "PRIMARY_CHARACTERIZATION",
            "Fisher_OR",
            num(oh2_ge2["odds_ratio_a_vs_b"]),
            math.nan,
            math.nan,
            num(oh2_ge2["p_value"]),
            num(nh2_ge2["odds_ratio_a_vs_b"]),
            math.nan,
            math.nan,
            num(nh2_ge2["p_value"]),
            "ROBUST",
            "RETAIN: difference remains concentrated at phenotypic-or-deeper resolution.",
        ),
        old_new_result_row(
            "H2",
            "Clinical utility Enterobacterales vs Gram-positive",
            "CHARACTERIZATION",
            "Fisher_OR",
            num(oh2_util["odds_ratio_a_vs_b"]),
            math.nan,
            math.nan,
            num(oh2_util["p_value"]),
            num(nh2_util["odds_ratio_a_vs_b"]),
            math.nan,
            math.nan,
            num(nh2_util["p_value"]),
            "ROBUST_NULL",
            "RETAIN: no material organism-group difference in registered clinical utility.",
        ),
        old_new_result_row(
            "H2",
            "Proportional-odds threshold interaction",
            "DIAGNOSTIC",
            "interaction_OR_ratio",
            num(old_h2d["interaction_or_ratio"]),
            math.nan,
            math.nan,
            num(old_h2d["p_value"]),
            num(new_h2d["interaction_or_ratio"]),
            math.nan,
            math.nan,
            num(new_h2d["p_value"]),
            "ROBUST_HETEROGENEITY",
            "RETAIN threshold-specific/permutation H2 hierarchy.",
        ),
        old_new_result_row(
            "H3",
            "Depth-4 quantitative AMR-mechanism diagnostic prevalence",
            "PRIMARY",
            "depth4_n",
            num(old_h3["depth4_n"]),
            math.nan,
            num(old_h3["one_sided_exact_95_upper_percent"]),
            math.nan,
            num(new_h3["depth4_n"]),
            math.nan,
            num(new_h3["one_sided_exact_95_upper_percent"]),
            math.nan,
            "STRENGTHENED",
            "RETAIN and strengthen precision: 0/573 with narrower one-sided upper bound.",
        ),
        old_new_result_row(
            "H4",
            "Categorical diagnostic depth x clinical utility",
            "PRIMARY",
            "permutation_chi_square",
            num(old_h4p["chi_square"]),
            math.nan,
            math.nan,
            num(old_h4p["permutation_p_value"]),
            num(new_h4p["chi_square"]),
            math.nan,
            math.nan,
            num(new_h4p["permutation_p_value"]),
            "ROBUST",
            "RETAIN categorical depth-utility association.",
        ),
        old_new_result_row(
            "H4",
            "Depth 2 vs depth 1 clinical utility",
            "PRIMARY_CHARACTERIZATION",
            "Fisher_OR",
            num(oh4_21["odds_ratio_a_vs_b"]),
            math.nan,
            math.nan,
            num(oh4_21["p_value"]),
            num(nh4_21["odds_ratio_a_vs_b"]),
            math.nan,
            math.nan,
            num(nh4_21["p_value"]),
            "ROBUST_NULL",
            "RETAIN: no evidence of further utility increase from depth 1 to depth 2.",
        ),
        old_new_result_row(
            "H4",
            "CORE_AMR vs BROAD clinical utility",
            "STRUCTURAL",
            "Fisher_OR",
            num(oh4_core_broad["odds_ratio_a_vs_b"]),
            math.nan,
            math.nan,
            num(oh4_core_broad["p_value"]),
            num(nh4_core_broad["odds_ratio_a_vs_b"]),
            math.nan,
            math.nan,
            num(nh4_core_broad["p_value"]),
            "ROBUST_ATTENUATED",
            "RETAIN: utility remains substantially more common in CORE_AMR studies.",
        ),
        old_new_result_row(
            "H4",
            "Within CORE_AMR: depth 2 vs depth 1 clinical utility",
            "STRUCTURAL",
            "Fisher_OR",
            num(oh4_core21["odds_ratio_a_vs_b"]),
            math.nan,
            math.nan,
            num(oh4_core21["p_value"]),
            num(nh4_core21["odds_ratio_a_vs_b"]),
            math.nan,
            math.nan,
            num(nh4_core21["p_value"]),
            "ROBUST_NULL",
            "RETAIN: within CORE_AMR, no further utility increase from depth 1 to 2.",
        ),
    ]

    main_compare = pd.DataFrame(comparisons)
    write_tsv(
        main_compare,
        out / f"main_H1_H4_comparison_v3_2_8_vs_v3_3_1.tsv",
    )

    # ------------------------------------------------------------------
    # H4 depth-specific utility comparison.
    # ------------------------------------------------------------------
    h4_depth_rows = []
    for depth in [0, 1, 2, 3]:
        o = select_row(old_h4d, "depth", str(depth))
        n = select_row(new_h4d, "depth", str(depth))
        h4_depth_rows.append({
            "depth": depth,
            "v3_2_8_n": int(float(o["n"])),
            "v3_2_8_utility_yes": int(float(o["utility_yes"])),
            "v3_2_8_utility_percent": num(o["utility_percent"]),
            "v3_3_1_n": int(float(n["n"])),
            "v3_3_1_utility_yes": int(float(n["utility_yes"])),
            "v3_3_1_utility_percent": num(n["utility_percent"]),
            "percentage_point_change": (
                num(n["utility_percent"]) - num(o["utility_percent"])
            ),
        })
    h4_depth_compare = pd.DataFrame(h4_depth_rows)
    write_tsv(
        h4_depth_compare,
        out / f"H4_utility_by_depth_comparison_{COMPARISON_TAG}.tsv",
    )

    # ------------------------------------------------------------------
    # H2 modality/FDR comparison for the key phenotypic AST token.
    # ------------------------------------------------------------------
    old_modtok = read_tsv(paths["old_h2_modtok"])
    new_modtok = read_tsv(paths["new_h2_modtok"])
    old_ast = select_row(old_modtok, "token", "PHENOTYPIC_AST_OR_MIC")
    new_ast = select_row(new_modtok, "token", "PHENOTYPIC_AST_OR_MIC")

    h2_modality_compare = pd.DataFrame([{
        "analysis": "Overall primary-modality distribution",
        "v3_2_8_p_value": num(old_h2mp["permutation_p_value"]),
        "v3_3_1_p_value": num(new_h2mp["permutation_p_value"]),
        "assessment": "ROBUST",
    }, {
        "analysis": "PHENOTYPIC_AST_OR_MIC token",
        "v3_2_8_p_value": num(old_ast["p_value"]),
        "v3_2_8_fdr_bh": num(old_ast["fdr_bh"]),
        "v3_3_1_p_value": num(new_ast["p_value"]),
        "v3_3_1_fdr_bh": num(new_ast["fdr_bh"]),
        "v3_3_1_entero_percent": num(new_ast["a_percent"]),
        "v3_3_1_gram_positive_percent": num(new_ast["b_percent"]),
        "assessment": "ROBUST_AFTER_FDR",
    }])
    write_tsv(
        h2_modality_compare,
        out / f"H2_modality_comparison_{COMPARISON_TAG}.tsv",
    )

    # ------------------------------------------------------------------
    # P35 robustness comparison within v3.3.1.
    # ------------------------------------------------------------------
    p35 = read_tsv(paths["p35_main"])

    def p35_row(hyp: str, phrase: str) -> pd.Series:
        sub = p35[p35["hypothesis"].map(text) == hyp]
        return select_contains(sub, "analysis", phrase)

    p35_specs = [
        ("H1", "Diagnostic depth per 5-year increase", nh1["odds_ratio"], nh1["p_value"]),
        (
            "H2",
            "Enterobacterales vs Gram-positive ordinal depth permutation",
            new_h2p["mean_depth_difference_a_minus_b"],
            new_h2p["permutation_p_value_two_sided"],
        ),
        ("H2", "Depth >=2 Enterobacterales vs Gram-positive", nh2_ge2["odds_ratio_a_vs_b"], nh2_ge2["p_value"]),
        ("H3", "Depth-4 quantitative AMR mechanism evaluation", new_h3["depth4_n"], ""),
        ("H4", "Categorical diagnostic depth x clinical utility", new_h4p["chi_square"], new_h4p["permutation_p_value"]),
        ("H4", "CORE_AMR vs BROAD clinical utility", nh4_core_broad["odds_ratio_a_vs_b"], nh4_core_broad["p_value"]),
        ("H4", "Within CORE_AMR: depth 2 vs depth 1", nh4_core21["odds_ratio_a_vs_b"], nh4_core21["p_value"]),
    ]

    p35_rows = []
    for hyp, phrase, primary_effect, primary_p in p35_specs:
        pr = p35_row(hyp, phrase)
        p35_rows.append({
            "hypothesis": hyp,
            "analysis": text(pr["analysis"]),
            "primary_573_effect": num(primary_effect),
            "primary_573_p_value": num(primary_p),
            "P35_565_effect": num(pr["effect"]),
            "P35_565_ci_low": num(pr["ci_low"]),
            "P35_565_ci_high": num(pr["ci_high"]),
            "P35_565_p_value": num(pr["p_value"]),
            "direction_preserved": (
                "YES"
                if (
                    math.isnan(num(primary_effect))
                    or math.isnan(num(pr["effect"]))
                    or num(primary_effect) == 0
                    or num(pr["effect"]) == 0
                    or math.copysign(1, num(primary_effect))
                    == math.copysign(1, num(pr["effect"]))
                )
                else "NO"
            ),
        })
    p35_compare = pd.DataFrame(p35_rows)
    write_tsv(
        p35_compare,
        out / f"P35_robustness_comparison_{COMPARISON_TAG}.tsv",
    )

    # ------------------------------------------------------------------
    # Explicit claim-transition memo.
    # ------------------------------------------------------------------
    old_h1_p = num(oh1["p_value"])
    new_h1_p = num(nh1["p_value"])
    old_h1_diag_p = num(oh1d["p_value"])
    new_h1_diag_p = num(nh1d["p_value"])

    memo = f"""# Frozen analysis comparison and manuscript claim transition

**Project:** ClinicalTrials.gov bacterial/AMR diagnostic landscape  
**Comparison:** v3.2.8 (historic 360-study analysis) vs v3.3.1 (final 573-study analysis)  
**Comparison package:** {COMPARISON_VERSION}  
**Generated:** {now_utc()}

## Purpose

This document records how the final screening-coverage rescue and cohort expansion changed the already-frozen H1-H4 results. It is an interpretation/provenance artifact only. No model was re-fit by this comparison package.

The v3.2.8 and v3.3.1 analysis directories passed their own SHA-256 manifests before extraction. The historic 360-study release is a strict subset of the final 573-study primary release; {len(added)} studies were added and no historic eligible record was removed.

## Cohort expansion

The final eligible cohort increased from **{len(old_release)} to {len(new_release)} studies**.

The {len(added)} newly eligible studies contributed:

- **{int((added["final_amr_depth"]=="0").sum())} depth-0 studies**
- **{int((added["final_amr_depth"]=="1").sum())} depth-1 studies**
- **{int((added["final_amr_depth"]=="2").sum())} depth-2 studies**
- **{int((added["final_amr_depth"]=="3").sum())} depth-3 studies**
- **{int((added["final_amr_depth"]=="4").sum())} depth-4 studies**

Among newly added studies beginning in 2020 or later, **{recent_added_d0}/{recent_added_n} ({recent_added_d0_pct:.1f}%) were depth 0**. This is the key compositional reason the final temporal result is more conservative than the historic 360-study result.

## H1 — temporal progression: CLAIM REQUIRES REVISION

Historic v3.2.8:
- proportional-odds OR per 5 years = **{num(oh1["odds_ratio"]):.2f}**
- 95% CI **{num(oh1["ci_low"]):.2f}–{num(oh1["ci_high"]):.2f}**
- p = **{old_h1_p:.4g}**
- cumulative-threshold interaction p = **{old_h1_diag_p:.4g}**

Final v3.3.1:
- proportional-odds OR per 5 years = **{num(nh1["odds_ratio"]):.2f}**
- 95% CI **{num(nh1["ci_low"]):.2f}–{num(nh1["ci_high"]):.2f}**
- p = **{new_h1_p:.4g}**
- depth >=1 OR = **{num(nh1_ge1["odds_ratio"]):.2f}**, p = **{num(nh1_ge1["p_value"]):.4g}**
- depth >=2 OR = **{num(nh1_ge2["odds_ratio"]):.2f}**, p = **{num(nh1_ge2["p_value"]):.4g}**
- cumulative-threshold interaction p = **{new_h1_diag_p:.4g}**

**Final interpretation:** the expanded cohort does not demonstrate a statistically supported monotonic increase in diagnostic depth over calendar time. Threshold-specific results suggest, at most, a modest increase in progression beyond organism-only testing, without evidence that phenotypic-or-deeper AMR resolution increased over time.

**Manuscript action:** remove the historic claim that the overall landscape shifted significantly toward greater resistance resolution over time. Explain transparently that improved screening coverage added many recent BROAD/depth-0 studies.

## H2 — organism-group differences: RETAIN

Final v3.3.1:
- Enterobacterales mean depth = **{num(new_h2p["a_mean_depth"]):.3f}**
- Gram-positive mean depth = **{num(new_h2p["b_mean_depth"]):.3f}**
- difference = **{num(new_h2p["mean_depth_difference_a_minus_b"]):.3f}**
- permutation p = **{num(new_h2p["permutation_p_value_two_sided"]):.4g}**
- depth >=2 exact OR = **{num(nh2_ge2["odds_ratio_a_vs_b"]):.2f}**, p = **{num(nh2_ge2["p_value"]):.4g}**
- year-adjusted depth >=2 OR = **{num(nh2_adj2["odds_ratio"]):.2f}**, 95% CI **{num(nh2_adj2["ci_low"]):.2f}–{num(nh2_adj2["ci_high"]):.2f}**, p = **{num(nh2_adj2["p_value"]):.4g}**
- H2 threshold-interaction diagnostic p = **{num(new_h2d["p_value"]):.4g}**
- clinical-utility exact OR = **{num(nh2_util["odds_ratio_a_vs_b"]):.2f}**, p = **{num(nh2_util["p_value"]):.4g}**

Rare-pathogen sensitivity:
- mean-depth difference = **{num(new_h2rp["mean_depth_difference_a_minus_b"]):.3f}**
- permutation p = **{num(new_h2rp["permutation_p_value_two_sided"]):.4g}**
- depth >=2 OR = **{num(select_row(new_h2re, "outcome", "depth_ge2")["odds_ratio_a_vs_b"]):.2f}**
- depth >=2 p = **{num(select_row(new_h2re, "outcome", "depth_ge2")["p_value"]):.4g}**

**Final interpretation:** Enterobacterales-focused studies disproportionately reach phenotypic AST-level resolution compared with Gram-positive studies. The effect is threshold-specific rather than a uniform proportional-odds shift. Clinical-utility registration does not materially differ between the organism groups.

## H3 — quantitative mechanistic resolution: RETAIN AND STRENGTHEN

Historic v3.2.8: **0/{int(float(old_h3["n"]))}** depth-4 studies; one-sided exact 95% upper bound **{num(old_h3["one_sided_exact_95_upper_percent"]):.2f}%**.

Final v3.3.1: **0/{int(float(new_h3["n"]))}** depth-4 studies; one-sided exact 95% upper bound **{num(new_h3["one_sided_exact_95_upper_percent"]):.2f}%**. Depth 3 remains **{int(float(new_h3["depth3_n"]))} studies**.

**Final interpretation:** no eligible registered study evaluated a quantitative AMR-mechanism measurement as the diagnostic resistance output. Integrated multimechanism interpretation remains exceptionally rare.

The separate 15-study v3.3.0 keyword queue remains nonbinding QC and must not be promoted to 15 formal near misses or used to reopen frozen depth.

## H4 — clinical translation: RETAIN WITH SHARPER FRAMING

Final v3.3.1 utility prevalence:
"""

    for _, r in new_h4d.sort_values("depth").iterrows():
        memo += (
            f"- depth {text(r['depth'])}: **{int(float(r['utility_yes']))}/"
            f"{int(float(r['n']))} ({num(r['utility_percent']):.1f}%)**\n"
        )

    memo += f"""
Primary categorical depth × utility permutation p = **{num(new_h4p["permutation_p_value"]):.4g}**.

Depth 2 versus depth 1:
- OR = **{num(nh4_21["odds_ratio_a_vs_b"]):.2f}**
- p = **{num(nh4_21["p_value"]):.4g}**

CORE_AMR versus BROAD:
- exact OR = **{num(nh4_core_broad["odds_ratio_a_vs_b"]):.2f}**
- p = **{num(nh4_core_broad["p_value"]):.4g}**

Within CORE_AMR, depth 2 versus depth 1:
- OR = **{num(nh4_core21["odds_ratio_a_vs_b"]):.2f}**
- p = **{num(nh4_core21["p_value"]):.4g}**

**Final interpretation:** clinical-utility evaluation is substantially more common after studies enter the AMR-focused diagnostic space and at depths 1–2 than among organism-only studies, but there is no evidence of an additional increase from categorical resistance detection to phenotypic susceptibility testing. Depth 3 remains descriptive.

## P35 imaging sensitivity: ROBUSTNESS CONFIRMED

The frozen P35 sensitivity excludes eight already-eligible pathogen-directed imaging studies, yielding n={EXPECTED_P35_N}. It does not alter the primary cohort.

The P35 results preserve the substantive final conclusions:
- H1 remains nonsignificant.
- H2 ordinal and depth >=2 conclusions remain supported.
- H3 remains 0 depth-4 studies.
- H4 categorical and CORE-vs-BROAD conclusions remain supported.
- Within-CORE depth 2 versus depth 1 remains unsupported.

Use `P35_robustness_comparison_{COMPARISON_TAG}.tsv` for exact values.

## Final claim hierarchy for manuscript v2

1. **The registered bacterial diagnostic landscape remains dominated by organism-level evaluation, and the final expanded cohort does not show convincing temporal progression toward deeper AMR resolution.**
2. **Enterobacterales studies disproportionately reach phenotypic AST-level resolution relative to Gram-positive studies.**
3. **Quantitative AMR-mechanism diagnostic evaluation is absent in the eligible registry cohort, and integrated multimechanism interpretation is exceptionally rare.**
4. **Clinical-utility evaluation is more common once studies enter the AMR-resolving diagnostic space, but does not continue increasing from categorical resistance detection to phenotypic susceptibility testing.**
5. **Expansion of bacterial diagnostic research has not been matched by corresponding progression toward quantitatively resolved AMR mechanisms or progressively deeper clinical-utility evaluation.**

## Versioning decision

- v3.2.7 and v3.2.8 remain frozen historical provenance.
- v3.2.9 remains the definitive screening freeze.
- v3.3.0 remains the definitive 573-study analytic data release.
- v3.3.1 remains the definitive final H1-H4 analysis.
- This {COMPARISON_VERSION} package is an interpretation/provenance layer only and must not replace or modify any of those frozen artifacts.

## Next manuscript-development gate

After this comparison package is frozen, build the manuscript-v2 source bundle from:
1. final v3.2.9 screening-freeze artifacts;
2. frozen v3.3.0 analytic release;
3. frozen v3.3.1 statistical analysis;
4. frozen P35 sensitivity artifacts;
5. this comparison package;
6. manuscript v1 as editorial source only.

The manuscript-v2 study-selection flow must use the final v3.2.9 screening universe rather than the superseded 2,097-reviewed / 360-eligible historical flow.
"""

    memo_path = out / f"MANUSCRIPT_CLAIM_TRANSITION_{COMPARISON_TAG}.md"
    memo_path.write_text(memo, encoding="utf-8")

    # Machine-readable summary.
    summary = {
        "created_at": now_utc(),
        "comparison_version": COMPARISON_VERSION,
        "comparison_pass": True,
        "analysis_versions": {
            "historic": "v3.2.8",
            "final": "v3.3.1",
        },
        "release_versions": {
            "historic": "v3.2.7",
            "final": "v3.3.0",
        },
        "input_integrity": {
            "historic_analysis_manifest": old_manifest,
            "final_analysis_manifest": new_manifest,
            "historic_release_sha256": sha256(args.old_release),
            "final_release_sha256": sha256(args.new_release),
        },
        "cohort": {
            "historic_n": len(old_release),
            "final_n": len(new_release),
            "added_n": len(added),
            "lost_n": len(lost_ids),
            "added_depth_distribution": depth_counts(added),
            "added_2020_or_later_n": recent_added_n,
            "added_2020_or_later_depth0_n": recent_added_d0,
            "added_2020_or_later_depth0_percent": recent_added_d0_pct,
        },
        "claim_transition": {
            "H1": "REVISE",
            "H2": "RETAIN",
            "H3": "RETAIN_AND_STRENGTHEN",
            "H4": "RETAIN_WITH_SHARPER_FRAMING",
            "P35": "ROBUSTNESS_CONFIRMED",
        },
        "next_gate": (
            "Freeze this comparison package, then assemble manuscript-v2 source "
            "materials and rebuild the final v3.2.9 study-selection flow before "
            "editing manuscript text."
        ),
    }
    summary_path = out / f"comparison_summary_{COMPARISON_TAG}.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    readme = f"""ClinicalTrials.gov bacterial/AMR diagnostic landscape
Frozen analysis comparison package {COMPARISON_VERSION}

Purpose
-------
Compare the frozen historic v3.2.8 analysis with the frozen final v3.3.1
analysis after the v3.2.9 screening-coverage rescue and v3.3.0 cohort rebuild.

This package does not refit models.

Authoritative statistical results remain in their source analysis directories.
This comparison package only extracts and juxtaposes those frozen results and
records the manuscript claim transition.

Primary manuscript consequence
------------------------------
H1 changes materially and must be revised. H2 remains robust. H3 is
strengthened by the larger zero-event denominator. H4 remains robust with
sharper emphasis on entry into the AMR-focused space rather than monotonic
deepening of utility.

After this package passes integrity checks, make it read-only. The next task is
to assemble manuscript-v2 sources and rebuild the final screening flow from
the v3.2.9 screening freeze.
"""
    readme_path = out / f"README_COMPARISON_{COMPARISON_TAG}.txt"
    readme_path.write_text(readme, encoding="utf-8")

    # File manifest/checksums.
    outputs = [
        out / f"cohort_expansion_summary_{COMPARISON_TAG}.tsv",
        out / f"added_213_by_era_and_depth_{COMPARISON_TAG}.tsv",
        out / f"main_H1_H4_comparison_v3_2_8_vs_v3_3_1.tsv",
        out / f"H4_utility_by_depth_comparison_{COMPARISON_TAG}.tsv",
        out / f"H2_modality_comparison_{COMPARISON_TAG}.tsv",
        out / f"P35_robustness_comparison_{COMPARISON_TAG}.tsv",
        memo_path,
        summary_path,
        readme_path,
    ]

    manifest_rows = []
    for p in outputs:
        parsed = ""
        if p.suffix == ".tsv":
            parsed = str(len(read_tsv(p)))
        manifest_rows.append({
            "filename": p.name,
            "sha256": sha256(p),
            "bytes": p.stat().st_size,
            "parsed_rows_if_tsv": parsed,
        })

    manifest_path = out / f"comparison_manifest_{COMPARISON_TAG}.tsv"
    write_tsv(pd.DataFrame(manifest_rows), manifest_path)
    outputs.append(manifest_path)

    sums_path = out / "SHA256SUMS.txt"
    with sums_path.open("w", encoding="utf-8") as h:
        for p in sorted(outputs, key=lambda x: x.name):
            h.write(f"{sha256(p)}  {p.name}\n")

    print("V3.2.8 -> V3.3.1 COMPARISON PACKAGE: PASS")
    print(f"Historic eligible n: {len(old_release)}")
    print(f"Final eligible n: {len(new_release)}")
    print(f"Added eligible n: {len(added)}")
    print(f"Lost historic eligible n: {len(lost_ids)}")
    print(f"Added depth counts: {depth_counts(added)}")
    print(
        "Added >=2020 depth-0: "
        f"{recent_added_d0}/{recent_added_n} ({recent_added_d0_pct:.1f}%)"
    )
    print("Claim transition: H1=REVISE; H2=RETAIN; H3=STRENGTHEN; H4=RETAIN")
    print(f"Output directory: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
