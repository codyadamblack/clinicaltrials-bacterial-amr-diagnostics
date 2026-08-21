#!/usr/bin/env python3
"""
Final refined H1-H4 analysis v3.3.1 for the ClinicalTrials.gov
bacterial/AMR diagnostic landscape.

Frozen data release
-------------------
Primary analytic release: v3.3.0, 573 eligible studies.
P35 imaging sensitivity: v3.3.0, 565 eligible studies after the eight
prespecified pathogen-directed imaging exclusions.

Inferential hierarchy
----------------------
This script carries forward the already-established v3.2.8 statistical
hierarchy. It does NOT change eligibility, descriptors, diagnostic depth,
outcomes, hypotheses, or primary-vs-supporting model status.

H1
  Primary: proportional-odds model of diagnostic depth on study start year,
  OR per 5-year increase.
  Supporting: logistic threshold models at depth >=1 and >=2.
  Diagnostic: stacked cumulative-logit interaction comparing the two
  cumulative thresholds.

H2
  Primary omnibus: permutation test of mean ordinal depth difference between
  Enterobacterales and Gram-positive studies.
  Primary characterization: exact Fisher tests at depth >=1 and depth >=2.
  Additional registered H2 outcome: exact clinical-utility comparison.
  Adjusted characterization: threshold logistic models adjusted for year.
  Diagnostic/supporting: threshold-interaction diagnostic and proportional-
  odds models.
  Prespecified sensitivity: exclude H2_EXCLUDE_TYPHOID_PLAGUE records and
  repeat the revised depth hierarchy.
  Modality/utility-category characterization: permutation/Fisher with
  Benjamini-Hochberg FDR, as in v3.2.8.

H3
  Primary: depth-4 prevalence and one-sided exact 95% upper bound.
  Descriptive: depth-3 records and previously frozen quantitative-mechanism
  near-miss flags. The separate 15-study v3.3.0 heuristic queue is copied
  only for QC/provenance when supplied and never changes depth or flags.

H4
  Primary omnibus: permutation chi-square for utility across depth categories.
  Primary characterization: depth-specific Wilson 95% CIs and exact pairwise
  comparisons.
  Structural: CORE_AMR_DIAGNOSTIC versus BROAD_BACTERIAL_DIAGNOSTIC and
  within-CORE depth 2 versus depth 1.
  Supporting: numeric-depth logistic trend models and their previously
  specified sensitivity sets.
  Depth 3 remains descriptive because n=3.

Bias-control boundary
---------------------
- Primary analysis is run on the frozen 573-study release.
- The already-frozen P35 definition is run as a sensitivity analysis using
  the same hierarchy.
- No manuscript conclusion is auto-written by this script.
- No v3.2.8 result is hard-coded into interpretation text.
- Do not alter the frozen v3.3.0 data or this analysis hierarchy after viewing
  v3.3.1 results. Genuine corrections require a new version.

All Monte Carlo permutation tests use fixed seeds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.miscmodels.ordinal_model import OrderedModel
from statsmodels.stats.proportion import proportion_confint


ANALYSIS_VERSION = "v3.3.1"
ANALYSIS_TAG = "v3_3_1"
RELEASE_VERSION = "v3.3.0"
PREDECESSOR_HIERARCHY = "v3.2.8"

EXPECTED_PRIMARY_SHA256 = (
    "686a99d7e33b78822c7e402478589cfec086d8a58b7156468b43e7aaa609c4b2"
)
EXPECTED_P35_SHA256 = (
    "4444890e3fe414add7a9d88cb40eb2cc98621dc409f69ce04e91b04912cb91d2"
)
EXPECTED_H3_QC_SHA256 = (
    "778a05f261375f3fe3006f7a59f77bd40b6074a193f7d5c06e9334930d84a4c0"
)

EXPECTED_PRIMARY_N = 573
EXPECTED_PRIMARY_DEPTHS = {0: 435, 1: 83, 2: 52, 3: 3, 4: 0}
EXPECTED_PRIMARY_H2 = {
    "GRAM_POSITIVE": 62,
    "ENTEROBACTERALES": 25,
    "OTHER_EXCLUDED": 486,
}
EXPECTED_PRIMARY_UTILITY = {"NO": 312, "YES": 261}
EXPECTED_PRIMARY_STATUS = {"FINAL": 521, "FINAL_WITH_UNCERTAINTY": 52}
EXPECTED_PRIMARY_STRATUM = {
    "BROAD_BACTERIAL_DIAGNOSTIC": 433,
    "CORE_AMR_DIAGNOSTIC": 140,
}

EXPECTED_P35_N = 565
EXPECTED_P35_EXCLUDED_IDS = {
    "NCT01378728",
    "NCT02450942",
    "NCT02491164",
    "NCT02558062",
    "NCT03091361",
    "NCT03290690",
    "NCT05285072",
    "NCT06986512",
}

EXPECTED_H2_RARE_IDS = {
    "NCT00128466",
    "NCT02689193",
    "NCT04673487",
    "NCT04688996",
    "NCT04801602",
}

RARE_H2_FLAG = "H2_EXCLUDE_TYPHOID_PLAGUE"
H3_NEAR_MISS_FLAG = "H3_QUANTITATIVE_MECHANISM_NEAR_MISS"
BETALACTA_NCT = "NCT03147807"

REQUIRED_COLUMNS = {
    "nct_id",
    "final_primary_eligible",
    "final_amr_depth",
    "final_stratum",
    "start_year",
    "study_type",
    "final_primary_diagnostic_modality",
    "final_all_diagnostic_modalities",
    "final_organism_group",
    "final_h2_comparison_group",
    "final_clinical_utility_endpoint_categories",
    "final_analytical_endpoint_categories",
    "final_clinical_utility_any",
    "final_preanalytical_flag",
    "descriptor_adjudication_status",
    "special_sensitivity_flags",
    "p35_imaging_sensitivity_exclude",
}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def text(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
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


def parse_multiselect(value: Any) -> set[str]:
    return {x.strip() for x in text(value).split("|") if x.strip()}


def safe_exp(x: float) -> float:
    try:
        return math.exp(float(x))
    except OverflowError:
        return math.inf


def bh_fdr(pvalues: Iterable[float]) -> list[float]:
    p = np.asarray(list(pvalues), dtype=float)
    if len(p) == 0:
        return []
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * len(p) / np.arange(1, len(p) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    result = np.empty(len(p), dtype=float)
    result[order] = adj
    return result.tolist()


def normalized_counter(series: pd.Series) -> dict[str, int]:
    c = Counter(series.map(text))
    return dict(sorted(c.items()))


def depth_counter(series: pd.Series) -> dict[int, int]:
    vals = pd.to_numeric(series, errors="raise").astype(int)
    c = Counter(vals)
    return {d: c.get(d, 0) for d in [0, 1, 2, 3, 4]}


def require_equal(observed: Any, expected: Any, label: str) -> None:
    if observed != expected:
        raise SystemExit(
            f"{label} mismatch.\nObserved: {observed}\nExpected: {expected}"
        )


def validate_primary(df: pd.DataFrame, path: Path) -> None:
    require_equal(
        sha256(path),
        EXPECTED_PRIMARY_SHA256,
        "Primary input SHA-256",
    )
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise SystemExit(f"Primary release missing columns: {sorted(missing)}")
    require_equal(len(df), EXPECTED_PRIMARY_N, "Primary row count")
    require_equal(df["nct_id"].nunique(), EXPECTED_PRIMARY_N, "Primary unique NCT count")
    require_equal(set(df["final_primary_eligible"].map(text)), {"YES"}, "Primary eligibility values")
    require_equal(depth_counter(df["final_amr_depth"]), EXPECTED_PRIMARY_DEPTHS, "Primary depth")
    require_equal(
        normalized_counter(df["final_h2_comparison_group"]),
        EXPECTED_PRIMARY_H2,
        "Primary H2",
    )
    require_equal(
        normalized_counter(df["final_clinical_utility_any"]),
        EXPECTED_PRIMARY_UTILITY,
        "Primary utility",
    )
    require_equal(
        normalized_counter(df["descriptor_adjudication_status"]),
        EXPECTED_PRIMARY_STATUS,
        "Primary descriptor status",
    )
    require_equal(
        normalized_counter(df["final_stratum"]),
        EXPECTED_PRIMARY_STRATUM,
        "Primary stratum",
    )
    p35 = set(
        df.loc[
            df["p35_imaging_sensitivity_exclude"].map(text) == "YES",
            "nct_id",
        ].map(text)
    )
    require_equal(p35, EXPECTED_P35_EXCLUDED_IDS, "Primary P35 flagged NCT set")

    rare = set(
        df.loc[
            df["special_sensitivity_flags"].map(
                lambda x: RARE_H2_FLAG in parse_multiselect(x)
            ),
            "nct_id",
        ].map(text)
    )
    require_equal(rare, EXPECTED_H2_RARE_IDS, "Primary H2 rare-pathogen flag set")

    if BETALACTA_NCT not in set(df["nct_id"].map(text)):
        raise SystemExit(f"Required BetaLACTA sensitivity record missing: {BETALACTA_NCT}")


def validate_p35(primary: pd.DataFrame, p35: pd.DataFrame, path: Path) -> None:
    require_equal(sha256(path), EXPECTED_P35_SHA256, "P35 input SHA-256")
    missing = REQUIRED_COLUMNS - set(p35.columns)
    if missing:
        raise SystemExit(f"P35 release missing columns: {sorted(missing)}")
    require_equal(len(p35), EXPECTED_P35_N, "P35 row count")
    require_equal(p35["nct_id"].nunique(), EXPECTED_P35_N, "P35 unique NCT count")

    primary_ids = set(primary["nct_id"].map(text))
    p35_ids = set(p35["nct_id"].map(text))
    require_equal(primary_ids - p35_ids, EXPECTED_P35_EXCLUDED_IDS, "P35 removed NCT set")
    require_equal(p35_ids - primary_ids, set(), "P35 unexpected NCT set")
    require_equal(
        set(p35["p35_imaging_sensitivity_exclude"].map(text)),
        {"NO"},
        "P35 retained flag values",
    )

    # All shared columns must be byte-equivalent at parsed cell level for
    # retained NCTs; the sensitivity file is a strict row subset.
    shared = [c for c in primary.columns if c in p35.columns]
    px = primary.set_index("nct_id", drop=False)
    sx = p35.set_index("nct_id", drop=False)
    changed = []
    for nct in sorted(p35_ids):
        for c in shared:
            if text(px.at[nct, c]) != text(sx.at[nct, c]):
                changed.append((nct, c))
                if len(changed) >= 20:
                    break
        if len(changed) >= 20:
            break
    if changed:
        raise SystemExit(
            "P35 retained-row values differ from primary release. "
            f"First mismatches: {changed}"
        )


def prepare(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    d["depth"] = pd.to_numeric(d["final_amr_depth"], errors="raise").astype(int)
    d["start_year_num"] = pd.to_numeric(d["start_year"], errors="coerce")
    d["year5"] = (d["start_year_num"] - 2020.0) / 5.0
    d["utility_yes"] = (
        d["final_clinical_utility_any"].map(text) == "YES"
    ).astype(int)
    d["depth_ge1"] = (d["depth"] >= 1).astype(int)
    d["depth_ge2"] = (d["depth"] >= 2).astype(int)
    d["h2_group"] = d["final_h2_comparison_group"].map(text)
    d["is_entero"] = (d["h2_group"] == "ENTEROBACTERALES").astype(int)
    d["study_type_clean"] = d["study_type"].map(text)
    d["preanalytical_yes"] = (
        d["final_preanalytical_flag"].map(text) == "YES"
    )
    d["fully_final"] = (
        d["descriptor_adjudication_status"].map(text) == "FINAL"
    )
    d["core_amr"] = (
        d["final_stratum"].map(text) == "CORE_AMR_DIAGNOSTIC"
    ).astype(int)
    d["rare_h2_sensitivity"] = d["special_sensitivity_flags"].map(
        lambda x: RARE_H2_FLAG in parse_multiselect(x)
    )
    d["h3_near_miss"] = d["special_sensitivity_flags"].map(
        lambda x: H3_NEAR_MISS_FLAG in parse_multiselect(x)
    )

    def era(y: float) -> str:
        if pd.isna(y):
            return "MISSING"
        if y <= 2009:
            return "<=2009"
        if y <= 2014:
            return "2010-2014"
        if y <= 2019:
            return "2015-2019"
        return ">=2020"

    d["start_era"] = d["start_year_num"].map(era)
    return d


def ordinal_model(
    data: pd.DataFrame,
    predictors: list[str],
    label: str,
) -> dict[str, Any]:
    dd = data.dropna(subset=["depth", *predictors]).copy()
    out = {
        "analysis": label,
        "n": len(dd),
        "status": "FAIL",
        "focal_term": predictors[0],
        "estimate": np.nan,
        "odds_ratio": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "p_value": np.nan,
        "note": "",
    }
    try:
        model = OrderedModel(
            dd["depth"].astype(int),
            dd[predictors].astype(float),
            distr="logit",
        )
        fit = model.fit(method="bfgs", disp=False, maxiter=2000)
        term = predictors[0]
        b = float(fit.params[term])
        se = float(fit.bse[term])
        out.update(
            status="PASS",
            estimate=b,
            odds_ratio=safe_exp(b),
            ci_low=safe_exp(b - 1.96 * se),
            ci_high=safe_exp(b + 1.96 * se),
            p_value=float(fit.pvalues[term]),
        )
    except Exception as exc:
        out["note"] = repr(exc)
    return out


def logistic_model(
    formula: str,
    data: pd.DataFrame,
    focal: str,
    label: str,
) -> dict[str, Any]:
    out = {
        "analysis": label,
        "formula": formula,
        "n": len(data),
        "status": "FAIL",
        "focal_term": focal,
        "estimate": np.nan,
        "odds_ratio": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "p_value": np.nan,
        "note": "",
    }
    try:
        fit = smf.glm(
            formula,
            data=data,
            family=sm.families.Binomial(),
        ).fit()
        b = float(fit.params[focal])
        se = float(fit.bse[focal])
        out.update(
            n=int(fit.nobs),
            status="PASS",
            estimate=b,
            odds_ratio=safe_exp(b),
            ci_low=safe_exp(b - 1.96 * se),
            ci_high=safe_exp(b + 1.96 * se),
            p_value=float(fit.pvalues[focal]),
        )
    except Exception as exc:
        out["note"] = repr(exc)
    return out


def fisher_2x2(
    a_yes: int,
    a_no: int,
    b_yes: int,
    b_no: int,
    label: str,
) -> dict[str, Any]:
    odds, p = fisher_exact(
        [[a_yes, a_no], [b_yes, b_no]],
        alternative="two-sided",
    )
    return {
        "analysis": label,
        "a_yes": a_yes,
        "a_no": a_no,
        "b_yes": b_yes,
        "b_no": b_no,
        "odds_ratio_a_vs_b": odds,
        "p_value": p,
    }


def fisher_group_outcome(
    data: pd.DataFrame,
    group_col: str,
    a: Any,
    b: Any,
    outcome_col: str,
    label: str,
) -> dict[str, Any]:
    aa = data[group_col] == a
    bb = data[group_col] == b
    a_yes = int((aa & (data[outcome_col] == 1)).sum())
    a_no = int((aa & (data[outcome_col] == 0)).sum())
    b_yes = int((bb & (data[outcome_col] == 1)).sum())
    b_no = int((bb & (data[outcome_col] == 0)).sum())
    out = fisher_2x2(a_yes, a_no, b_yes, b_no, label)
    out.update(
        group_a=a,
        group_b=b,
        outcome=outcome_col,
    )
    return out


def stacked_threshold_diagnostic(
    data: pd.DataFrame,
    focal: str,
    adjust_terms: list[str],
    label: str,
) -> dict[str, Any]:
    rows = []
    for _, r in data.iterrows():
        for threshold, threshold_label in [(1, "ge1"), (2, "ge2")]:
            item = {
                "nct_id": r["nct_id"],
                "threshold": threshold_label,
                "y": int(r["depth"] >= threshold),
                focal: r[focal],
            }
            for term in adjust_terms:
                item[term] = r[term]
            rows.append(item)

    long = pd.DataFrame(rows).dropna(subset=[focal, *adjust_terms])
    formula = f"y ~ {focal} + C(threshold) + {focal}:C(threshold)"
    if adjust_terms:
        formula += " + " + " + ".join(adjust_terms)

    interaction = f"{focal}:C(threshold)[T.ge2]"
    out = {
        "analysis": label,
        "n_studies": long["nct_id"].nunique(),
        "n_stacked_rows": len(long),
        "formula": formula,
        "interaction_term": interaction,
        "interaction_estimate": np.nan,
        "interaction_or_ratio": np.nan,
        "p_value": np.nan,
        "status": "FAIL",
        "note": "",
    }
    try:
        fit = smf.glm(
            formula,
            data=long,
            family=sm.families.Binomial(),
        ).fit(
            cov_type="cluster",
            cov_kwds={"groups": long["nct_id"]},
        )
        b = float(fit.params[interaction])
        out.update(
            status="PASS",
            interaction_estimate=b,
            interaction_or_ratio=safe_exp(b),
            p_value=float(fit.pvalues[interaction]),
        )
    except Exception as exc:
        out["note"] = repr(exc)
    return out


def permutation_mean_depth(
    data: pd.DataFrame,
    group_col: str,
    group_a: str,
    group_b: str,
    n_perm: int,
    seed: int,
    label: str,
) -> dict[str, Any]:
    sub = data[data[group_col].isin([group_a, group_b])].copy()
    groups = sub[group_col].to_numpy()
    depth = sub["depth"].to_numpy(dtype=float)
    mask_a = groups == group_a
    n_a = int(mask_a.sum())
    obs = depth[mask_a].mean() - depth[~mask_a].mean()

    rng = np.random.default_rng(seed)
    exceed = 0
    indices = np.arange(len(depth))
    for _ in range(n_perm):
        selected = rng.choice(indices, size=n_a, replace=False)
        perm_a = np.zeros(len(depth), dtype=bool)
        perm_a[selected] = True
        stat = depth[perm_a].mean() - depth[~perm_a].mean()
        if abs(stat) >= abs(obs) - 1e-12:
            exceed += 1

    return {
        "analysis": label,
        "group_a": group_a,
        "group_b": group_b,
        "a_n": n_a,
        "b_n": len(depth) - n_a,
        "a_mean_depth": float(depth[mask_a].mean()),
        "b_mean_depth": float(depth[~mask_a].mean()),
        "mean_depth_difference_a_minus_b": float(obs),
        "permutations": n_perm,
        "permutation_p_value_two_sided": (exceed + 1) / (n_perm + 1),
    }


def permutation_chisq(
    row_group: np.ndarray,
    binary_outcome: np.ndarray,
    row_levels: list[Any],
    n_perm: int,
    seed: int,
    label: str,
) -> dict[str, Any]:
    level_to_code = {level: i for i, level in enumerate(row_levels)}
    codes = np.asarray([level_to_code[x] for x in row_group], dtype=int)
    y = np.asarray(binary_outcome, dtype=int)
    k = len(row_levels)
    row_sizes = np.bincount(codes, minlength=k).astype(float)
    yes_total = int(y.sum())
    no_total = len(y) - yes_total

    yes_obs = np.bincount(codes[y == 1], minlength=k).astype(float)
    no_obs = row_sizes - yes_obs
    observed = np.column_stack([no_obs, yes_obs])
    expected = np.column_stack(
        [
            row_sizes * no_total / len(y),
            row_sizes * yes_total / len(y),
        ]
    )
    if np.any(expected == 0):
        raise SystemExit(f"{label}: zero expected cell in permutation chi-square")
    chi_obs = float(np.sum((observed - expected) ** 2 / expected))

    rng = np.random.default_rng(seed)
    exceed = 0
    indices = np.arange(len(y))
    for _ in range(n_perm):
        selected = rng.choice(indices, size=yes_total, replace=False)
        yes_counts = np.bincount(codes[selected], minlength=k).astype(float)
        no_counts = row_sizes - yes_counts
        table = np.column_stack([no_counts, yes_counts])
        stat = float(np.sum((table - expected) ** 2 / expected))
        if stat >= chi_obs - 1e-12:
            exceed += 1

    return {
        "analysis": label,
        "chi_square": chi_obs,
        "permutations": n_perm,
        "permutation_p_value": (exceed + 1) / (n_perm + 1),
        "row_levels": "|".join(map(str, row_levels)),
    }


def permutation_binary_group_by_category(
    groups: np.ndarray,
    categories: np.ndarray,
    group_a: str,
    group_b: str,
    category_levels: list[str],
    n_perm: int,
    seed: int,
    label: str,
) -> dict[str, Any]:
    groups = np.asarray(groups)
    categories = np.asarray(categories)
    cat_to_code = {cat: i for i, cat in enumerate(category_levels)}
    codes = np.asarray([cat_to_code[x] for x in categories], dtype=int)
    k = len(category_levels)

    n_a = int((groups == group_a).sum())
    n_b = int((groups == group_b).sum())
    totals = np.bincount(codes, minlength=k).astype(float)
    a_obs = np.bincount(codes[groups == group_a], minlength=k).astype(float)
    b_obs = totals - a_obs
    observed = np.vstack([a_obs, b_obs])
    expected = np.vstack(
        [
            totals * n_a / len(groups),
            totals * n_b / len(groups),
        ]
    )
    if np.any(expected == 0):
        raise SystemExit(f"{label}: zero expected cell in modality permutation")
    chi_obs = float(np.sum((observed - expected) ** 2 / expected))

    rng = np.random.default_rng(seed)
    exceed = 0
    indices = np.arange(len(groups))
    for _ in range(n_perm):
        selected = rng.choice(indices, size=n_a, replace=False)
        a_counts = np.bincount(codes[selected], minlength=k).astype(float)
        b_counts = totals - a_counts
        table = np.vstack([a_counts, b_counts])
        stat = float(np.sum((table - expected) ** 2 / expected))
        if stat >= chi_obs - 1e-12:
            exceed += 1

    return {
        "analysis": label,
        "chi_square": chi_obs,
        "permutations": n_perm,
        "permutation_p_value": (exceed + 1) / (n_perm + 1),
        "groups": f"{group_a}|{group_b}",
        "categories": "|".join(category_levels),
    }


def utility_by_depth(data: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for depth in sorted(data["depth"].unique()):
        g = data[data["depth"] == depth]
        n = len(g)
        yes = int(g["utility_yes"].sum())
        low, high = proportion_confint(
            yes,
            n,
            alpha=0.05,
            method="wilson",
        )
        rows.append(
            {
                "analysis_set": label,
                "depth": int(depth),
                "n": n,
                "utility_yes": yes,
                "utility_no": n - yes,
                "utility_percent": 100 * yes / n,
                "wilson_95ci_low_percent": 100 * low,
                "wilson_95ci_high_percent": 100 * high,
            }
        )
    return pd.DataFrame(rows)


def token_level_fisher(
    data: pd.DataFrame,
    multiselect_col: str,
    group_col: str,
    group_a: str,
    group_b: str,
    family: str,
) -> pd.DataFrame:
    sub = data[data[group_col].isin([group_a, group_b])].copy()
    sets = sub[multiselect_col].map(parse_multiselect)
    tokens = sorted(set().union(*sets.tolist())) if len(sets) else []

    rows = []
    for token in tokens:
        present = sets.map(lambda s: token in s)
        a = sub[group_col] == group_a
        b = sub[group_col] == group_b
        a_yes = int((a & present).sum())
        a_no = int((a & ~present).sum())
        b_yes = int((b & present).sum())
        b_no = int((b & ~present).sum())
        odds, p = fisher_exact([[a_yes, a_no], [b_yes, b_no]])
        rows.append(
            {
                "family": family,
                "token": token,
                "group_a": group_a,
                "group_b": group_b,
                "a_n": int(a.sum()),
                "a_yes": a_yes,
                "a_percent": 100 * a_yes / int(a.sum()),
                "b_n": int(b.sum()),
                "b_yes": b_yes,
                "b_percent": 100 * b_yes / int(b.sum()),
                "odds_ratio_group_a_vs_b": odds,
                "p_value": p,
            }
        )

    out = pd.DataFrame(rows)
    if len(out):
        out["fdr_bh"] = bh_fdr(out["p_value"].tolist())
    else:
        out["fdr_bh"] = pd.Series(dtype=float)
    return out


def depth_distribution(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for depth in [0, 1, 2, 3, 4]:
        n = int((data["depth"] == depth).sum())
        rows.append(
            {
                "depth": depth,
                "n": n,
                "percent": 100 * n / len(data),
            }
        )
    return pd.DataFrame(rows)


def figure_depth_by_era(d: pd.DataFrame, out: Path) -> None:
    era_order = ["<=2009", "2010-2014", "2015-2019", ">=2020"]
    plot = (
        d[d["start_era"] != "MISSING"]
        .groupby(["start_era", "depth"])
        .size()
        .reset_index(name="n")
    )
    totals = plot.groupby("start_era")["n"].sum().to_dict()
    plot["percent"] = 100 * plot["n"] / plot["start_era"].map(totals)
    p = plot.pivot(index="start_era", columns="depth", values="percent").fillna(0)
    p = p.reindex(era_order).fillna(0)
    ax = p.plot(kind="bar", stacked=True, figsize=(8.0, 5.4))
    ax.set_xlabel("")
    ax.set_ylabel("Studies within start-year era (%)")
    ax.set_title("Diagnostic depth by study start-year era")
    ax.legend(title="Diagnostic depth", frameon=False)
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()


def figure_h2_depth(d: pd.DataFrame, out: Path) -> None:
    h2 = d[d["h2_group"].isin(["GRAM_POSITIVE", "ENTEROBACTERALES"])]
    dist = h2.groupby(["h2_group", "depth"]).size().reset_index(name="n")
    totals = h2.groupby("h2_group").size().to_dict()
    dist["percent"] = 100 * dist["n"] / dist["h2_group"].map(totals)
    p = dist.pivot(index="h2_group", columns="depth", values="percent").fillna(0)
    p = p.reindex(["GRAM_POSITIVE", "ENTEROBACTERALES"])
    p.index = ["Gram-positive", "Enterobacterales"]
    ax = p.plot(kind="bar", stacked=True, figsize=(7.5, 5.3))
    ax.set_xlabel("")
    ax.set_ylabel("Studies within organism group (%)")
    ax.set_title("Diagnostic depth by organism group")
    ax.legend(title="Diagnostic depth", frameon=False)
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()


def figure_h4_depth(h4_depth: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.4))
    x = h4_depth["depth"].to_numpy()
    y = h4_depth["utility_percent"].to_numpy()
    lower = y - h4_depth["wilson_95ci_low_percent"].to_numpy()
    upper = h4_depth["wilson_95ci_high_percent"].to_numpy() - y
    ax.errorbar(x, y, yerr=np.vstack([lower, upper]), fmt="o", capsize=4)
    for _, row in h4_depth.iterrows():
        ax.annotate(
            f"{int(row['utility_yes'])}/{int(row['n'])}",
            (row["depth"], row["utility_percent"]),
            xytext=(0, 9),
            textcoords="offset points",
            ha="center",
        )
    ax.set_xlabel("Diagnostic depth")
    ax.set_ylabel("Studies with registered clinical utility (%)")
    ax.set_title("Clinical-utility evaluation by diagnostic depth")
    ax.set_xticks([0, 1, 2, 3])
    ymax = min(100, max(80, float(h4_depth["wilson_95ci_high_percent"].max()) + 12))
    ax.set_ylim(0, ymax)
    d3 = h4_depth[h4_depth["depth"] == 3]
    if len(d3):
        ax.text(
            3,
            min(6, ymax * 0.08),
            f"depth 3: n={int(d3.iloc[0]['n'])}",
            ha="center",
            va="bottom",
        )
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()


def figure_h4_stratum(summary: pd.DataFrame, out: Path) -> None:
    sp = summary.copy()
    sp["label"] = sp["final_stratum"].map(
        {
            "BROAD_BACTERIAL_DIAGNOSTIC": "Broad bacterial\ndiagnostic",
            "CORE_AMR_DIAGNOSTIC": "Core AMR\ndiagnostic",
        }
    )
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    ax.bar(sp["label"], sp["utility_percent"])
    ymax = min(100, max(75, float(sp["utility_percent"].max()) + 15))
    for i, row in sp.reset_index(drop=True).iterrows():
        ax.text(
            i,
            row["utility_percent"] + 2,
            f"{int(row['utility_yes'])}/{int(row['n'])}",
            ha="center",
        )
    ax.set_ylabel("Studies with registered clinical utility (%)")
    ax.set_ylim(0, ymax)
    ax.set_title("Clinical-utility evaluation by diagnostic stratum")
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()


def run_analysis_set(
    source: pd.DataFrame,
    outdir: Path,
    label: str,
    n_perm: int,
    seed: int,
    make_figures: bool,
) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)
    d = prepare(source)

    # Descriptive anchors.
    write_tsv(depth_distribution(d), outdir / f"table_depth_distribution_{ANALYSIS_TAG}.tsv")
    write_tsv(
        d.groupby("study_type_clean").size().reset_index(name="n").sort_values("n", ascending=False),
        outdir / f"table_study_type_{ANALYSIS_TAG}.tsv",
    )
    write_tsv(
        d.groupby("final_primary_diagnostic_modality").size().reset_index(name="n").sort_values("n", ascending=False),
        outdir / f"table_primary_modality_{ANALYSIS_TAG}.tsv",
    )
    write_tsv(
        d[d["descriptor_adjudication_status"].map(text) == "FINAL_WITH_UNCERTAINTY"][
            [
                c for c in [
                    "nct_id",
                    "brief_title",
                    "registry_brief_title",
                    "final_amr_depth",
                    "final_stratum",
                    "final_primary_diagnostic_modality",
                    "descriptor_adjudication_notes",
                ]
                if c in d.columns
            ]
        ].copy(),
        outdir / f"table_final_with_uncertainty_records_{ANALYSIS_TAG}.tsv",
    )

    # H1
    h1 = d[d["start_year_num"].notna()].copy()
    h1_primary = pd.DataFrame(
        [
            ordinal_model(
                h1,
                ["year5"],
                "H1 primary proportional odds: depth per 5-year increase",
            ),
            logistic_model(
                "depth_ge1 ~ year5",
                h1,
                "year5",
                "H1 support: depth>=1 per 5 years",
            ),
            logistic_model(
                "depth_ge2 ~ year5",
                h1,
                "year5",
                "H1 support: depth>=2 per 5 years",
            ),
        ]
    )
    write_tsv(h1_primary, outdir / f"H1_primary_and_threshold_models_{ANALYSIS_TAG}.tsv")

    h1_diag = pd.DataFrame(
        [
            stacked_threshold_diagnostic(
                h1,
                "year5",
                [],
                "H1 cumulative-threshold slope diagnostic: >=1 vs >=2",
            )
        ]
    )
    write_tsv(h1_diag, outdir / f"H1_proportional_odds_diagnostic_{ANALYSIS_TAG}.tsv")

    era_order = ["<=2009", "2010-2014", "2015-2019", ">=2020", "MISSING"]
    era_depth = d.groupby(["start_era", "depth"]).size().reset_index(name="n")
    era_n = d.groupby("start_era").size().to_dict()
    era_depth["era_n"] = era_depth["start_era"].map(era_n)
    era_depth["percent_within_era"] = 100 * era_depth["n"] / era_depth["era_n"]
    era_depth["start_era"] = pd.Categorical(
        era_depth["start_era"],
        categories=era_order,
        ordered=True,
    )
    era_depth = era_depth.sort_values(["start_era", "depth"])
    write_tsv(era_depth, outdir / f"H1_depth_by_era_{ANALYSIS_TAG}.tsv")

    era_stratum = (
        d[d["start_era"] != "MISSING"]
        .groupby(["start_era", "final_stratum"])
        .size()
        .reset_index(name="n")
    )
    era_totals = (
        d[d["start_era"] != "MISSING"]
        .groupby("start_era")
        .size()
        .to_dict()
    )
    era_stratum["era_n"] = era_stratum["start_era"].map(era_totals)
    era_stratum["percent_within_era"] = (
        100 * era_stratum["n"] / era_stratum["era_n"]
    )
    write_tsv(era_stratum, outdir / f"H1_stratum_by_era_{ANALYSIS_TAG}.tsv")

    h1_sens = []
    for sens_label, sub in {
        "Exclude preanalytical": h1[~h1["preanalytical_yes"]],
        "Exclude FINAL_WITH_UNCERTAINTY": h1[h1["fully_final"]],
        "CORE_AMR_DIAGNOSTIC only": h1[h1["core_amr"] == 1],
    }.items():
        h1_sens.append(
            ordinal_model(
                sub,
                ["year5"],
                f"H1 sensitivity: {sens_label}",
            )
        )
    beta = h1.copy()
    beta.loc[beta["nct_id"] == BETALACTA_NCT, "depth"] = 1
    h1_sens.append(
        ordinal_model(
            beta,
            ["year5"],
            "H1 sensitivity: BetaLACTA depth 2->1",
        )
    )
    write_tsv(pd.DataFrame(h1_sens), outdir / f"H1_sensitivity_models_{ANALYSIS_TAG}.tsv")

    # H2
    h2 = d[d["h2_group"].isin(["GRAM_POSITIVE", "ENTEROBACTERALES"])].copy()

    h2_dist = h2.groupby(["h2_group", "depth"]).size().reset_index(name="n")
    h2_totals = h2.groupby("h2_group").size().to_dict()
    h2_dist["group_n"] = h2_dist["h2_group"].map(h2_totals)
    h2_dist["percent_within_group"] = 100 * h2_dist["n"] / h2_dist["group_n"]
    write_tsv(h2_dist, outdir / f"H2_depth_distribution_{ANALYSIS_TAG}.tsv")

    h2_perm = pd.DataFrame(
        [
            permutation_mean_depth(
                h2,
                "h2_group",
                "ENTEROBACTERALES",
                "GRAM_POSITIVE",
                n_perm,
                seed,
                "H2 primary omnibus ordinal permutation test",
            )
        ]
    )
    write_tsv(h2_perm, outdir / f"H2_primary_ordinal_permutation_{ANALYSIS_TAG}.tsv")

    h2_exact = pd.DataFrame(
        [
            fisher_group_outcome(
                h2,
                "h2_group",
                "ENTEROBACTERALES",
                "GRAM_POSITIVE",
                "depth_ge1",
                "H2 exact: depth>=1",
            ),
            fisher_group_outcome(
                h2,
                "h2_group",
                "ENTEROBACTERALES",
                "GRAM_POSITIVE",
                "depth_ge2",
                "H2 exact: depth>=2",
            ),
            fisher_group_outcome(
                h2,
                "h2_group",
                "ENTEROBACTERALES",
                "GRAM_POSITIVE",
                "utility_yes",
                "H2 exact: clinical utility",
            ),
        ]
    )
    write_tsv(h2_exact, outdir / f"H2_primary_threshold_fisher_{ANALYSIS_TAG}.tsv")

    h2_adj = pd.DataFrame(
        [
            logistic_model(
                "depth_ge1 ~ is_entero + year5",
                h2,
                "is_entero",
                "H2 adjusted: depth>=1, Enterobacterales vs Gram-positive",
            ),
            logistic_model(
                "depth_ge2 ~ is_entero + year5",
                h2,
                "is_entero",
                "H2 adjusted: depth>=2, Enterobacterales vs Gram-positive",
            ),
        ]
    )
    write_tsv(h2_adj, outdir / f"H2_year_adjusted_threshold_models_{ANALYSIS_TAG}.tsv")

    h2_po_diag = pd.DataFrame(
        [
            stacked_threshold_diagnostic(
                h2,
                "is_entero",
                [],
                "H2 proportional-odds diagnostic, unadjusted",
            ),
            stacked_threshold_diagnostic(
                h2,
                "is_entero",
                ["year5"],
                "H2 proportional-odds diagnostic, adjusted for start year",
            ),
        ]
    )
    write_tsv(h2_po_diag, outdir / f"H2_proportional_odds_diagnostic_{ANALYSIS_TAG}.tsv")

    h2_support_po = pd.DataFrame(
        [
            ordinal_model(
                h2,
                ["is_entero"],
                "H2 supporting proportional-odds model",
            ),
            ordinal_model(
                h2,
                ["is_entero", "year5"],
                "H2 supporting proportional-odds model adjusted for year",
            ),
        ]
    )
    write_tsv(h2_support_po, outdir / f"H2_supporting_proportional_odds_models_{ANALYSIS_TAG}.tsv")

    h2r = h2[~h2["rare_h2_sensitivity"]].copy()
    h2r_perm = pd.DataFrame(
        [
            permutation_mean_depth(
                h2r,
                "h2_group",
                "ENTEROBACTERALES",
                "GRAM_POSITIVE",
                n_perm,
                seed + 11,
                "H2 sensitivity omnibus excluding typhoid/plague",
            )
        ]
    )
    write_tsv(h2r_perm, outdir / f"H2_rare_pathogen_ordinal_permutation_{ANALYSIS_TAG}.tsv")

    h2r_exact = pd.DataFrame(
        [
            fisher_group_outcome(
                h2r,
                "h2_group",
                "ENTEROBACTERALES",
                "GRAM_POSITIVE",
                "depth_ge1",
                "H2 sensitivity exact depth>=1",
            ),
            fisher_group_outcome(
                h2r,
                "h2_group",
                "ENTEROBACTERALES",
                "GRAM_POSITIVE",
                "depth_ge2",
                "H2 sensitivity exact depth>=2",
            ),
            fisher_group_outcome(
                h2r,
                "h2_group",
                "ENTEROBACTERALES",
                "GRAM_POSITIVE",
                "utility_yes",
                "H2 sensitivity exact clinical utility",
            ),
        ]
    )
    write_tsv(h2r_exact, outdir / f"H2_rare_pathogen_threshold_fisher_{ANALYSIS_TAG}.tsv")

    h2r_adj = pd.DataFrame(
        [
            logistic_model(
                "depth_ge1 ~ is_entero + year5",
                h2r,
                "is_entero",
                "H2 sensitivity adjusted depth>=1",
            ),
            logistic_model(
                "depth_ge2 ~ is_entero + year5",
                h2r,
                "is_entero",
                "H2 sensitivity adjusted depth>=2",
            ),
        ]
    )
    write_tsv(
        h2r_adj,
        outdir / f"H2_rare_pathogen_year_adjusted_threshold_models_{ANALYSIS_TAG}.tsv",
    )

    modalities = sorted(h2["final_primary_diagnostic_modality"].map(text).unique())
    h2_modality = pd.DataFrame(
        [
            permutation_binary_group_by_category(
                h2["h2_group"].to_numpy(),
                h2["final_primary_diagnostic_modality"].map(text).to_numpy(),
                "ENTEROBACTERALES",
                "GRAM_POSITIVE",
                modalities,
                n_perm,
                seed + 21,
                "H2 overall primary-modality distribution",
            )
        ]
    )
    write_tsv(
        h2_modality,
        outdir / f"H2_primary_modality_permutation_{ANALYSIS_TAG}.tsv",
    )

    write_tsv(
        token_level_fisher(
            h2,
            "final_all_diagnostic_modalities",
            "h2_group",
            "ENTEROBACTERALES",
            "GRAM_POSITIVE",
            "all_diagnostic_modalities",
        ),
        outdir / f"H2_modality_token_fisher_FDR_{ANALYSIS_TAG}.tsv",
    )
    write_tsv(
        token_level_fisher(
            h2,
            "final_clinical_utility_endpoint_categories",
            "h2_group",
            "ENTEROBACTERALES",
            "GRAM_POSITIVE",
            "clinical_utility_endpoint_categories",
        ),
        outdir / f"H2_utility_token_fisher_FDR_{ANALYSIS_TAG}.tsv",
    )

    # H3
    depth4 = int((d["depth"] == 4).sum())
    upper = (
        1 - 0.05 ** (1 / len(d))
        if depth4 == 0
        else np.nan
    )
    h3 = pd.DataFrame(
        [
            {
                "n": len(d),
                "depth4_n": depth4,
                "depth4_percent": 100 * depth4 / len(d),
                "one_sided_exact_95_upper_proportion": upper,
                "one_sided_exact_95_upper_percent": 100 * upper,
                "depth3_n": int((d["depth"] == 3).sum()),
            }
        ]
    )
    write_tsv(h3, outdir / f"H3_depth4_exact_bound_{ANALYSIS_TAG}.tsv")

    title_col = (
        "registry_brief_title"
        if "registry_brief_title" in d.columns
        else "brief_title"
        if "brief_title" in d.columns
        else "nct_id"
    )
    h3_cols = [
        c for c in [
            "nct_id",
            title_col,
            "depth",
            "final_primary_diagnostic_modality",
            "final_clinical_utility_any",
            "final_analytical_endpoint_categories",
            "special_sensitivity_flags",
            "descriptor_adjudication_notes",
        ]
        if c in d.columns
    ]
    h3_records = d[(d["depth"] >= 3) | d["h3_near_miss"]][h3_cols].copy()
    write_tsv(
        h3_records,
        outdir / f"H3_depth3_and_frozen_quantitative_near_miss_records_{ANALYSIS_TAG}.tsv",
    )

    # H4
    h4_depth = utility_by_depth(d, label)
    write_tsv(h4_depth, outdir / f"H4_utility_by_depth_{ANALYSIS_TAG}.tsv")

    h4_perm = pd.DataFrame(
        [
            permutation_chisq(
                d["depth"].to_numpy(),
                d["utility_yes"].to_numpy(),
                [0, 1, 2, 3],
                n_perm,
                seed + 31,
                "H4 primary omnibus categorical depth x utility",
            )
        ]
    )
    write_tsv(h4_perm, outdir / f"H4_primary_categorical_permutation_{ANALYSIS_TAG}.tsv")

    h4_pair = pd.DataFrame(
        [
            fisher_group_outcome(
                d,
                "depth",
                1,
                0,
                "utility_yes",
                "H4 exact: depth 1 vs depth 0",
            ),
            fisher_group_outcome(
                d,
                "depth",
                2,
                0,
                "utility_yes",
                "H4 exact: depth 2 vs depth 0",
            ),
            fisher_group_outcome(
                d,
                "depth",
                2,
                1,
                "utility_yes",
                "H4 exact: depth 2 vs depth 1",
            ),
            fisher_group_outcome(
                d,
                "depth",
                3,
                2,
                "utility_yes",
                "H4 descriptive exact: depth 3 vs depth 2",
            ),
        ]
    )
    write_tsv(h4_pair, outdir / f"H4_pairwise_depth_fisher_{ANALYSIS_TAG}.tsv")

    d02 = d[d["depth"] <= 2].copy()
    shape_row = {
        "analysis": "H4 linear-vs-categorical depth shape test restricted to depths 0-2",
        "n": len(d02),
        "likelihood_ratio": np.nan,
        "df": np.nan,
        "p_value": np.nan,
        "status": "FAIL",
        "note": (
            "Depth 3 excluded from this likelihood-ratio shape check because "
            "n=3; depth 3 remains included descriptively and in the omnibus "
            "permutation test."
        ),
    }
    try:
        lin02 = smf.glm(
            "utility_yes ~ depth",
            data=d02,
            family=sm.families.Binomial(),
        ).fit()
        cat02 = smf.glm(
            "utility_yes ~ C(depth)",
            data=d02,
            family=sm.families.Binomial(),
        ).fit()
        lr = 2 * (cat02.llf - lin02.llf)
        df_lr = int(cat02.df_model - lin02.df_model)
        shape_row.update(
            likelihood_ratio=lr,
            df=df_lr,
            p_value=float(chi2.sf(lr, df_lr)),
            status="PASS",
        )
    except Exception as exc:
        shape_row["note"] += f" Error: {exc!r}"
    shape = pd.DataFrame([shape_row])
    write_tsv(shape, outdir / f"H4_depth_shape_diagnostic_{ANALYSIS_TAG}.tsv")

    stratum_summary = (
        d.groupby("final_stratum")
        .agg(n=("nct_id", "size"), utility_yes=("utility_yes", "sum"))
        .reset_index()
    )
    stratum_summary["utility_no"] = (
        stratum_summary["n"] - stratum_summary["utility_yes"]
    )
    stratum_summary["utility_percent"] = (
        100 * stratum_summary["utility_yes"] / stratum_summary["n"]
    )
    write_tsv(
        stratum_summary,
        outdir / f"H4_utility_by_stratum_{ANALYSIS_TAG}.tsv",
    )

    h4_stratum = pd.DataFrame(
        [
            fisher_group_outcome(
                d,
                "final_stratum",
                "CORE_AMR_DIAGNOSTIC",
                "BROAD_BACTERIAL_DIAGNOSTIC",
                "utility_yes",
                "H4 structural exact: CORE_AMR vs BROAD bacterial diagnostic",
            ),
            logistic_model(
                "utility_yes ~ core_amr",
                d,
                "core_amr",
                "H4 structural logistic: CORE_AMR vs BROAD",
            ),
            logistic_model(
                "utility_yes ~ core_amr + year5 + C(study_type_clean)",
                d.dropna(subset=["start_year_num"]),
                "core_amr",
                "H4 structural adjusted: CORE_AMR vs BROAD + year + study type",
            ),
        ]
    )
    write_tsv(h4_stratum, outdir / f"H4_stratum_models_{ANALYSIS_TAG}.tsv")

    core = d[d["core_amr"] == 1].copy()
    core_depth = utility_by_depth(core, f"{label}: CORE_AMR_DIAGNOSTIC")
    write_tsv(
        core_depth,
        outdir / f"H4_core_AMR_utility_by_depth_{ANALYSIS_TAG}.tsv",
    )
    core_tests = pd.DataFrame(
        [
            fisher_group_outcome(
                core,
                "depth",
                2,
                1,
                "utility_yes",
                "H4 CORE_AMR exact: depth 2 vs depth 1",
            ),
            fisher_group_outcome(
                core,
                "depth",
                3,
                2,
                "utility_yes",
                "H4 CORE_AMR descriptive exact: depth 3 vs depth 2",
            ),
        ]
    )
    write_tsv(
        core_tests,
        outdir / f"H4_core_AMR_depth_fisher_{ANALYSIS_TAG}.tsv",
    )

    h4_trend = pd.DataFrame(
        [
            logistic_model(
                "utility_yes ~ depth",
                d,
                "depth",
                "H4 secondary trend: utility per depth level",
            ),
            logistic_model(
                "utility_yes ~ depth + year5 + C(study_type_clean)",
                d.dropna(subset=["start_year_num"]),
                "depth",
                "H4 secondary adjusted trend: depth + year + study type",
            ),
        ]
    )
    write_tsv(
        h4_trend,
        outdir / f"H4_secondary_linear_trend_models_{ANALYSIS_TAG}.tsv",
    )

    h4_sens = []
    for sens_label, sub in {
        "Exclude preanalytical": d[~d["preanalytical_yes"]],
        "Exclude FINAL_WITH_UNCERTAINTY": d[d["fully_final"]],
        "CORE_AMR_DIAGNOSTIC only": d[d["core_amr"] == 1],
    }.items():
        h4_sens.append(
            logistic_model(
                "utility_yes ~ depth",
                sub,
                "depth",
                f"H4 secondary trend sensitivity: {sens_label}",
            )
        )
    beta_h4 = d.copy()
    beta_h4.loc[beta_h4["nct_id"] == BETALACTA_NCT, "depth"] = 1
    h4_sens.append(
        logistic_model(
            "utility_yes ~ depth",
            beta_h4,
            "depth",
            "H4 secondary trend sensitivity: BetaLACTA depth 2->1",
        )
    )
    write_tsv(
        pd.DataFrame(h4_sens),
        outdir / f"H4_secondary_trend_sensitivity_models_{ANALYSIS_TAG}.tsv",
    )

    if make_figures:
        figure_depth_by_era(
            d,
            outdir / f"figure_H1_depth_by_era_{ANALYSIS_TAG}.png",
        )
        figure_h2_depth(
            d,
            outdir / f"figure_H2_depth_by_group_{ANALYSIS_TAG}.png",
        )
        figure_h4_depth(
            h4_depth,
            outdir / f"figure_H4_utility_by_depth_{ANALYSIS_TAG}.png",
        )
        figure_h4_stratum(
            stratum_summary,
            outdir / f"figure_H4_utility_by_stratum_{ANALYSIS_TAG}.png",
        )

    # Neutral numeric extraction table. No conclusion text.
    h1r = h1_primary.iloc[0]
    h2p = h2_perm.iloc[0]
    h2ge1 = h2_exact[h2_exact["outcome"] == "depth_ge1"].iloc[0]
    h2ge2 = h2_exact[h2_exact["outcome"] == "depth_ge2"].iloc[0]
    h2u = h2_exact[h2_exact["outcome"] == "utility_yes"].iloc[0]
    h4om = h4_perm.iloc[0]
    core12 = core_tests[
        core_tests["analysis"] == "H4 CORE_AMR exact: depth 2 vs depth 1"
    ].iloc[0]
    core_vs_broad = h4_stratum[
        h4_stratum["analysis"]
        == "H4 structural exact: CORE_AMR vs BROAD bacterial diagnostic"
    ].iloc[0]

    numeric = pd.DataFrame(
        [
            {
                "analysis_set": label,
                "hypothesis": "H1",
                "analysis": "Diagnostic depth per 5-year increase in start year",
                "hierarchy": "PRIMARY",
                "effect_type": "OR",
                "effect": h1r["odds_ratio"],
                "ci_low": h1r["ci_low"],
                "ci_high": h1r["ci_high"],
                "p_value": h1r["p_value"],
            },
            {
                "analysis_set": label,
                "hypothesis": "H2",
                "analysis": "Enterobacterales vs Gram-positive ordinal depth permutation",
                "hierarchy": "PRIMARY",
                "effect_type": "mean_depth_difference",
                "effect": h2p["mean_depth_difference_a_minus_b"],
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p_value": h2p["permutation_p_value_two_sided"],
            },
            {
                "analysis_set": label,
                "hypothesis": "H2",
                "analysis": "Depth >=1 Enterobacterales vs Gram-positive",
                "hierarchy": "PRIMARY_CHARACTERIZATION",
                "effect_type": "Fisher_OR",
                "effect": h2ge1["odds_ratio_a_vs_b"],
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p_value": h2ge1["p_value"],
            },
            {
                "analysis_set": label,
                "hypothesis": "H2",
                "analysis": "Depth >=2 Enterobacterales vs Gram-positive",
                "hierarchy": "PRIMARY_CHARACTERIZATION",
                "effect_type": "Fisher_OR",
                "effect": h2ge2["odds_ratio_a_vs_b"],
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p_value": h2ge2["p_value"],
            },
            {
                "analysis_set": label,
                "hypothesis": "H2",
                "analysis": "Clinical utility Enterobacterales vs Gram-positive",
                "hierarchy": "CHARACTERIZATION",
                "effect_type": "Fisher_OR",
                "effect": h2u["odds_ratio_a_vs_b"],
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p_value": h2u["p_value"],
            },
            {
                "analysis_set": label,
                "hypothesis": "H3",
                "analysis": "Depth-4 quantitative AMR mechanism evaluation",
                "hierarchy": "PRIMARY",
                "effect_type": "prevalence_n",
                "effect": depth4,
                "ci_low": np.nan,
                "ci_high": float(h3.iloc[0]["one_sided_exact_95_upper_percent"]),
                "p_value": np.nan,
            },
            {
                "analysis_set": label,
                "hypothesis": "H4",
                "analysis": "Categorical diagnostic depth x clinical utility",
                "hierarchy": "PRIMARY",
                "effect_type": "permutation_chi_square",
                "effect": h4om["chi_square"],
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p_value": h4om["permutation_p_value"],
            },
            {
                "analysis_set": label,
                "hypothesis": "H4",
                "analysis": "CORE_AMR vs BROAD clinical utility",
                "hierarchy": "STRUCTURAL",
                "effect_type": "Fisher_OR",
                "effect": core_vs_broad["odds_ratio_a_vs_b"],
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p_value": core_vs_broad["p_value"],
            },
            {
                "analysis_set": label,
                "hypothesis": "H4",
                "analysis": "Within CORE_AMR: depth 2 vs depth 1 clinical utility",
                "hierarchy": "STRUCTURAL",
                "effect_type": "Fisher_OR",
                "effect": core12["odds_ratio_a_vs_b"],
                "ci_low": np.nan,
                "ci_high": np.nan,
                "p_value": core12["p_value"],
            },
        ]
    )
    write_tsv(
        numeric,
        outdir / f"TABLE_MAIN_HYPOTHESIS_NUMERIC_RESULTS_{ANALYSIS_TAG}.tsv",
    )

    required_status = {
        "H1_primary_model": text(h1_primary.iloc[0]["status"]),
        "H1_PO_diagnostic": text(h1_diag.iloc[0]["status"]),
        "H2_PO_diagnostic_unadjusted": text(h2_po_diag.iloc[0]["status"]),
        "H2_PO_diagnostic_adjusted": text(h2_po_diag.iloc[1]["status"]),
        "H4_shape_diagnostic": text(shape.iloc[0]["status"]),
    }
    failed = [k for k, v in required_status.items() if v != "PASS"]
    if failed:
        raise SystemExit(
            f"{label}: required model/diagnostic failure: {', '.join(failed)}"
        )

    return {
        "analysis_set": label,
        "n": len(d),
        "start_year_complete_case_n": int(d["start_year_num"].notna().sum()),
        "h2_primary_n": len(h2),
        "h2_rare_pathogen_sensitivity_n": len(h2r),
        "rare_h2_flagged_n": int(d["rare_h2_sensitivity"].sum()),
        "depth_distribution": depth_counter(d["final_amr_depth"]),
        "h2_distribution": normalized_counter(d["final_h2_comparison_group"]),
        "utility_distribution": normalized_counter(d["final_clinical_utility_any"]),
        "stratum_distribution": normalized_counter(d["final_stratum"]),
        "depth4_n": depth4,
        "depth3_n": int((d["depth"] == 3).sum()),
        "frozen_h3_near_miss_flag_n": int(d["h3_near_miss"].sum()),
        "required_status": required_status,
    }


def build_manifest(root: Path, path: Path) -> None:
    rows = []
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        if p == path or p.name == "SHA256SUMS.txt":
            continue
        rel = p.relative_to(root)
        parsed_rows = ""
        if p.suffix == ".tsv":
            try:
                parsed_rows = str(len(read_tsv(p)))
            except Exception:
                parsed_rows = "PARSE_ERROR"
        rows.append(
            {
                "relative_path": str(rel),
                "sha256": sha256(p),
                "bytes": p.stat().st_size,
                "parsed_rows_if_tsv": parsed_rows,
            }
        )
    write_tsv(pd.DataFrame(rows), path)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--primary", required=True, type=Path)
    ap.add_argument("--p35", required=True, type=Path)
    ap.add_argument("--h3-new-qc", type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--permutations", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=3312026)
    args = ap.parse_args()

    for p in [args.primary, args.p35]:
        if not p.exists():
            raise SystemExit(f"Required input missing: {p}")
    if args.h3_new_qc and not args.h3_new_qc.exists():
        raise SystemExit(f"H3 QC file not found: {args.h3_new_qc}")
    if args.permutations < 1000:
        raise SystemExit("--permutations must be >=1000")

    out = args.output_dir.expanduser().resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(
            f"Output directory is not empty: {out}\n"
            "Refusing to overwrite a prior analysis."
        )
    out.mkdir(parents=True, exist_ok=True)

    warnings.filterwarnings("ignore")

    primary = read_tsv(args.primary)
    p35 = read_tsv(args.p35)
    validate_primary(primary, args.primary)
    validate_p35(primary, p35, args.p35)

    h3_qc_summary = None
    if args.h3_new_qc:
        require_equal(
            sha256(args.h3_new_qc),
            EXPECTED_H3_QC_SHA256,
            "New H3 heuristic QC SHA-256",
        )
        q = read_tsv(args.h3_new_qc)
        required_q = {
            "nct_id",
            "final_amr_depth",
            "qc_status",
            "analytic_effect",
        }
        if not required_q <= set(q.columns):
            raise SystemExit(
                f"H3 QC queue missing columns: {sorted(required_q-set(q.columns))}"
            )
        if not (q["qc_status"].map(text) == "HEURISTIC_CANDIDATE_ONLY").all():
            raise SystemExit("Unexpected H3 QC status")
        if not (q["analytic_effect"].map(text) == "NONE").all():
            raise SystemExit("H3 QC queue has non-NONE analytic effect")
        shutil.copy2(
            args.h3_new_qc,
            out / f"H3_NEW_HEURISTIC_CANDIDATES_NONBINDING_{ANALYSIS_TAG}.tsv",
        )
        h3_qc_summary = {
            "rows": len(q),
            "unique_nct_ids": q["nct_id"].nunique(),
            "input_sha256": sha256(args.h3_new_qc),
            "analytic_effect": "NONE",
        }

    carryforward = f"""ClinicalTrials.gov bacterial/AMR diagnostic landscape
Analysis methods carryforward {ANALYSIS_VERSION}

Frozen data release: {RELEASE_VERSION}
Predecessor inferential hierarchy: {PREDECESSOR_HIERARCHY}
Created: {now_utc()}

This is not a new post-result statistical amendment.

The v3.2.8 hierarchy is carried forward unchanged:
- H1 proportional odds primary; cumulative-threshold models supporting.
- H2 ordinal permutation omnibus and threshold Fisher tests primary; adjusted
  threshold models characterize year adjustment; proportional odds supporting;
  typhoid/plague sensitivity retained.
- H3 depth-4 prevalence/exact one-sided upper bound unchanged.
- H4 categorical/permutation primary with depth-specific and structural
  analyses; numeric-depth trend models remain secondary.

The only planned additions are:
1. analysis of the frozen expanded 573-study v3.3.0 cohort, and
2. repetition under the separately frozen P35 imaging sensitivity cohort.

The 15-study new H3 keyword queue, when supplied, is nonbinding QC provenance
only. It does not alter depth, near-miss flags, or the H3 denominator.

No manuscript interpretation is generated automatically. Interpretation and
comparison with v3.2.8 occur only after this v3.3.1 output is complete,
checksummed, and frozen.
"""
    (out / f"ANALYSIS_METHODS_CARRYFORWARD_{ANALYSIS_TAG}.txt").write_text(
        carryforward,
        encoding="utf-8",
    )

    primary_dir = out / "primary_573"
    p35_dir = out / "P35_sensitivity_565"

    primary_summary = run_analysis_set(
        primary,
        primary_dir,
        "PRIMARY_573",
        args.permutations,
        args.seed,
        make_figures=True,
    )
    p35_summary = run_analysis_set(
        p35,
        p35_dir,
        "P35_SENSITIVITY_565",
        args.permutations,
        args.seed + 1000,
        make_figures=False,
    )

    # Combine neutral numeric summaries for later human comparison.
    pmain = read_tsv(
        primary_dir / f"TABLE_MAIN_HYPOTHESIS_NUMERIC_RESULTS_{ANALYSIS_TAG}.tsv"
    )
    psens = read_tsv(
        p35_dir / f"TABLE_MAIN_HYPOTHESIS_NUMERIC_RESULTS_{ANALYSIS_TAG}.tsv"
    )
    write_tsv(
        pd.concat([pmain, psens], ignore_index=True),
        out / f"TABLE_PRIMARY_AND_P35_NUMERIC_RESULTS_{ANALYSIS_TAG}.tsv",
    )

    summary = {
        "created_at": now_utc(),
        "analysis_version": ANALYSIS_VERSION,
        "frozen_release_version": RELEASE_VERSION,
        "predecessor_statistical_hierarchy": PREDECESSOR_HIERARCHY,
        "analysis_pass": True,
        "permutations": args.permutations,
        "primary_seed": args.seed,
        "p35_seed": args.seed + 1000,
        "input_sha256": {
            "primary": sha256(args.primary),
            "p35": sha256(args.p35),
            "h3_new_qc": (
                sha256(args.h3_new_qc)
                if args.h3_new_qc
                else None
            ),
            "analysis_script": sha256(Path(__file__).resolve()),
        },
        "primary": primary_summary,
        "p35_sensitivity": p35_summary,
        "p35_removed_nct_ids": sorted(EXPECTED_P35_EXCLUDED_IDS),
        "new_h3_heuristic_qc": h3_qc_summary,
        "interpretation_boundary": (
            "Numerical analysis only. Freeze/checksum output before comparison "
            "with v3.2.8 or manuscript revision."
        ),
    }
    summary_path = out / f"FINAL_ANALYSIS_RUN_SUMMARY_{ANALYSIS_TAG}.json"
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    readme = f"""ClinicalTrials.gov bacterial/AMR diagnostic landscape
Final refined H1-H4 analysis {ANALYSIS_VERSION}

Inputs
------
Primary: frozen {RELEASE_VERSION} 573-study cohort
P35: frozen {RELEASE_VERSION} 565-study sensitivity cohort

Method
------
The inferential hierarchy is carried forward from v3.2.8 without post-result
re-ranking. Primary results are under primary_573/. The P35 sensitivity run is
under P35_sensitivity_565/.

Do not use the P35 directory as a replacement primary analysis. It is a
prespecified sensitivity analysis.

Do not change study coding, depth, descriptors, P35 membership, or statistical
hierarchy in response to these results. Any genuine correction requires a new
version and explicit provenance.

The root numeric table is intentionally interpretation-free. Compare with
v3.2.8 only after this directory is checksummed and made read-only.
"""
    (out / f"README_ANALYSIS_{ANALYSIS_TAG}.txt").write_text(
        readme,
        encoding="utf-8",
    )

    manifest = out / f"analysis_manifest_{ANALYSIS_TAG}.tsv"
    build_manifest(out, manifest)

    sums = out / "SHA256SUMS.txt"
    files = [
        p for p in sorted(out.rglob("*"))
        if p.is_file() and p != sums
    ]
    with sums.open("w", encoding="utf-8") as handle:
        for p in files:
            handle.write(f"{sha256(p)}  {p.relative_to(out)}\n")

    # Console output intentionally contains no inferential estimates or p-values.
    print(f"REFINED H1-H4 ANALYSIS {ANALYSIS_VERSION}: PASS")
    print(f"Primary n: {primary_summary['n']}")
    print(f"P35 sensitivity n: {p35_summary['n']}")
    print(f"Primary H2 n: {primary_summary['h2_primary_n']}")
    print(
        "Primary H2 rare-pathogen sensitivity n: "
        f"{primary_summary['h2_rare_pathogen_sensitivity_n']}"
    )
    print(f"Primary depth 4 n: {primary_summary['depth4_n']}")
    print(f"Primary depth 3 n: {primary_summary['depth3_n']}")
    print(f"Output directory: {out}")
    print("Freeze/checksum this directory before inspecting inferential tables.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
