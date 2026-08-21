#!/usr/bin/env python3
"""
Prepare the post-screening descriptor stage for the ClinicalTrials.gov
bacterial AMR diagnostic landscape.

Version family
--------------
Screening/adjudication freeze: v3.2.9
Descriptor-stage preparation: v3.3.0-prep

This script DOES NOT:
- reopen screening;
- change any frozen v3.2.9 eligibility/stratum/depth decision;
- change any historic v3.2.7 descriptor;
- run H1-H4 statistics.

It DOES:
1. Validate the historic and v3.2.9 screening inputs.
2. Construct a full 4,236-record review-provenance master.
3. Preserve the existing 291-row v3.2.9 final-decision ledger as an
   adjudication/override layer.
4. Identify exactly 213 newly eligible studies.
5. Recover frozen ClinicalTrials.gov metadata/evidence for those 213 NCTs.
6. Generate two independently shuffled blinded canonical-descriptor packets,
   one for SB and one for ZB.
7. Create a private descriptor-review key and carry-forward QC flags.
8. Create a full projected-cohort P35 imaging-sensitivity candidate queue.
9. Create deterministic reviewer tar.gz packages and SHA-256 manifests.

Expected final screening counts
-------------------------------
Full reviewed universe: 4,236
Eligible: 573 = 360 historic + 213 newly eligible
Excluded: 3,647
Uncertain: 16 = 8 historic + 8 newly uncertain

Newly eligible depth:
d0 = 177
d1 = 23
d2 = 13
d3 = 0
d4 = 0

Projected final cohort depth before descriptor coding:
d0 = 435
d1 = 83
d2 = 52
d3 = 3
d4 = 0
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import random
import re
import shutil
import sys
import tarfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


VERSION = "v3.3.0-prep"
SCREENING_VERSION = "v3.2.9"

EXPECTED = {
    "historic_all": 2097,
    "historic_eligible": 360,
    "ib_rescue": 1332,
    "ib_nonflag": 400,
    "ib_registry_negative": 300,
    "ib_total": 2032,
    "tb_expansion": 107,
    "new_total": 2139,
    "full_review_master": 4236,
    "ledger": 291,
    "new_eligible": 213,
    "new_uncertain": 8,
    "new_excluded": 1918,
    "projected_eligible": 573,
    "projected_excluded": 3647,
    "projected_uncertain": 16,
}

EXPECTED_NEW_DEPTH = {"0": 177, "1": 23, "2": 13, "3": 0, "4": 0}
EXPECTED_PROJECTED_DEPTH = {"0": 435, "1": 83, "2": 52, "3": 3, "4": 0}

SB_SEED = 3302026
ZB_SEED = 3312026

PRIMARY_STRATA = {"CORE_AMR_DIAGNOSTIC", "BROAD_BACTERIAL_DIAGNOSTIC"}

SCREENING_STRATA = {
    "CORE_AMR_DIAGNOSTIC",
    "BROAD_BACTERIAL_DIAGNOSTIC",
    "HOST_RESPONSE_DIAGNOSTIC",
    "CLINICAL_SYNDROMIC_SUPPORT",
    "MECHANISM_SUPPORT",
    "SPECIAL_PATHOGEN_DIAGNOSTIC",
    "SURVEILLANCE_SUPPORT",
    "THERAPEUTIC_SUPPORT",
    "NONINFECTIOUS_OR_UNRELATED",
    "OTHER",
}

# Historic v3.2.7 used HOST_RESPONSE_SUPPORT; v3.2.9 uses
# HOST_RESPONSE_DIAGNOSTIC. Preserve the historic label in the provenance master
# without treating that versioned vocabulary difference as an error.
MASTER_ALLOWED_STRATA = SCREENING_STRATA | {"HOST_RESPONSE_SUPPORT"}

SCREENING_EXCLUSION_REASONS = {
    "NONBACTERIAL_OR_UNRELATED",
    "THERAPEUTIC_ONLY",
    "SURVEILLANCE_ONLY",
    "MECHANISM_ONLY",
    "HOST_RESPONSE_ONLY",
    "CLINICAL_SYNDROMIC_ONLY",
    "SPECIAL_PATHOGEN_SEPARATE",
    "PREVENTION_OR_VACCINE_ONLY",
    "COLONIZATION_WITHOUT_DIAGNOSTIC_EVALUATION",
    "NO_DIRECT_DIAGNOSTIC_EVALUATION",
    "INSUFFICIENT_INFORMATION",
}

DESCRIPTOR_REVIEW_FIELDS = [
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

DESCRIPTOR_PACKET_FIELDS = [
    "descriptor_review_id",
    "nct_id",
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
    *DESCRIPTOR_REVIEW_FIELDS,
]

P35_TERMS = {
    "PET": r"\bPET(?:/CT)?\b",
    "SPECT": r"\bSPECT\b",
    "MRI": r"\bMRI\b|\bmagnetic resonance\b",
    "ULTRASOUND": r"\bultrasound\b|\bultrason",
    "IMAGING": r"\bimaging\b|\bimage[- ]guided\b",
    "MICROSCOPY": r"\bmicroscop",
    "OPTICAL": r"\boptical\b|\bendomicroscop",
    "FLUORESCENCE": r"\bfluorescen",
    "SPECTROSCOPY": r"\bspectroscop",
    "TOMOGRAPHY": r"\btomograph",
    "RADIOGRAPHY": r"\bradiograph|\bx[- ]?ray\b",
}

HISTORIC_EXPECTED_HASH = (
    "a59a9ec30d188533c2e4508bb8044150fb86e563208273cf10c4563b1543bda6"
)
FINAL_LEDGER_EXPECTED_HASH = (
    "0369a72cb625aa07f276d98a30923ab3f6fff4f72b91e29d8fa23060eb53d76c"
)
SCREEN_CODEBOOK_EXPECTED_HASH = (
    "dd08d32f497717503a3135670e565b7499304fad27246f77b4872cbe827e9fa4"
)
DESCRIPTOR_CODEBOOK_EXPECTED_HASH = (
    "a984f470979f0b914de67a4d989573dd84725f6f74fe9a3801530487a1387d6c"
)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def text(v: Any) -> str:
    return str(v or "").strip()


def upper(v: Any) -> str:
    return text(v).upper()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def set_csv_limit() -> None:
    n = sys.maxsize
    while True:
        try:
            csv.field_size_limit(n)
            return
        except OverflowError:
            n //= 10


def read_tsv(path: Path) -> pd.DataFrame:
    set_csv_limit()
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        encoding="utf-8-sig",
    )


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, lineterminator="\n")


def nested(obj: Any, *keys: str, default: Any = "") -> Any:
    cur = obj
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def join_unique(values: Iterable[Any]) -> str:
    out = []
    for x in values or []:
        s = text(x)
        if s and s not in out:
            out.append(s)
    return " | ".join(out)


def date_year(v: Any) -> str:
    m = re.match(r"^(\d{4})", text(v))
    return m.group(1) if m else ""


def outcome_text(outcomes: Any) -> str:
    rows = []
    for x in outcomes or []:
        if not isinstance(x, dict):
            continue
        parts = [
            text(x.get("measure")),
            text(x.get("description")),
            text(x.get("timeFrame")),
        ]
        val = " | ".join(p for p in parts if p)
        if val:
            rows.append(val)
    return " || ".join(rows)


def unwrap_study(obj: Any) -> Any:
    if isinstance(obj, dict) and "protocolSection" in obj:
        return obj
    if isinstance(obj, dict) and isinstance(obj.get("study"), dict):
        return obj["study"]
    return obj


def extract_registry_descriptor_evidence(study: dict, source_shard: str) -> dict:
    study = unwrap_study(study)
    p = study.get("protocolSection", {}) or {}
    ident = p.get("identificationModule", {}) or {}
    status = p.get("statusModule", {}) or {}
    design = p.get("designModule", {}) or {}
    sponsor = p.get("sponsorCollaboratorsModule", {}) or {}
    cond = p.get("conditionsModule", {}) or {}
    arms = p.get("armsInterventionsModule", {}) or {}
    outcomes = p.get("outcomesModule", {}) or {}
    desc = p.get("descriptionModule", {}) or {}
    contacts = p.get("contactsLocationsModule", {}) or {}

    interventions = arms.get("interventions", []) or []
    locations = contacts.get("locations", []) or []
    lead = sponsor.get("leadSponsor", {}) or {}
    enrollment = design.get("enrollmentInfo", {}) or {}

    brief_summary = text(desc.get("briefSummary"))
    detailed = text(desc.get("detailedDescription"))
    summary_parts = []
    for piece in [brief_summary, detailed]:
        if piece and piece not in summary_parts:
            summary_parts.append(piece)

    countries = sorted(
        {
            text(x.get("country"))
            for x in locations
            if isinstance(x, dict) and text(x.get("country"))
        }
    )

    start_date = text(nested(status, "startDateStruct", "date", default=""))
    nct = upper(ident.get("nctId"))

    return {
        "nct_id": nct,
        "brief_title": text(ident.get("briefTitle")),
        "official_title": text(ident.get("officialTitle")),
        "conditions": join_unique(cond.get("conditions", [])),
        "keywords": join_unique(cond.get("keywords", [])),
        "intervention_names": join_unique(
            x.get("name") for x in interventions if isinstance(x, dict)
        ),
        "intervention_types": join_unique(
            x.get("type") for x in interventions if isinstance(x, dict)
        ),
        "primary_outcomes": outcome_text(outcomes.get("primaryOutcomes", [])),
        "secondary_outcomes": outcome_text(outcomes.get("secondaryOutcomes", [])),
        "other_outcomes": outcome_text(outcomes.get("otherOutcomes", [])),
        "summary": " | ".join(summary_parts),
        "clinicaltrials_url": f"https://clinicaltrials.gov/study/{nct}",
        "study_type": text(design.get("studyType")),
        "overall_status": text(status.get("overallStatus")),
        "start_date": start_date,
        "start_year": date_year(start_date),
        "has_results": (
            "YES"
            if bool(study.get("hasResults"))
            or isinstance(study.get("resultsSection"), dict)
            else "NO"
        ),
        "enrollment_count": text(enrollment.get("count")),
        "enrollment_type": text(enrollment.get("type")),
        "lead_sponsor_name": text(lead.get("name")),
        "lead_sponsor_class": text(lead.get("class")),
        "countries": " | ".join(countries),
        "source_registry_shard": source_shard,
    }


def scan_registry_targets(parts_dir: Path, targets: set[str]) -> dict[str, dict]:
    shards = sorted(parts_dir.glob("part-*.jsonl.gz"))
    if not shards:
        raise SystemExit(f"No part-*.jsonl.gz found in {parts_dir}")

    found: dict[str, dict] = {}
    for i, shard in enumerate(shards, start=1):
        with gzip.open(shard, "rt", encoding="utf-8") as f:
            for line_no, line in enumerate(f, start=1):
                if not line.strip():
                    continue
                try:
                    obj = json.loads(line)
                except json.JSONDecodeError as exc:
                    raise SystemExit(f"JSON error {shard}:{line_no}: {exc}")
                study = unwrap_study(obj)
                nct = upper(
                    nested(
                        study,
                        "protocolSection",
                        "identificationModule",
                        "nctId",
                        default="",
                    )
                )
                if nct in targets and nct not in found:
                    found[nct] = extract_registry_descriptor_evidence(
                        study, shard.name
                    )
        print(
            f"[registry {i}/{len(shards)}] found {len(found)}/{len(targets)}",
            file=sys.stderr,
        )
        if len(found) == len(targets):
            break

    missing = sorted(targets - set(found))
    if missing:
        raise SystemExit(
            "Frozen registry recovery incomplete. Missing NCTs:\n"
            + "\n".join(missing)
        )
    return found


def validate_hash(path: Path, expected: str, label: str) -> None:
    observed = sha256(path)
    if observed != expected:
        raise SystemExit(
            f"{label} SHA-256 mismatch\nexpected {expected}\nobserved {observed}\n"
            f"path {path}"
        )


def validate_screen_codebook(path: Path) -> None:
    validate_hash(path, SCREEN_CODEBOOK_EXPECTED_HASH, "Screening codebook")
    cb = read_tsv(path)
    if not {"field", "allowed_values", "instruction"} <= set(cb.columns):
        raise SystemExit("Unexpected screening codebook columns.")
    row = cb.loc[cb["field"].eq("manual_final_stratum")]
    if len(row) != 1:
        raise SystemExit("Screening codebook missing manual_final_stratum.")
    allowed_strata = {x.strip() for x in row.iloc[0]["allowed_values"].split("|")}
    if allowed_strata != SCREENING_STRATA:
        raise SystemExit(
            "Screening stratum vocabulary differs from expected v3.2.9 set."
        )
    row = cb.loc[cb["field"].eq("manual_exclusion_reason")]
    if len(row) != 1:
        raise SystemExit("Screening codebook missing manual_exclusion_reason.")
    allowed_reasons = {x.strip() for x in row.iloc[0]["allowed_values"].split("|")}
    if allowed_reasons != SCREENING_EXCLUSION_REASONS:
        raise SystemExit(
            "Screening exclusion vocabulary differs from expected v3.2.9 set."
        )


def validate_descriptor_codebook(path: Path) -> None:
    validate_hash(path, DESCRIPTOR_CODEBOOK_EXPECTED_HASH, "Descriptor codebook")
    cb = read_tsv(path)
    required = {
        "field",
        "position",
        "allowed_value_or_instruction",
        "multi_select",
    }
    if not required <= set(cb.columns):
        raise SystemExit("Unexpected descriptor codebook columns.")


def normalize_lowrisk_exclusion(stratum: str, source_reason: str) -> str:
    if source_reason in SCREENING_EXCLUSION_REASONS:
        return source_reason

    mapping = {
        "HOST_RESPONSE_DIAGNOSTIC": "HOST_RESPONSE_ONLY",
        "CLINICAL_SYNDROMIC_SUPPORT": "CLINICAL_SYNDROMIC_ONLY",
        "MECHANISM_SUPPORT": "MECHANISM_ONLY",
        "SPECIAL_PATHOGEN_DIAGNOSTIC": "SPECIAL_PATHOGEN_SEPARATE",
        "SURVEILLANCE_SUPPORT": "SURVEILLANCE_ONLY",
        "THERAPEUTIC_SUPPORT": "THERAPEUTIC_ONLY",
        "NONINFECTIOUS_OR_UNRELATED": "NONBACTERIAL_OR_UNRELATED",
        "OTHER": "NO_DIRECT_DIAGNOSTIC_EVALUATION",
    }
    if stratum not in mapping:
        raise SystemExit(
            f"Cannot normalize low-risk exclusion reason: "
            f"stratum={stratum!r} reason={source_reason!r}"
        )
    return mapping[stratum]


def validate_final_logic(row: dict) -> None:
    nct = row["nct_id"]
    elig = row["final_primary_eligible"]
    stratum = row["final_stratum"]
    depth = row["final_amr_depth"]
    reason = row["final_exclusion_reason"]

    if elig == "YES":
        if stratum not in PRIMARY_STRATA:
            raise SystemExit(f"{nct}: YES with nonprimary stratum {stratum}")
        if depth not in {"0", "1", "2", "3", "4"}:
            raise SystemExit(f"{nct}: YES with invalid depth {depth}")
        if reason:
            raise SystemExit(f"{nct}: YES with exclusion reason {reason}")
    elif elig == "NO":
        if stratum not in MASTER_ALLOWED_STRATA - PRIMARY_STRATA:
            raise SystemExit(f"{nct}: NO with invalid stratum {stratum}")
        if depth != "NA":
            raise SystemExit(f"{nct}: NO with depth {depth}")
        if reason not in SCREENING_EXCLUSION_REASONS:
            raise SystemExit(f"{nct}: NO with invalid reason {reason}")
    elif elig == "UNCERTAIN":
        if stratum in PRIMARY_STRATA or stratum not in MASTER_ALLOWED_STRATA:
            raise SystemExit(f"{nct}: UNCERTAIN with invalid stratum {stratum}")
        if depth != "NA":
            raise SystemExit(f"{nct}: UNCERTAIN with depth {depth}")
        if reason != "INSUFFICIENT_INFORMATION":
            raise SystemExit(
                f"{nct}: UNCERTAIN must use INSUFFICIENT_INFORMATION"
            )
    else:
        raise SystemExit(f"{nct}: invalid eligibility {elig!r}")


def historic_master_rows(hist: pd.DataFrame) -> list[dict]:
    rows = []
    for _, r in hist.iterrows():
        row = {
            "nct_id": upper(r["nct_id"]),
            "review_universe_segment": "HISTORIC_V3_2_7",
            "review_source_file": "all_reviewed_final_enriched_v3_2_7.tsv",
            "review_source_id": text(r.get("source_packet_id")),
            "brief_title": text(r.get("brief_title")),
            "clinicaltrials_url": text(r.get("clinicaltrials_url")),
            "source_primary_eligible": text(r.get("final_primary_eligible")),
            "source_stratum": text(r.get("final_stratum")),
            "source_amr_depth": text(r.get("final_amr_depth")),
            "source_exclusion_reason": text(r.get("final_exclusion_reason")),
            "source_reviewer_initials": text(r.get("final_reviewer_initials")),
            "source_review_status": text(r.get("final_review_status")),
            "source_notes": text(r.get("final_notes")),
            "in_v3_2_9_final_291_ledger": "NO",
            "final_decision_source": "HISTORIC_V3_2_7_FROZEN",
            "final_primary_eligible": text(r.get("final_primary_eligible")),
            "final_stratum": text(r.get("final_stratum")),
            "final_amr_depth": text(r.get("final_amr_depth")),
            "final_exclusion_reason": text(r.get("final_exclusion_reason")),
            "final_reviewer_initials": text(r.get("final_reviewer_initials")),
            "final_review_status": text(r.get("final_review_status")),
            "final_decision_basis": "Carried forward unchanged from frozen analytic release v3.2.7.",
            "historic_360_eligible": (
                "YES" if text(r.get("final_primary_eligible")) == "YES" else "NO"
            ),
            "new_eligible_for_descriptor_coding": "NO",
        }
        validate_final_logic(row)
        rows.append(row)
    return rows


def new_source_rows(
    df: pd.DataFrame,
    segment: str,
    source_filename: str,
    reviewer: str,
) -> list[dict]:
    rows = []
    for _, r in df.iterrows():
        source_reason = text(r["manual_exclusion_reason"])
        elig = upper(r["manual_primary_eligible"])
        stratum = upper(r["manual_final_stratum"])
        depth = upper(r["manual_amr_depth"])
        if elig == "NO":
            final_reason = normalize_lowrisk_exclusion(stratum, source_reason)
        elif elig == "UNCERTAIN":
            final_reason = "INSUFFICIENT_INFORMATION"
        else:
            final_reason = ""

        row = {
            "nct_id": upper(r["nct_id"]),
            "review_universe_segment": segment,
            "review_source_file": source_filename,
            "review_source_id": text(r["coverage_audit_id"]),
            "brief_title": text(r["brief_title"]),
            "clinicaltrials_url": text(r["clinicaltrials_url"]),
            "source_primary_eligible": elig,
            "source_stratum": stratum,
            "source_amr_depth": depth,
            "source_exclusion_reason": source_reason,
            "source_reviewer_initials": text(r["reviewer_initials"]),
            "source_review_status": text(r["review_status"]),
            "source_notes": text(r["reviewer_notes"]),
            "in_v3_2_9_final_291_ledger": "NO",
            "final_decision_source": (
                "IB_V3_2_9_LOW_RISK_CLOSED_BY_STOPPING_RULE"
                if reviewer == "IB"
                else "TB_V3_2_9_EXHAUSTIVE_EXPANSION_LOW_RISK"
            ),
            "final_primary_eligible": elig,
            "final_stratum": stratum,
            "final_amr_depth": depth,
            "final_exclusion_reason": final_reason,
            "final_reviewer_initials": reviewer,
            "final_review_status": "FINAL",
            "final_decision_basis": (
                "Low-risk record carried from the completed v3.2.9 reviewer "
                "packet because it did not meet the prespecified high-stakes "
                "confirmation trigger; the source stratum stopping rule was satisfied."
                if reviewer == "IB"
                else
                "Low-risk record carried from the completed exhaustive failed-strata "
                "expansion; it did not meet the prespecified Packet B confirmation trigger."
            ),
            "historic_360_eligible": "NO",
            "new_eligible_for_descriptor_coding": "NO",
        }
        validate_final_logic(row)
        rows.append(row)
    return rows


def apply_ledger_overrides(
    rows_by_nct: dict[str, dict], ledger: pd.DataFrame
) -> None:
    for _, r in ledger.iterrows():
        nct = upper(r["nct_id"])
        if nct not in rows_by_nct:
            raise SystemExit(f"Ledger NCT absent from new review sources: {nct}")
        row = rows_by_nct[nct]
        row["in_v3_2_9_final_291_ledger"] = "YES"
        row["final_decision_source"] = text(r["decision_source"])
        row["final_primary_eligible"] = upper(r["final_primary_eligible"])
        row["final_stratum"] = upper(r["final_stratum"])
        row["final_amr_depth"] = upper(r["final_amr_depth"])
        row["final_exclusion_reason"] = upper(r["final_exclusion_reason"])
        row["final_reviewer_initials"] = text(r["final_adjudicator_initials"])
        row["final_review_status"] = text(r["final_review_status"])
        row["final_decision_basis"] = text(r["final_decision_basis"])
        row["new_eligible_for_descriptor_coding"] = text(
            r["new_eligible_for_descriptor_coding"]
        )
        validate_final_logic(row)


def make_descriptor_packet(
    source: pd.DataFrame, reviewer: str, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    order = list(range(len(source)))
    random.Random(seed).shuffle(order)
    rows = []
    keys = []

    prefix = f"{reviewer}NEWDR"
    depth_output = {
        "0": "ORGANISM_ONLY",
        "1": "BINARY_OR_CATEGORICAL_RESISTANCE",
        "2": "PHENOTYPIC_AST_MIC_ZONE",
        "3": "INTEGRATED_MULTIMECHANISM",
        "4": "QUANTITATIVE_AMR_MECHANISM",
    }

    for seq, idx in enumerate(order, start=1):
        src = source.iloc[idx]
        rid = f"{prefix}-{seq:04d}"
        row = {c: text(src.get(c, "")) for c in DESCRIPTOR_PACKET_FIELDS}
        row["descriptor_review_id"] = rid
        for field in DESCRIPTOR_REVIEW_FIELDS:
            row[field] = ""
        row["index_test_output_type"] = depth_output[src["final_amr_depth"]]
        rows.append(row)
        keys.append(
            {
                "reviewer": reviewer,
                "descriptor_review_id": rid,
                "nct_id": src["nct_id"],
                "original_new_eligible_row_number": str(idx + 2),
            }
        )

    return (
        pd.DataFrame(rows, columns=DESCRIPTOR_PACKET_FIELDS),
        pd.DataFrame(
            keys,
            columns=[
                "reviewer",
                "descriptor_review_id",
                "nct_id",
                "original_new_eligible_row_number",
            ],
        ),
    )


def descriptor_instructions(reviewer: str, n: int) -> str:
    return f"""\
Canonical descriptor dual review for newly eligible v3.2.9 rescue studies
Descriptor-stage version: {VERSION}
Reviewer: {reviewer}
Rows: {n}

Purpose
-------
Independently code the canonical analytic descriptors for the newly confirmed
eligible studies added by screening/adjudication v3.2.9.

This stage DOES NOT reopen:
- primary eligibility;
- BROAD versus CORE final stratum;
- diagnostic depth.

Those fields are frozen and displayed only as context.

Evidence hierarchy
------------------
1. Use the supplied frozen ClinicalTrials.gov evidence first.
2. Open the official ClinicalTrials.gov registration only when a descriptor
   cannot be resolved from the supplied evidence.
3. Do not use publications, manufacturer webpages, package inserts, product
   specifications, or general platform knowledge to infer a feature absent
   from the registration.
4. Set registry_lookup_performed = YES only if you actually open the official
   registration during this descriptor pass.
5. Preserve UNCERTAIN when the registration remains insufficient.

Required fields
---------------
Complete every reviewer field:
- registry_lookup_performed
- primary_diagnostic_modality
- all_diagnostic_modalities
- organism_group
- gram_group
- h2_comparison_group
- analytical_endpoint_categories
- clinical_utility_endpoint_categories
- clinical_utility_any
- preanalytical_flag
- amr_reporting_intervention_flag
- mixed_viral_bacterial_panel_flag
- direct_patient_specimen_flag
- index_test_output_type
- descriptor_notes
- reviewer_initials
- review_status

Binding descriptor rules
------------------------
- Use the supplied canonical_descriptor_codebook_v3_2_5.tsv.
- Pipe-separate multi-select values with | and no surrounding spaces.
- Taxonomic Enterobacterales includes Salmonella and Yersinia for the H2 rule.
- direct_patient_specimen_flag = YES when at least one registered index-test
  stream is performed directly on a patient specimen.
- Positive blood-culture broth is not a direct patient specimen.
- clinical_utility_any = YES only when a registered endpoint evaluates
  patient management, treatment, clinical outcome, infection control, or
  resource use.
- Do not infer assay chemistry, organism coverage, AMR capability, specimen
  matrix, or endpoint meaning from product knowledge.
- index_test_output_type has been prefilled from frozen depth. Do not use
  descriptor review to change frozen depth. If you believe the prefill is
  inconsistent with the registration, explain the concern in descriptor_notes
  and mark NEEDS_DISCUSSION rather than editing eligibility/depth.

Quality control
---------------
- Do not add, remove, reorder, or sort rows.
- Do not edit evidence or frozen final_stratum/final_amr_depth columns.
- Use reviewer_initials = {reviewer} on every row.
- Use review_status = COMPLETE unless genuine descriptor ambiguity remains;
  then use NEEDS_DISCUSSION and explain it.
- Reopen the saved TSV before return and confirm row count and headers.

Return filename
---------------
{reviewer}_new_eligible_canonical_descriptor_review_v3_3_0_completed.tsv

The TSV is the authoritative return. A concise summary memo is welcome but
not required.
"""


def deterministic_targz(source_dir: Path, out_path: Path) -> None:
    """
    Create deterministic tar.gz: fixed member metadata and gzip mtime=0.
    """
    import gzip as _gzip

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as raw:
        with _gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, mtime=0
        ) as gz:
            with tarfile.open(fileobj=gz, mode="w") as tf:
                for p in sorted(source_dir.rglob("*")):
                    if not p.is_file():
                        continue
                    arcname = str(Path(source_dir.name) / p.relative_to(source_dir))
                    data = p.read_bytes()
                    info = tarfile.TarInfo(name=arcname)
                    info.size = len(data)
                    info.mtime = 0
                    info.uid = 0
                    info.gid = 0
                    info.uname = ""
                    info.gname = ""
                    info.mode = 0o644
                    tf.addfile(info, io.BytesIO(data))


def copy_if_supplied(path: Path | None, dest: Path) -> None:
    if path is None:
        return
    if not path.exists():
        raise SystemExit(f"Optional policy file not found: {path}")
    shutil.copy2(path, dest / path.name)


def build_p35_candidate_queue(
    historic_eligible: pd.DataFrame,
    new_source: pd.DataFrame,
    imaging_audit: pd.DataFrame | None,
) -> pd.DataFrame:
    rows = []

    hist_cols = [
        "brief_title",
        "summary",
        "intervention_names",
        "primary_outcomes",
        "secondary_outcomes",
    ]
    for _, r in historic_eligible.iterrows():
        hay = " ".join(text(r.get(c)) for c in hist_cols)
        matches = [name for name, pat in P35_TERMS.items() if re.search(pat, hay, re.I)]
        if not matches:
            continue
        rows.append(
            {
                "nct_id": upper(r["nct_id"]),
                "cohort_source": "HISTORIC_V3_2_7",
                "brief_title": text(r.get("brief_title")),
                "final_stratum": text(r.get("final_stratum")),
                "final_amr_depth": text(r.get("final_amr_depth")),
                "historic_primary_modality": text(
                    r.get("final_primary_diagnostic_modality")
                ),
                "matched_p35_terms": "|".join(matches),
                "prior_v3_2_9_p35_ruling": "",
                "prior_v3_2_9_p35_override_status": "",
                "sensitivity_classification": "REQUIRES_CENTRAL_CLASSIFICATION",
                "sensitivity_rationale": "",
                "central_adjudicator_initials": "",
                "audit_status": "PENDING",
            }
        )

    audit_map = {}
    if imaging_audit is not None:
        for _, r in imaging_audit.iterrows():
            audit_map[upper(r["nct_id"])] = r.to_dict()

    new_cols = [
        "brief_title",
        "summary",
        "intervention_names",
        "primary_outcomes",
        "secondary_outcomes",
    ]
    for _, r in new_source.iterrows():
        hay = " ".join(text(r.get(c)) for c in new_cols)
        matches = [name for name, pat in P35_TERMS.items() if re.search(pat, hay, re.I)]
        if not matches:
            continue
        nct = upper(r["nct_id"])
        a = audit_map.get(nct, {})
        ruling = text(a.get("p35_ruling") or a.get("final_p35_ruling"))
        override = text(
            a.get("override_status")
            or a.get("p35_override_status")
            or a.get("override")
        )

        classification = "REQUIRES_CENTRAL_CLASSIFICATION"
        rationale = ""

        if ruling == "P35_QUALIFIES_REGISTERED_PATHOGEN_DERIVED_SIGNAL":
            classification = "INCLUDE_P35_IMAGING_SENSITIVITY"
            rationale = (
                "v3.2.9 P35 audit identifies a registered pathogen-derived "
                "imaging signal."
            )
        elif ruling == "P35_QUALIFIES_DIRECT_OBSERVATION_OF_ORGANISM":
            classification = "EXCLUDE_DIRECT_MICROSCOPY_NOT_P35_IMAGING_SENSITIVITY"
            rationale = (
                "Direct microscopy/visualization of the organism is not the "
                "registration-wording imaging boundary targeted by this sensitivity."
            )
        elif ruling == "P35_NOT_APPLICABLE_PATHOGEN_DERIVED_SPECTROSCOPY":
            classification = "EXCLUDE_SPECTROSCOPY_P35_NOT_APPLICABLE"
            rationale = (
                "The v3.2.9 audit classified this as pathogen-derived "
                "spectroscopy rather than P35 imaging."
            )

        rows.append(
            {
                "nct_id": nct,
                "cohort_source": "NEW_V3_2_9",
                "brief_title": text(r.get("brief_title")),
                "final_stratum": text(r.get("final_stratum")),
                "final_amr_depth": text(r.get("final_amr_depth")),
                "historic_primary_modality": "",
                "matched_p35_terms": "|".join(matches),
                "prior_v3_2_9_p35_ruling": ruling,
                "prior_v3_2_9_p35_override_status": override,
                "sensitivity_classification": classification,
                "sensitivity_rationale": rationale,
                "central_adjudicator_initials": "",
                "audit_status": (
                    "PREFILLED_FROM_V3_2_9_P35"
                    if classification != "REQUIRES_CENTRAL_CLASSIFICATION"
                    else "PENDING"
                ),
            }
        )

    out = pd.DataFrame(rows)
    if len(out):
        out = out.sort_values(
            ["cohort_source", "nct_id"], kind="stable"
        ).reset_index(drop=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()

    ap.add_argument("--historic-all-reviewed", required=True, type=Path)
    ap.add_argument("--historic-eligible", required=True, type=Path)

    ap.add_argument("--ib-rescue", required=True, type=Path)
    ap.add_argument("--ib-nonflag", required=True, type=Path)
    ap.add_argument("--ib-registry-negative", required=True, type=Path)
    ap.add_argument("--tb-expansion", required=True, type=Path)

    ap.add_argument("--final-ledger", required=True, type=Path)
    ap.add_argument("--freeze-summary", required=True, type=Path)
    ap.add_argument("--screening-codebook", required=True, type=Path)
    ap.add_argument("--descriptor-codebook", required=True, type=Path)

    ap.add_argument("--registry-parts-dir", required=True, type=Path)

    ap.add_argument("--imaging-audit", type=Path, default=None)

    ap.add_argument("--h2-policy", type=Path, default=None)
    ap.add_argument("--direct-specimen-policy", type=Path, default=None)
    ap.add_argument("--descriptor-lookup-note", type=Path, default=None)

    ap.add_argument("--output-dir", required=True, type=Path)

    args = ap.parse_args()

    for p in [
        args.historic_all_reviewed,
        args.historic_eligible,
        args.ib_rescue,
        args.ib_nonflag,
        args.ib_registry_negative,
        args.tb_expansion,
        args.final_ledger,
        args.freeze_summary,
        args.screening_codebook,
        args.descriptor_codebook,
    ]:
        if not p.exists():
            raise SystemExit(f"Required input missing: {p}")

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise SystemExit(f"Output directory is not empty: {args.output_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    validate_hash(
        args.historic_eligible,
        HISTORIC_EXPECTED_HASH,
        "Historic eligible cohort",
    )
    validate_hash(
        args.final_ledger,
        FINAL_LEDGER_EXPECTED_HASH,
        "v3.2.9 final decision ledger",
    )
    validate_screen_codebook(args.screening_codebook)
    validate_descriptor_codebook(args.descriptor_codebook)

    freeze = json.loads(args.freeze_summary.read_text(encoding="utf-8"))
    if not freeze.get("freeze_criteria_all_pass"):
        raise SystemExit("Freeze summary does not report all criteria PASS.")
    if freeze.get("newly_eligible_for_descriptor_coding", {}).get("n") != 213:
        raise SystemExit("Freeze summary does not report 213 newly eligible.")
    if freeze.get("projected_v3_3_0_cohort_size") not in (None, 573):
        # defensive for alternate nesting
        raise SystemExit("Unexpected projected cohort size in freeze summary.")

    hist_all = read_tsv(args.historic_all_reviewed)
    hist_elig = read_tsv(args.historic_eligible)
    ib_rescue = read_tsv(args.ib_rescue)
    ib_nonflag = read_tsv(args.ib_nonflag)
    ib_rnga = read_tsv(args.ib_registry_negative)
    tb_exp = read_tsv(args.tb_expansion)
    ledger = read_tsv(args.final_ledger)

    if len(hist_all) != EXPECTED["historic_all"]:
        raise SystemExit(f"Historic all: expected 2097, got {len(hist_all)}")
    if len(hist_elig) != EXPECTED["historic_eligible"]:
        raise SystemExit(f"Historic eligible: expected 360, got {len(hist_elig)}")
    if len(ib_rescue) != EXPECTED["ib_rescue"]:
        raise SystemExit(f"IB rescue: expected 1332, got {len(ib_rescue)}")
    if len(ib_nonflag) != EXPECTED["ib_nonflag"]:
        raise SystemExit(f"IB nonflag: expected 400, got {len(ib_nonflag)}")
    if len(ib_rnga) != EXPECTED["ib_registry_negative"]:
        raise SystemExit(f"IB registry-negative: expected 300, got {len(ib_rnga)}")
    if len(tb_exp) != EXPECTED["tb_expansion"]:
        raise SystemExit(f"TB expansion: expected 107, got {len(tb_exp)}")
    if len(ledger) != EXPECTED["ledger"]:
        raise SystemExit(f"Final ledger: expected 291, got {len(ledger)}")

    hist_ids = set(hist_all["nct_id"].map(upper))
    hist_elig_ids = set(hist_elig["nct_id"].map(upper))
    hist_yes_ids = set(
        hist_all.loc[
            hist_all["final_primary_eligible"].map(upper).eq("YES"), "nct_id"
        ].map(upper)
    )
    if hist_elig_ids != hist_yes_ids:
        raise SystemExit("Historic eligible file does not equal YES subset of historic all.")
    if len(hist_ids) != len(hist_all):
        raise SystemExit("Duplicate NCT IDs in historic all-reviewed file.")

    ib_all = pd.concat([ib_rescue, ib_nonflag, ib_rnga], ignore_index=True)
    if len(ib_all) != EXPECTED["ib_total"]:
        raise SystemExit("IB packet total is not 2,032.")
    if set(ib_all["reviewer_initials"]) != {"IB"}:
        raise SystemExit("IB packets contain reviewer initials other than IB.")
    if set(tb_exp["reviewer_initials"]) != {"TB"}:
        raise SystemExit("TB expansion contains reviewer initials other than TB.")

    ib_ids = set(ib_all["nct_id"].map(upper))
    tb_ids = set(tb_exp["nct_id"].map(upper))
    if len(ib_ids) != len(ib_all):
        raise SystemExit("Duplicate NCT ID across IB v3.2.9 packets.")
    if len(tb_ids) != len(tb_exp):
        raise SystemExit("Duplicate NCT ID in TB expansion.")
    if hist_ids & ib_ids:
        raise SystemExit("Historic and IB v3.2.9 universes overlap.")
    if hist_ids & tb_ids:
        raise SystemExit("Historic and TB expansion universes overlap.")
    if ib_ids & tb_ids:
        raise SystemExit("IB and TB expansion universes overlap.")

    ledger_ids = set(ledger["nct_id"].map(upper))
    if not ledger_ids <= (ib_ids | tb_ids):
        raise SystemExit("Final ledger contains NCT outside new review universe.")

    # Build full master.
    master_rows = historic_master_rows(hist_all)

    for df, segment, filename in [
        (
            ib_rescue,
            "V3_2_9_SUPPORT_RESCUE",
            args.ib_rescue.name,
        ),
        (
            ib_nonflag,
            "V3_2_9_SUPPORT_NONFLAG_AUDIT",
            args.ib_nonflag.name,
        ),
        (
            ib_rnga,
            "V3_2_9_REGISTRY_NEGATIVE_AUDIT",
            args.ib_registry_negative.name,
        ),
        (
            tb_exp,
            "V3_2_9_FAILED_STRATA_EXHAUSTIVE_EXPANSION",
            args.tb_expansion.name,
        ),
    ]:
        reviewer = "IB" if "IB" in set(df["reviewer_initials"]) else "TB"
        master_rows.extend(
            new_source_rows(df, segment, filename, reviewer)
        )

    rows_by_nct = {r["nct_id"]: r for r in master_rows}
    if len(rows_by_nct) != len(master_rows):
        raise SystemExit("Full master has duplicate NCT IDs before ledger overrides.")

    apply_ledger_overrides(rows_by_nct, ledger)

    master = pd.DataFrame(list(rows_by_nct.values()))
    master = master.sort_values("nct_id", kind="stable").reset_index(drop=True)

    if len(master) != EXPECTED["full_review_master"]:
        raise SystemExit(f"Full master expected 4236, got {len(master)}")

    counts = Counter(master["final_primary_eligible"])
    expected_counts = {
        "YES": EXPECTED["projected_eligible"],
        "NO": EXPECTED["projected_excluded"],
        "UNCERTAIN": EXPECTED["projected_uncertain"],
    }
    if dict(counts) != expected_counts:
        raise SystemExit(
            f"Full master eligibility mismatch: observed {dict(counts)}, "
            f"expected {expected_counts}"
        )

    master_path = (
        args.output_dir
        / "provenance"
        / "screening_v3_2_9_full_review_master.tsv"
    )
    write_tsv(master, master_path)

    # Scope clarification without changing frozen JSON.
    clarification = f"""\
Screening v3.2.9 provenance clarification
Created: {now_utc()}

The frozen file screening_v3_2_9_final_decision_ledger.tsv contains 291 rows:
282 Packet A confirmation/adjudication records and 9 Packet B records.

Its freeze-summary criterion C8 should therefore be read as:
"The final-decision ledger covers every Packet A and Packet B
confirmation/adjudication NCT exactly once."

It was not intended to enumerate every record ever reviewed.

This companion file:
screening_v3_2_9_full_review_master.tsv

provides one row per unique NCT in the complete reviewed universe:
- historic v3.2.7 reviewed universe: 2,097
- v3.2.9 rescue/audit universe: 2,032
- v3.2.9 failed-strata exhaustive expansion: 107
- total: 4,236

The 291-row ledger is applied as the authoritative v3.2.9 override layer.
This clarification does not reopen or modify any frozen screening decision.
"""
    clarification_path = (
        args.output_dir
        / "provenance"
        / "PROVENANCE_CLARIFICATION_LEDGER_SCOPE_v3_2_9.txt"
    )
    clarification_path.write_text(clarification, encoding="utf-8")

    # Identify new eligible.
    new_eligible = master[
        master["new_eligible_for_descriptor_coding"].eq("YES")
    ].copy()

    if len(new_eligible) != EXPECTED["new_eligible"]:
        raise SystemExit(
            f"Expected 213 new eligible from ledger, got {len(new_eligible)}"
        )
    if set(new_eligible["nct_id"]) & hist_elig_ids:
        raise SystemExit("New eligible set overlaps historic eligible cohort.")

    new_depth = Counter(new_eligible["final_amr_depth"])
    for d, n in EXPECTED_NEW_DEPTH.items():
        if new_depth.get(d, 0) != n:
            raise SystemExit(
                f"New depth {d}: observed {new_depth.get(d,0)}, expected {n}"
            )

    nct_list_path = (
        args.output_dir
        / "descriptor_source"
        / "newly_eligible_nct_ids_v3_3_0.txt"
    )
    nct_list_path.parent.mkdir(parents=True, exist_ok=True)
    nct_list_path.write_text(
        "\n".join(sorted(new_eligible["nct_id"])) + "\n",
        encoding="utf-8",
    )

    # Frozen registry recovery for 213.
    registry = scan_registry_targets(
        args.registry_parts_dir,
        set(new_eligible["nct_id"]),
    )

    descriptor_rows = []
    for _, r in new_eligible.sort_values("nct_id").iterrows():
        nct = r["nct_id"]
        reg = dict(registry[nct])
        reg["final_stratum"] = r["final_stratum"]
        reg["final_amr_depth"] = r["final_amr_depth"]
        reg["screening_final_decision_source"] = r["final_decision_source"]
        reg["screening_final_decision_basis"] = r["final_decision_basis"]
        descriptor_rows.append(reg)

    descriptor_source = pd.DataFrame(descriptor_rows)

    source_cols = [
        "nct_id",
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
        "final_stratum",
        "final_amr_depth",
        "study_type",
        "overall_status",
        "start_date",
        "start_year",
        "has_results",
        "enrollment_count",
        "enrollment_type",
        "lead_sponsor_name",
        "lead_sponsor_class",
        "countries",
        "source_registry_shard",
        "screening_final_decision_source",
        "screening_final_decision_basis",
    ]
    descriptor_source = descriptor_source[source_cols]

    descriptor_source_path = (
        args.output_dir
        / "descriptor_source"
        / "newly_eligible_descriptor_source_enriched_v3_3_0.tsv"
    )
    write_tsv(descriptor_source, descriptor_source_path)

    metadata_summary = {
        "created_at": now_utc(),
        "version": VERSION,
        "newly_eligible_rows": len(descriptor_source),
        "unique_nct_ids": descriptor_source["nct_id"].nunique(),
        "registry_records_recovered": len(registry),
        "start_year_nonblank": int(
            descriptor_source["start_year"].ne("").sum()
        ),
        "study_type_distribution": dict(
            Counter(
                x if x else "MISSING"
                for x in descriptor_source["study_type"]
            )
        ),
        "overall_status_distribution": dict(
            Counter(
                x if x else "MISSING"
                for x in descriptor_source["overall_status"]
            )
        ),
        "has_results_distribution": dict(
            Counter(descriptor_source["has_results"])
        ),
        "new_depth_distribution": dict(
            sorted(Counter(descriptor_source["final_amr_depth"]).items())
        ),
    }
    metadata_summary_path = (
        args.output_dir
        / "descriptor_source"
        / "newly_eligible_registry_recovery_summary_v3_3_0.json"
    )
    metadata_summary_path.write_text(
        json.dumps(metadata_summary, indent=2) + "\n",
        encoding="utf-8",
    )

    # Projected 573 depth check.
    historic_depth = Counter(
        hist_elig["final_amr_depth"].map(upper)
    )
    projected_depth = {
        d: historic_depth.get(d, 0) + new_depth.get(d, 0)
        for d in ["0", "1", "2", "3", "4"]
    }
    if projected_depth != EXPECTED_PROJECTED_DEPTH:
        raise SystemExit(
            f"Projected depth mismatch: {projected_depth} "
            f"!= {EXPECTED_PROJECTED_DEPTH}"
        )

    # Build descriptor packets.
    packet_source = descriptor_source.copy()
    sb_packet, sb_key = make_descriptor_packet(
        packet_source, "SB", SB_SEED
    )
    zb_packet, zb_key = make_descriptor_packet(
        packet_source, "ZB", ZB_SEED
    )

    private_dir = args.output_dir / "private"
    private_dir.mkdir(parents=True, exist_ok=True)
    key = pd.concat([sb_key, zb_key], ignore_index=True)
    key_path = private_dir / "canonical_descriptor_review_private_key_v3_3_0.tsv"
    write_tsv(key, key_path)

    # Carry-forward flags are private QC, not reviewer-facing.
    carry = pd.DataFrame(
        freeze.get("carry_forward_flags_for_descriptor_coding", [])
    )
    if len(carry):
        carry_path = private_dir / "descriptor_carry_forward_flags_v3_3_0.tsv"
        write_tsv(carry, carry_path)
    else:
        carry_path = private_dir / "descriptor_carry_forward_flags_v3_3_0.tsv"
        pd.DataFrame(columns=["nct_id", "flag"]).to_csv(
            carry_path, sep="\t", index=False
        )

    # P35 projected-cohort sensitivity candidate queue.
    imaging_audit = (
        read_tsv(args.imaging_audit)
        if args.imaging_audit is not None and args.imaging_audit.exists()
        else None
    )
    p35_queue = build_p35_candidate_queue(
        hist_elig, descriptor_source, imaging_audit
    )
    p35_path = (
        args.output_dir
        / "sensitivity"
        / "P35_full_cohort_imaging_sensitivity_candidate_queue_v3_3_0.tsv"
    )
    write_tsv(p35_queue, p35_path)

    p35_note = f"""\
P35 full-cohort sensitivity-definition note
Created: {now_utc()}

Purpose
-------
Freeze the exact set of eligible studies to exclude in the P35
pathogen-directed-imaging sensitivity BEFORE running v3.3.1.

Primary cohort
--------------
The primary projected cohort remains unchanged at 573 eligible studies.

Sensitivity target
------------------
The intended sensitivity is narrow:
exclude eligible studies whose PRIMARY qualifying diagnostic is an imaging
approach and whose eligibility depends on the registration explicitly stating
that the imaging signal is bacteria/pathogen-derived.

Do NOT automatically exclude:
- direct conventional microscopy merely because it is optical;
- pathogen-derived spectroscopy that the P35 audit classified as not imaging;
- studies containing imaging where eligibility rests on a separate in-scope
  bacterial assay.

Important historical check
--------------------------
The v3.2.9 P35 audit covered the new rescue/expansion records. The historic
360-study cohort predates that audit. Therefore this preparation generates
P35_full_cohort_imaging_sensitivity_candidate_queue_v3_3_0.tsv so the
sensitivity definition is symmetric across the entire projected 573-study
cohort.

One historic record that clearly warrants inspection is NCT02558062
(Microdosing of BAC ONE to the Distal Lung), whose frozen registration
describes an imaging agent intended to label bacteria in the human lung.

Complete the candidate queue centrally and freeze the final sensitivity NCT
set before inspecting v3.3.1 hypothesis results. This sensitivity classification
does NOT alter primary eligibility, stratum, or depth.
"""
    p35_note_path = (
        args.output_dir
        / "sensitivity"
        / "README_P35_FULL_COHORT_SENSITIVITY_v3_3_0.txt"
    )
    p35_note_path.write_text(p35_note, encoding="utf-8")

    # Reviewer dispatches.
    dispatch_root = args.output_dir / "reviewer_dispatch"
    package_paths = {}

    for reviewer, packet in [("SB", sb_packet), ("ZB", zb_packet)]:
        d = dispatch_root / f"{reviewer}_new_eligible_descriptor_review_v3_3_0"
        d.mkdir(parents=True, exist_ok=True)

        packet_name = (
            f"{reviewer}_new_eligible_canonical_descriptor_review_blinded_"
            "v3_3_0.tsv"
        )
        packet_path = d / packet_name
        write_tsv(packet, packet_path)

        shutil.copy2(
            args.descriptor_codebook,
            d / args.descriptor_codebook.name,
        )
        (d / f"REVIEW_INSTRUCTIONS_{reviewer}_v3_3_0.txt").write_text(
            descriptor_instructions(reviewer, len(packet)),
            encoding="utf-8",
        )

        copy_if_supplied(args.h2_policy, d)
        copy_if_supplied(args.direct_specimen_policy, d)
        copy_if_supplied(args.descriptor_lookup_note, d)

        manifest = d / "SHA256SUMS.txt"
        with manifest.open("w", encoding="utf-8") as f:
            for p in sorted(
                x for x in d.iterdir()
                if x.is_file() and x.name != "SHA256SUMS.txt"
            ):
                f.write(f"{sha256(p)}  {p.name}\n")

        archive = (
            args.output_dir
            / f"{reviewer}_New_Eligible_Descriptor_Reviewer_Package_v3_3_0.tar.gz"
        )
        deterministic_targz(d, archive)
        package_paths[reviewer] = {
            "packet": str(packet_path),
            "dispatch_dir": str(d),
            "archive": str(archive),
            "archive_sha256": sha256(archive),
        }

    # Validate blinded dispatch content.
    for reviewer in ["SB", "ZB"]:
        d = Path(package_paths[reviewer]["dispatch_dir"])
        bad = [
            p.name
            for p in d.iterdir()
            if "private" in p.name.lower()
            or "carry_forward" in p.name.lower()
            or "final_decision_ledger" in p.name.lower()
        ]
        if bad:
            raise SystemExit(
                f"{reviewer} dispatch contains prohibited private files: {bad}"
            )

    # Overall summary.
    summary = {
        "created_at": now_utc(),
        "version": VERSION,
        "screening_version": SCREENING_VERSION,
        "screening_remains_frozen": True,
        "full_review_master": {
            "rows": len(master),
            "eligibility_counts": dict(Counter(master["final_primary_eligible"])),
            "segment_counts": dict(Counter(master["review_universe_segment"])),
            "v3_2_9_ledger_override_rows": int(
                master["in_v3_2_9_final_291_ledger"].eq("YES").sum()
            ),
        },
        "historic_cohort": {
            "eligible_n": len(hist_elig),
            "sha256": sha256(args.historic_eligible),
        },
        "newly_eligible": {
            "n": len(new_eligible),
            "depth_distribution": dict(sorted(new_depth.items())),
            "registry_recovered_n": len(registry),
        },
        "projected_v3_3_0": {
            "eligible_n": EXPECTED["projected_eligible"],
            "depth_distribution": projected_depth,
        },
        "descriptor_review": {
            "SB_rows": len(sb_packet),
            "ZB_rows": len(zb_packet),
            "SB_seed": SB_SEED,
            "ZB_seed": ZB_SEED,
            "review_design": (
                "Independent dual coding of only the 213 newly eligible studies; "
                "historic 360 descriptors remain frozen."
            ),
            "packages": package_paths,
        },
        "p35_sensitivity": {
            "candidate_queue_rows": len(p35_queue),
            "must_be_frozen_before_v3_3_1": True,
            "primary_cohort_unchanged": True,
        },
        "next_gate": (
            "Send the two independent descriptor packets to SB and ZB. "
            "In parallel, centrally finalize the P35 full-cohort imaging "
            "sensitivity candidate queue. After both descriptor returns are "
            "received, reconcile descriptor fields and perform neutral final "
            "descriptor adjudication before building analytic release v3.3.0."
        ),
        "input_sha256": {
            "historic_all_reviewed": sha256(args.historic_all_reviewed),
            "historic_eligible": sha256(args.historic_eligible),
            "ib_rescue": sha256(args.ib_rescue),
            "ib_nonflag": sha256(args.ib_nonflag),
            "ib_registry_negative": sha256(args.ib_registry_negative),
            "tb_expansion": sha256(args.tb_expansion),
            "final_ledger": sha256(args.final_ledger),
            "freeze_summary": sha256(args.freeze_summary),
            "screening_codebook": sha256(args.screening_codebook),
            "descriptor_codebook": sha256(args.descriptor_codebook),
        },
    }

    summary_path = args.output_dir / "PREP_SUMMARY_v3_3_0.json"
    summary_path.write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    readme = f"""\
Post-screening descriptor-stage preparation
Version: {VERSION}

PASS conditions
---------------
- Full reviewed universe = 4,236 unique NCTs.
- Final eligible = 573.
- Historic eligible = 360.
- Newly eligible requiring descriptor coding = 213.
- New depth = 177 d0, 23 d1, 13 d2, 0 d3, 0 d4.
- Projected depth = 435 d0, 83 d1, 52 d2, 3 d3, 0 d4.
- SB and ZB each receive 213 independently shuffled descriptor rows.
- Historic descriptors are not copied into reviewer packets and are not reopened.
- P35 full-cohort sensitivity definition is generated for central finalization
  before v3.3.1 results are inspected.

Do not run H1-H4 yet.
"""
    (args.output_dir / "README_NEXT_STAGE_v3_3_0.txt").write_text(
        readme, encoding="utf-8"
    )

    # Master SHA256 manifest, excluding reviewer archives' internal redundancy only
    # by hashing every file present at completion except the manifest itself.
    manifest_path = args.output_dir / "SHA256SUMS.txt"
    with manifest_path.open("w", encoding="utf-8") as f:
        for p in sorted(
            x for x in args.output_dir.rglob("*")
            if x.is_file() and x.name != "SHA256SUMS.txt"
        ):
            f.write(f"{sha256(p)}  {p.relative_to(args.output_dir)}\n")

    print("V3.3.0 DESCRIPTOR-STAGE PREPARATION: PASS")
    print(json.dumps(summary, indent=2))
    print(f"\nOutput directory: {args.output_dir}")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
