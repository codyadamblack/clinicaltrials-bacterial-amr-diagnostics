#!/usr/bin/env python3
"""Prepare independent confirmation and required support-strata expansion after IB v3.2.9 review.

This script does NOT change the frozen v3.2.7 analytic release or v3.2.8 analysis.
It:
  1) validates the three IB-completed v3.2.9 packets against their blinded originals;
  2) prepares an independent second-review packet for every high-stakes IB row
     (YES, UNCERTAIN, NEEDS_DISCUSSION, or depth >=1);
  3) identifies the source classifier strata of eligible SUPPORT_NONFLAG_AUDIT rows;
  4) exhaustively prepares the remaining nonflag population of each failed source stratum;
  5) writes private keys, summary/provenance, instructions, and SHA-256 checksums.

Expected result for the reviewed data described in the project:
  - 282 rows in independent second-review packet;
  - failed nonflag source strata: CLINICAL_SYNDROMIC_SUPPORT and HOST_RESPONSE_DIAGNOSTIC;
  - 107 remaining nonflag records to review (67 + 40).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

VERSION = "3.2.9-confirmation-expansion"

EVIDENCE_FIELDS = [
    "brief_title",
    "conditions_keywords",
    "intervention_names",
    "primary_outcomes",
    "secondary_outcomes",
    "summary",
    "clinicaltrials_url",
]
REVIEW_FIELDS = [
    "registry_lookup_performed",
    "manual_primary_eligible",
    "manual_final_stratum",
    "manual_amr_depth",
    "manual_exclusion_reason",
    "reviewer_notes",
    "reviewer_initials",
    "review_status",
]
PACKET_FIELDS = ["coverage_audit_id", "nct_id", *EVIDENCE_FIELDS, *REVIEW_FIELDS]
SUPPORT_STRATA = [
    "CLINICAL_SYNDROMIC_SUPPORT",
    "SURVEILLANCE_SUPPORT",
    "HOST_RESPONSE_DIAGNOSTIC",
    "THERAPEUTIC_SUPPORT",
]
RESCUE_FLAGS = [
    "has_diagnostic_test_intervention",
    "structured_diagnostic",
    "diagnostic_evaluation_high",
    "direct_microbial_analyte_high",
    "direct_pathogen_diagnostic_intent",
    "antibiotic_action_high",
]
ALLOWED_ELIG = {"YES", "NO", "UNCERTAIN"}
ALLOWED_STATUS = {"COMPLETE", "NEEDS_DISCUSSION"}
ALLOWED_LOOKUP = {"YES", "NO"}
ALLOWED_DEPTH = {"", "0", "1", "2", "3", "4", "NA"}


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def text(x) -> str:
    return "" if x is None else str(x).strip()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = [{k: ("" if v is None else v) for k, v in row.items()} for row in r]
        return rows, list(r.fieldnames or [])


def write_tsv(path: Path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def key_by(rows, key):
    out = {}
    for r in rows:
        v = text(r.get(key))
        if not v:
            raise SystemExit(f"Blank {key} encountered")
        if v in out:
            raise SystemExit(f"Duplicate {key}: {v}")
        out[v] = r
    return out


def validate_completed(name, completed, original, expected_n, expected_prefix):
    rows, fields = read_tsv(completed)
    orig, orig_fields = read_tsv(original)
    if fields != PACKET_FIELDS:
        raise SystemExit(f"{name}: completed columns differ from expected 17-field schema")
    if orig_fields != PACKET_FIELDS:
        raise SystemExit(f"{name}: original columns differ from expected 17-field schema")
    if len(rows) != expected_n or len(orig) != expected_n:
        raise SystemExit(f"{name}: expected {expected_n} rows; completed={len(rows)}, original={len(orig)}")

    for i, (a, b) in enumerate(zip(rows, orig), 1):
        if text(a["coverage_audit_id"]) != text(b["coverage_audit_id"]):
            raise SystemExit(f"{name}: row-order/id mismatch at row {i}")
        if not text(a["coverage_audit_id"]).startswith(expected_prefix + "-"):
            raise SystemExit(f"{name}: unexpected packet id {a['coverage_audit_id']}")
        if text(a["nct_id"]) != text(b["nct_id"]):
            raise SystemExit(f"{name}: nct mismatch at row {i}")
        for fld in EVIDENCE_FIELDS:
            if a.get(fld, "") != b.get(fld, ""):
                raise SystemExit(f"{name}: evidence field changed: {a['coverage_audit_id']} {fld}")

        elig = text(a.get("manual_primary_eligible")).upper()
        status = text(a.get("review_status")).upper()
        lookup = text(a.get("registry_lookup_performed")).upper()
        depth = text(a.get("manual_amr_depth")).upper()
        initials = text(a.get("reviewer_initials")).upper()
        if elig not in ALLOWED_ELIG:
            raise SystemExit(f"{name}: invalid eligibility at {a['coverage_audit_id']}: {elig}")
        if status not in ALLOWED_STATUS:
            raise SystemExit(f"{name}: invalid status at {a['coverage_audit_id']}: {status}")
        if lookup not in ALLOWED_LOOKUP:
            raise SystemExit(f"{name}: invalid registry lookup at {a['coverage_audit_id']}: {lookup}")
        if depth not in ALLOWED_DEPTH:
            raise SystemExit(f"{name}: invalid depth at {a['coverage_audit_id']}: {depth}")
        if initials != "IB":
            raise SystemExit(f"{name}: reviewer initials not IB at {a['coverage_audit_id']}: {initials}")
        if not text(a.get("manual_final_stratum")):
            raise SystemExit(f"{name}: blank final stratum at {a['coverage_audit_id']}")
        if not text(a.get("reviewer_notes")):
            raise SystemExit(f"{name}: blank reviewer notes at {a['coverage_audit_id']}")
        if elig == "YES" and depth not in {"0", "1", "2", "3", "4"}:
            raise SystemExit(f"{name}: eligible row lacks valid depth at {a['coverage_audit_id']}")
    return rows


def blank_review_row(new_id, source_row):
    out = {"coverage_audit_id": new_id, "nct_id": text(source_row.get("nct_id"))}
    for fld in EVIDENCE_FIELDS:
        out[fld] = source_row.get(fld, "")
    for fld in REVIEW_FIELDS:
        out[fld] = ""
    return out


def is_high_stakes(row):
    elig = text(row.get("manual_primary_eligible")).upper()
    status = text(row.get("review_status")).upper()
    depth = text(row.get("manual_amr_depth")).upper()
    return elig in {"YES", "UNCERTAIN"} or status == "NEEDS_DISCUSSION" or depth in {"1", "2", "3", "4"}


def rescue_flagged(row):
    return any(text(row.get(f)).lower() in {"1", "true", "yes", "y"} for f in RESCUE_FLAGS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ib-rescue", required=True, type=Path)
    ap.add_argument("--ib-nonflag", required=True, type=Path)
    ap.add_argument("--ib-registry-negative", required=True, type=Path)
    ap.add_argument("--original-rescue", required=True, type=Path)
    ap.add_argument("--original-nonflag", required=True, type=Path)
    ap.add_argument("--original-registry-negative", required=True, type=Path)
    ap.add_argument("--rescue-private-key", required=True, type=Path)
    ap.add_argument("--nonflag-private-key", required=True, type=Path)
    ap.add_argument("--registry-private-key", required=True, type=Path)
    ap.add_argument("--retained-all", required=True, type=Path)
    ap.add_argument("--review-provenance", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    for p in [
        args.ib_rescue, args.ib_nonflag, args.ib_registry_negative,
        args.original_rescue, args.original_nonflag, args.original_registry_negative,
        args.rescue_private_key, args.nonflag_private_key, args.registry_private_key,
        args.retained_all, args.review_provenance,
    ]:
        if not p.exists():
            raise SystemExit(f"Missing input: {p}")

    out = args.output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"Output directory is not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)

    resc = validate_completed("SUPPORT_RESCUE", args.ib_rescue, args.original_rescue, 1332, "RESC")
    saud = validate_completed("SUPPORT_NONFLAG_AUDIT", args.ib_nonflag, args.original_nonflag, 400, "SAUD")
    rnga = validate_completed("REGISTRY_NEGATIVE_AUDIT", args.ib_registry_negative, args.original_registry_negative, 300, "RNGA")

    # Validate expected headline counts.
    all_rows = resc + saud + rnga
    elig_counts = Counter(text(r["manual_primary_eligible"]).upper() for r in all_rows)
    status_counts = Counter(text(r["review_status"]).upper() for r in all_rows)
    depth_counts = Counter(text(r["manual_amr_depth"]).upper() for r in all_rows if text(r["manual_amr_depth"]).upper() in {"0","1","2","3","4"})
    if elig_counts != Counter({"NO": 1817, "YES": 214, "UNCERTAIN": 1}):
        raise SystemExit(f"Unexpected combined eligibility counts: {dict(elig_counts)}")
    if status_counts.get("NEEDS_DISCUSSION", 0) != 160:
        raise SystemExit(f"Expected 160 NEEDS_DISCUSSION; found {status_counts.get('NEEDS_DISCUSSION',0)}")
    if depth_counts != Counter({"0": 177, "1": 21, "2": 16}):
        raise SystemExit(f"Unexpected eligible depth counts: {dict(depth_counts)}")

    # Private keys for source classifier stratum and provenance.
    resc_key_rows, _ = read_tsv(args.rescue_private_key)
    saud_key_rows, _ = read_tsv(args.nonflag_private_key)
    rnga_key_rows, _ = read_tsv(args.registry_private_key)
    resc_key = key_by(resc_key_rows, "coverage_audit_id")
    saud_key = key_by(saud_key_rows, "coverage_audit_id")
    rnga_key = key_by(rnga_key_rows, "coverage_audit_id")

    packet_meta = {
        "RESC": ("SUPPORT_RESCUE", resc_key),
        "SAUD": ("SUPPORT_NONFLAG_AUDIT", saud_key),
        "RNGA": ("REGISTRY_NEGATIVE_AUDIT", rnga_key),
    }

    # Independent confirmation packet: union, not additive count.
    high = [r for r in all_rows if is_high_stakes(r)]
    high.sort(key=lambda r: (text(r["nct_id"]), text(r["coverage_audit_id"])))
    if len(high) != 282:
        raise SystemExit(f"Expected 282 high-stakes rows for second review; found {len(high)}")

    second_blinded = []
    second_private = []
    for i, r in enumerate(high, 1):
        sid = f"SR2-{i:04d}"
        original_id = text(r["coverage_audit_id"])
        prefix = original_id.split("-", 1)[0]
        source_name, pkey = packet_meta[prefix]
        meta = pkey.get(original_id, {})
        second_blinded.append(blank_review_row(sid, r))
        second_private.append({
            "second_review_id": sid,
            "nct_id": text(r["nct_id"]),
            "original_coverage_audit_id": original_id,
            "original_packet": source_name,
            "predicted_stratum": text(meta.get("predicted_stratum")),
            "ib_registry_lookup_performed": text(r.get("registry_lookup_performed")),
            "ib_primary_eligible": text(r.get("manual_primary_eligible")),
            "ib_final_stratum": text(r.get("manual_final_stratum")),
            "ib_amr_depth": text(r.get("manual_amr_depth")),
            "ib_exclusion_reason": text(r.get("manual_exclusion_reason")),
            "ib_review_status": text(r.get("review_status")),
            "ib_notes": text(r.get("reviewer_notes")),
            "trigger_yes_or_uncertain": "YES" if text(r.get("manual_primary_eligible")).upper() in {"YES","UNCERTAIN"} else "NO",
            "trigger_needs_discussion": "YES" if text(r.get("review_status")).upper() == "NEEDS_DISCUSSION" else "NO",
            "trigger_depth_ge1": "YES" if text(r.get("manual_amr_depth")) in {"1","2","3","4"} else "NO",
        })

    second_blinded_path = out / "v3_2_9_high_stakes_second_review_blinded.tsv"
    second_private_path = out / "v3_2_9_high_stakes_second_review_private_key.tsv"
    write_tsv(second_blinded_path, second_blinded, PACKET_FIELDS)
    write_tsv(second_private_path, second_private, list(second_private[0].keys()))

    # Identify nonflag audit source strata with IB eligible findings.
    saud_yes = [r for r in saud if text(r["manual_primary_eligible"]).upper() == "YES"]
    if len(saud_yes) != 4:
        raise SystemExit(f"Expected 4 nonflag-audit eligible rows; found {len(saud_yes)}")
    failed_map = []
    failed_strata = set()
    for r in saud_yes:
        m = saud_key[text(r["coverage_audit_id"])]
        ps = text(m.get("predicted_stratum"))
        if ps not in SUPPORT_STRATA:
            raise SystemExit(f"Nonflag positive has unexpected source stratum: {r['coverage_audit_id']} {ps}")
        failed_strata.add(ps)
        failed_map.append({
            "coverage_audit_id": text(r["coverage_audit_id"]),
            "nct_id": text(r["nct_id"]),
            "predicted_source_stratum": ps,
            "ib_final_stratum": text(r["manual_final_stratum"]),
            "ib_amr_depth": text(r["manual_amr_depth"]),
            "ib_review_status": text(r["review_status"]),
            "brief_title": text(r["brief_title"]),
        })

    expected_failed = {"CLINICAL_SYNDROMIC_SUPPORT", "HOST_RESPONSE_DIAGNOSTIC"}
    if failed_strata != expected_failed:
        raise SystemExit(f"Expected failed strata {sorted(expected_failed)}; found {sorted(failed_strata)}")

    failed_map_path = out / "v3_2_9_nonflag_positive_source_strata.tsv"
    write_tsv(failed_map_path, failed_map, list(failed_map[0].keys()))

    # Reconstruct remaining nonflag pool in failed strata.
    retained, retained_fields = read_tsv(args.retained_all)
    req = {"nct_id", "predicted_stratum", *EVIDENCE_FIELDS, *RESCUE_FLAGS}
    missing = req - set(retained_fields)
    if missing:
        raise SystemExit(f"retained_all missing required fields: {sorted(missing)}")
    provenance, _ = read_tsv(args.review_provenance)
    previously_reviewed = {text(r.get("nct_id")) for r in provenance if text(r.get("nct_id"))}
    if len(previously_reviewed) != 2097:
        raise SystemExit(f"Expected 2,097 prior reviewed NCTs; found {len(previously_reviewed)}")

    rescue_ids = {text(r["nct_id"]) for r in resc_key_rows}
    audited_nonflag_ids = {text(r["nct_id"]) for r in saud_key_rows}

    support_unreviewed = [
        r for r in retained
        if text(r.get("predicted_stratum")) in SUPPORT_STRATA
        and text(r.get("nct_id")) not in previously_reviewed
    ]
    nonflag_pool = [r for r in support_unreviewed if text(r["nct_id"]) not in rescue_ids and not rescue_flagged(r)]

    expansion = [
        r for r in nonflag_pool
        if text(r.get("predicted_stratum")) in failed_strata
        and text(r.get("nct_id")) not in audited_nonflag_ids
    ]
    expansion.sort(key=lambda r: (SUPPORT_STRATA.index(text(r["predicted_stratum"])), text(r["nct_id"])))

    exp_counts = Counter(text(r["predicted_stratum"]) for r in expansion)
    expected_exp = Counter({"CLINICAL_SYNDROMIC_SUPPORT": 67, "HOST_RESPONSE_DIAGNOSTIC": 40})
    if exp_counts != expected_exp or len(expansion) != 107:
        raise SystemExit(f"Unexpected nonflag expansion: n={len(expansion)}, counts={dict(exp_counts)}")

    exp_blinded = []
    exp_private = []
    for i, r in enumerate(expansion, 1):
        eid = f"NFEX-{i:04d}"
        exp_blinded.append(blank_review_row(eid, r))
        exp_private.append({
            "expansion_review_id": eid,
            "nct_id": text(r["nct_id"]),
            "predicted_stratum": text(r["predicted_stratum"]),
            "classification_reason": text(r.get("classification_reason")),
            "infection_score": text(r.get("infection_score")),
            "diagnostic_score": text(r.get("diagnostic_score")),
            "amr_score": text(r.get("amr_score")),
            "near_miss_score": text(r.get("near_miss_score")),
            **{f: text(r.get(f)) for f in RESCUE_FLAGS},
        })

    exp_blinded_path = out / "v3_2_9_failed_strata_nonflag_expansion_blinded.tsv"
    exp_private_path = out / "v3_2_9_failed_strata_nonflag_expansion_private_key.tsv"
    write_tsv(exp_blinded_path, exp_blinded, PACKET_FIELDS)
    write_tsv(exp_private_path, exp_private, list(exp_private[0].keys()))

    # Summary.
    packet_counts = {
        "SUPPORT_RESCUE": Counter(text(r["manual_primary_eligible"]).upper() for r in resc),
        "SUPPORT_NONFLAG_AUDIT": Counter(text(r["manual_primary_eligible"]).upper() for r in saud),
        "REGISTRY_NEGATIVE_AUDIT": Counter(text(r["manual_primary_eligible"]).upper() for r in rnga),
    }
    summary = {
        "created_at": now_utc(),
        "version": VERSION,
        "ib_attestation_treatment": "Human-attested independent review per PI instruction",
        "ib_total_rows": len(all_rows),
        "ib_eligibility_counts": dict(elig_counts),
        "ib_needs_discussion_n": status_counts.get("NEEDS_DISCUSSION", 0),
        "ib_depth_counts_nonblank": dict(depth_counts),
        "ib_packet_eligibility_counts": {k: dict(v) for k,v in packet_counts.items()},
        "independent_second_review_union_n": len(high),
        "second_review_selection_rule": "manual_primary_eligible in {YES,UNCERTAIN} OR review_status=NEEDS_DISCUSSION OR manual_amr_depth>=1",
        "nonflag_eligible_n": len(saud_yes),
        "nonflag_eligible_source_strata": sorted(failed_strata),
        "remaining_nonflag_expansion_n": len(expansion),
        "remaining_nonflag_expansion_counts": dict(exp_counts),
        "registry_negative_audit_yes_n": sum(text(r["manual_primary_eligible"]).upper()=="YES" for r in rnga),
        "next_gate": (
            "Do not alter the frozen analytic cohort yet. Complete independent second review of 282 high-stakes rows and exhaustive review of 107 remaining nonflag records in the two failed source strata. "
            "Any expansion YES/UNCERTAIN/NEEDS_DISCUSSION/depth>=1 requires independent confirmation before final screening freeze."
        ),
        "input_sha256": {k: sha256(v) for k,v in {
            "ib_rescue": args.ib_rescue,
            "ib_nonflag": args.ib_nonflag,
            "ib_registry_negative": args.ib_registry_negative,
            "original_rescue": args.original_rescue,
            "original_nonflag": args.original_nonflag,
            "original_registry_negative": args.original_registry_negative,
            "rescue_private_key": args.rescue_private_key,
            "nonflag_private_key": args.nonflag_private_key,
            "registry_private_key": args.registry_private_key,
            "retained_all": args.retained_all,
            "review_provenance": args.review_provenance,
        }.items()},
    }
    summary_path = out / "v3_2_9_confirmation_expansion_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    instructions = f"""FINAL SCREENING CONFIRMATION + EXPANSION v3.2.9\n\nPURPOSE\nThis is the final coverage closure step before rebuilding the analytic cohort. The frozen v3.2.7 release and v3.2.8 statistics remain historical and must not be edited.\n\nPACKET A: v3_2_9_high_stakes_second_review_blinded.tsv\nRows: {len(high)}\nReview every row independently, without seeing IB decisions. Use the binding v3.2.9 eligibility/depth protocol. This packet contains the union of all IB rows coded YES or UNCERTAIN, NEEDS_DISCUSSION, or depth >=1.\n\nPACKET B: v3_2_9_failed_strata_nonflag_expansion_blinded.tsv\nRows: {len(expansion)}\nThis is exhaustive, not sampled. It contains every still-unreviewed nonflag record in the two classifier source strata that produced eligible misses in the 400-row audit:\n- CLINICAL_SYNDROMIC_SUPPORT: {exp_counts['CLINICAL_SYNDROMIC_SUPPORT']}\n- HOST_RESPONSE_DIAGNOSTIC: {exp_counts['HOST_RESPONSE_DIAGNOSTIC']}\nReview every row using the identical eligibility threshold.\n\nREVIEWER FIELDS\nComplete all eight review fields: registry_lookup_performed, manual_primary_eligible, manual_final_stratum, manual_amr_depth, manual_exclusion_reason, reviewer_notes, reviewer_initials, review_status. Eligible rows require depth 0-4. Use NEEDS_DISCUSSION for substantive ambiguity.\n\nESCALATION\n1. Any disagreement with IB in Packet A requires adjudication.\n2. Any YES, UNCERTAIN, NEEDS_DISCUSSION, or depth >=1 in Packet B requires independent confirmation before inclusion/exclusion is frozen.\n3. Registry-negative expansion is NOT required because the new audit found 0/300 eligible records; retain that as a completed sensitivity safeguard.\n4. Do not run final H1-H4 analyses until screening adjudication is frozen and descriptors are completed for every newly eligible study.\n"""
    instructions_path = out / "REVIEW_INSTRUCTIONS_CONFIRMATION_EXPANSION_v3_2_9.txt"
    instructions_path.write_text(instructions, encoding="utf-8")

    # Checksum all outputs except checksum file itself.
    outputs = [
        second_blinded_path, second_private_path, failed_map_path,
        exp_blinded_path, exp_private_path, summary_path, instructions_path,
    ]
    checksum_path = out / "SHA256SUMS.txt"
    with checksum_path.open("w", encoding="utf-8") as f:
        for p in sorted(outputs, key=lambda x: x.name):
            f.write(f"{sha256(p)}  {p.name}\n")

    print("V3.2.9 CONFIRMATION + EXPANSION PREPARATION: PASS")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
