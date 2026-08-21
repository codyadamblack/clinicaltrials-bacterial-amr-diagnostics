#!/usr/bin/env python3
"""
ClinicalTrials.gov bacterial AMR diagnostic landscape screen — version 3.2.5.

Version 3.2.5 is a staged classifier designed after review of the full v2.4
registry output. It deliberately separates:

1. High-value bacterial/infectious-disease relevance
2. Direct pathogen diagnostic intent
3. Host-response diagnostic intent
4. AMR diagnostic depth
5. Mechanism-only, surveillance-only, therapeutic-only, and special-pathogen
   supporting strata
6. Near-miss and random registry-negative audit pools
7. A blinded held-out validation sample that excludes all development controls

Only title, conditions/keywords, intervention names, and registered outcomes can
establish primary bacterial or diagnostic relevance. Summary text can support a
classification but cannot create a primary diagnostic candidate by itself.
Eligibility criteria never establish eligibility.

No third-party packages are required.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import json
import os
import random
import re
import shutil
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VERSION = "3.2.5"
DEFAULT_INPUT = Path("/mnt/d/clinicaltrials/ctgov_api_v2")
DEFAULT_OUTPUT = Path(
    "/mnt/d/clinicaltrials/ctgov_api_v2/projects/"
    "bacterial_amr_diagnostics/screen_v3_2_5"
)

DEVELOPMENT_CONTROL_IDS = {
    "NCT04479189",
    "NCT04283422",
    "NCT03475472",
    "NCT05642767",
    "NCT07469436",
    "NCT00258869",
    "NCT03782454",
    "NCT02176122",
    "NCT03477422",
    "NCT00400946",
    "NCT07141771",
    "NCT02246647",
    "NCT03018925",
    # Version 3.2.5 regression/development controls.
    "NCT00591240",
    "NCT01640886",
    "NCT01981993",
    "NCT03096405",
    "NCT03841162",
    "NCT06765135",
    # Version 3.2.5 regression/development controls.
    "NCT00331019",
    "NCT03199287",
    "NCT02908399",
    "NCT03846921",
    # Version 3.2.5 regression/development controls.
    "NCT05187871",
    "NCT03932942",
    "NCT04190303",
    "NCT05060679",
    "NCT00701948",
    "NCT04323553",
    # Version 3.2.5 host-response regression/development controls.
    "NCT06045416",
    "NCT07211997",
}

HIGH_VALUE_FIELDS = {
    "title",
    "conditions_keywords",
    "intervention_names",
    "primary_outcomes",
    "secondary_outcomes",
}

PRIMARY_RELEVANCE_FIELDS = {
    "title",
    "conditions_keywords",
    "intervention_names",
    "primary_outcomes",
}

MID_VALUE_FIELDS = {
    "intervention_descriptions",
    "arm_descriptions",
    "summary",
}

LOW_VALUE_FIELDS = {"eligibility"}

# Fields searched for supporting evidence after the primary relevance gate.
# Eligibility and arm narrative are deliberately excluded because they caused
# substantial context leakage in v2.4 and add major runtime cost.
EVIDENCE_SCAN_FIELDS = HIGH_VALUE_FIELDS | {
    "summary",
    "intervention_descriptions",
}
MECHANISM_SCAN_FIELDS = EVIDENCE_SCAN_FIELDS

DEPTH_LABELS = {
    0: "organism_identification_only",
    1: "binary_resistance_marker",
    2: "phenotypic_ast_or_mic",
    3: "multimechanism_or_resistance_prediction",
    4: "quantitative_resistance_mechanism",
}

PREDICTED_STRATA = [
    "CORE_AMR_DIAGNOSTIC",
    "BROAD_BACTERIAL_DIAGNOSTIC",
    "HOST_RESPONSE_DIAGNOSTIC",
    "CLINICAL_SYNDROMIC_SUPPORT",
    "MECHANISM_SUPPORT",
    "SPECIAL_PATHOGEN_DIAGNOSTIC",
    "SURVEILLANCE_SUPPORT",
    "THERAPEUTIC_SUPPORT",
]

PRIMARY_DIAGNOSTIC_STRATA = {
    "CORE_AMR_DIAGNOSTIC",
    "BROAD_BACTERIAL_DIAGNOSTIC",
}

VALIDATION_TARGETS = {
    "CORE_AMR_DIAGNOSTIC": 70,
    "BROAD_BACTERIAL_DIAGNOSTIC": 80,
    "HOST_RESPONSE_DIAGNOSTIC": 30,
    "CLINICAL_SYNDROMIC_SUPPORT": 20,
    "MECHANISM_SUPPORT": 40,
    "SPECIAL_PATHOGEN_DIAGNOSTIC": 35,
    "SURVEILLANCE_SUPPORT": 20,
    "THERAPEUTIC_SUPPORT": 20,
    "NEAR_MISS": 70,
    "RANDOM_REGISTRY_NEGATIVE": 50,
}


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def get_nested(obj: dict[str, Any], *keys: str, default: Any = None) -> Any:
    current: Any = obj
    for key in keys:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return current


def as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def clean_text(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        # ClinicalTrials.gov titles occasionally contain Unicode hyphens or
        # non-breaking spaces. Normalize them before regex classification so
        # visually identical phrases such as "point-of-care" behave identically.
        value = value.translate(
            {
                0x00AD: ord("-"),  # soft hyphen
                0x2010: ord("-"),  # hyphen
                0x2011: ord("-"),  # non-breaking hyphen
                0x2012: ord("-"),  # figure dash
                0x2013: ord("-"),  # en dash
                0x2014: ord("-"),  # em dash
                0x2015: ord("-"),  # horizontal bar
                0x2212: ord("-"),  # minus sign
                0x00A0: ord(" "),  # non-breaking space
                0x202F: ord(" "),  # narrow non-breaking space
            }
        )
        return re.sub(r"\s+", " ", value).strip()
    if isinstance(value, (int, float, bool)):
        return str(value)
    if isinstance(value, list):
        return " | ".join(filter(None, (clean_text(v) for v in value)))
    if isinstance(value, dict):
        return " | ".join(filter(None, (clean_text(v) for v in value.values())))
    return str(value)


def join_unique(values: Iterable[Any]) -> str:
    output: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = clean_text(value)
        if text and text not in seen:
            seen.add(text)
            output.append(text)
    return " | ".join(output)


def compile_group(mapping: dict[str, list[str]]) -> dict[str, re.Pattern[str]]:
    return {
        key: re.compile("|".join(patterns), re.IGNORECASE)
        for key, patterns in mapping.items()
    }


def stable_hash_int(value: str) -> int:
    return int(hashlib.sha256(value.encode("utf-8")).hexdigest()[:16], 16)


def stable_fraction(value: str) -> float:
    return stable_hash_int(value) / float(16**16)


def matching_terms(text: str, pattern: re.Pattern[str], limit: int = 20) -> list[str]:
    terms: list[str] = []
    seen: set[str] = set()
    for match in pattern.finditer(text):
        term = match.group(0).strip()
        key = term.casefold()
        if key not in seen:
            seen.add(key)
            terms.append(term)
        if len(terms) >= limit:
            break
    return terms


def snippets_for_pattern(
    text: str,
    pattern: re.Pattern[str],
    *,
    radius: int = 170,
    limit: int = 3,
) -> list[str]:
    output: list[str] = []
    for match in pattern.finditer(text):
        start = max(0, match.start() - radius)
        end = min(len(text), match.end() + radius)
        snippet = re.sub(r"\s+", " ", text[start:end]).strip()
        if snippet not in output:
            output.append(snippet)
        if len(output) >= limit:
            break
    return output


def near_patterns(
    text: str,
    left: re.Pattern[str],
    right: re.Pattern[str],
    *,
    window: int,
) -> bool:
    left_matches = list(left.finditer(text))
    right_matches = list(right.finditer(text))
    for a in left_matches:
        for b in right_matches:
            if a.start() <= b.start():
                distance = b.start() - a.end()
            else:
                distance = a.start() - b.end()
            if distance <= window:
                return True
    return False


def pattern_evidence_by_field(
    fields: dict[str, str],
    pattern: re.Pattern[str],
    *,
    allowed_fields: set[str] | None = None,
) -> tuple[list[str], list[str], list[str]]:
    matched_fields: list[str] = []
    terms: list[str] = []
    snippets: list[str] = []
    for field_name, field_text in fields.items():
        if allowed_fields is not None and field_name not in allowed_fields:
            continue
        if pattern.search(field_text):
            matched_fields.append(field_name)
            terms.extend(matching_terms(field_text, pattern))
            snippets.extend(
                f"{field_name}: {snippet}"
                for snippet in snippets_for_pattern(field_text, pattern)
            )
    return (
        sorted(set(matched_fields)),
        list(dict.fromkeys(terms))[:30],
        list(dict.fromkeys(snippets))[:12],
    )


def atomic_write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temp, path)


def write_delimited(
    path: Path,
    rows: Iterable[dict[str, Any]],
    fieldnames: list[str],
    delimiter: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=delimiter,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temp, path)


def write_csv_and_tsv(
    output_dir: Path,
    stem: str,
    rows: list[dict[str, Any]],
    fieldnames: list[str],
) -> None:
    write_delimited(output_dir / f"{stem}.csv", rows, fieldnames, ",")
    write_delimited(output_dir / f"{stem}.tsv", rows, fieldnames, "\t")


class TopKRows:
    """Keep the highest-scoring rows using bounded memory."""

    def __init__(self, limit: int):
        self.limit = limit
        self.heap: list[tuple[float, int, dict[str, Any]]] = []

    def add(self, score: float, unique_key: str, row: dict[str, Any]) -> None:
        tie = -stable_hash_int(unique_key)
        item = (score, tie, row)
        if len(self.heap) < self.limit:
            heapq.heappush(self.heap, item)
        elif item[:2] > self.heap[0][:2]:
            heapq.heapreplace(self.heap, item)

    def rows(self) -> list[dict[str, Any]]:
        return [
            item[2]
            for item in sorted(
                self.heap,
                key=lambda item: (-item[0], -item[1]),
            )
        ]


class LowestHashRows:
    """Keep a deterministic registry-wide random-like sample."""

    def __init__(self, limit: int):
        self.limit = limit
        self.heap: list[tuple[int, str, dict[str, Any]]] = []

    def add(self, unique_key: str, row: dict[str, Any]) -> None:
        # Negative value creates a max-heap behavior on Python's min-heap.
        value = -stable_hash_int(unique_key)
        item = (value, unique_key, row)
        if len(self.heap) < self.limit:
            heapq.heappush(self.heap, item)
        elif item[:2] > self.heap[0][:2]:
            heapq.heapreplace(self.heap, item)

    def rows(self) -> list[dict[str, Any]]:
        return [
            item[2]
            for item in sorted(
                self.heap,
                key=lambda item: (-item[0], item[1]),
            )
        ]


# ---------------------------------------------------------------------------
# Curated terminology
# ---------------------------------------------------------------------------

ORGANISM_PATTERNS = compile_group(
    {
        "enterobacterales": [
            r"\benterobacterales\b",
            r"\benterobacteriaceae\b",
            r"\bescherichia coli\b",
            r"\be\.?\s*coli\b",
            r"\bklebsiella pneumoniae\b",
            r"\bk\.?\s*pneumoniae\b",
            r"\bklebsiella\b",
            r"\benterobacter(?: cloacae)?\b",
            r"\bcitrobacter\b",
            r"\bserratia\b",
            r"\bproteus(?: mirabilis)?\b",
            r"\bmorganella\b",
            r"\bprovidencia\b",
            r"\bsalmonella\b",
            r"\bshigella\b",
            r"\byersinia\b",
        ],
        "s_aureus": [
            r"\bstaphylococcus aureus\b",
            r"\bs\.?\s*aureus\b",
            r"\bmrsa\b",
            r"\bmssa\b",
        ],
        "other_gram_positive": [
            r"\benterococcus(?: faecalis| faecium)?\b",
            r"\bvre\b",
            r"\bstreptococcus(?: pneumoniae| pyogenes| agalactiae)?\b",
            r"\bgroup [ab] strep(?:tococcus|tococci)?\b",
            r"\bgroup b streptococc(?:us|i)?\b",
            r"\bpneumococcus\b",
            r"\bcoagulase[- ]negative staphylococc",
            r"\bstaphylococcus epidermidis\b",
            r"\blisteria monocytogenes\b",
            r"\bgram[- ]positive bacter",
        ],
        "nonfermenter": [
            r"\bpseudomonas aeruginosa\b",
            r"\bp\.?\s*aeruginosa\b",
            r"\bacinetobacter baumannii\b",
            r"\ba\.?\s*baumannii\b",
            r"\bstenotrophomonas maltophilia\b",
            r"\bburkholderia cepacia\b",
            r"\bnon[- ]ferment",
        ],
        "other_bacterial": [
            r"\bhaemophilus influenzae\b",
            r"\bmoraxella catarrhalis\b",
            r"\blegionella pneumophila\b",
            r"\bcampylobacter\b",
            r"\bvibrio\b",
            r"\bbacteroides\b",
            r"\bneisseria meningitidis\b",
            r"\bmeningococcus\b",
            r"\bborrelia\b",
            r"\brickettsia\b",
            r"\bbartonella\b",
            r"\bbrucella\b",
            r"\bclostridium perfringens\b",
            r"\banaerobic bacter",
        ],
        "mixed_pan_bacterial": [
            r"\bbacterial infection",
            r"\bbacterial pathogen",
            r"\bbacterial isolate",
            r"\bclinical isolate",
            r"\bgram[- ]negative bacter",
            r"\bgram[- ]positive bacter",
            r"\beskape\b",
            r"\bpan[- ]bacterial\b",
            r"\bpolymicrobial infection",
        ],
    }
)

SPECIAL_PATHOGEN_PATTERNS = compile_group(
    {
        "mycobacteria": [
            r"\bmycobacter(?:ium|ia|ial)\b",
            r"\btuberculosis\b",
            r"\bmtb\b",
            r"\bnon[- ]tuberculous mycobacter",
        ],
        "h_pylori": [r"\bhelicobacter pylori\b", r"\bh\.?\s*pylori\b"],
        "c_difficile": [
            r"\bclostridioides difficile\b",
            r"\bclostridium difficile\b",
            r"\bc\.?\s*difficile\b",
        ],
        "bacterial_sti": [
            r"\bneisseria gonorrhoeae\b",
            r"\bgonorrh(?:ea|oea)\b",
            r"\bchlamydia trachomatis\b",
            r"\bsyphilis\b",
            r"\btreponema pallidum\b",
        ],
    }
)

BACTERIAL_SYNDROME_PATTERNS = compile_group(
    {
        "bloodstream_sepsis": [
            r"\bbloodstream infection",
            r"\bbacteremia\b",
            r"\bbacteraemia\b",
            r"\bpositive blood culture",
            r"\bsepticemia\b",
            r"\bsepsis\b",
        ],
        "urinary": [
            r"\burinary tract infection",
            r"\bcomplicated uti\b",
            r"\bpyelonephritis\b",
            r"\burosepsis\b",
        ],
        "respiratory": [
            r"\bbacterial pneumonia\b",
            r"\bcommunity[- ]acquired pneumonia\b",
            r"\bhospital[- ]acquired pneumonia\b",
            r"\bventilator[- ]associated pneumonia\b",
            r"\blower respiratory tract infection\b",
            r"\bupper respiratory tract infection\b",
            r"\brespiratory tract infection(?:s)?\b",
            r"\bacute respiratory infection(?:s)?\b",
        ],
        "skin_soft_tissue": [
            r"\bskin and soft tissue infection",
            r"\bssti\b",
            r"\bwound infection",
            r"\bsurgical site infection",
            r"\bcellulitis\b",
            r"\babscess\b",
        ],
        "bone_joint": [
            r"\bosteomyelitis\b",
            r"\bprosthetic joint infection",
            r"\bbone and joint infection",
        ],
        "endocarditis": [r"\binfective endocarditis\b", r"\bendocarditis\b"],
        "intra_abdominal": [
            r"\bintra[- ]abdominal infection",
            r"\bpostoperative (?:intra[- ]?)?abdominal infection",
            r"\babdominal infection",
            r"\bperitonitis\b",
            r"\bcholangitis\b",
        ],
        "cns": [
            r"\bbacterial meningitis\b",
            r"\bmeningitis\b",
            r"\bencephalitis\b",
            r"\bbrain infection(?:s)?\b",
            r"\bbrain abscess\b",
            r"\babscess brain\b",
            r"\bventriculitis\b",
        ],
        "colonization_surveillance": [
            r"\bcolonization\b",
            r"\bcolonisation\b",
            r"\bcarriage\b",
            r"\brectal swab\b",
            r"\bsurveillance culture\b",
        ],
    }
)

GENERAL_INFECTION_PATTERN = re.compile(
    r"\b(infectious disease|infection(?:s)?|infection diagnosis|"
    r"infection diagnostic|suspected infection|confirmed infection|"
    r"serious infection|acute infection|microbiolog(?:y|ical)|"
    r"positive culture|blood culture|clinical microbiology)\b",
    re.IGNORECASE,
)

NONBACTERIAL_PATHOGEN_PATTERN = re.compile(
    r"\b(viral|virus|influenza|respiratory syncytial virus|\brsv\b|"
    r"sars[- ]?cov[- ]?2|covid|hiv|hepatitis|fungal|fungus|candida|"
    r"aspergillus|parasit(?:e|ic)|malaria|plasmodium)\b",
    re.IGNORECASE,
)

DIRECT_DIAGNOSTIC_PATTERNS = compile_group(
    {
        "molecular_pcr": [
            r"\bpcr\b",
            r"\bpolymerase chain reaction\b",
            r"\bqpcr\b",
            r"\brt[- ]?pcr\b",
            r"\bdigital pcr\b",
            r"\bddpcr\b",
            r"\blamp(?: assay)?\b",
            r"\bnucleic acid amplification",
            r"\bnaat\b",
            r"\bmultiplex (?:pcr|panel|assay)",
            r"\bmolecular diagnostic",
            r"\bmolecular detection",
            r"\bmolecular identification",
            r"\bmolecular characteri[sz]ation",
            r"\bfilmarray\b",
            r"\bbiofire\b",
            r"\bverigene\b",
            r"\bgenexpert\b",
            r"\bxpert(?: mtb/rif| xdr)?\b",
            r"\bqiastat\b",
            r"\bt2mr\b",
            r"\bbcid\b",
        ],
        "sequencing": [
            r"\bwhole[- ]genome sequencing\b",
            r"\bwgs\b",
            r"\bnext[- ]generation sequencing\b",
            r"\bmetagenomic sequencing\b",
            r"\bshotgun metagenom",
            r"\bnanopore\b",
            r"\bminion\b",
            r"\b16s rrna sequencing\b",
        ],
        "phenotypic_ast": [
            r"\bantimicrobial susceptibility test(?:ing)?\b",
            r"\bantibiotic susceptibility test(?:ing)?\b",
            r"\brapid ast\b",
            r"\bphenotypic ast\b",
            r"\bphenotypic susceptibility\b",
            r"\bgrowth[- ]based susceptibility\b",
            r"\bminimum inhibitory concentration\b",
            r"\bmic testing\b",
            r"\bbroth microdilution\b",
            r"\bdisk diffusion\b",
            r"\bdisc diffusion\b",
            r"\baccelerate pheno\b",
            r"\bresistell\b",
        ],
        "mass_spectrometry": [
            r"\bmaldi[- ]?tof\b",
            r"\bmass spectrometr(?:y|ic)\b",
        ],
        "culture_microscopy": [
            r"\bbacterial culture\b",
            r"\bblood culture identification\b",
            r"\brapid culture\b",
            r"\bchromogenic (?:agar|medium)\b",
            r"\bmicroscopy[- ]based\b",
            r"\bgram stain\b",
        ],
        "antigen_immunoassay": [
            r"\bantigen test(?:ing)?\b",
            r"\burinary antigen\b",
            r"\blateral flow\b",
            r"\benzyme immunoassay\b",
            r"\bantigen immunoassay\b",
            r"\belisa\b",
            r"\bimmunochromatograph",
        ],
        "biosensor_other": [
            r"\bbiosensor\b",
            r"\bmicrofluidic\b",
            r"\bpoint[-\u2010-\u2015\u2212 ]of[-\u2010-\u2015\u2212 ]care (?:test(?:ing)?|assay|diagnostic(?: testing)?)",
            r"\bdiagnostic platform\b",
            r"\bautomated identification system\b",
        ],
        "bacterial_dna_quantification": [
            r"\bquantification of (?:bacterial|mycobacterial|mtb) dna\b",
            r"\b(?:bacterial|mycobacterial|mtb) dna quantification\b",
            r"\b(?:microbial|bacterial|mycobacterial|mtb) cell[- ]free dna\b",
            r"\bbacterial load(?: measurement| quantification)?\b",
            r"\bquantitative bacterial dna\b",
        ],
    }
)

# Host-derived immune analytes are distinct from direct microbial targets.
# These patterns capture patient antibody/cell/protein responses used to infer
# infection etiology. They intentionally avoid generic "antibody" wording so
# antibody-based capture of a bacterial antigen remains a direct assay.
HOST_IMMUNE_ANALYTE_PATTERN = re.compile(
    r"\b(?:"
    r"antibody[- ]secreting cell(?:s)?|(?:asc|plasmablast)\s+elispot|"
    r"b[- ]cell(?:ular)?\s+(?:response|activation|assay|diagnostic(?:s)?)|"
    r"t[- ]cell(?:ular)?\s+(?:response|activation|assay|diagnostic(?:s)?)|"
    r"host immune response|host[- ]response biomarker(?:s)?|"
    r"serologic(?:al)?\s+(?:response|assay|testing|diagnos(?:is|tic))|"
    r"antibody response|pathogen[- ]specific antibodies|"
    r"borrelia[- ]specific asc(?:s)?|borrelia[- ]specific antibody[- ]secreting cell(?:s)?|"
    r"febridx|myxovirus resistance protein a|\bmxa\b"
    r")\b",
    re.IGNORECASE,
)

HOST_IMMUNE_DIAGNOSTIC_PATTERN = re.compile(
    r"\b(?:"
    r"febridx|"
    r"(?:host immune response|host[- ]response biomarker(?:s)?|"
    r"antibody[- ]secreting cell(?:s)?|(?:asc|plasmablast)\s+elispot|"
    r"b[- ]cell(?:ular)?\s+(?:response|assay|diagnostic(?:s)?)|"
    r"t[- ]cell(?:ular)?\s+(?:response|assay|diagnostic(?:s)?)|"
    r"serologic(?:al)?\s+(?:response|assay|testing)|antibody response|"
    r"myxovirus resistance protein a|\bmxa\b)"
    r".{0,100}\b(?:diagnos(?:is|tic)|test(?:ing)?|assay|validation|"
    r"sensitivity|specificity|differentiat(?:e|ion|ing)|distinguish)"
    r"|(?:diagnos(?:is|tic)|test(?:ing)?|assay|validation|sensitivity|"
    r"specificity|differentiat(?:e|ion|ing)|distinguish)"
    r".{0,100}\b(?:host immune response|host[- ]response biomarker(?:s)?|"
    r"antibody[- ]secreting cell(?:s)?|(?:asc|plasmablast)\s+elispot|"
    r"b[- ]cell(?:ular)?\s+(?:response|assay|diagnostic(?:s)?)|"
    r"t[- ]cell(?:ular)?\s+(?:response|assay|diagnostic(?:s)?)|"
    r"serologic(?:al)?\s+(?:response|assay|testing)|antibody response|"
    r"myxovirus resistance protein a|\bmxa\b)"
    r")",
    re.IGNORECASE,
)

# Positive-control pattern for direct pathogen-derived analytes. This flag is
# used for auditability and self-tests; host-response routing is based on the
# measured patient response, not merely the presence of an immunoassay.
DIRECT_MICROBIAL_ANALYTE_PATTERN = re.compile(
    r"\b(?:bacterial|microbial|pathogen|pneumococcal|streptococcal|"
    r"staphylococcal)\s+(?:antigen|toxin|capsular antigen|cell[- ]wall antigen|"
    r"protein|dna|rna|nucleic acid)\b|"
    r"\b(?:bacterial|microbial|pathogen)\s+nucleic[- ]acid\b",
    re.IGNORECASE,
)

# FebriDx and the paired MxA/CRP signature are host-response tests designed to
# distinguish bacterial from viral infection; flag them as mixed even when the
# registration omits the literal phrase "bacterial versus viral."
FEBRIDX_PATTERN = re.compile(r"\bfebridx(?:®)?\b", re.IGNORECASE)
MXA_PATTERN = re.compile(
    r"\b(?:myxovirus resistance protein a|mx[- ]?a)\b", re.IGNORECASE
)
CRP_PATTERN = re.compile(
    r"\b(?:c[- ]reactive protein|crp)\b", re.IGNORECASE
)

HOST_RESPONSE_PATTERN = re.compile(
    r"\b(procalcitonin|c[- ]reactive protein|\bcrp\b|host[- ]response|"
    r"host gene expression|host transcriptom|immune signature|"
    r"cytokine signature|inflammatory marker(?:s)?|trail[- /]?ip[- ]?10|memed bv|"
    r"alpha[- ]defensin|volatile organic compounds?|breath biomarker|"
    r"breathomics|metabolomic signature|proteomic signature|"
    r"bacterial versus viral|bacterial vs viral|presepsin|"
    r"febridx|myxovirus resistance protein a|mx[- ]?a|"
    r"antibody[- ]secreting cell(?:s)?|(?:asc|plasmablast)\s+elispot|"
    r"b[- ]cell(?:ular)?\s+(?:response|assay|diagnostic(?:s)?)|"
    r"t[- ]cell(?:ular)?\s+(?:response|assay|diagnostic(?:s)?)|"
    r"serologic(?:al)?\s+(?:response|assay|testing)|antibody response|"
    r"soluble triggering receptor expressed on myeloid cells|supar)\b",
    re.IGNORECASE,
)

HOST_DIAGNOSTIC_PATTERN = re.compile(
    r"\b(bacterial versus viral|bacterial vs viral|"
    r"distinguish(?:ing)? between "
    r"(?:(?:acute|suspected|confirmed|community[- ]acquired|systemic)\s+){0,3}"
    r"bacterial and viral(?: infections?)?|"
    r"differentiat(?:e|ing|ion) .* bacterial|"
    r"diagnos(?:is|tic) of "
    r"(?:(?:acute|suspected|confirmed|nosocomial|community[- ]acquired|"
    r"hospital[- ]acquired|neonatal|late[- ]onset)\s+){0,3}"
    r"(?:bacterial )?infection(?:s)?|diagnostic accuracy|"
    r"identify bacterial infection|detect bacterial infection|"
    r"etiolog(?:y|ic) diagnosis|memed bv|triverity|febridx|"
    r"b[- ]cell diagnostics?|t[- ]cell diagnostics?|"
    r"serologic(?:al)? diagnos(?:is|tic)|host immune response)\b",
    re.IGNORECASE,
)

PROGNOSTIC_PATTERN = re.compile(
    r"\b(prognos(?:is|tic)|predict(?:ion|or|ive) of (?:mortality|outcome|"
    r"severity|organ failure)|risk stratification|disease severity|"
    r"mortality prediction)\b",
    re.IGNORECASE,
)

DIAGNOSTIC_DEVELOPMENT_PATTERN = re.compile(
    r"\b(develop(?:ment|ing)? (?:of )?(?:a |new |novel )?diagnostic tests?|"
    r"diagnostic test development|infection diagnostic toolkit|"
    r"diagnostic stewardship|rapid diagnostic technolog(?:y|ies)|"
    r"pathogen identification assay|bacterial diagnostic assay|"
    r"molecular diagnosis of (?:bacterial )?infection)\b",
    re.IGNORECASE,
)

DIAGNOSTIC_EVALUATION_PATTERN = re.compile(
    r"\b(diagnostic accuracy|diagnostic performance|assay performance|"
    r"test performance|analytical validation|clinical validation|"
    r"diagnostic sensitivity|diagnostic specificity|"
    r"sensitivity and specificity|positive predictive value|"
    r"negative predictive value|\bppv\b|\bnpv\b|"
    r"receiver operating characteristic|area under the receiver operating|"
    r"\bauroc\b|agreement with|concordance|"
    r"limit of detection|analytical sensitivity|reference standard|"
    r"gold standard|turnaround time|time to result|time to identification|"
    r"rapid identification|rapid detection|clinical utility|"
    r"diagnostic yield|diagnostic test development|develop diagnostic|"
    r"validation of (?:the )?(?:assay|test|panel|platform)|"
    r"compare(?:d|s|ing)? (?:the )?(?:assay|test|panel|platform))\b",
    re.IGNORECASE,
)

DIAGNOSTIC_ROLE_PATTERN = re.compile(
    r"\b(diagnos(?:is|tic)|detect(?:ion|ing)|identif(?:ication|y|ying)|"
    r"characteri[sz]ation of bacteria|"
    r"quantification of (?:bacterial|mycobacterial|mtb) dna|"
    r"pathogen panel|molecular panel|assay|reference method)\b",
    re.IGNORECASE,
)

MICROBIAL_ASSAY_FOCUS_PATTERN = re.compile(
    r"\b(assay for pathogen identification|pathogen identification assay|"
    r"rapid (?:test|assay) (?:to |for )?(?:detect|identify) (?:a )?"
    r"(?:bacterium|bacteria|bacterial pathogen|pathogen)|"
    r"bacterial detection assay|bacterial identification test|"
    r"microbial identification test|pathogen profile assay|"
    r"blood culture identification panel)\b",
    re.IGNORECASE,
)

# Version 3.2: unmistakable named bacterial/AMR assay language in a title.
# This pathway recovers registrations such as “MRSA/MSSA Blood Culture Test”
# and “MRSA-PCR” even when older registry records lack detailed structured
# diagnostic metadata.
NAMED_BACTERIAL_ASSAY_TITLE_PATTERN = re.compile(
    r"\b(?:"
    r"(?:keypath\s+)?mrsa\s*(?:/|and|-)\s*mssa(?:\s+blood\s+culture)?"
    r"\s+(?:test|assay|panel)|"
    r"(?:keypath\s+)?mssa\s*(?:/|and|-)\s*mrsa(?:\s+blood\s+culture)?"
    r"\s+(?:test|assay|panel)|"
    r"mrsa\s*[-/]?\s*pcr|mssa\s*[-/]?\s*pcr|"
    r"(?:carbapenemase|esbl|cre|cpe|methicillin resistance|"
    r"vancomycin resistance|resistance gene)\s+"
    r"(?:test|assay|pcr|panel)|"
    r"(?:blood culture|bacterial|microbial|pathogen)\s+"
    r"(?:identification|detection)\s+(?:test|assay|panel)|"
    r"rapid\s+(?:pathogen|bacterial|microbial)\s+"
    r"(?:identification|detection|test|assay)|"
    r"(?:antimicrobial|antibiotic)\s+susceptibility\s+"
    r"(?:test|testing|assay)|"
    r"rapid\s+ast|phenotypic\s+ast"
    r")\b",
    re.IGNORECASE,
)


# Version 3.2: a title can establish direct bacterial diagnostic intent when it
# explicitly describes detection, identification, screening, or diagnosis of a
# curated named bacterial organism. The organism match is evaluated separately
# from this wording pattern to prevent generic "rapid detection" studies from
# entering the primary cohort.
NAMED_ORGANISM_DETECTION_WORDING_PATTERN = re.compile(
    r"\b(?:rapid\s+(?:detection|identification|screening|diagnosis)|"
    r"(?:detection|identification|screening|diagnosis)\s+of|"
    r"(?:test|assay|screen)\s+for)\b",
    re.IGNORECASE,
)

# Clinical metagenomics in a bacterial infectious syndrome is inherently a
# direct pathogen-diagnostic strategy even when the registry title omits the
# words diagnostic, detection, or identification. Transmission genomics and
# generic microbiome studies do not match this title pattern.
CLINICAL_METAGENOMICS_TITLE_PATTERN = re.compile(
    r"\b(?:clinical\s+metagenomics?|"
    r"metagenomic(?:\s+sequencing)?\s+(?:diagnos(?:is|tic)|testing)|"
    r"diagnostic\s+metagenomics?)\b",
    re.IGNORECASE,
)


# Version 3.2.5: rapid or molecular diagnosis of a defined infectious syndrome
# can establish direct microbial diagnostic intent when a direct laboratory
# modality is evaluated. The modality requirement prevents clinical scores,
# imaging, and host biomarkers from entering the primary cohort.
SYNDROME_DIRECT_DIAGNOSTIC_TITLE_PATTERN = re.compile(
    r"\b(?:rapid|molecular|microbiologic(?:al)?|metagenomic|etiologic(?:al)?)\s+"
    r"diagnos(?:is|tic)(?:\s+(?:test|testing|study|strategy))?\s+(?:of|for)\s+"
    r"(?:suspected\s+|postoperative\s+|acute\s+|complicated\s+)?"
    r"(?:abdominal|intra[- ]abdominal|bloodstream|urinary tract|respiratory tract|"
    r"lung|wound|surgical site|bone and joint|prosthetic joint|brain|"
    r"central nervous system)?\s*infection(?:s)?\b",
    re.IGNORECASE,
)

# Point-of-care pathogen testing is an explicit direct diagnostic intervention
# even when an older registration omits words such as accuracy or detection.
POINT_OF_CARE_PATHOGEN_TESTING_TITLE_PATTERN = re.compile(
    r"\b(?:point[-\u2010-\u2015\u2212 ]of[-\u2010-\u2015\u2212 ]care|bedside)\s+"
    r"(?:test(?:ing)?|assay|diagnostic(?: testing)?)\s+"
    r"(?:of|for)?\s*(?:respiratory\s+|bacterial\s+|microbial\s+)?"
    r"pathogen(?:s)?\b|"
    r"\b(?:respiratory\s+|bacterial\s+|microbial\s+)?pathogen(?:s)?\s+"
    r"(?:point[-\u2010-\u2015\u2212 ]of[-\u2010-\u2015\u2212 ]care|bedside)\s+(?:test(?:ing)?|assay)\b",
    re.IGNORECASE,
)

# Implementation packages that improve syndromic or microbiological diagnosis
# without evaluating a defined direct pathogen assay remain clinical/syndromic
# support rather than broad bacterial diagnostics.
SYNDROMIC_DIAGNOSTIC_IMPLEMENTATION_TITLE_PATTERN = re.compile(
    r"\b(?:improv(?:e|ing)|implementation|intervention|care pathway|"
    r"diagnostic stewardship).*\bdiagnos(?:is|tic)\b.*\b"
    r"(?:infection(?:s)?|meningitis|encephalitis|sepsis|pneumonia)\b|"
    r"\bdiagnos(?:is|tic)\b.*\b(?:management|care pathway|implementation)\b.*"
    r"\b(?:infection(?:s)?|meningitis|encephalitis|sepsis|pneumonia)\b",
    re.IGNORECASE,
)

# A direct primary diagnostic must target a microbe, microbial analyte, culture,
# resistance determinant, or susceptibility phenotype. A biomarker measured in
# patients with sepsis is not enough by itself. Conditions/keywords alone do
# not establish this target.
MICROBIAL_TARGET_PATTERN = re.compile(
    r"\b(?:"
    r"pathogen(?:s)?|bacterium|bacteria|bacterial|microbial|microbiologic(?:al)?|"
    r"organism identification|blood culture|positive culture|clinical isolate|"
    r"bacterial dna|microbial dna|bacterial rna|microbial rna|16s rrna|"
    r"bacterial antigen|microbial antigen|mrsa|mssa|cre|cpe|esbl|vre|"
    r"carbapenemase|beta[- ]lactamase|β[- ]lactamase|resistance gene|"
    r"resistance marker|antimicrobial susceptibility|antibiotic susceptibility|"
    r"susceptibility testing|minimum inhibitory concentration|mic testing|"
    r"rapid ast|phenotypic ast|culture identification"
    r")\b",
    re.IGNORECASE,
)

# Organ-injury and prognostic targets are retained as clinical/syndromic
# support unless a separate direct microbial target is actually evaluated.
ORGAN_INJURY_PROGNOSTIC_TARGET_PATTERN = re.compile(
    r"\b(?:acute kidney injury|\baki\b|renal injury|kidney injury|"
    r"organ failure|multiple organ dysfunction|cardiac injury|myocardial injury|"
    r"hepatic injury|liver injury|acute lung injury|mortality prediction|"
    r"severity prediction|prognostic biomarker|risk stratification|"
    r"disease severity|predict(?:ion|ive) of mortality)\b",
    re.IGNORECASE,
)

# Mixed respiratory/syndromic panels remain broad direct diagnostics but are
# flagged so bacterial-only, mixed viral/bacterial, and host-response studies
# can be analyzed separately.
MIXED_VIRAL_BACTERIAL_PANEL_PATTERN = re.compile(
    r"\b(?:bacterial and viral|viral and bacterial|bacterial versus viral|"
    r"bacterial vs viral|respiratory pathogen panel|multiplex respiratory panel|"
    r"syndromic respiratory panel|biofire spotfire|spotfire|"
    r"filmarray respiratory|respiratory pathogen testing|"
    r"testing of respiratory pathogens?)\b",
    re.IGNORECASE,
)

REFERENCE_STANDARD_PATTERN = re.compile(
    r"\b(culture as (?:a |the )?reference|compared with culture|"
    r"compared to culture|reference culture|standard microbiological|"
    r"routine microbiology|conventional microbiology|"
    r"phenotypic susceptibility as (?:a |the )?reference|"
    r"broth microdilution as (?:a |the )?reference)\b",
    re.IGNORECASE,
)

ANTIBIOTIC_ACTION_PATTERN = re.compile(
    r"\b(time to active therap|time to effective therap|"
    r"time to optimal therap|time to targeted therap|"
    r"appropriate antibiotic|appropriate antimicrobial|"
    r"antibiotic selection|antimicrobial selection|"
    r"antibiotic exposure|antimicrobial exposure|"
    r"antibiotic duration|days of therapy|de[- ]escalation|"
    r"antimicrobial stewardship|antibiotic stewardship|"
    r"change in antibiotic|guide antibiotic|targeted treatment|"
    r"definitive therapy|stop empiric vancomycin|discontinue vancomycin)\b",
    re.IGNORECASE,
)

SURVEILLANCE_PATTERN = re.compile(
    r"\b(prevalence|incidence|epidemiology|molecular epidemiology|"
    r"surveillance|colonization|colonisation|carriage|screening program|"
    r"outbreak investigation|resistance surveillance|antibiogram|"
    r"transmission(?: rate| analysis| dynamics)?|genomic epidemiology|"
    r"clonal relatedness|phylogenetic transmission)\b",
    re.IGNORECASE,
)

THERAPEUTIC_PATTERN = re.compile(
    r"\b(randomi[sz]ed treatment|antibiotic treatment|antimicrobial treatment|"
    r"compare(?:d|s|ing)? .* versus .* antibiotic|definitive treatment|"
    r"efficacy of .* antibiotic|noninferiority.*antibiotic|"
    r"therapy for .* infection|treatment of .* infection)\b",
    re.IGNORECASE,
)

GENERAL_AMR_PATTERN = re.compile(
    r"\b(antimicrobial resistance|antibiotic resistance|"
    r"drug[- ]resistant bacter|multidrug[- ]resistant|"
    r"multi[- ]drug resistant|extensively drug[- ]resistant|"
    r"resistant isolate|resistance mechanism|resistance determinant|"
    r"resistance gene|non[- ]susceptib|nonsusceptib|\bmdro\b|\bxdr\b)\b",
    re.IGNORECASE,
)

AST_PATTERN = re.compile(
    r"\b(antimicrobial susceptibility test(?:ing)?|"
    r"antibiotic susceptibility test(?:ing)?|susceptibility testing|"
    r"rapid ast|phenotypic ast|minimum inhibitory concentration|"
    r"\bmic testing\b|broth microdilution|disk diffusion|disc diffusion|"
    r"resistance phenotype|susceptibility phenotype)\b",
    re.IGNORECASE,
)

ESBL_CRE_PATTERN = re.compile(
    r"\b(esbl|extended[- ]spectrum (?:beta|β)[- ]lactamase|"
    r"carbapenem[- ]resistant|carbapenemase[- ]producing|"
    r"carbapenemase detection|\bcre\b|\bcpe\b|"
    r"ertapenem[- ]resistant|ceftriaxone non[- ]susceptible|"
    r"methicillin[- ]resistant|vancomycin[- ]resistant|\bmrsa\b|\bvre\b)\b",
    re.IGNORECASE,
)

# Curated resistance genes. This intentionally does not accept arbitrary words
# beginning with "bla" and therefore cannot match bladder, blast, black,
# Bland-Altman, blinatumomab, or similar non-AMR words.
BETA_LACTAMASE_GENE_PATTERN = re.compile(
    r"\bbla(?:"
    r"CTX[- ]?M(?:[- ]?\d+)?|KPC(?:[- ]?\d+)?|NDM(?:[- ]?\d+)?|"
    r"OXA[- ]?\d+|VIM(?:[- ]?\d+)?|IMP(?:[- ]?\d+)?|"
    r"CMY(?:[- ]?\d+)?|SHV(?:[- ]?\d+)?|TEM(?:[- ]?\d+)?|"
    r"GES(?:[- ]?\d+)?|PER(?:[- ]?\d+)?|VEB(?:[- ]?\d+)?|"
    r"DHA(?:[- ]?\d+)?|ADC(?:[- ]?\d+)?|PDC(?:[- ]?\d+)?|"
    r"ACC(?:[- ]?\d+)?|ACT(?:[- ]?\d+)?|MIR(?:[- ]?\d+)?|"
    r"FOX(?:[- ]?\d+)?|MOX(?:[- ]?\d+)?|LAT(?:[- ]?\d+)?|"
    r"BIL(?:[- ]?\d+)?|ROB(?:[- ]?\d+)?|CARB(?:[- ]?\d+)?|"
    r"SFO(?:[- ]?\d+)?"
    r")\b",
    re.IGNORECASE,
)

OTHER_RESISTANCE_GENE_PATTERN = re.compile(
    r"\b(mecA|mecC|vanA|vanB|vanC|mcr[- ]?[1-9]|"
    r"qnrA|qnrB|qnrS|aac\(6['’]?-Ib-cr\)|"
    r"tet\([A-Z]\)|erm\([A-Z]\)|ermA|ermB|ermC|"
    r"cfr|optrA|poxtA|dfrA\d+|sul[123]|fosA\d*|"
    r"rpoB mutation|gyrA mutation|parC mutation|"
    r"16S rRNA methylase|armA|rmtA|rmtB|rmtC|rmtD|rmtF)\b",
    re.IGNORECASE,
)

RESISTANCE_MARKER_PATTERN = re.compile(
    rf"(?:{BETA_LACTAMASE_GENE_PATTERN.pattern}|"
    rf"{OTHER_RESISTANCE_GENE_PATTERN.pattern})",
    re.IGNORECASE,
)

RESISTANCE_MARKER_GENERIC_PATTERN = re.compile(
    r"\b(carbapenemase gene|beta[- ]lactamase gene|β[- ]lactamase gene|"
    r"methicillin resistance gene|vancomycin resistance gene|"
    r"resistance marker|resistance mutation)\b",
    re.IGNORECASE,
)

RESISTANCE_PREDICTION_PATTERN = re.compile(
    r"\b(genotypic susceptibility|genotypic resistance|"
    r"resistance prediction|predict(?:ion|ing) of antimicrobial resistance|"
    r"whole[- ]genome sequencing.*resistance|"
    r"resistome|genotype[- ]phenotype concordance|"
    r"multiple resistance mechanisms|multimechanism|"
    r"combined resistance mechanisms|mechanism[- ]informed)\b",
    re.IGNORECASE,
)

COPY_NUMBER_PATTERN = re.compile(
    r"\b(copy number|gene dosage|copy[- ]number variation|\bcnv\b|"
    r"multiple copies|increased copies|gene amplification)\b",
    re.IGNORECASE,
)

EXPRESSION_PATTERN = re.compile(
    r"\b(gene expression|transcript abundance|transcription level|"
    r"expression level|overexpression|mrna abundance)\b",
    re.IGNORECASE,
)

PORIN_PATTERN = re.compile(
    r"\b(porin(?: loss| deficiency| expression)?|ompC|ompF|"
    r"ompK35|ompK36|outer membrane porin|outer membrane protein loss)\b",
    re.IGNORECASE,
)

EFFLUX_PATTERN = re.compile(
    r"\b(efflux pump|acrA|acrB|tolC|mexA|mexB|mexC|mexD|mexE|mexF|"
    r"adeA|adeB|adeC|oqxA|oqxB)\b",
    re.IGNORECASE,
)

ENZYME_ACTIVITY_PATTERN = re.compile(
    r"\b(beta[- ]lactamase activity|β[- ]lactamase activity|"
    r"carbapenemase activity|carbapenem hydrolysis|"
    r"beta[- ]lactam hydrolysis|β[- ]lactam hydrolysis)\b",
    re.IGNORECASE,
)

PROTEIN_ABUNDANCE_PATTERN = re.compile(
    r"\b(protein abundance|protein level|enzyme abundance|"
    r"beta[- ]lactamase concentration|β[- ]lactamase concentration)\b",
    re.IGNORECASE,
)

QUANTITATIVE_MEASUREMENT_PATTERN = re.compile(
    r"\b(quantif(?:y|ication)|quantitative|absolute abundance|"
    r"relative abundance|copies per|copy number|expression level|"
    r"activity level|protein level|concentration)\b",
    re.IGNORECASE,
)

ONCOLOGY_PATTERN = re.compile(
    r"\b(cancer|carcinoma|tumou?r|leukemia|leukaemia|lymphoma|myeloma|"
    r"chemotherapy|neoadjuvant|metastatic|solid tumou?r|her2|egfr|"
    r"pd[- ]?1|pd[- ]?l1|blinatumomab|blast count)\b",
    re.IGNORECASE,
)

CARDIO_METABOLIC_PATTERN = re.compile(
    r"\b(heart failure|coronary artery|myocardial|diabetes|glucose monitor|"
    r"hypertension|preeclampsia|erectile dysfunction|frailty|alzheimer|"
    r"multiple sclerosis|sarcopenia)\b",
    re.IGNORECASE,
)

IMAGING_ONLY_PATTERN = re.compile(
    r"\b(mri|magnetic resonance|computed tomography|\bct scan\b|"
    r"ultrasound|echocardiograph|pet scan|radiograph|x[- ]ray|imaging)\b",
    re.IGNORECASE,
)

MICROBIOME_PATTERN = re.compile(
    r"\b(gut microbiome|gut microbiota|fecal microbiota|faecal microbiota|"
    r"microbiome composition|microbiota composition|dysbiosis|probiotic|"
    r"prebiotic|fecal microbiota transplantation)\b",
    re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Registry extraction
# ---------------------------------------------------------------------------


def extract_study_fields(study: dict[str, Any]) -> tuple[dict[str, str], dict[str, Any]]:
    protocol = get_nested(study, "protocolSection", default={}) or {}
    identification = protocol.get("identificationModule", {}) or {}
    status = protocol.get("statusModule", {}) or {}
    description = protocol.get("descriptionModule", {}) or {}
    conditions_module = protocol.get("conditionsModule", {}) or {}
    design = protocol.get("designModule", {}) or {}
    arms_module = protocol.get("armsInterventionsModule", {}) or {}
    outcomes_module = protocol.get("outcomesModule", {}) or {}
    eligibility_module = protocol.get("eligibilityModule", {}) or {}
    sponsor_module = protocol.get("sponsorCollaboratorsModule", {}) or {}
    lead_sponsor = sponsor_module.get("leadSponsor", {}) or {}
    enrollment = design.get("enrollmentInfo", {}) or {}

    title = join_unique(
        [
            identification.get("briefTitle"),
            identification.get("officialTitle"),
            identification.get("acronym"),
        ]
    )
    conditions_keywords = join_unique(
        as_list(conditions_module.get("conditions"))
        + as_list(conditions_module.get("keywords"))
    )

    intervention_names: list[Any] = []
    intervention_descriptions: list[Any] = []
    has_diagnostic_test_intervention = False
    for intervention in as_list(arms_module.get("interventions")):
        if not isinstance(intervention, dict):
            continue
        intervention_type = clean_text(intervention.get("type")).upper()
        if intervention_type == "DIAGNOSTIC_TEST":
            has_diagnostic_test_intervention = True
        intervention_names.extend(
            [
                intervention.get("type"),
                intervention.get("name"),
                *as_list(intervention.get("otherNames")),
            ]
        )
        intervention_descriptions.append(intervention.get("description"))

    arm_values: list[Any] = []
    for arm in as_list(arms_module.get("armGroups")):
        if isinstance(arm, dict):
            arm_values.extend(
                [
                    arm.get("label"),
                    arm.get("type"),
                    arm.get("description"),
                    *as_list(arm.get("interventionNames")),
                ]
            )

    def outcome_text(outcomes: list[Any]) -> str:
        values: list[Any] = []
        for outcome in outcomes:
            if isinstance(outcome, dict):
                values.extend(
                    [
                        outcome.get("measure"),
                        outcome.get("description"),
                        outcome.get("timeFrame"),
                    ]
                )
        return join_unique(values)

    primary_outcomes = outcome_text(as_list(outcomes_module.get("primaryOutcomes")))
    secondary_outcomes = outcome_text(
        as_list(outcomes_module.get("secondaryOutcomes"))
        + as_list(outcomes_module.get("otherOutcomes"))
    )

    summary = join_unique(
        [
            description.get("briefSummary"),
            description.get("detailedDescription"),
            eligibility_module.get("studyPopulation"),
        ]
    )

    fields = {
        "title": title,
        "conditions_keywords": conditions_keywords,
        "intervention_names": join_unique(intervention_names),
        "intervention_descriptions": join_unique(intervention_descriptions),
        "arm_descriptions": join_unique(arm_values),
        "primary_outcomes": primary_outcomes,
        "secondary_outcomes": secondary_outcomes,
        "summary": summary,
        "eligibility": clean_text(eligibility_module.get("eligibilityCriteria")),
    }

    metadata = {
        "nct_id": clean_text(identification.get("nctId")),
        "brief_title": clean_text(identification.get("briefTitle")),
        "official_title": clean_text(identification.get("officialTitle")),
        "overall_status": clean_text(status.get("overallStatus")),
        "study_type": clean_text(design.get("studyType")),
        "primary_purpose": clean_text(
            get_nested(design, "designInfo", "primaryPurpose")
        ).upper(),
        "has_diagnostic_test_intervention": int(has_diagnostic_test_intervention),
        "phases": join_unique(as_list(design.get("phases"))),
        "enrollment_count": enrollment.get("count", ""),
        "enrollment_type": clean_text(enrollment.get("type")),
        "lead_sponsor_name": clean_text(lead_sponsor.get("name")),
        "lead_sponsor_class": clean_text(lead_sponsor.get("class")),
        "start_date": clean_text(get_nested(status, "startDateStruct", "date")),
        "completion_date": clean_text(
            get_nested(status, "completionDateStruct", "date")
        ),
        "study_first_post_date": clean_text(
            get_nested(status, "studyFirstPostDateStruct", "date")
        ),
        "results_first_post_date": clean_text(
            get_nested(status, "resultsFirstPostDateStruct", "date")
        ),
        "has_results": int(bool(study.get("hasResults", False))),
    }
    return fields, metadata


# ---------------------------------------------------------------------------
# Staged evidence evaluation
# ---------------------------------------------------------------------------


def group_field_evidence(
    fields: dict[str, str],
    patterns: dict[str, re.Pattern[str]],
    *,
    allowed_fields: set[str] | None = None,
) -> tuple[dict[str, int], dict[str, list[str]], dict[str, list[str]]]:
    flags: dict[str, int] = {}
    evidence_fields: dict[str, list[str]] = {}
    evidence_terms: dict[str, list[str]] = {}
    for category, pattern in patterns.items():
        matched_fields, terms, _ = pattern_evidence_by_field(
            fields,
            pattern,
            allowed_fields=allowed_fields,
        )
        flags[category] = int(bool(matched_fields))
        evidence_fields[category] = matched_fields
        evidence_terms[category] = terms
    return flags, evidence_fields, evidence_terms


def fields_intersect(
    evidence_fields: dict[str, list[str]],
    allowed: set[str],
) -> bool:
    return any(
        field_name in allowed
        for matched_fields in evidence_fields.values()
        for field_name in matched_fields
    )


def categories_in_fields(
    evidence_fields: dict[str, list[str]],
    allowed: set[str],
) -> list[str]:
    return sorted(
        category
        for category, matched_fields in evidence_fields.items()
        if any(field_name in allowed for field_name in matched_fields)
    )


def detect_mechanisms(
    fields: dict[str, str],
) -> tuple[dict[str, int], dict[str, list[str]], dict[str, list[str]]]:
    flags = {
        "copy_number": 0,
        "gene_expression": 0,
        "porin": 0,
        "efflux": 0,
        "enzyme_activity": 0,
        "protein_abundance": 0,
        "multimechanism": 0,
    }
    evidence_fields: dict[str, list[str]] = defaultdict(list)
    snippets: dict[str, list[str]] = defaultdict(list)

    amr_anchor = re.compile(
        rf"(?:{RESISTANCE_MARKER_PATTERN.pattern}|"
        rf"{RESISTANCE_MARKER_GENERIC_PATTERN.pattern}|"
        rf"{PORIN_PATTERN.pattern}|{EFFLUX_PATTERN.pattern}|"
        rf"beta[- ]lactamase|β[- ]lactamase|carbapenemase|"
        rf"resistance gene|resistance mechanism)",
        re.IGNORECASE,
    )

    for field_name, text in fields.items():
        if field_name not in MECHANISM_SCAN_FIELDS:
            continue
        if near_patterns(text, COPY_NUMBER_PATTERN, amr_anchor, window=180):
            flags["copy_number"] = 1
            evidence_fields["copy_number"].append(field_name)
            snippets["copy_number"].extend(
                f"{field_name}: {snippet}"
                for snippet in snippets_for_pattern(text, COPY_NUMBER_PATTERN)
            )

        if near_patterns(text, EXPRESSION_PATTERN, amr_anchor, window=180):
            flags["gene_expression"] = 1
            evidence_fields["gene_expression"].append(field_name)
            snippets["gene_expression"].extend(
                f"{field_name}: {snippet}"
                for snippet in snippets_for_pattern(text, EXPRESSION_PATTERN)
            )

        if PORIN_PATTERN.search(text):
            flags["porin"] = 1
            evidence_fields["porin"].append(field_name)
            snippets["porin"].extend(
                f"{field_name}: {snippet}"
                for snippet in snippets_for_pattern(text, PORIN_PATTERN)
            )

        if EFFLUX_PATTERN.search(text):
            flags["efflux"] = 1
            evidence_fields["efflux"].append(field_name)
            snippets["efflux"].extend(
                f"{field_name}: {snippet}"
                for snippet in snippets_for_pattern(text, EFFLUX_PATTERN)
            )

        if ENZYME_ACTIVITY_PATTERN.search(text):
            flags["enzyme_activity"] = 1
            evidence_fields["enzyme_activity"].append(field_name)
            snippets["enzyme_activity"].extend(
                f"{field_name}: {snippet}"
                for snippet in snippets_for_pattern(text, ENZYME_ACTIVITY_PATTERN)
            )

        if near_patterns(text, PROTEIN_ABUNDANCE_PATTERN, amr_anchor, window=180):
            flags["protein_abundance"] = 1
            evidence_fields["protein_abundance"].append(field_name)
            snippets["protein_abundance"].extend(
                f"{field_name}: {snippet}"
                for snippet in snippets_for_pattern(text, PROTEIN_ABUNDANCE_PATTERN)
            )

        if RESISTANCE_PREDICTION_PATTERN.search(text):
            flags["multimechanism"] = 1
            evidence_fields["multimechanism"].append(field_name)
            snippets["multimechanism"].extend(
                f"{field_name}: {snippet}"
                for snippet in snippets_for_pattern(text, RESISTANCE_PREDICTION_PATTERN)
            )

    if sum(flags.values()) >= 2:
        flags["multimechanism"] = 1

    return (
        flags,
        {key: sorted(set(value)) for key, value in evidence_fields.items()},
        {
            key: list(dict.fromkeys(value))[:10]
            for key, value in snippets.items()
        },
    )


def evaluate_study(study: dict[str, Any]) -> dict[str, Any]:
    fields, metadata = extract_study_fields(study)
    nct_id = metadata["nct_id"]
    if not nct_id:
        return {"nct_id": "", "predicted_stratum": "EXCLUDED"}

    # Cheap stage-0 prefilter. Most registry records are unrelated to
    # infectious diseases. Avoid running the expensive diagnostic, AMR, and
    # mechanism regex families unless a bacterial/infection/special-pathogen
    # term occurs in a primary relevance field.
    prefilter_text = " | ".join(
        fields[name] for name in PRIMARY_RELEVANCE_FIELDS
    )
    prefilter_match = bool(GENERAL_INFECTION_PATTERN.search(prefilter_text))
    if not prefilter_match:
        prefilter_match = any(
            pattern.search(prefilter_text)
            for pattern in ORGANISM_PATTERNS.values()
        )
    if not prefilter_match:
        prefilter_match = any(
            pattern.search(prefilter_text)
            for pattern in SPECIAL_PATHOGEN_PATTERNS.values()
        )
    if not prefilter_match:
        prefilter_match = any(
            pattern.search(prefilter_text)
            for pattern in BACTERIAL_SYNDROME_PATTERNS.values()
        )

    if not prefilter_match:
        return {
            **metadata,
            "predicted_stratum": "EXCLUDED",
            "classification_reason": (
                "no bacterial, infectious-disease, or special-pathogen evidence "
                "in primary relevance fields"
            ),
            "diagnostic_depth_level": "",
            "diagnostic_depth_label": "",
            "infection_score": 0,
            "diagnostic_score": 0,
            "amr_score": 0,
            "near_miss_score": 0,
            "development_control": int(nct_id in DEVELOPMENT_CONTROL_IDS),
            "conditions_keywords": fields["conditions_keywords"][:5000],
            "intervention_names": fields["intervention_names"][:5000],
            "primary_outcomes": fields["primary_outcomes"][:7000],
            "secondary_outcomes": fields["secondary_outcomes"][:7000],
            "summary": fields["summary"][:9000],
            "clinicaltrials_url": f"https://clinicaltrials.gov/study/{nct_id}",
        }

    organism_flags, organism_fields, organism_terms = group_field_evidence(
        fields, ORGANISM_PATTERNS, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    special_flags, special_fields, special_terms = group_field_evidence(
        fields, SPECIAL_PATHOGEN_PATTERNS, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    syndrome_flags, syndrome_fields, syndrome_terms = group_field_evidence(
        fields, BACTERIAL_SYNDROME_PATTERNS, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    modality_flags, modality_fields, modality_terms = group_field_evidence(
        fields, DIRECT_DIAGNOSTIC_PATTERNS, allowed_fields=EVIDENCE_SCAN_FIELDS
    )

    direct_modality_high = fields_intersect(modality_fields, HIGH_VALUE_FIELDS)
    direct_modality_primary = fields_intersect(
        modality_fields, PRIMARY_RELEVANCE_FIELDS
    )

    organism_primary = fields_intersect(organism_fields, PRIMARY_RELEVANCE_FIELDS)
    organism_high = fields_intersect(organism_fields, HIGH_VALUE_FIELDS)
    syndrome_primary = fields_intersect(syndrome_fields, PRIMARY_RELEVANCE_FIELDS)
    syndrome_high = fields_intersect(syndrome_fields, HIGH_VALUE_FIELDS)
    special_primary = fields_intersect(special_fields, PRIMARY_RELEVANCE_FIELDS)
    special_high = fields_intersect(special_fields, HIGH_VALUE_FIELDS)

    general_infection_fields, general_infection_terms, general_infection_snips = (
        pattern_evidence_by_field(
            fields, GENERAL_INFECTION_PATTERN, allowed_fields=EVIDENCE_SCAN_FIELDS
        )
    )
    general_infection_primary = any(
        field_name in PRIMARY_RELEVANCE_FIELDS
        for field_name in general_infection_fields
    )

    host_fields, host_terms, host_snips = pattern_evidence_by_field(
        fields, HOST_RESPONSE_PATTERN, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    host_immune_fields, host_immune_terms, host_immune_snips = (
        pattern_evidence_by_field(
            fields,
            HOST_IMMUNE_ANALYTE_PATTERN,
            allowed_fields=EVIDENCE_SCAN_FIELDS,
        )
    )
    direct_microbial_analyte_fields, direct_microbial_analyte_terms, (
        direct_microbial_analyte_snips
    ) = pattern_evidence_by_field(
        fields,
        DIRECT_MICROBIAL_ANALYTE_PATTERN,
        allowed_fields=EVIDENCE_SCAN_FIELDS,
    )
    host_high = any(
        field_name in HIGH_VALUE_FIELDS
        for field_name in set(host_fields) | set(host_immune_fields)
    )
    host_immune_analyte_high = any(
        field_name in HIGH_VALUE_FIELDS for field_name in host_immune_fields
    )
    direct_microbial_analyte_high = any(
        field_name in HIGH_VALUE_FIELDS
        for field_name in direct_microbial_analyte_fields
    )
    host_diagnostic_fields, host_diagnostic_terms, host_diagnostic_snips = (
        pattern_evidence_by_field(
            fields,
            HOST_DIAGNOSTIC_PATTERN,
            allowed_fields=EVIDENCE_SCAN_FIELDS,
        )
    )
    host_immune_diagnostic_fields, host_immune_diagnostic_terms, (
        host_immune_diagnostic_snips
    ) = pattern_evidence_by_field(
        fields,
        HOST_IMMUNE_DIAGNOSTIC_PATTERN,
        allowed_fields=EVIDENCE_SCAN_FIELDS,
    )
    host_diagnostic_high = any(
        field_name in HIGH_VALUE_FIELDS
        for field_name in set(host_diagnostic_fields)
        | set(host_immune_diagnostic_fields)
    )
    prognostic_fields, prognostic_terms, prognostic_snips = pattern_evidence_by_field(
        fields,
        PROGNOSTIC_PATTERN,
        allowed_fields=EVIDENCE_SCAN_FIELDS,
    )
    prognostic_high = any(
        field_name in HIGH_VALUE_FIELDS for field_name in prognostic_fields
    )
    development_fields, development_terms, development_snips = pattern_evidence_by_field(
        fields,
        DIAGNOSTIC_DEVELOPMENT_PATTERN,
        allowed_fields=EVIDENCE_SCAN_FIELDS,
    )
    development_high = any(
        field_name in HIGH_VALUE_FIELDS for field_name in development_fields
    )

    diagnostic_eval_fields, diagnostic_eval_terms, diagnostic_eval_snips = (
        pattern_evidence_by_field(
            fields,
            DIAGNOSTIC_EVALUATION_PATTERN,
            allowed_fields=EVIDENCE_SCAN_FIELDS,
        )
    )
    diagnostic_eval_high = any(
        field_name in HIGH_VALUE_FIELDS for field_name in diagnostic_eval_fields
    )

    diagnostic_role_fields, diagnostic_role_terms, diagnostic_role_snips = (
        pattern_evidence_by_field(
            fields,
            DIAGNOSTIC_ROLE_PATTERN,
            allowed_fields=EVIDENCE_SCAN_FIELDS,
        )
    )
    diagnostic_role_high = any(
        field_name in HIGH_VALUE_FIELDS for field_name in diagnostic_role_fields
    )
    microbial_focus_fields, microbial_focus_terms, microbial_focus_snips = (
        pattern_evidence_by_field(
            fields,
            MICROBIAL_ASSAY_FOCUS_PATTERN,
            allowed_fields={"title", "intervention_names"},
        )
    )

    named_bacterial_assay_title = bool(
        NAMED_BACTERIAL_ASSAY_TITLE_PATTERN.search(fields["title"])
    )
    named_organism_in_title = any(
        "title" in organism_fields.get(category, [])
        for category in [
            "enterobacterales",
            "s_aureus",
            "other_gram_positive",
            "nonfermenter",
            "other_bacterial",
        ]
    )
    named_organism_detection_title = bool(
        named_organism_in_title
        and NAMED_ORGANISM_DETECTION_WORDING_PATTERN.search(fields["title"])
    )
    clinical_metagenomics_title = bool(
        CLINICAL_METAGENOMICS_TITLE_PATTERN.search(fields["title"])
        and (
            any("title" in matched for matched in syndrome_fields.values())
            or any(
                field_name in {"title", "conditions_keywords"}
                for matched in syndrome_fields.values()
                for field_name in matched
            )
        )
    )
    syndrome_direct_diagnostic_title = bool(
        SYNDROME_DIRECT_DIAGNOSTIC_TITLE_PATTERN.search(fields["title"])
        and any(
            field_name in {"title", "conditions_keywords"}
            for matched in syndrome_fields.values()
            for field_name in matched
        )
    )
    point_of_care_pathogen_testing_title = bool(
        POINT_OF_CARE_PATHOGEN_TESTING_TITLE_PATTERN.search(fields["title"])
    )
    syndromic_diagnostic_implementation_title = bool(
        SYNDROMIC_DIAGNOSTIC_IMPLEMENTATION_TITLE_PATTERN.search(fields["title"])
    )
    microbial_target_fields, microbial_target_terms, microbial_target_snips = (
        pattern_evidence_by_field(
            fields,
            MICROBIAL_TARGET_PATTERN,
            allowed_fields={
                "title",
                "intervention_names",
                "primary_outcomes",
                "secondary_outcomes",
            },
        )
    )
    organ_injury_fields, organ_injury_terms, organ_injury_snips = (
        pattern_evidence_by_field(
            fields,
            ORGAN_INJURY_PROGNOSTIC_TARGET_PATTERN,
            allowed_fields=HIGH_VALUE_FIELDS,
        )
    )
    mixed_panel_fields, mixed_panel_terms, mixed_panel_snips = (
        pattern_evidence_by_field(
            fields,
            MIXED_VIRAL_BACTERIAL_PANEL_PATTERN,
            allowed_fields=EVIDENCE_SCAN_FIELDS,
        )
    )

    reference_fields, reference_terms, reference_snips = pattern_evidence_by_field(
        fields, REFERENCE_STANDARD_PATTERN, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    reference_high = any(field_name in HIGH_VALUE_FIELDS for field_name in reference_fields)

    action_fields, action_terms, action_snips = pattern_evidence_by_field(
        fields, ANTIBIOTIC_ACTION_PATTERN, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    action_high = any(field_name in HIGH_VALUE_FIELDS for field_name in action_fields)

    surveillance_fields, surveillance_terms, _ = pattern_evidence_by_field(
        fields, SURVEILLANCE_PATTERN, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    surveillance_high = any(
        field_name in HIGH_VALUE_FIELDS for field_name in surveillance_fields
    )

    therapeutic_fields, therapeutic_terms, _ = pattern_evidence_by_field(
        fields, THERAPEUTIC_PATTERN, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    treatment_primary_purpose = metadata["primary_purpose"] == "TREATMENT"
    treatment_or_prevention_purpose = metadata["primary_purpose"] in {
        "TREATMENT",
        "PREVENTION",
    }

    general_amr_fields, general_amr_terms, general_amr_snips = pattern_evidence_by_field(
        fields, GENERAL_AMR_PATTERN, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    ast_fields, ast_terms, ast_snips = pattern_evidence_by_field(
        fields, AST_PATTERN, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    esbl_fields, esbl_terms, esbl_snips = pattern_evidence_by_field(
        fields, ESBL_CRE_PATTERN, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    marker_fields, marker_terms, marker_snips = pattern_evidence_by_field(
        fields, RESISTANCE_MARKER_PATTERN, allowed_fields=EVIDENCE_SCAN_FIELDS
    )
    marker_generic_fields, marker_generic_terms, marker_generic_snips = (
        pattern_evidence_by_field(
            fields,
            RESISTANCE_MARKER_GENERIC_PATTERN,
            allowed_fields=EVIDENCE_SCAN_FIELDS,
        )
    )
    prediction_fields, prediction_terms, prediction_snips = pattern_evidence_by_field(
        fields, RESISTANCE_PREDICTION_PATTERN, allowed_fields=EVIDENCE_SCAN_FIELDS
    )

    mechanism_flags, mechanism_fields, mechanism_snips = detect_mechanisms(fields)

    bacterial_primary_context = bool(
        organism_primary
        or syndrome_primary
        or general_infection_primary
    )
    bacterial_high_context = bool(
        organism_high
        or syndrome_high
        or any(field_name in HIGH_VALUE_FIELDS for field_name in general_infection_fields)
    )

    # An outcome/intervention-only bacterial match must co-occur in the same
    # high-value field with diagnostic evidence. This blocks studies where
    # sepsis/infection is merely an adverse event and an unrelated diagnostic
    # test appears elsewhere in the record.
    bacterial_title_or_condition = bool(
        fields_intersect(organism_fields, {"title", "conditions_keywords"})
        or fields_intersect(syndrome_fields, {"title", "conditions_keywords"})
        or any(
            field_name in {"title", "conditions_keywords"}
            for field_name in general_infection_fields
        )
    )
    if not bacterial_title_or_condition:
        bacterial_evidence_fields = {
            field_name
            for matched_fields in organism_fields.values()
            for field_name in matched_fields
        } | {
            field_name
            for matched_fields in syndrome_fields.values()
            for field_name in matched_fields
        } | set(general_infection_fields)
        diagnostic_evidence_field_set = {
            field_name
            for matched_fields in modality_fields.values()
            for field_name in matched_fields
        } | set(diagnostic_eval_fields) | set(reference_fields) | set(host_fields)
        same_field_diagnostic_infection = bool(
            bacterial_evidence_fields
            & diagnostic_evidence_field_set
            & {"intervention_names", "primary_outcomes"}
        )
        bacterial_primary_context = bool(
            bacterial_primary_context and same_field_diagnostic_infection
        )

    microbial_target_field_set = set(microbial_target_fields)
    specific_organism_target_high = any(
        field_name in {
            "title",
            "intervention_names",
            "primary_outcomes",
            "secondary_outcomes",
        }
        for category in [
            "enterobacterales",
            "s_aureus",
            "other_gram_positive",
            "nonfermenter",
            "other_bacterial",
        ]
        for field_name in organism_fields.get(category, [])
    )
    microbial_target_high = bool(
        named_bacterial_assay_title
        or named_organism_detection_title
        or clinical_metagenomics_title
        or syndrome_direct_diagnostic_title
        or point_of_care_pathogen_testing_title
        or direct_microbial_analyte_high
        or microbial_target_field_set
        or specific_organism_target_high
    )
    organ_injury_target_high = bool(organ_injury_fields)
    evidence_text_for_host_panel = " | ".join(
        fields[field_name] for field_name in EVIDENCE_SCAN_FIELDS
    )
    febridx_host_panel = bool(FEBRIDX_PATTERN.search(evidence_text_for_host_panel))
    mxa_crp_host_panel = bool(
        MXA_PATTERN.search(evidence_text_for_host_panel)
        and CRP_PATTERN.search(evidence_text_for_host_panel)
    )
    host_mixed_bacterial_viral_test = bool(
        febridx_host_panel or mxa_crp_host_panel
    )
    mixed_viral_bacterial_panel = bool(
        mixed_panel_fields or host_mixed_bacterial_viral_test
    )

    structured_diagnostic = bool(
        metadata["primary_purpose"] == "DIAGNOSTIC"
        or metadata["has_diagnostic_test_intervention"]
    )

    # Version 3.2.5: a mixed viral/bacterial point-of-care pathogen panel is
    # not a nonbacterial-only study merely because viral targets are named in
    # the intervention or outcomes. Require a strong direct-testing anchor so
    # generic viral diagnostics cannot enter the bacterial primary cohort.
    mixed_direct_pathogen_panel = bool(
        point_of_care_pathogen_testing_title
        and mixed_viral_bacterial_panel
        and (structured_diagnostic or direct_modality_high or action_high)
    )
    modality_title = any(
        "title" in matched_fields for matched_fields in modality_fields.values()
    )
    diagnostic_title_focus = bool(
        named_bacterial_assay_title
        or named_organism_detection_title
        or clinical_metagenomics_title
        or syndrome_direct_diagnostic_title
        or point_of_care_pathogen_testing_title
        or microbial_focus_fields
        or (
            modality_title
            and (
                "title" in diagnostic_role_fields
                or "title" in diagnostic_eval_fields
            )
        )
    )

    explicit_diagnostic_development = bool(
        development_high or syndromic_diagnostic_implementation_title
    )

    infectious_primary_context = bool(
        bacterial_primary_context or special_primary
    )

    direct_pathogen_diagnostic_intent = bool(
        infectious_primary_context
        and microbial_target_high
        and not (host_high and not direct_modality_high)
        and not IMAGING_ONLY_PATTERN.search(
            " | ".join(
                [
                    fields["title"],
                    fields["intervention_names"],
                    fields["primary_outcomes"],
                ]
            )
        )
        and (
            named_bacterial_assay_title
            or named_organism_detection_title
            or clinical_metagenomics_title
            or (
                syndrome_direct_diagnostic_title
                and direct_modality_high
                and not host_high
            )
            or point_of_care_pathogen_testing_title
            or diagnostic_title_focus
            or (
                reference_high
                and direct_modality_high
            )
            or (
                diagnostic_eval_high
                and diagnostic_role_high
                and direct_modality_high
            )
            or (
                structured_diagnostic
                and diagnostic_role_high
                and direct_modality_primary
                and (diagnostic_eval_high or action_high)
            )
            or (
                action_high
                and diagnostic_role_high
                and direct_modality_primary
            )
        )
    )

    imaging_syndromic_diagnostic_intent = bool(
        bacterial_primary_context
        and not host_high
        and IMAGING_ONLY_PATTERN.search(
            " | ".join(
                [
                    fields["title"],
                    fields["intervention_names"],
                    fields["primary_outcomes"],
                ]
            )
        )
        and (
            diagnostic_eval_high
            or structured_diagnostic
            or diagnostic_role_high
        )
    )

    clinical_syndromic_diagnostic_intent = bool(
        imaging_syndromic_diagnostic_intent
        or (
            bacterial_primary_context
            and syndromic_diagnostic_implementation_title
            and not direct_pathogen_diagnostic_intent
            and not host_high
        )
        or (
            bacterial_primary_context
            and not microbial_target_high
            and not host_high
            and (
                explicit_diagnostic_development
                or (
                    diagnostic_eval_high
                    and diagnostic_role_high
                )
                or (
                    structured_diagnostic
                    and diagnostic_role_high
                )
            )
        )
    )

    host_response_diagnostic_intent = bool(
        infectious_primary_context
        and host_high
        and host_diagnostic_high
        and (
            host_immune_analyte_high
            or bool(host_fields)
            or host_mixed_bacterial_viral_test
        )
        and not (
            prognostic_high
            and not diagnostic_eval_high
            and not host_diagnostic_high
        )
    )

    # Direct and host-response assays can coexist. Direct pathogen testing takes
    # priority for primary classification, while the host flag remains available.
    any_diagnostic_intent = bool(
        direct_pathogen_diagnostic_intent
        or host_response_diagnostic_intent
        or clinical_syndromic_diagnostic_intent
    )

    general_amr_high = any(
        field_name in HIGH_VALUE_FIELDS for field_name in general_amr_fields
    )
    ast_high = any(field_name in HIGH_VALUE_FIELDS for field_name in ast_fields)
    esbl_high = any(field_name in HIGH_VALUE_FIELDS for field_name in esbl_fields)
    marker_high = any(field_name in HIGH_VALUE_FIELDS for field_name in marker_fields)
    marker_generic_high = any(
        field_name in HIGH_VALUE_FIELDS for field_name in marker_generic_fields
    )
    prediction_high = any(
        field_name in HIGH_VALUE_FIELDS for field_name in prediction_fields
    )

    explicit_amr_high = bool(
        general_amr_high
        or ast_high
        or esbl_high
        or marker_high
        or marker_generic_high
        or prediction_high
    )

    any_amr_evidence = bool(
        general_amr_fields
        or ast_fields
        or esbl_fields
        or marker_fields
        or marker_generic_fields
        or prediction_fields
        or any(mechanism_flags.values())
    )

    mechanism_high = any(
        field_name in HIGH_VALUE_FIELDS
        for matched_fields in mechanism_fields.values()
        for field_name in matched_fields
    )

    quantitative_high = bool(
        mechanism_high
        and any(
            QUANTITATIVE_MEASUREMENT_PATTERN.search(fields[field_name])
            for matched_fields in mechanism_fields.values()
            for field_name in matched_fields
            if field_name in HIGH_VALUE_FIELDS
        )
    )

    # Diagnostic depth is assigned only to direct pathogen diagnostics.
    depth_level: int | str = ""
    depth_label = ""
    if direct_pathogen_diagnostic_intent:
        if quantitative_high and any(mechanism_flags.values()):
            depth_level = 4
        elif prediction_high or mechanism_flags["multimechanism"]:
            depth_level = 3
        elif ast_high or modality_flags.get("phenotypic_ast", 0):
            depth_level = 2
        elif marker_high or marker_generic_high or esbl_high:
            depth_level = 1
        else:
            depth_level = 0
        depth_label = DEPTH_LABELS[int(depth_level)]

    # Background-risk flags cannot exclude a study if high-value bacterial and
    # diagnostic evidence is present.
    high_text = " | ".join(fields[name] for name in HIGH_VALUE_FIELDS)
    background_text = " | ".join(
        fields[name] for name in MID_VALUE_FIELDS | LOW_VALUE_FIELDS
    )
    oncology_high = bool(ONCOLOGY_PATTERN.search(high_text))
    oncology_background_only = bool(
        ONCOLOGY_PATTERN.search(background_text) and not oncology_high
    )
    cardiometabolic_high = bool(CARDIO_METABOLIC_PATTERN.search(high_text))
    microbiome_high = bool(MICROBIOME_PATTERN.search(high_text))
    imaging_only_high = bool(
        IMAGING_ONLY_PATTERN.search(high_text)
        and not direct_modality_high
        and not host_high
    )

    special_title_condition = bool(
        fields_intersect(special_fields, {"title", "conditions_keywords"})
    )
    special_diagnostic_same_field = bool(
        {
            field_name
            for matched_fields in special_fields.values()
            for field_name in matched_fields
        }
        & (
            {
                field_name
                for matched_fields in modality_fields.values()
                for field_name in matched_fields
            }
            | set(diagnostic_eval_fields)
            | set(reference_fields)
        )
        & {"intervention_names", "primary_outcomes"}
    )
    specific_non_special_organism_primary = any(
        any(field_name in PRIMARY_RELEVANCE_FIELDS for field_name in organism_fields.get(category, []))
        for category in [
            "enterobacterales",
            "s_aureus",
            "other_gram_positive",
            "nonfermenter",
            "other_bacterial",
        ]
    )
    special_only_context = bool(
        special_title_condition
        and not specific_non_special_organism_primary
    )

    strict_mechanism = any(mechanism_flags.values())
    clinical_diagnostic_evaluation = bool(
        diagnostic_eval_high
        or reference_high
        or action_high
        or structured_diagnostic
    )
    surveillance_title_condition = any(
        field_name in {"title", "conditions_keywords"}
        for field_name in surveillance_fields
    )
    surveillance_dominant = bool(
        surveillance_title_condition
        and not diagnostic_eval_high
        and not reference_high
        and not action_high
    )
    therapeutic_dominant = bool(
        treatment_or_prevention_purpose
        and not metadata["has_diagnostic_test_intervention"]
        and not action_high
        and not explicit_diagnostic_development
        and not diagnostic_title_focus
        and not named_bacterial_assay_title
        and not named_organism_detection_title
        and not clinical_metagenomics_title
        and not syndrome_direct_diagnostic_title
        and not point_of_care_pathogen_testing_title
    )

    explicit_bacterial_primary = bool(
        organism_primary
        or fields_intersect(
            {"mixed_pan_bacterial": organism_fields.get("mixed_pan_bacterial", [])},
            PRIMARY_RELEVANCE_FIELDS,
        )
        or explicit_amr_high
        or mixed_direct_pathogen_panel
    )
    nonbacterial_only_high_context = bool(
        NONBACTERIAL_PATHOGEN_PATTERN.search(high_text)
        and not explicit_bacterial_primary
        and not host_response_diagnostic_intent
        and not special_title_condition
    )
    special_diagnostic_focus = bool(
        special_only_context
        and any_diagnostic_intent
        and (
            diagnostic_eval_high
            or reference_high
            or action_high
            or explicit_diagnostic_development
            or diagnostic_title_focus
        )
        and not therapeutic_dominant
    )

    predicted_stratum = "EXCLUDED"
    classification_reason = "did not meet staged bacterial diagnostic criteria"

    mechanism_measurement_only = bool(
        bacterial_primary_context
        and strict_mechanism
        and any_amr_evidence
        and not clinical_diagnostic_evaluation
    )

    if special_diagnostic_focus:
        predicted_stratum = "SPECIAL_PATHOGEN_DIAGNOSTIC"
        classification_reason = (
            "special pathogen in a title/condition or same-field diagnostic context"
        )
    elif nonbacterial_only_high_context:
        predicted_stratum = "EXCLUDED"
        classification_reason = (
            "viral, fungal, parasitic, or other nonbacterial high-value context "
            "without independent bacterial evidence"
        )
    elif mechanism_measurement_only:
        predicted_stratum = "MECHANISM_SUPPORT"
        classification_reason = (
            "AMR mechanism study using an assay as a measurement tool without "
            "registered diagnostic evaluation or clinical utility"
        )
    elif host_response_diagnostic_intent:
        predicted_stratum = "HOST_RESPONSE_DIAGNOSTIC"
        classification_reason = (
            "host-response or nonpathogen biomarker diagnostic for suspected infection"
        )
    elif surveillance_dominant:
        predicted_stratum = "SURVEILLANCE_SUPPORT"
        classification_reason = (
            "surveillance or molecular epidemiology is the registered primary purpose"
        )
    elif therapeutic_dominant:
        predicted_stratum = "THERAPEUTIC_SUPPORT"
        classification_reason = (
            "treatment/prevention study without a diagnostic-guided intervention"
        )
    elif direct_pathogen_diagnostic_intent:
        if explicit_amr_high or int(depth_level or 0) >= 1:
            predicted_stratum = "CORE_AMR_DIAGNOSTIC"
            classification_reason = (
                "direct bacterial diagnostic intent plus high-value AMR/AST evidence"
            )
        else:
            predicted_stratum = "BROAD_BACTERIAL_DIAGNOSTIC"
            classification_reason = (
                "direct bacterial diagnostic intent without high-value AMR/AST evidence"
            )
    elif clinical_syndromic_diagnostic_intent:
        predicted_stratum = "CLINICAL_SYNDROMIC_SUPPORT"
        classification_reason = (
            "infection-context diagnostic development or validation without a "
            "direct pathogen assay or host-response analyte"
        )
    elif (
        bacterial_primary_context
        and organ_injury_target_high
        and (host_high or prognostic_high)
        and not direct_pathogen_diagnostic_intent
    ):
        predicted_stratum = "CLINICAL_SYNDROMIC_SUPPORT"
        classification_reason = (
            "infection-context organ-injury or prognostic biomarker support "
            "without direct microbial diagnostic intent"
        )
    elif bacterial_primary_context and strict_mechanism and any_amr_evidence:
        predicted_stratum = "MECHANISM_SUPPORT"
        classification_reason = (
            "bacterial AMR mechanism characterization without direct clinical diagnostic intent"
        )
    elif bacterial_primary_context and surveillance_high and not any_diagnostic_intent:
        predicted_stratum = "SURVEILLANCE_SUPPORT"
        classification_reason = (
            "bacterial surveillance/epidemiology without evaluated diagnostic intent"
        )
    elif (
        bacterial_primary_context
        and (treatment_primary_purpose or therapeutic_fields)
        and not any_diagnostic_intent
    ):
        predicted_stratum = "THERAPEUTIC_SUPPORT"
        classification_reason = (
            "bacterial therapeutic study without evaluated diagnostic intent"
        )

    # Clear unrelated high-value contexts are excluded unless direct bacterial
    # diagnostic evidence independently qualifies them.
    if predicted_stratum not in PRIMARY_DIAGNOSTIC_STRATA and not any_diagnostic_intent:
        if oncology_high or cardiometabolic_high or imaging_only_high:
            if predicted_stratum not in {
                "MECHANISM_SUPPORT",
                "SURVEILLANCE_SUPPORT",
                "THERAPEUTIC_SUPPORT",
            }:
                predicted_stratum = "EXCLUDED"
                classification_reason = "noninfectious high-value clinical context"

    # A study of organ injury, prognosis, or severity in an infection cohort is
    # not a direct microbial diagnostic unless an independent microbial target
    # was identified in title/intervention/outcome fields.
    if (
        predicted_stratum in PRIMARY_DIAGNOSTIC_STRATA
        and organ_injury_target_high
        and not microbial_target_high
    ):
        predicted_stratum = "CLINICAL_SYNDROMIC_SUPPORT"
        classification_reason = (
            "organ-injury/prognostic diagnostic in an infection population "
            "without an independently evaluated microbial target"
        )

    if predicted_stratum not in PRIMARY_DIAGNOSTIC_STRATA:
        depth_level = ""
        depth_label = ""

    diagnostic_score = 0
    diagnostic_score += 5 * int(direct_modality_primary)
    diagnostic_score += 4 * int(diagnostic_eval_high)
    diagnostic_score += 3 * int(reference_high)
    diagnostic_score += 3 * int(structured_diagnostic)
    diagnostic_score += 3 * int(action_high)
    diagnostic_score += 2 * int(diagnostic_role_high)
    diagnostic_score += 2 * int(host_high)
    diagnostic_score += 4 * int(named_bacterial_assay_title)
    diagnostic_score += 4 * int(named_organism_detection_title)
    diagnostic_score += 4 * int(clinical_metagenomics_title)
    diagnostic_score += 4 * int(syndrome_direct_diagnostic_title)
    diagnostic_score += 4 * int(point_of_care_pathogen_testing_title)
    diagnostic_score += 3 * int(syndromic_diagnostic_implementation_title)
    diagnostic_score += 2 * int(microbial_target_high)

    infection_score = 0
    infection_score += 6 * int(
        fields_intersect(organism_fields, {"title", "conditions_keywords"})
    )
    infection_score += 6 * int(
        fields_intersect(syndrome_fields, {"title", "conditions_keywords"})
    )
    infection_score += 3 * int(organism_primary)
    infection_score += 3 * int(syndrome_primary)
    infection_score += 2 * int(general_infection_primary)

    amr_score = 0
    amr_score += 4 * int(ast_high)
    amr_score += 4 * int(marker_high)
    amr_score += 3 * int(esbl_high)
    amr_score += 2 * int(general_amr_high)
    amr_score += 3 * int(prediction_high)
    amr_score += 2 * sum(mechanism_flags.values())

    near_miss_score = infection_score + diagnostic_score + min(amr_score, 5)

    nonbacterial_high = bool(NONBACTERIAL_PATHOGEN_PATTERN.search(high_text))
    if nonbacterial_high and not bacterial_primary_context:
        near_miss_score = max(0, near_miss_score - 6)

    row: dict[str, Any] = {
        **metadata,
        "predicted_stratum": predicted_stratum,
        "classification_reason": classification_reason,
        "diagnostic_depth_level": depth_level,
        "diagnostic_depth_label": depth_label,
        "infection_score": infection_score,
        "diagnostic_score": diagnostic_score,
        "amr_score": amr_score,
        "near_miss_score": near_miss_score,
        "bacterial_primary_context": int(bacterial_primary_context),
        "bacterial_high_context": int(bacterial_high_context),
        "special_primary_context": int(special_primary),
        "direct_pathogen_diagnostic_intent": int(direct_pathogen_diagnostic_intent),
        "named_bacterial_assay_title": int(named_bacterial_assay_title),
        "named_organism_detection_title": int(named_organism_detection_title),
        "clinical_metagenomics_title": int(clinical_metagenomics_title),
        "syndrome_direct_diagnostic_title": int(
            syndrome_direct_diagnostic_title
        ),
        "point_of_care_pathogen_testing_title": int(
            point_of_care_pathogen_testing_title
        ),
        "syndromic_diagnostic_implementation_title": int(
            syndromic_diagnostic_implementation_title
        ),
        "microbial_target_high": int(microbial_target_high),
        "organ_injury_target_high": int(organ_injury_target_high),
        "mixed_viral_bacterial_panel": int(mixed_viral_bacterial_panel),
        "host_mixed_bacterial_viral_test": int(
            host_mixed_bacterial_viral_test
        ),
        "febridx_host_panel": int(febridx_host_panel),
        "mxa_crp_host_panel": int(mxa_crp_host_panel),
        "mixed_direct_pathogen_panel": int(mixed_direct_pathogen_panel),
        "host_immune_analyte_high": int(host_immune_analyte_high),
        "direct_microbial_analyte_high": int(
            direct_microbial_analyte_high
        ),
        "host_response_diagnostic_intent": int(host_response_diagnostic_intent),
        "clinical_syndromic_diagnostic_intent": int(
            clinical_syndromic_diagnostic_intent
        ),
        "imaging_syndromic_diagnostic_intent": int(
            imaging_syndromic_diagnostic_intent
        ),
        "structured_diagnostic": int(structured_diagnostic),
        "diagnostic_evaluation_high": int(diagnostic_eval_high),
        "reference_standard_high": int(reference_high),
        "antibiotic_action_high": int(action_high),
        "explicit_amr_high": int(explicit_amr_high),
        "general_amr_high": int(general_amr_high),
        "ast_high": int(ast_high),
        "esbl_cre_high": int(esbl_high),
        "resistance_marker_high": int(marker_high or marker_generic_high),
        "resistance_prediction_high": int(prediction_high),
        "strict_mechanism": int(strict_mechanism),
        "oncology_high": int(oncology_high),
        "oncology_background_only": int(oncology_background_only),
        "cardiometabolic_high": int(cardiometabolic_high),
        "microbiome_high": int(microbiome_high),
        "imaging_only_high": int(imaging_only_high),
        "development_control": int(nct_id in DEVELOPMENT_CONTROL_IDS),
        "organism_categories": "|".join(
            categories_in_fields(organism_fields, HIGH_VALUE_FIELDS)
        ),
        "special_pathogen_categories": "|".join(
            categories_in_fields(special_fields, HIGH_VALUE_FIELDS)
        ),
        "syndrome_categories": "|".join(
            categories_in_fields(syndrome_fields, HIGH_VALUE_FIELDS)
        ),
        "diagnostic_modalities": "|".join(
            categories_in_fields(modality_fields, HIGH_VALUE_FIELDS)
        ),
        "mechanism_categories": "|".join(
            sorted(key for key, value in mechanism_flags.items() if value)
        ),
        "organism_evidence_fields_json": json.dumps(
            organism_fields, ensure_ascii=False, separators=(",", ":")
        ),
        "organism_evidence_terms_json": json.dumps(
            organism_terms, ensure_ascii=False, separators=(",", ":")
        ),
        "special_evidence_fields_json": json.dumps(
            special_fields, ensure_ascii=False, separators=(",", ":")
        ),
        "syndrome_evidence_fields_json": json.dumps(
            syndrome_fields, ensure_ascii=False, separators=(",", ":")
        ),
        "modality_evidence_fields_json": json.dumps(
            modality_fields, ensure_ascii=False, separators=(",", ":")
        ),
        "modality_evidence_terms_json": json.dumps(
            modality_terms, ensure_ascii=False, separators=(",", ":")
        ),
        "mechanism_evidence_fields_json": json.dumps(
            mechanism_fields, ensure_ascii=False, separators=(",", ":")
        ),
        "mechanism_evidence_snippets_json": json.dumps(
            mechanism_snips, ensure_ascii=False, separators=(",", ":")
        ),
        "diagnostic_evidence_fields": "|".join(diagnostic_eval_fields),
        "diagnostic_evidence_terms": "|".join(diagnostic_eval_terms),
        "diagnostic_evidence_snippets": " || ".join(diagnostic_eval_snips),
        "diagnostic_role_fields": "|".join(diagnostic_role_fields),
        "diagnostic_role_terms": "|".join(diagnostic_role_terms),
        "diagnostic_role_snippets": " || ".join(diagnostic_role_snips),
        "microbial_focus_fields": "|".join(microbial_focus_fields),
        "microbial_focus_terms": "|".join(microbial_focus_terms),
        "microbial_focus_snippets": " || ".join(microbial_focus_snips),
        "microbial_target_fields": "|".join(microbial_target_fields),
        "microbial_target_terms": "|".join(microbial_target_terms),
        "microbial_target_snippets": " || ".join(microbial_target_snips),
        "organ_injury_target_fields": "|".join(organ_injury_fields),
        "organ_injury_target_terms": "|".join(organ_injury_terms),
        "organ_injury_target_snippets": " || ".join(organ_injury_snips),
        "mixed_panel_fields": "|".join(mixed_panel_fields),
        "mixed_panel_terms": "|".join(mixed_panel_terms),
        "mixed_panel_snippets": " || ".join(mixed_panel_snips),
        "diagnostic_title_focus": int(diagnostic_title_focus),
        "reference_standard_fields": "|".join(reference_fields),
        "reference_standard_terms": "|".join(reference_terms),
        "reference_standard_snippets": " || ".join(reference_snips),
        "antibiotic_action_fields": "|".join(action_fields),
        "antibiotic_action_terms": "|".join(action_terms),
        "antibiotic_action_snippets": " || ".join(action_snips),
        "general_amr_fields": "|".join(general_amr_fields),
        "general_amr_terms": "|".join(general_amr_terms),
        "general_amr_snippets": " || ".join(general_amr_snips),
        "ast_fields": "|".join(ast_fields),
        "ast_terms": "|".join(ast_terms),
        "ast_snippets": " || ".join(ast_snips),
        "esbl_cre_fields": "|".join(esbl_fields),
        "esbl_cre_terms": "|".join(esbl_terms),
        "esbl_cre_snippets": " || ".join(esbl_snips),
        "resistance_marker_fields": "|".join(marker_fields + marker_generic_fields),
        "resistance_marker_terms": "|".join(
            list(dict.fromkeys(marker_terms + marker_generic_terms))
        ),
        "resistance_marker_snippets": " || ".join(
            list(dict.fromkeys(marker_snips + marker_generic_snips))
        ),
        "resistance_prediction_fields": "|".join(prediction_fields),
        "resistance_prediction_terms": "|".join(prediction_terms),
        "resistance_prediction_snippets": " || ".join(prediction_snips),
        "host_response_fields": "|".join(host_fields),
        "host_response_terms": "|".join(host_terms),
        "host_response_snippets": " || ".join(host_snips),
        "host_immune_analyte_fields": "|".join(host_immune_fields),
        "host_immune_analyte_terms": "|".join(host_immune_terms),
        "host_immune_analyte_snippets": " || ".join(host_immune_snips),
        "direct_microbial_analyte_fields": "|".join(
            direct_microbial_analyte_fields
        ),
        "direct_microbial_analyte_terms": "|".join(
            direct_microbial_analyte_terms
        ),
        "direct_microbial_analyte_snippets": " || ".join(
            direct_microbial_analyte_snips
        ),
        "host_diagnostic_fields": "|".join(host_diagnostic_fields),
        "host_diagnostic_terms": "|".join(host_diagnostic_terms),
        "host_diagnostic_snippets": " || ".join(host_diagnostic_snips),
        "host_immune_diagnostic_fields": "|".join(
            host_immune_diagnostic_fields
        ),
        "host_immune_diagnostic_terms": "|".join(
            host_immune_diagnostic_terms
        ),
        "host_immune_diagnostic_snippets": " || ".join(
            host_immune_diagnostic_snips
        ),
        "prognostic_fields": "|".join(prognostic_fields),
        "prognostic_terms": "|".join(prognostic_terms),
        "prognostic_snippets": " || ".join(prognostic_snips),
        "diagnostic_development_fields": "|".join(development_fields),
        "diagnostic_development_terms": "|".join(development_terms),
        "diagnostic_development_snippets": " || ".join(development_snips),
        "general_infection_fields": "|".join(general_infection_fields),
        "general_infection_terms": "|".join(general_infection_terms),
        "general_infection_snippets": " || ".join(general_infection_snips),
        "surveillance_fields": "|".join(surveillance_fields),
        "surveillance_terms": "|".join(surveillance_terms),
        "therapeutic_fields": "|".join(therapeutic_fields),
        "therapeutic_terms": "|".join(therapeutic_terms),
        "conditions_keywords": fields["conditions_keywords"][:5000],
        "intervention_names": fields["intervention_names"][:5000],
        "intervention_descriptions": fields["intervention_descriptions"][:7000],
        "primary_outcomes": fields["primary_outcomes"][:7000],
        "secondary_outcomes": fields["secondary_outcomes"][:7000],
        "summary": fields["summary"][:9000],
        "clinicaltrials_url": f"https://clinicaltrials.gov/study/{nct_id}",
    }

    for category, value in organism_flags.items():
        row[f"org_{category}"] = value
    for category, value in special_flags.items():
        row[f"special_{category}"] = value
    for category, value in syndrome_flags.items():
        row[f"syndrome_{category}"] = value
    for category, value in modality_flags.items():
        row[f"diag_{category}"] = value
    for category, value in mechanism_flags.items():
        row[f"mech_{category}"] = value

    return row


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------


def iter_studies(parts: list[Path]) -> Iterator[dict[str, Any]]:
    for part_number, part_path in enumerate(parts, start=1):
        with gzip.open(part_path, "rt", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                study = json.loads(line)
                if not isinstance(study, dict):
                    raise ValueError(
                        f"Expected JSON object in {part_path}, line {line_number}"
                    )
                yield study
        if part_number == 1 or part_number % 25 == 0 or part_number == len(parts):
            print(f"  scanned shards {part_number:,}/{len(parts):,}", flush=True)


def fieldnames_for_rows(rows: list[dict[str, Any]]) -> list[str]:
    preferred = [
        "nct_id",
        "predicted_stratum",
        "classification_reason",
        "diagnostic_depth_level",
        "diagnostic_depth_label",
        "infection_score",
        "diagnostic_score",
        "amr_score",
        "near_miss_score",
        "brief_title",
        "official_title",
        "overall_status",
        "study_type",
        "primary_purpose",
        "has_diagnostic_test_intervention",
        "phases",
        "enrollment_count",
        "enrollment_type",
        "lead_sponsor_name",
        "lead_sponsor_class",
        "start_date",
        "completion_date",
        "study_first_post_date",
        "results_first_post_date",
        "has_results",
        "organism_categories",
        "special_pathogen_categories",
        "syndrome_categories",
        "diagnostic_modalities",
        "mechanism_categories",
        "clinicaltrials_url",
    ]
    all_fields = {key for row in rows for key in row.keys()}
    return preferred + sorted(all_fields - set(preferred))


def compact_audit_row(row: dict[str, Any], source: str) -> dict[str, Any]:
    return {
        "nct_id": row.get("nct_id", ""),
        "audit_source": source,
        "near_miss_score": row.get("near_miss_score", ""),
        "infection_score": row.get("infection_score", ""),
        "diagnostic_score": row.get("diagnostic_score", ""),
        "amr_score": row.get("amr_score", ""),
        "brief_title": row.get("brief_title", ""),
        "conditions_keywords": row.get("conditions_keywords", ""),
        "intervention_names": row.get("intervention_names", ""),
        "primary_outcomes": row.get("primary_outcomes", ""),
        "secondary_outcomes": row.get("secondary_outcomes", ""),
        "summary": row.get("summary", ""),
        "clinicaltrials_url": row.get("clinicaltrials_url", ""),
        "development_control": row.get("development_control", 0),
    }


def deterministic_sample(
    rows: list[dict[str, Any]],
    n: int,
    *,
    salt: str,
) -> list[dict[str, Any]]:
    eligible = [
        row
        for row in rows
        if row.get("nct_id") not in DEVELOPMENT_CONTROL_IDS
    ]
    return sorted(
        eligible,
        key=lambda row: stable_hash_int(f"{salt}|{row['nct_id']}"),
    )[: min(n, len(eligible))]


def create_validation_set(
    retained: list[dict[str, Any]],
    near_misses: list[dict[str, Any]],
    random_negatives: list[dict[str, Any]],
    *,
    seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    by_stratum: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in retained:
        by_stratum[str(row["predicted_stratum"])].append(row)

    selected: list[tuple[str, dict[str, Any]]] = []
    for stratum, target in VALIDATION_TARGETS.items():
        if stratum == "NEAR_MISS":
            source_rows = near_misses
        elif stratum == "RANDOM_REGISTRY_NEGATIVE":
            source_rows = random_negatives
        else:
            source_rows = by_stratum.get(stratum, [])
        chosen = deterministic_sample(
            source_rows,
            target,
            salt=f"v3_2_5-validation-{seed}-{stratum}",
        )
        selected.extend((stratum, row) for row in chosen)

    rng = random.Random(seed)
    rng.shuffle(selected)

    blinded_rows: list[dict[str, Any]] = []
    key_rows: list[dict[str, Any]] = []
    for index, (source_stratum, row) in enumerate(selected, start=1):
        validation_id = f"V3-{index:04d}"
        blinded_rows.append(
            {
                "validation_id": validation_id,
                "nct_id": row.get("nct_id", ""),
                "brief_title": row.get("brief_title", ""),
                "conditions_keywords": row.get("conditions_keywords", ""),
                "intervention_names": row.get("intervention_names", ""),
                "primary_outcomes": row.get("primary_outcomes", ""),
                "secondary_outcomes": row.get("secondary_outcomes", ""),
                "summary": row.get("summary", ""),
                "clinicaltrials_url": row.get("clinicaltrials_url", ""),
                "manual_primary_eligible": "",
                "manual_final_stratum": "",
                "manual_amr_depth": "",
                "manual_exclusion_reason": "",
                "manual_notes": "",
                "reviewer_1": "",
                "reviewer_2_primary_eligible": "",
                "reviewer_2_final_stratum": "",
                "reviewer_2_amr_depth": "",
                "reviewer_2_notes": "",
            }
        )
        key_rows.append(
            {
                "validation_id": validation_id,
                "nct_id": row.get("nct_id", ""),
                "validation_source": source_stratum,
                "predicted_stratum": row.get("predicted_stratum", "EXCLUDED"),
                "predicted_amr_depth": row.get("diagnostic_depth_level", ""),
                "near_miss_score": row.get("near_miss_score", ""),
                "development_control": row.get("development_control", 0),
            }
        )
    return blinded_rows, key_rows


def count_binary_columns(
    rows: list[dict[str, Any]],
    columns: list[str],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    denominator = len(rows)
    for column in columns:
        count = sum(int(row.get(column, 0) or 0) for row in rows)
        output.append(
            {
                "category": column,
                "count": count,
                "percent_of_retained": round(100 * count / denominator, 4)
                if denominator
                else 0,
            }
        )
    return sorted(output, key=lambda item: (-item["count"], item["category"]))


def make_inventory(root: Path) -> None:
    rows: list[dict[str, Any]] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name == "file_inventory.tsv":
            continue
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        rows.append(
            {
                "relative_path": path.relative_to(root).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": digest.hexdigest(),
            }
        )
    write_delimited(
        root / "file_inventory.tsv",
        rows,
        ["relative_path", "size_bytes", "sha256"],
        "\t",
    )


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------


def run_self_tests() -> None:
    positives = [
        "blaNDM-1",
        "blaKPC",
        "blaCTX-M-15",
        "blaOXA-48",
        "mecA",
        "vanA",
        "mcr-1",
    ]
    negatives = [
        "bladder cancer",
        "blast count",
        "black participants",
        "Bland-Altman analysis",
        "blinatumomab",
        "blastocarb",
    ]
    for value in positives:
        if not RESISTANCE_MARKER_PATTERN.search(value):
            raise AssertionError(f"Expected resistance marker match: {value}")
    for value in negatives:
        if RESISTANCE_MARKER_PATTERN.search(value):
            raise AssertionError(f"Unexpected resistance marker match: {value}")

    if not near_patterns(
        "blaCTX-M copy number was quantified",
        COPY_NUMBER_PATTERN,
        BETA_LACTAMASE_GENE_PATTERN,
        window=180,
    ):
        raise AssertionError("Expected copy-number mechanism match")

    if near_patterns(
        "PCR amplification was performed for bladder tissue",
        COPY_NUMBER_PATTERN,
        BETA_LACTAMASE_GENE_PATTERN,
        window=180,
    ):
        raise AssertionError("Unexpected copy-number mechanism match")

    if not HOST_IMMUNE_ANALYTE_PATTERN.search(
        "Borrelia-specific antibody-secreting cell ELISpot"
    ):
        raise AssertionError("Expected host immune analyte match")
    if not HOST_IMMUNE_DIAGNOSTIC_PATTERN.search(
        "Borrelia B-cell Diagnostics using antibody-secreting cells"
    ):
        raise AssertionError("Expected host immune diagnostic match")
    if not FEBRIDX_PATTERN.search("FebriDx Pediatric Validation Study"):
        raise AssertionError("Expected FebriDx match")
    if not MXA_PATTERN.search("Myxovirus resistance protein A (MxA)"):
        raise AssertionError("Expected MxA match")
    if not DIRECT_MICROBIAL_ANALYTE_PATTERN.search(
        "pneumococcal antigen in urine"
    ):
        raise AssertionError("Expected direct microbial antigen match")
    if HOST_IMMUNE_ANALYTE_PATTERN.search(
        "direct pneumococcal antigen detection"
    ):
        raise AssertionError("Direct microbial antigen misread as host response")

    title_positives = [
        "Study of the Performance of the KeyPath MRSA/MSSA Blood Culture Test",
        "Impact MRSA-PCR on Patient Management",
        "Rapid Carbapenemase Test in Blood Cultures",
    ]
    title_negatives = [
        "Validation of a Urinary Biomarker for AKI in Sepsis",
        "MRI Diagnosis of Cardiac Injury",
    ]
    for value in title_positives:
        if not NAMED_BACTERIAL_ASSAY_TITLE_PATTERN.search(value):
            raise AssertionError(f"Expected named bacterial assay title: {value}")
    for value in title_negatives:
        if NAMED_BACTERIAL_ASSAY_TITLE_PATTERN.search(value):
            raise AssertionError(f"Unexpected named bacterial assay title: {value}")

    if not MICROBIAL_TARGET_PATTERN.search("MRSA blood culture test"):
        raise AssertionError("Expected direct microbial target")
    if MICROBIAL_TARGET_PATTERN.search("urinary biomarker for acute kidney injury"):
        raise AssertionError("Unexpected microbial target in organ-injury biomarker")
    if not ORGAN_INJURY_PROGNOSTIC_TARGET_PATTERN.search(
        "acute kidney injury in sepsis"
    ):
        raise AssertionError("Expected organ-injury target")
    if not MIXED_VIRAL_BACTERIAL_PANEL_PATTERN.search(
        "point of care respiratory pathogen testing"
    ):
        raise AssertionError("Expected mixed respiratory panel flag")

    if not NAMED_ORGANISM_DETECTION_WORDING_PATTERN.search(
        "Rapid Detection of Group B Strep"
    ):
        raise AssertionError("Expected named-organism detection wording")
    if NAMED_ORGANISM_DETECTION_WORDING_PATTERN.search(
        "Thermal Images in Bacterial Pneumonia"
    ):
        raise AssertionError("Unexpected named-organism detection wording")
    if not CLINICAL_METAGENOMICS_TITLE_PATTERN.search(
        "Clinical Metagenomics of Infective Endocarditis"
    ):
        raise AssertionError("Expected clinical metagenomics title")
    if CLINICAL_METAGENOMICS_TITLE_PATTERN.search(
        "Whole Genome Sequencing to Determine Transmission Rate"
    ):
        raise AssertionError("Unexpected clinical metagenomics title")
    if not SYNDROME_DIRECT_DIAGNOSTIC_TITLE_PATTERN.search(
        "Study for Rapid Diagnosis of Postoperative Abdominal Infection"
    ):
        raise AssertionError("Expected syndrome direct diagnostic title")
    if SYNDROME_DIRECT_DIAGNOSTIC_TITLE_PATTERN.search(
        "Improving Diagnosis and Management of Suspected Brain Infections"
    ):
        raise AssertionError("Unexpected direct-assay title for implementation study")
    if not POINT_OF_CARE_PATHOGEN_TESTING_TITLE_PATTERN.search(
        "Point-of-care Testing of Respiratory Pathogens at Pediatric Emergency Room"
    ):
        raise AssertionError("Expected point-of-care pathogen testing title")
    if not POINT_OF_CARE_PATHOGEN_TESTING_TITLE_PATTERN.search(
        "Point\u2011of\u2011care Testing of Respiratory Pathogens at Pediatric Emergency Room"
    ):
        raise AssertionError("Expected Unicode-hyphen point-of-care title")
    if clean_text("Point\u2011of\u2011care") != "Point-of-care":
        raise AssertionError("Unicode hyphen normalization failed")
    if not SYNDROMIC_DIAGNOSTIC_IMPLEMENTATION_TITLE_PATTERN.search(
        "Improving Diagnosis and Management of Suspected Brain Infections Globally"
    ):
        raise AssertionError("Expected syndromic diagnostic implementation title")

    def synthetic_study(
        nct_id: str,
        title: str,
        conditions: list[str],
        intervention_name: str,
        intervention_type: str,
        outcome: str,
        purpose: str = "DIAGNOSTIC",
        summary: str = "",
        study_type: str = "INTERVENTIONAL",
    ) -> dict[str, Any]:
        return {
            "protocolSection": {
                "identificationModule": {
                    "nctId": nct_id,
                    "briefTitle": title,
                },
                "conditionsModule": {
                    "conditions": conditions,
                    "keywords": [],
                },
                "designModule": {
                    "studyType": study_type,
                    "designInfo": {"primaryPurpose": purpose},
                },
                "armsInterventionsModule": {
                    "interventions": [
                        {
                            "type": intervention_type,
                            "name": intervention_name,
                            "description": summary,
                        }
                    ]
                },
                "outcomesModule": {
                    "primaryOutcomes": [{"measure": outcome}]
                },
                "descriptionModule": {"briefSummary": summary},
            },
            "hasResults": False,
        }

    regression_cases = [
        (
            synthetic_study(
                "NCTSYN001",
                "Performance of the KeyPath MRSA/MSSA Blood Culture Test",
                ["Staphylococcus aureus bacteremia"],
                "KeyPath MRSA/MSSA blood culture test",
                "DIAGNOSTIC_TEST",
                "Sensitivity and specificity for MRSA and MSSA",
            ),
            "CORE_AMR_DIAGNOSTIC",
            1,
        ),
        (
            synthetic_study(
                "NCTSYN002",
                "Impact MRSA-PCR on Patient Management",
                ["Staphylococcus aureus bacteremia"],
                "MRSA PCR",
                "DIAGNOSTIC_TEST",
                "Time to targeted antibiotic therapy",
                purpose="TREATMENT",
            ),
            "CORE_AMR_DIAGNOSTIC",
            1,
        ),
        (
            synthetic_study(
                "NCTSYN003",
                "Validation of a Urinary Biomarker for AKI in Sepsis",
                ["Sepsis"],
                "Mass spectrometry urinary AKI biomarker",
                "DIAGNOSTIC_TEST",
                "Diagnostic accuracy for acute kidney injury",
            ),
            "CLINICAL_SYNDROMIC_SUPPORT",
            "",
        ),
        (
            synthetic_study(
                "NCTSYN004",
                "Point of Care Respiratory Pathogen Testing for Antibiotic Stewardship",
                ["Upper respiratory infection"],
                "BioFire SpotFire respiratory pathogen panel",
                "DIAGNOSTIC_TEST",
                "Antibiotic prescribing",
            ),
            "BROAD_BACTERIAL_DIAGNOSTIC",
            0,
        ),
        (
            synthetic_study(
                "NCTSYN005",
                "Rapid Detection of Group B Strep",
                ["Group B Streptococcus"],
                "Rapid Group B Streptococcus test",
                "DIAGNOSTIC_TEST",
                "Detection of Group B Streptococcus colonization",
            ),
            "BROAD_BACTERIAL_DIAGNOSTIC",
            0,
        ),
        (
            synthetic_study(
                "NCTSYN006",
                "Clinical Metagenomics of Infective Endocarditis",
                ["Infective Endocarditis"],
                "Metagenomic sequencing",
                "DIAGNOSTIC_TEST",
                "Diagnostic yield of metagenomic sequencing",
            ),
            "BROAD_BACTERIAL_DIAGNOSTIC",
            0,
        ),
        (
            synthetic_study(
                "NCTSYN007",
                "Thermal Images on Smartphones to Diagnose Bacterial Pneumonia",
                ["Bacterial Pneumonia"],
                "Smartphone thermal imaging",
                "DIAGNOSTIC_TEST",
                "Diagnostic accuracy for bacterial pneumonia",
            ),
            "CLINICAL_SYNDROMIC_SUPPORT",
            "",
        ),
        (
            synthetic_study(
                "NCTSYN008",
                "Alpha-defensin as a Diagnostic Means to Distinguish Between Acute Bacterial and Viral Infections",
                ["Bacterial Infections", "Viral Infection"],
                "Alpha-defensin assay",
                "DIAGNOSTIC_TEST",
                "Diagnostic accuracy for bacterial versus viral infection",
            ),
            "HOST_RESPONSE_DIAGNOSTIC",
            "",
        ),
        (
            synthetic_study(
                "NCTSYN009",
                "Study for Rapid Diagnosis of Postoperative Abdominal Infection",
                ["Postoperative Abdominal Infection"],
                "Metagenomic next-generation sequencing of drainage fluid",
                "DIAGNOSTIC_TEST",
                "Diagnostic accuracy compared with conventional culture",
            ),
            "BROAD_BACTERIAL_DIAGNOSTIC",
            0,
        ),
        (
            synthetic_study(
                "NCTSYN010",
                "Point-of-care Testing of Respiratory Pathogens at Pediatric Emergency Room",
                ["Respiratory Tract Infections"],
                "QIAstat-Dx respiratory SARS-CoV-2 panel",
                "DIAGNOSTIC_TEST",
                "Antibiotic consumption and hospital admission",
                summary=(
                    "Point-of-care multiplex testing for influenza, RSV, "
                    "SARS-CoV-2, and bacterial respiratory pathogens to "
                    "support antimicrobial stewardship."
                ),
            ),
            "BROAD_BACTERIAL_DIAGNOSTIC",
            0,
        ),
        (
            synthetic_study(
                "NCTSYN011",
                "BIGlobal Intervention Study: Improving Diagnosis and Management of Suspected Brain Infections Globally",
                ["Meningitis", "Encephalitis", "Brain Abscess"],
                "Pragmatic multi-component package",
                "OTHER",
                "Percentage achieving a microbiological diagnosis",
                purpose="",
                study_type="OBSERVATIONAL",
            ),
            "CLINICAL_SYNDROMIC_SUPPORT",
            "",
        ),
        (
            synthetic_study(
                "NCTSYN012",
                "Presepsin:Gelsolin Ratio in Sepsis-related Organ Dysfunction",
                ["Sepsis", "Organ Dysfunction", "Prognosis"],
                "Presepsin and gelsolin biomarker measurements",
                "DIAGNOSTIC_TEST",
                "Prediction of mortality and organ failure",
                purpose="",
                study_type="OBSERVATIONAL",
            ),
            "CLINICAL_SYNDROMIC_SUPPORT",
            "",
        ),
        (
            synthetic_study(
                "NCTSYN013",
                "New Dosages of Inflammatory Markers for the Early Diagnosis of Nosocomial Bacterial Infections of the Newborn",
                ["Nosocomial Bacterial Infection", "Late Onset Neonatal Sepsis"],
                "Inflammatory markers including interleukins and CRP",
                "DIAGNOSTIC_TEST",
                "Diagnosis of nosocomial bacterial infection",
            ),
            "HOST_RESPONSE_DIAGNOSTIC",
            "",
        ),
        (
            synthetic_study(
                "NCTSYN014",
                "Comparison of Pulsed-field Gel Electrophoresis and Whole Genome Sequencing to Determine Transmission Rate of ESBL-producing E. coli",
                ["ESBL-producing Escherichia coli"],
                "Whole genome sequencing and PFGE",
                "DIAGNOSTIC_TEST",
                "Transmission rate and clonal relatedness",
                purpose="",
                study_type="OBSERVATIONAL",
            ),
            "SURVEILLANCE_SUPPORT",
            "",
        ),
        (
            synthetic_study(
                "NCTSYN015",
                "Borrelia B-cell Diagnostics",
                ["Lyme Disease", "Borrelia burgdorferi Infection"],
                "Borrelia-specific antibody-secreting cell ELISpot assay",
                "DIAGNOSTIC_TEST",
                "Sensitivity and specificity for early Lyme disease",
                summary=(
                    "The test measures the patient's Borrelia-specific "
                    "B-cell and antibody-secreting-cell response."
                ),
            ),
            "HOST_RESPONSE_DIAGNOSTIC",
            "",
        ),
        (
            synthetic_study(
                "NCTSYN016",
                "FebriDx Pediatric Validation Study",
                ["Acute Respiratory Infections"],
                "FebriDx point-of-care host immune response assay",
                "DIAGNOSTIC_TEST",
                "Diagnostic accuracy for bacterial versus viral infection",
                summary=(
                    "The host-response test measures myxovirus resistance "
                    "protein A (MxA) and C-reactive protein (CRP)."
                ),
            ),
            "HOST_RESPONSE_DIAGNOSTIC",
            "",
        ),
        (
            synthetic_study(
                "NCTSYN017",
                "Rapid Pneumococcal Antigen Detection in Pneumonia",
                ["Bacterial Pneumonia", "Streptococcus pneumoniae"],
                "Urinary pneumococcal antigen immunoassay",
                "DIAGNOSTIC_TEST",
                "Sensitivity and specificity for pneumococcal antigen",
            ),
            "BROAD_BACTERIAL_DIAGNOSTIC",
            0,
        ),
    ]
    for study, expected_stratum, expected_depth in regression_cases:
        row = evaluate_study(study)
        if row.get("predicted_stratum") != expected_stratum:
            raise AssertionError(
                f"Synthetic regression failed: expected {expected_stratum}, "
                f"observed {row.get('predicted_stratum')} for "
                f"{row.get('brief_title')}"
            )
        if row.get("diagnostic_depth_level") != expected_depth:
            raise AssertionError(
                f"Synthetic depth failed: expected {expected_depth}, "
                f"observed {row.get('diagnostic_depth_level')} for "
                f"{row.get('brief_title')}"
            )
    for index in (3, 7, 9):
        mixed_row = evaluate_study(regression_cases[index][0])
        if int(mixed_row.get("mixed_viral_bacterial_panel", 0) or 0) != 1:
            raise AssertionError(
                f"Expected mixed viral/bacterial panel flag for "
                f"{mixed_row.get('brief_title')}"
            )

    febridx_row = evaluate_study(regression_cases[15][0])
    if int(febridx_row.get("mixed_viral_bacterial_panel", 0) or 0) != 1:
        raise AssertionError("Expected FebriDx mixed viral/bacterial flag")
    if int(febridx_row.get("host_mixed_bacterial_viral_test", 0) or 0) != 1:
        raise AssertionError("Expected FebriDx host mixed-test flag")

    borrelia_host_row = evaluate_study(regression_cases[14][0])
    if int(borrelia_host_row.get("host_immune_analyte_high", 0) or 0) != 1:
        raise AssertionError("Expected Borrelia host immune analyte flag")

    antigen_row = evaluate_study(regression_cases[16][0])
    if int(antigen_row.get("direct_microbial_analyte_high", 0) or 0) != 1:
        raise AssertionError("Expected direct microbial analyte flag")
    if int(antigen_row.get("host_response_diagnostic_intent", 0) or 0) != 0:
        raise AssertionError("Direct bacterial antigen misrouted to host response")

    poc_mixed_row = evaluate_study(regression_cases[9][0])
    if int(poc_mixed_row.get("mixed_direct_pathogen_panel", 0) or 0) != 1:
        raise AssertionError(
            "Expected mixed direct pathogen panel to bypass nonbacterial-only guard"
        )

    print("Version 3.2.5 self-tests: PASS")


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Staged ClinicalTrials.gov bacterial AMR diagnostic screen v3.2.5."
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-parts", type=int, default=None)
    parser.add_argument("--near-miss-pool-size", type=int, default=1500)
    parser.add_argument("--random-negative-pool-size", type=int, default=1000)
    parser.add_argument("--validation-seed", type=int, default=20260712)
    parser.add_argument("--no-zip", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run_self_tests()
    if args.self_test:
        return 0

    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    parts = sorted((input_dir / "parts").glob("part-*.jsonl.gz"))
    if not parts:
        raise SystemExit(f"No part-*.jsonl.gz files found under {input_dir / 'parts'}")
    if args.max_parts is not None:
        parts = parts[: args.max_parts]

    output_dir.mkdir(parents=True, exist_ok=True)

    retained: list[dict[str, Any]] = []
    near_pool = TopKRows(args.near_miss_pool_size)
    random_negative_pool = LowestHashRows(args.random_negative_pool_size)
    registry_stratum_counts: Counter[str] = Counter()
    studies_scanned = 0

    print(f"Scanning {len(parts):,} ClinicalTrials.gov shard(s)...", flush=True)
    for study in iter_studies(parts):
        studies_scanned += 1
        row = evaluate_study(study)
        nct_id = str(row.get("nct_id", ""))
        if not nct_id:
            continue
        stratum = str(row["predicted_stratum"])
        registry_stratum_counts[stratum] += 1
        if stratum in PREDICTED_STRATA:
            retained.append(row)
        else:
            compact = compact_audit_row(row, "NEAR_MISS")
            if int(row.get("near_miss_score", 0) or 0) >= 8:
                near_pool.add(
                    float(row.get("near_miss_score", 0) or 0),
                    nct_id,
                    compact,
                )
            random_compact = compact_audit_row(row, "RANDOM_REGISTRY_NEGATIVE")
            random_negative_pool.add(nct_id, random_compact)

    retained.sort(
        key=lambda row: (
            str(row["predicted_stratum"]),
            -int(row.get("diagnostic_score", 0) or 0),
            str(row["nct_id"]),
        )
    )
    near_misses = near_pool.rows()
    random_negatives = random_negative_pool.rows()

    if not retained:
        raise SystemExit("No retained records were generated.")

    retained_fields = fieldnames_for_rows(retained)
    write_csv_and_tsv(
        output_dir,
        "retained_all_v3_2_5",
        retained,
        retained_fields,
    )

    for stratum in PREDICTED_STRATA:
        rows = [row for row in retained if row["predicted_stratum"] == stratum]
        write_csv_and_tsv(
            output_dir,
            stratum.lower(),
            rows,
            retained_fields,
        )

    audit_fields = [
        "nct_id",
        "audit_source",
        "near_miss_score",
        "infection_score",
        "diagnostic_score",
        "amr_score",
        "brief_title",
        "conditions_keywords",
        "intervention_names",
        "primary_outcomes",
        "secondary_outcomes",
        "summary",
        "clinicaltrials_url",
        "development_control",
    ]
    write_csv_and_tsv(
        output_dir,
        "near_miss_audit_pool_v3_2_5",
        near_misses,
        audit_fields,
    )
    write_csv_and_tsv(
        output_dir,
        "random_registry_negative_pool_v3_2_5",
        random_negatives,
        audit_fields,
    )

    validation_blinded, validation_key = create_validation_set(
        retained,
        near_misses,
        random_negatives,
        seed=args.validation_seed,
    )
    blinded_fields = list(validation_blinded[0].keys()) if validation_blinded else []
    key_fields = list(validation_key[0].keys()) if validation_key else []
    write_csv_and_tsv(
        output_dir,
        "validation_set_v3_2_5_blinded",
        validation_blinded,
        blinded_fields,
    )
    write_csv_and_tsv(
        output_dir,
        "validation_key_v3_2_5_keep_blinded",
        validation_key,
        key_fields,
    )

    tier_counts = Counter(row["predicted_stratum"] for row in retained)
    depth_counts = Counter(
        (
            "NA"
            if row.get("diagnostic_depth_level", "") in {"", None}
            else str(row.get("diagnostic_depth_level"))
        )
        for row in retained
        if row["predicted_stratum"] in PRIMARY_DIAGNOSTIC_STRATA
    )

    write_csv_and_tsv(
        output_dir,
        "counts_by_stratum",
        [
            {
                "predicted_stratum": key,
                "count": value,
                "percent_of_retained": round(100 * value / len(retained), 4),
            }
            for key, value in tier_counts.most_common()
        ],
        ["predicted_stratum", "count", "percent_of_retained"],
    )
    write_csv_and_tsv(
        output_dir,
        "counts_by_diagnostic_depth",
        [
            {
                "diagnostic_depth_level": key,
                "diagnostic_depth_label": DEPTH_LABELS.get(int(key), "NA")
                if key.isdigit()
                else "NA",
                "count": value,
            }
            for key, value in sorted(depth_counts.items())
        ],
        ["diagnostic_depth_level", "diagnostic_depth_label", "count"],
    )

    organism_columns = [f"org_{key}" for key in ORGANISM_PATTERNS]
    special_columns = [f"special_{key}" for key in SPECIAL_PATHOGEN_PATTERNS]
    syndrome_columns = [f"syndrome_{key}" for key in BACTERIAL_SYNDROME_PATTERNS]
    modality_columns = [f"diag_{key}" for key in DIRECT_DIAGNOSTIC_PATTERNS]
    mechanism_columns = [
        "mech_copy_number",
        "mech_gene_expression",
        "mech_porin",
        "mech_efflux",
        "mech_enzyme_activity",
        "mech_protein_abundance",
        "mech_multimechanism",
    ]

    for stem, columns in [
        ("counts_by_organism", organism_columns),
        ("counts_by_special_pathogen", special_columns),
        ("counts_by_syndrome", syndrome_columns),
        ("counts_by_modality", modality_columns),
        ("counts_by_mechanism", mechanism_columns),
        ("counts_by_panel_type", ["mixed_viral_bacterial_panel"]),
    ]:
        rows = count_binary_columns(retained, columns)
        write_csv_and_tsv(
            output_dir,
            stem,
            rows,
            ["category", "count", "percent_of_retained"],
        )

    for stem, field in [
        ("counts_by_status", "overall_status"),
        ("counts_by_study_type", "study_type"),
        ("counts_by_sponsor_class", "lead_sponsor_class"),
    ]:
        counter = Counter(str(row.get(field, "") or "MISSING") for row in retained)
        rows = [
            {
                field: value,
                "count": count,
                "percent_of_retained": round(100 * count / len(retained), 4),
            }
            for value, count in counter.most_common()
        ]
        write_csv_and_tsv(
            output_dir,
            stem,
            rows,
            [field, "count", "percent_of_retained"],
        )

    summary = {
        "version": VERSION,
        "created_at": utc_now(),
        "input_directory": str(input_dir),
        "output_directory": str(output_dir),
        "shards_scanned": len(parts),
        "studies_scanned": studies_scanned,
        "retained_count": len(retained),
        "retained_percent_of_registry": round(
            100 * len(retained) / studies_scanned, 5
        ),
        "registry_stratum_counts_including_excluded": dict(registry_stratum_counts),
        "retained_stratum_counts": dict(tier_counts),
        "primary_diagnostic_count": sum(
            tier_counts.get(key, 0) for key in PRIMARY_DIAGNOSTIC_STRATA
        ),
        "near_miss_pool_count": len(near_misses),
        "random_negative_pool_count": len(random_negatives),
        "validation_set_count": len(validation_blinded),
        "validation_targets": VALIDATION_TARGETS,
        "development_controls_excluded_from_validation": sorted(
            DEVELOPMENT_CONTROL_IDS
        ),
        "diagnostic_depth_counts": dict(depth_counts),
        "strict_copy_number_count": sum(
            int(row.get("mech_copy_number", 0) or 0) for row in retained
        ),
        "strict_gene_expression_count": sum(
            int(row.get("mech_gene_expression", 0) or 0) for row in retained
        ),
        "strict_porin_count": sum(
            int(row.get("mech_porin", 0) or 0) for row in retained
        ),
        "strict_enzyme_activity_count": sum(
            int(row.get("mech_enzyme_activity", 0) or 0) for row in retained
        ),
        "mixed_viral_bacterial_panel_count": sum(
            int(row.get("mixed_viral_bacterial_panel", 0) or 0)
            for row in retained
        ),
        "notes": [
            "This is a screened cohort, not a final analytic cohort.",
            "Do not open validation_key_v3_2_5_keep_blinded.tsv while adjudicating.",
            "All 31 development/regression controls were excluded from the held-out validation set.",
            "Run score_bacterial_amr_v3_2_5_validation.py after manual adjudication.",
        ],
    }
    atomic_write_json(output_dir / "summary_v3_2_5.json", summary)

    report = f"""# ClinicalTrials.gov bacterial AMR diagnostic screen v3.2.5

- Registry studies scanned: **{studies_scanned:,}**
- Retained records across all strata: **{len(retained):,}**
- Primary direct diagnostic records: **{summary['primary_diagnostic_count']:,}**
- Core AMR diagnostic: **{tier_counts.get('CORE_AMR_DIAGNOSTIC', 0):,}**
- Broad bacterial diagnostic: **{tier_counts.get('BROAD_BACTERIAL_DIAGNOSTIC', 0):,}**
- Host-response diagnostic: **{tier_counts.get('HOST_RESPONSE_DIAGNOSTIC', 0):,}**
- Clinical/syndromic diagnostic support: **{tier_counts.get('CLINICAL_SYNDROMIC_SUPPORT', 0):,}**
- Mechanism support: **{tier_counts.get('MECHANISM_SUPPORT', 0):,}**
- Special-pathogen diagnostic: **{tier_counts.get('SPECIAL_PATHOGEN_DIAGNOSTIC', 0):,}**
- Surveillance support: **{tier_counts.get('SURVEILLANCE_SUPPORT', 0):,}**
- Therapeutic support: **{tier_counts.get('THERAPEUTIC_SUPPORT', 0):,}**
- Mixed viral/bacterial panel records: **{summary['mixed_viral_bacterial_panel_count']:,}**
- Held-out blinded validation records: **{len(validation_blinded):,}**

## Required next step

Manually adjudicate `validation_set_v3_2_5_blinded.tsv` without opening the key.
Then run `score_bacterial_amr_v3_2_5_validation.py`. The retained cohort must not be
used for prevalence estimates, organism comparisons, or novelty conclusions
until the independent validation gate passes.
"""
    (output_dir / "screen_v3_2_5_report.md").write_text(report, encoding="utf-8")

    validation_readme = """# V3.1 blinded validation instructions

Do not open `validation_key_v3_2_5_keep_blinded.tsv` before adjudication.

Complete the following fields in `validation_set_v3_2_5_blinded.tsv` or `.csv`:

- manual_primary_eligible: YES / NO / UNCERTAIN
- manual_final_stratum: CORE_AMR_DIAGNOSTIC / BROAD_BACTERIAL_DIAGNOSTIC /
  HOST_RESPONSE_DIAGNOSTIC / CLINICAL_SYNDROMIC_SUPPORT / MECHANISM_SUPPORT /
  SPECIAL_PATHOGEN_DIAGNOSTIC / SURVEILLANCE_SUPPORT /
  THERAPEUTIC_SUPPORT / NONINFECTIOUS_OR_UNRELATED / OTHER
- manual_amr_depth: 0 / 1 / 2 / 3 / 4 / NA
- manual_exclusion_reason
- manual_notes
- reviewer_1

Primary eligibility means the record is an evaluated direct bacterial diagnostic
appropriate for the primary landscape analysis. Host-response, mechanism,
surveillance, therapeutic, and special-pathogen records are not primary even
when scientifically relevant.

After adjudication, run `score_bacterial_amr_v3_2_5_validation.py` using the edited
blinded file and the untouched key.
"""
    (output_dir / "README_VALIDATION_V3_2_3.md").write_text(
        validation_readme, encoding="utf-8"
    )

    pattern_manifest = {
        "version": VERSION,
        "curated_beta_lactamase_gene_regex": BETA_LACTAMASE_GENE_PATTERN.pattern,
        "other_resistance_gene_regex": OTHER_RESISTANCE_GENE_PATTERN.pattern,
        "named_bacterial_assay_title_regex": (
            NAMED_BACTERIAL_ASSAY_TITLE_PATTERN.pattern
        ),
        "syndrome_direct_diagnostic_title_regex": (
            SYNDROME_DIRECT_DIAGNOSTIC_TITLE_PATTERN.pattern
        ),
        "point_of_care_pathogen_testing_title_regex": (
            POINT_OF_CARE_PATHOGEN_TESTING_TITLE_PATTERN.pattern
        ),
        "syndromic_diagnostic_implementation_title_regex": (
            SYNDROMIC_DIAGNOSTIC_IMPLEMENTATION_TITLE_PATTERN.pattern
        ),
        "host_immune_analyte_regex": HOST_IMMUNE_ANALYTE_PATTERN.pattern,
        "host_immune_diagnostic_regex": HOST_IMMUNE_DIAGNOSTIC_PATTERN.pattern,
        "direct_microbial_analyte_regex": DIRECT_MICROBIAL_ANALYTE_PATTERN.pattern,
        "febridx_regex": FEBRIDX_PATTERN.pattern,
        "mxa_regex": MXA_PATTERN.pattern,
        "crp_regex": CRP_PATTERN.pattern,
        "microbial_target_regex": MICROBIAL_TARGET_PATTERN.pattern,
        "organ_injury_prognostic_target_regex": (
            ORGAN_INJURY_PROGNOSTIC_TARGET_PATTERN.pattern
        ),
        "mixed_viral_bacterial_panel_regex": (
            MIXED_VIRAL_BACTERIAL_PANEL_PATTERN.pattern
        ),
        "high_value_fields": sorted(HIGH_VALUE_FIELDS),
        "primary_relevance_fields": sorted(PRIMARY_RELEVANCE_FIELDS),
        "predicted_strata": PREDICTED_STRATA,
        "diagnostic_depth_labels": DEPTH_LABELS,
    }
    atomic_write_json(output_dir / "pattern_manifest_v3_2_5.json", pattern_manifest)

    make_inventory(output_dir)

    if not args.no_zip:
        zip_path = shutil.make_archive(
            str(output_dir),
            "zip",
            root_dir=output_dir.parent,
            base_dir=output_dir.name,
        )
        print(f"ZIP archive: {zip_path}", flush=True)

    print(f"\nScreen v3.2.5 complete: {output_dir}", flush=True)
    print(f"Retained records: {len(retained):,}", flush=True)
    print(f"Validation records: {len(validation_blinded):,}", flush=True)
    print(f"Open first: {output_dir / 'summary_v3_2_5.json'}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        raise SystemExit(130)
