# Reproducibility workflow

The workflow separates deterministic registry retrieval, human adjudication,
analytic-release construction, statistical analysis, and manuscript comparison.

## 1. Registry-wide deterministic retrieval and prioritization

Primary script

`screening/screen_bacterial_amr_diagnostics_v3_2_5.py`

Frozen pattern record

`screening/pattern_manifest_v3_2_5.json`

The script processes the complete ClinicalTrials.gov API v2 census and routes
records into diagnostic/support categories or provisional exclusion. It also
constructs the near-miss and deterministic registry-negative retrieval pools.

## 2. Coverage rescue and v3.2.9 screening closure

Supporting script

`screening/prepare_v3_2_9_confirmation_and_expansion.py`

Final eligibility decisions were established by human review and adjudication.
The frozen screening outputs are deposited with the companion Zenodo dataset.

## 3. Descriptor-stage preparation and adjudication

Scripts

`adjudication/prepare_v3_3_0_descriptor_stage.py`

`adjudication/prepare_v3_3_0_neutral_descriptor_adjudication.py`

`adjudication/validate_v3_3_0_descriptor_dispatch.py`

`adjudication/validate_and_freeze_neutral_descriptor_adjudication_v3_3_0.py`

Human adjudication occurs between preparation and final validation. The frozen
input and final output files are included in the data deposit.

## 4. P35 imaging sensitivity finalization

Script

`sensitivity/finalize_P35_full_cohort_sensitivity_v3_3_0.py`

## 5. Frozen analytic release

Script

`release/build_analytic_release_v3_3_0.py`

Primary release

`eligible_primary_cohort_final_v3_3_0.tsv`

Primary n = 573.

## 6. Final H1-H4 statistical analysis

Script

`analysis/run_refined_H1_H4_analysis_v3_3_1.py`

This script generates the final hypothesis tests, sensitivity analyses,
diagnostics, figures, and numerical result tables.

## 7. Historical versus final comparison

Script

`provenance/build_v3_2_8_v3_3_1_comparison_package.py`

This step documents how coverage rescue changed the analytic cohort and
manuscript claims without refitting the historical analysis.
