# ClinicalTrials.gov bacterial antimicrobial resistance diagnostic landscape

This repository contains the executable code supporting the manuscript on
diagnostic depth and clinical utility in registered bacterial antimicrobial
resistance diagnostic studies.

## Scientific workflow

The repository preserves the code used for

1. deterministic ClinicalTrials.gov registry retrieval and prioritization
2. v3.2.9 coverage rescue preparation
3. descriptor-stage preparation and freeze validation
4. P35 imaging sensitivity finalization
5. construction of the frozen v3.3.0 analytic release
6. final v3.3.1 H1-H4 statistical analysis
7. comparison of the historical and final analytic results

Human screening and descriptor adjudication are represented by frozen input,
decision, and output datasets in the companion Zenodo data deposit rather than
being represented as automated classifications.

## Data source

The original source was the ClinicalTrials.gov API v2 snapshot with a data
timestamp of 2026-07-10T09:00:05 and 593,334 public study records.

The companion dataset contains derived screening and analytic products. It is
not a current mirror of ClinicalTrials.gov.

## Repository identifiers

GitHub: https://github.com/codyadamblack/clinicaltrials-bacterial-amr-diagnostics

Dataset DOI: https://doi.org/10.5281/zenodo.22046951

Software DOI: https://doi.org/10.5281/zenodo.22046953

## Environment

`environment.yml` provides the full conda environment without a machine-specific
prefix.

`environment_minimal.yml` records the explicitly requested conda dependencies.

`docs/conda_explicit_linux-64.txt` preserves the exact conda package build set
used in the analysis environment.

`requirements.txt` provides a portable package and version listing.

## Reproduction

See `docs/execution_order.md`.

A small public test fixture for the registry-screening stage will be added before
the versioned software release is frozen.

## License

MIT License for the authors' code.

Underlying ClinicalTrials.gov records remain subject to ClinicalTrials.gov terms
and conditions.
