# Screening regression test data

The frozen v3.2.5 screening program contains an embedded synthetic regression
suite invoked with `--self-test`. The suite constructs synthetic
ClinicalTrials.gov-style records spanning primary bacterial diagnostics,
AMR-focused diagnostics, host-response diagnostics, surveillance,
clinical-syndromic support, and mixed viral-bacterial panels.

Because the synthetic fixtures are embedded in the authoritative frozen Python
source, they remain version-locked to the exact screening implementation.
Run them with

```bash
bash tests/test_screening_self_test.sh
```

The wrapper requires the exact message

`Version 3.2.5 self-tests: PASS`

and exits nonzero if the regression suite fails.
