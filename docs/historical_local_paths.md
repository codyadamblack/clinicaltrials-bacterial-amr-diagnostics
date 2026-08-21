# Historical local path defaults

The frozen v3.2.5 registry-screening script is preserved byte-for-byte as used
for the study. Its `DEFAULT_INPUT` and `DEFAULT_OUTPUT` values therefore retain
the original WSL development paths under `/mnt/d/clinicaltrials`.

Those defaults are provenance, not required installation locations.

For a portable run, pass explicit paths on the command line:

```bash
bash screening/run_screen_v3_2_5.sh \
    /path/to/ctgov_api_v2 \
    /path/to/output \
    --no-zip
```

The underlying Python program exposes `--input-dir` and `--output-dir` for this
purpose. The frozen Python file should not be edited merely to change local
defaults because preserving it allows verification against the authoritative
study-time SHA-256.
