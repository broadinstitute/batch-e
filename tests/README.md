# Tests

## `synthetic_v1.json` — WDL smoke test

Small synthetic dataset that runs the full WDL (stage_interval scatter,
hailrunner analysis, generate_report) in ~10 min on a 2-worker Dataproc.

### Regenerate the underlying dataset

The generator is a hailrunner-compatible Hail script (same pattern as
`batch_e/batch_e.py`): plain `argparse` CLI, `hl.init()` up front, all GCS I/O
through `hl.hadoop_open`. It can be run three ways.

**Via hailrunner on Dataproc (matches production path):**

```
hailrunner run \
    --staging-bucket gs://fc-secure-b337ea24-f011-4291-9725-f553a53d6e94 \
    --script https://raw.githubusercontent.com/broadinstitute/batch-e/refs/heads/main/scripts/make_synthetic_test.py \
    --workers 2 --worker-type n1-standard-4 --driver-type n1-standard-4 \
    -- --gcs-prefix gs://fc-secure-b337ea24-f011-4291-9725-f553a53d6e94/batch_e_synth_v1
```

**Local script against a live Dataproc cluster:**

```
hailrunner run \
    --staging-bucket gs://fc-secure-.../ \
    --script scripts/make_synthetic_test.py \
    --workers 2 --worker-type n1-standard-4 --driver-type n1-standard-4 \
    -- --gcs-prefix gs://fc-secure-.../batch_e_synth_v1
```

**Locally with hail installed (fastest for iteration):**

```
python scripts/make_synthetic_test.py \
    --gcs-prefix gs://fc-secure-.../batch_e_synth_v1
```

All three produce the same layout under `--gcs-prefix`:

- `synthetic.mt/` — 300 samples x ~5000 variants on chr20:1-2 Mb
- `ancestry.tsv`, `comparison.tsv`
- `intervals/{ACMG59,Low_Mappability,GC_gt_85,GC_lt_25,HighConf_Genome}.bed.gz`

### Submit the WDL

```
cromwell submit wdl/batch_e.wdl -i tests/synthetic_v1.json
```

Outputs land at
`gs://fc-secure-.../batch_effect_results/synthetic_v1/` plus a rendered
`report.html`.
