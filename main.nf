#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * QC Test Pipeline
 *
 * Simulates a typical NGS QC workflow and produces qc_metrics.json
 * with per-sample metrics: q30_pct, total_reads, dup_rate, gc_pct,
 * mean_coverage, insert_size_median.
 *
 * params.mode:
 *   "normal"  — all samples within expected range (use to build baseline)
 *   "outlier" — 3 samples pushed beyond 3σ on 2 metrics each (triggers fail/warn)
 *
 * Run normal first to build history, then run outlier to see anomaly detection.
 */

params.mode        = "normal"   // "normal" | "outlier"
params.num_samples = 100
params.read_length = 150
params.num_reads   = 200000
params.run_seed    = (new Date().format('yyyyMMdd')).toInteger()

process FASTQC_SIM {
    tag "${sample_id}"
    memory '128 MB'
    cpus 1

    input:
    tuple val(sample_id), val(mode), val(run_seed)

    output:
    tuple val(sample_id), path("${sample_id}_fastqc.json")

    script:
    """
    python3 -c "
import random, hashlib, json

sample_seed = int(hashlib.md5('${sample_id}'.encode()).hexdigest(), 16) % (2**31)
rng = random.Random(sample_seed)

# Stable per-sample base values
base_q30   = rng.uniform(83.0, 87.0)
base_reads = rng.randint(4_600_000, 5_400_000)
base_dup   = rng.uniform(0.09, 0.17)
base_gc    = rng.uniform(47.0, 51.0)

# Small daily drift
drift_seed = sample_seed ^ ${run_seed}
drng = random.Random(drift_seed)

mode = '${mode}'
outlier_samples = ['S003', 'S027', 'S081']  # fixed set that will be outliers

if mode == 'outlier' and '${sample_id}' in outlier_samples:
    # Push q30 and dup_rate well beyond 3sigma
    q30  = round(base_q30  - drng.gauss(6.0, 0.3), 2)   # ~4 sigma low
    dup  = round(base_dup  + drng.gauss(0.12, 0.01), 4)  # ~4 sigma high
else:
    q30  = round(base_q30  + drng.gauss(0, 0.8), 2)
    dup  = round(base_dup  + drng.gauss(0, 0.008), 4)

reads = int(base_reads + drng.gauss(0, 80_000))
gc    = round(base_gc  + drng.gauss(0, 0.5), 2)

metrics = {
    'q30_pct':     max(50.0, min(99.0, q30)),
    'total_reads': max(1_000_000, reads),
    'dup_rate':    max(0.0, min(0.99, dup)),
    'gc_pct':      max(30.0, min(70.0, gc)),
}
with open('${sample_id}_fastqc.json', 'w') as f:
    json.dump(metrics, f)
"
    """
}

process ALIGNMENT_QC {
    tag "${sample_id}"
    memory '128 MB'
    cpus 1

    input:
    tuple val(sample_id), val(mode), val(run_seed)

    output:
    tuple val(sample_id), path("${sample_id}_align.json")

    script:
    """
    python3 -c "
import random, hashlib, json

sample_seed = int(hashlib.md5('${sample_id}'.encode()).hexdigest(), 16) % (2**31)
rng = random.Random(sample_seed)

base_cov    = rng.uniform(29.0, 39.0)
base_insert = rng.randint(285, 335)

drift_seed = sample_seed ^ ${run_seed}
drng = random.Random(drift_seed)

mode = '${mode}'
outlier_samples = ['S003', 'S027', 'S081']

if mode == 'outlier' and '${sample_id}' in outlier_samples:
    cov    = round(base_cov    - drng.gauss(8.0, 0.5), 2)   # ~4 sigma low
    insert = int(base_insert   + drng.gauss(60, 3))           # ~5 sigma high
else:
    cov    = round(base_cov    + drng.gauss(0, 1.0), 2)
    insert = int(base_insert   + drng.gauss(0, 6))

metrics = {
    'mean_coverage':      max(1.0, cov),
    'insert_size_median': max(100, insert),
}
with open('${sample_id}_align.json', 'w') as f:
    json.dump(metrics, f)
"
    """
}

process AGGREGATE_QC {
    memory '256 MB'
    cpus 1

    publishDir "results", mode: 'copy'

    input:
    val(all_samples)
    path(fastqc_jsons)
    path(align_jsons)

    output:
    path("qc_metrics.json")

    script:
    def sample_list = all_samples.join(' ')
    """
    python3 -c "
import json, glob

samples = '${sample_list}'.split()
fastqc_files = sorted(glob.glob('*_fastqc.json'))
align_files  = sorted(glob.glob('*_align.json'))

def load_by_prefix(files):
    out = {}
    for f in files:
        sample = '_'.join(f.split('_')[:-1])  # handle S001, S002 etc.
        out[sample] = json.load(open(f))
    return out

fastqc = load_by_prefix(fastqc_files)
align  = load_by_prefix(align_files)

result = {'samples': {}}
for s in samples:
    merged = {}
    merged.update(fastqc.get(s, {}))
    merged.update(align.get(s, {}))
    result['samples'][s] = merged

with open('qc_metrics.json', 'w') as f:
    json.dump(result, f, indent=2)
print('Aggregated', len(result['samples']), 'samples, mode=${params.mode}')
"
    """
}

workflow {
    def sample_ids = (1..params.num_samples).collect { String.format("S%03d", it) }

    samples_ch = Channel.fromList(sample_ids)
        .map { sid -> tuple(sid, params.mode, params.run_seed) }

    fastqc_ch = FASTQC_SIM(samples_ch)
    align_ch  = ALIGNMENT_QC(samples_ch)

    AGGREGATE_QC(
        fastqc_ch.map { it[0] }.collect(),
        fastqc_ch.map { it[1] }.collect(),
        align_ch.map  { it[1] }.collect(),
    )
}
