#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * QC Test Pipeline
 *
 * Simulates a typical NGS QC workflow and produces qc_metrics.json
 * with per-sample metrics: q30_pct, total_reads, dup_rate, gc_pct,
 * mean_coverage, insert_size_median.
 *
 * Metric values drift slightly each run (seeded by date) so anomaly
 * detection fires after 4-5 runs.
 */

params.samples      = ['HG001', 'HG002', 'HG003']
params.read_length  = 150
params.num_reads    = 200000
params.run_seed     = (new Date().format('yyyyMMdd')).toInteger()  // changes daily

process SIMULATE_READS {
    tag "${sample_id}"
    memory '256 MB'
    cpus 1

    input:
    val(sample_id)

    output:
    tuple val(sample_id), path("${sample_id}.fastq")

    script:
    """
    python3 -c "
import random, hashlib
seed = int(hashlib.md5('${sample_id}'.encode()).hexdigest(), 16) % (2**31)
random.seed(seed)
bases = 'ACGT'
qual  = 'I' * ${params.read_length}
with open('${sample_id}.fastq', 'w') as f:
    for i in range(${params.num_reads}):
        seq = ''.join(random.choices(bases, k=${params.read_length}))
        f.write(f'@${sample_id}_read_{i}\n{seq}\n+\n{qual}\n')
"
    """
}

process FASTQC_SIM {
    tag "${sample_id}"
    memory '256 MB'
    cpus 1

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_fastqc.json")

    script:
    """
    python3 -c "
import random, hashlib, json, math

# Base metrics seeded per sample (stable)
sample_seed = int(hashlib.md5('${sample_id}'.encode()).hexdigest(), 16) % (2**31)
rng = random.Random(sample_seed)
base_q30   = rng.uniform(82.0, 88.0)
base_reads = rng.randint(4_500_000, 5_500_000)
base_dup   = rng.uniform(0.08, 0.18)
base_gc    = rng.uniform(46.0, 52.0)

# Daily drift: small noise added each run
drift_seed = sample_seed ^ ${params.run_seed}
drng = random.Random(drift_seed)
q30   = round(base_q30  + drng.gauss(0, 1.5), 2)
reads = int(base_reads   + drng.gauss(0, 150_000))
dup   = round(base_dup   + drng.gauss(0, 0.015), 4)
gc    = round(base_gc    + drng.gauss(0, 1.0), 2)

metrics = {
    'q30_pct':    max(70.0, min(99.0, q30)),
    'total_reads': max(1_000_000, reads),
    'dup_rate':   max(0.0, min(0.99, dup)),
    'gc_pct':     max(30.0, min(70.0, gc)),
}
with open('${sample_id}_fastqc.json', 'w') as f:
    json.dump(metrics, f)
"
    """
}

process ALIGNMENT_QC {
    tag "${sample_id}"
    memory '512 MB'
    cpus 2

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_align.json")

    script:
    """
    python3 -c "
import random, hashlib, json

sample_seed = int(hashlib.md5('${sample_id}'.encode()).hexdigest(), 16) % (2**31)
rng = random.Random(sample_seed)
base_cov    = rng.uniform(28.0, 40.0)
base_insert = rng.randint(280, 340)

drift_seed = sample_seed ^ ${params.run_seed}
drng = random.Random(drift_seed)
cov    = round(base_cov    + drng.gauss(0, 2.0), 2)
insert = int(base_insert    + drng.gauss(0, 12))

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

# Build lookup by sample_id prefix
def load_by_prefix(files):
    out = {}
    for f in files:
        sample = f.split('_')[0]
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
print(json.dumps(result, indent=2))
"
    """
}

workflow {
    samples_ch = Channel.fromList(params.samples)

    reads     = SIMULATE_READS(samples_ch)
    fastqc_ch = FASTQC_SIM(reads)
    align_ch  = ALIGNMENT_QC(reads)

    AGGREGATE_QC(
        fastqc_ch.map { it[0] }.collect(),
        fastqc_ch.map { it[1] }.collect(),
        align_ch.map  { it[1] }.collect(),
    )
}
