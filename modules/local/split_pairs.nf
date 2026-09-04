// Step 2 -- split the QC-passing pairs into balanced per-task chunks.
//
// Whole targets stay together: the simulation builds perturbation status and guide-level effect
// sizes once per target, so splitting a target across tasks would duplicate that setup.
//
// n_splits affects scheduling only, not results -- seeds derive from
// (seed, target, rep, effect_size), so output is invariant to the split layout. That invariance is
// also what lets the Nextflow and sbatch runners be compared split by split.
//
// --target-overhead is deliberately left at its default of 0. It would be ~2.0 pair-equivalents at
// the current cost model, on splits already holding 34-36 pairs, which moves the array's tail task
// by about a minute; see docs/status.md, "Open questions".

process SPLIT_PAIRS {
    tag "${meta.id}"

// NOT published. The splits are parallelisation bookkeeping, not a result: they are 1,000 files and
// 1,000 inodes per sample, and they are fully regenerable -- the LPT bin packing is deterministic
// given pairs.tsv and n_splits, which docs/status.md confirms reproduces the on-disk splits exactly.
// The pairs themselves are already published in prepared/pairs.tsv. They stay in the work directory
// where Nextflow needs them for the fan-out, and go when it is cleaned.

    input:
    tuple val(meta), path(pairs)

    output:
    tuple val(meta), path('split_*.tsv'), emit: splits

    script:
    """
    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        Rscript ${projectDir}/src/split_pairs.R \\
            --pairs ${pairs} \\
            --n-splits ${params.n_splits} \\
            --outdir .

    n_written=\$(ls split_*.tsv | wc -l)
    if [ "\${n_written}" -ne "${params.n_splits}" ]; then
        echo "ERROR: wrote \${n_written} splits but expected ${params.n_splits}." >&2
        exit 1
    fi

    # Every pair must appear exactly once across the splits, or power is computed from the wrong
    # number of replicates.
    pairs_in=\$(( \$(wc -l < ${pairs}) - 1 ))
    pairs_out=\$(cat split_*.tsv | grep -vc '^grna_target' || true)
    if [ "\${pairs_in}" -ne "\${pairs_out}" ]; then
        echo "ERROR: pair count changed during splitting (\${pairs_in} -> \${pairs_out})." >&2
        exit 1
    fi
    """

    stub:
    """
    for i in \$(seq -w 1 ${params.n_splits}); do touch split_\${i}.tsv; done
    """
}
