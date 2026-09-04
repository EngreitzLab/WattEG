// Step 4b -- one Parquet file per effect size, from that effect size's per-split output.
//
// The per-replicate output is the only thing power can be re-derived from, so it is kept -- but as
// 1,000 gzipped TSVs it cost 30-90 s of parsing per read, and the reduced-design study reads it
// hundreds of times. One Parquet file reads in a second or two and takes the published output from
// 6,000 files to 6.
//
// This is the ONLY published form of the per-replicate data. The per-split files stay in the work
// directory, which is what makes the run resumable and preemption-tolerant, and means a failed
// consolidation cannot lose anything.
//
// Memory: this holds one effect size's rows in a data frame before writing -- 34,886 pairs x 100
// replicates x 10 columns, measured at ~2 GB peak for compute_power.R doing the same read. 16 GB
// leaves ample headroom; raise it before n_splits if a dataset is much larger.

process CONSOLIDATE_REPLICATES {
    tag "${meta.id} es${effect_size}"

    publishDir "${params.outdir}/${meta.id}/per_replicate", mode: params.publish_mode

    input:
    tuple val(meta), val(effect_size), path(simulations, stageAs: 'sim/*')

    output:
    tuple val(meta), val(effect_size), path("replicates_es${effect_size}.parquet"), emit: parquet

    script:
    """
    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        Rscript ${projectDir}/src/consolidate_replicates.R \\
            --simulations sim \\
            --out replicates_es${effect_size}.parquet
    """

    stub:
    """
    touch replicates_es${effect_size}.parquet
    """
}
