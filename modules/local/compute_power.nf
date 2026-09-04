// Step 5 -- turn per-replicate simulation rows into a power estimate per pair.
//
// Power is the fraction of replicates whose p-value clears the discovery threshold, and the
// threshold comes from the sceptre object's own @discovery_result rather than from a nominal alpha,
// so it reflects the multiple-testing correction actually applied.
//
// Wilson intervals, not the normal approximation: at 100 replicates the normal interval gives
// [0, 0] for a pair that never cleared the threshold, which is exactly the case the analysis cares
// about. `power_ci_low` is what gets thresholded to certify a negative -- see
// docs/output.md#interpreting-negatives.

process COMPUTE_POWER {
    tag "${meta.id} es${effect_size}"

    publishDir "${params.outdir}/${meta.id}/power", mode: params.publish_mode

    input:
    tuple val(meta), val(effect_size), path(simulations, stageAs: 'sim/*')
    tuple val(meta2), path(threshold)

    output:
    tuple val(meta), val(effect_size), path("power_es${effect_size}.tsv"), emit: power

    script:
    """
    # A comma-separated list rather than a glob: compute_power.R accepts either, but an explicit
    # list is a single token and its ordering is stable.
    # Both forms matched: the per-replicate output is gzipped, but an older sweep being
    # re-aggregated is plain .tsv. Missing the .gz here would have found no files and produced an
    # empty power table rather than an error.
    sim_list=\$(ls sim/*.tsv.gz sim/*.tsv 2>/dev/null | paste -sd, -)
    echo "combining \$(ls sim/*.tsv.gz sim/*.tsv 2>/dev/null | wc -l) simulation file(s)"

    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        Rscript ${projectDir}/src/compute_power.R \\
            --simulations "\${sim_list}" \\
            --threshold-file ${threshold} \\
            --conf-level ${params.conf_level} \\
            --out power_es${effect_size}.tsv
    """

    stub:
    """
    touch power_es${effect_size}.tsv
    """
}
