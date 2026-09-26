// Step 5 -- a power estimate per pair, from the simulation tasks' per-pair counts or, when
// simulations were chunked across tasks, from the consolidated per-simulation rows. Both routes give
// the same table (watteg.power.power_from_counts; tests/test_power.py).
//
// Power is the fraction of simulations whose p-value clears the discovery threshold, and the
// threshold comes from the sceptre object's own @discovery_result rather than from a nominal alpha,
// so it reflects the multiple-testing correction actually applied.
//
// Wilson intervals, not the normal approximation: at 100 simulations the normal interval gives
// [0, 0] for a pair that never cleared the threshold, which is exactly the case the analysis cares
// about. `power_ci_low` is what gets thresholded to establish a negative -- see
// docs/output.md#interpreting-negatives.

process COMPUTE_POWER {
    tag "${meta.id} es${effect_size}"

    publishDir { "${params.outdir}/${meta.id}/power" }, mode: params.publish_mode

    input:
    // Either every simulation task's per-pair counts for one effect size (*.partial.tsv.gz), or
    // one consolidated Parquet from CONSOLIDATE_REPLICATES; the file list below is kept general so
    // an older sweep's TSVs still work.
    //
    // The threshold rides in the same tuple, joined on meta in main.nf. It used to arrive as a second
    // channel built with `.first()`, so in a multi-sample run every sample was scored against sample
    // 1's discovery threshold.
    tuple val(meta), val(effect_size), path(simulations, stageAs: 'sim/*'), path(threshold)

    output:
    tuple val(meta), val(effect_size), path("power_es${effect_size}.tsv"), emit: power

    script:
    """
    # An explicit list rather than a glob, so the ordering is stable and the count is checked.
    # All three forms matched: the per-replicate output is gzipped, the consolidated one is
    # Parquet, and an older sweep being re-aggregated is plain .tsv. Missing the .gz here would
    # find no files and produce an empty power table rather than an error.
    #
    # nullglob, NOT `ls`. Only one of the three patterns ever matches, so `ls` exits 2 for the
    # other two, and under Nextflow's `bash -e` that killed the task on this line with no output
    # at all -- every real Python run, first seen on moi5 cis 2026-09-25. R's module escaped it
    # only because it piped `ls` into `paste`, whose status is the one bash checks.
    shopt -s nullglob
    partial_list=(sim/*.partial.tsv.gz)
    sim_list=(sim/*.parquet sim/*.tsv.gz sim/*.tsv)
    shopt -u nullglob
    if [ \${#partial_list[@]} -gt 0 ]; then
        source_arg=--partials
        sim_list=("\${partial_list[@]}")
    else
        source_arg=--simulations
    fi
    echo "combining \${#sim_list[@]} file(s) (\${source_arg})"
    if [ \${#sim_list[@]} -eq 0 ]; then
        echo "ERROR: no simulation files staged under sim/" >&2
        exit 1
    fi

    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        watteg-compute-power \\
            \${source_arg} "\${sim_list[@]}" \\
            --threshold-file ${threshold} \\
            --conf-level ${params.conf_level} \\
            --out power_es${effect_size}.tsv
    """

    stub:
    """
    touch power_es${effect_size}.tsv
    """
}
