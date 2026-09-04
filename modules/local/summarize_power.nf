// Step 6 -- one row per pair across every effect size, and the deliverable column.
//
// The deliverable is `min_detectable_effect_size`, not six power columns: "we would have caught a
// knockdown of >= 25%, so the absence of a call rules out effects that large" is the sentence the
// analysis exists to produce. It is bracketed by `_ci_low` / `_ci_high`, and the interval inverts,
// because power rises with effect size.
//
// Qualifying requires clearing the threshold at an effect size AND at every larger one tested.
// Taking the first effect size that clears would let noise bias every pair towards looking more
// detectable than it is.

process SUMMARIZE_POWER {
    tag "${meta.id}"

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_mode

    input:
    // One tuple, joined on meta upstream. Two separate input channels would be paired positionally
    // and could match one sample's power tables to another sample's sim_input.
    //
    // sim_input is here so the summary can carry the gene's dispersion and normalised mean, sparing
    // every covariate model a 16 MB RDS read for the other half of SE^2 ~ (1/n)(1/mu + 1/theta).
    tuple val(meta), path(power_files, stageAs: 'power/*'), path(sim_input)

    output:
    tuple val(meta), path('power_summary.tsv'), emit: summary

    script:
    """
    power_list=\$(ls power/*.tsv | paste -sd, -)
    echo "summarizing \$(ls power/*.tsv | wc -l) effect size(s)"

    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        Rscript ${projectDir}/src/summarize_power.R \\
            --power "\${power_list}" \\
            --power-threshold ${params.power_threshold} \\
            --sim-input ${sim_input} \\
            --out power_summary.tsv
    """

    stub:
    """
    touch power_summary.tsv
    """
}
