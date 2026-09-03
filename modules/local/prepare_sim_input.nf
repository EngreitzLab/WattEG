// Step 1 -- reduce the sceptre object to what the simulation actually needs.
//
// One input file, not four: the other three the old pipeline demanded duplicate data already inside
// the object (@discovery_pairs_with_info, @grna_target_data_frame, @discovery_result). See
// docs/status.md, "Measured improvements".
//
// The emitted sceptre_template.rds has its @response_matrix emptied and, importantly, its inherited
// @response_precomputations cleared by slim_sceptre_object() -- that clearing is what makes an
// accidental as_is run impossible to reintroduce downstream.

process PREPARE_SIM_INPUT {
    tag "${meta.id}"

    publishDir "${params.outdir}/${meta.id}/prepared", mode: params.publish_mode

    input:
    tuple val(meta), path(sceptre_object), path(response_odm)

    output:
    tuple val(meta), path('sim_input.rds'),         emit: sim_input
    tuple val(meta), path('sceptre_template.rds'),  emit: template
    tuple val(meta), path('pairs.tsv'),             emit: pairs
    tuple val(meta), path('grna_targets.tsv'),      emit: grna_targets
    tuple val(meta), path('discovery_threshold.txt'), emit: threshold
    tuple val(meta), path('analysis_mode.tsv'),     emit: analysis_mode
    path 'versions.yml',                            emit: versions

    script:
    // response_odm is [] (no files staged) for every sample whose response matrix isn't odm-backed.
    def response_odm_flag = response_odm ? "--response-odm ${response_odm}" : ''
    def control_cells = params.n_control_cells ? "--n-control-cells ${params.n_control_cells}" : ''
    def batches       = params.cell_batches    ? "--cell-batches ${params.cell_batches}"       : ''
    def alpha         = params.alpha           ? "--alpha ${params.alpha}"                     : ''
    """
    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        Rscript ${projectDir}/src/prepare_sim_input.R \\
            --sceptre-object ${sceptre_object} \\
            --outdir . \\
            ${response_odm_flag} ${control_cells} ${batches} ${alpha}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(pixi run --frozen --manifest-path ${projectDir}/pixi.toml Rscript -e 'cat(strsplit(R.version.string, " ")[[1]][3])')
        sceptre: \$(pixi run --frozen --manifest-path ${projectDir}/pixi.toml Rscript -e 'cat(as.character(utils::packageVersion("sceptre")))')
    END_VERSIONS
    """

    stub:
    """
    touch sim_input.rds sceptre_template.rds pairs.tsv grna_targets.tsv discovery_threshold.txt analysis_mode.tsv
    echo '"${task.process}": {}' > versions.yml
    """
}
