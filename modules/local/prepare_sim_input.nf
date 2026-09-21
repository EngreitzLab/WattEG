// Step 1 -- reduce the dataset to what the simulation actually needs.
//
// The input is a .h5mu written by pysceptre's export, NOT a sceptre object: nothing in the Python
// path reads R, which is what lets the environment be one Python package. The export must have been
// written with --all-cells, because the poscounts size factors are a per-cell reduction over every
// cell in the object and reproducing R's needs every cell. See src/watteg/expression.py.
//
// There is no sceptre_template.rds any more. That was an R object carrying the covariate matrix and
// the analysis parameters; those now live in sim_input.h5 and analysis_mode.tsv.

process PREPARE_SIM_INPUT {
    tag "${meta.id}"

    publishDir { "${params.outdir}/${meta.id}/prepared" }, mode: params.publish_mode

    input:
    tuple val(meta), path(dataset)

    output:
    tuple val(meta), path('sim_input.h5'),            emit: sim_input
    tuple val(meta), path('pairs.tsv'),               emit: pairs
    tuple val(meta), path('pairs_with_info.tsv'),     emit: pairs_with_info, optional: true
    tuple val(meta), path('grna_targets.tsv'),        emit: grna_targets
    tuple val(meta), path('discovery_threshold.txt'), emit: threshold
    tuple val(meta), path('analysis_mode.tsv'),       emit: analysis_mode
    path 'versions.yml',                              emit: versions

    script:
    def alpha = params.alpha ? "--threshold ${params.alpha}" : ''
    """
    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        watteg-prepare-sim-input \\
            --dataset ${dataset} \\
            --outdir . \\
            --n-jobs ${task.cpus} \\
            ${alpha}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(pixi run --frozen --manifest-path ${projectDir}/pixi.toml python -c 'import platform; print(platform.python_version())')
        watteg: \$(pixi run --frozen --manifest-path ${projectDir}/pixi.toml python -c 'import watteg; print(watteg.__version__)')
        pysceptre: \$(pixi run --frozen --manifest-path ${projectDir}/pixi.toml python -c 'import pysceptre; print(getattr(pysceptre, "__version__", "unknown"))')
    END_VERSIONS
    """

    stub:
    """
    touch sim_input.h5 pairs.tsv pairs_with_info.tsv grna_targets.tsv discovery_threshold.txt analysis_mode.tsv
    echo '"${task.process}": {}' > versions.yml
    """
}
