// Step 2b -- each gene's null-model fit, once per simulation, for every simulation task to reuse.
//
// Runs only under --driver fast --null-fits reuse (the defaults). Without it each simulation refits
// every gene once per target it is tested with, which is the largest cost left in a simulated
// pair-test (~49 ms of it on cis); the real screen fits each gene once for ~135 targets. With it,
// each (gene, simulation) is fitted once, on an independent draw with no knockdown keyed
// rng_for(seed, "__null_fit__|" + gene, simulation, 0), and every target and effect size reuses the
// fit. R's FIT_NULL_MODELS approximation, measured to leave the moi5 cis call rate unchanged
// (docs/pysceptre-backend.md, section 13, item 3c; src/watteg/null_fits.py).
//
// One task per sample, covering every gene in pairs.tsv and every simulation 1..num_replicates, so
// every simulation task, whatever its split or simulation chunk, finds its fits in one file. It
// parallelises over simulations with task.cpus workers: 244 genes x 100 simulations is ~26 CPU
// minutes on moi5 cis. The file is ~2.5 MB and records the seed, the baseline and a digest of
// sim_input.h5; POWER_SIMULATION refuses one that does not match its own run.

process FIT_NULL_MODELS {
    tag "${meta.id}"

    publishDir { "${params.outdir}/${meta.id}/prepared" }, mode: params.publish_mode

    input:
    // sim_input.h5 and analysis_mode.tsv keep their names when staged, so --prepared . finds them.
    tuple val(meta), path(sim_input), path(pairs), path(analysis_mode)

    output:
    tuple val(meta), path('null_fits.h5'), emit: fits

    script:
    """
    # One thread per worker, as in POWER_SIMULATION.
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        watteg-fit-null-models \\
            --prepared . \\
            --pairs ${pairs} \\
            --reps ${params.num_replicates} \\
            --seed ${params.seed} \\
            --expression-model ${params.expression_model} \\
            --n-jobs ${task.cpus} \\
            --out null_fits.h5
    """

    stub:
    """
    touch null_fits.h5
    """
}
