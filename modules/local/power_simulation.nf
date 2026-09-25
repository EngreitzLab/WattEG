// Step 3 -- the simulation. One task per (split, effect size, simulation chunk).
//
// This is where essentially all of the compute goes.
//
// WHY THE FAN-OUT IS OVER (split x effect size) AND NOT A LOOP OVER EFFECT SIZES
//
// Looping inside the task would share the per-task setup across effect sizes. That setup is a
// second or so against a task measured in minutes, while the loop would multiply each task's
// duration by the number of effect sizes and divide the parallelism by the same factor.
//
// NULL-MODEL FITS. Under --driver fast --null-fits reuse each gene's null model comes from
// FIT_NULL_MODELS (one fit per gene per simulation, on an independent draw with no knockdown), staged
// here as null_fits.h5 and checked against this run's seed, baseline and sim_input before use. Under
// --null-fits refit the placeholder is empty and every gene is refitted on each simulation's own
// counts for every target, the exact configuration. R hoisted the fits for speed too; this is the
// same approximation, with its seed-matching guard.
//
// POWER INSIDE THE TASK. A task that holds all simulations of its pairs (reps_per_chunk ==
// num_replicates, the default) writes each pair's counts -- simulations called, simulations used,
// the fold-change and cell-count sums -- and COMPUTE_POWER adds them up into the same power table
// the per-simulation rows give, byte for byte (tests/test_power.py). The per-simulation rows, 74M
// of them for moi5 trans, are written only under keep_per_simulation, or when simulations are
// chunked across tasks, where a pair's counts would come from several tasks and the rows are what
// power is computed from.

process POWER_SIMULATION {
    tag "${meta.id} ${split.baseName} es${effect_size} reps ${rep_offset + 1}-${rep_offset + reps}"

    // NOT published. These per-split files are what makes the run resumable and
    // preemption-tolerant -- many tasks write concurrently and no single-file format supports that.
    // COMPUTE_POWER publishes the power table; CONSOLIDATE_REPLICATES, when it runs, turns the
    // per-simulation rows into one Parquet file per effect size and publishes that. They survive in
    // the work directory until Nextflow cleans it, so a failed later step cannot lose data.

    input:
    // analysis_mode.tsv is read by the CLI from --prepared (here, the task directory). It was not
    // staged until 2026-09-24, so every real run of this process died with FileNotFoundError and
    // every Python number so far came from direct CLI calls; the stub never reads it.
    //
    // null_fits is FIT_NULL_MODELS' file under --null-fits reuse, and an empty placeholder otherwise.
    // threshold is the sample's discovery threshold, which the per-pair counts are taken at.
    tuple val(meta), path(sim_input), path(grna_targets), path(analysis_mode), path(pairs_with_info),
          path(null_fits), path(threshold), path(split), val(effect_size), val(rep_offset), val(reps)

    output:
    tuple val(meta), val(effect_size), path(out_name), emit: sim, optional: true
    tuple val(meta), val(effect_size), path(counts_name), emit: partials, optional: true

    script:
    // The simulation range is in the filenames so chunks of one split cannot collide and a stray
    // file is attributable.
    def stem = "${split.baseName}_es${effect_size}_rep${rep_offset}"
    out_name = "${stem}.tsv.gz"
    counts_name = "${stem}.partial.tsv.gz"
    def whole = (params.reps_per_chunk as int) == (params.num_replicates as int)
    def write_rows = params.keep_per_simulation || !whole
    def fits_arg = null_fits ? "--null-fits-file ${null_fits}" : ''
    def rows_arg = write_rows ? "--out ${stem}.tsv" : ''
    def counts_arg = whole ? "--partials-out ${stem}.partial.tsv --threshold-file ${threshold}" : ''
    """
    # One thread per worker. The task's parallelism is its --n-jobs worker processes; a BLAS or
    # OpenMP thread pool in the parent when it forks them is oversubscription at best and a known
    # cause of crashed children at worst.
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        watteg-run-power-simulation \\
            --prepared . \\
            --pairs ${split} \\
            --effect-size ${effect_size} \\
            --reps ${reps} \\
            --rep-offset ${rep_offset} \\
            --guide-spread-c ${params.guide_spread_c} \\
            --estimand ${params.estimand} \\
            --seed ${params.seed} \\
            --n-jobs ${task.cpus} \\
            --expression-model ${params.expression_model} \\
            --permutations ${params.permutations} \\
            --nulls ${params.nulls} \\
            --driver ${params.driver} \\
            --null-fits ${params.null_fits} ${fits_arg} \\
            ${rows_arg} ${counts_arg}

    # A short file means rows were lost, which compute_power would otherwise absorb as a smaller
    # denominator for the affected pairs. A target that perturbs no cell is skipped and writes
    # nothing, so it shows here too.
    n_pairs=\$(( \$(wc -l < ${split}) - 1 ))
    check_lines() {
        local file=\$1 expected=\$2 what=\$3
        local actual=\$(wc -l < "\${file}")
        if [ "\${actual}" -ne "\${expected}" ]; then
            echo "ERROR: \${file} has \${actual} lines, expected \${expected} (\${what} + header)." >&2
            echo "A target that perturbs no cell is skipped and writes no rows; look for 'skipped' above." >&2
            exit 1
        fi
    }
    if [ -f ${stem}.tsv ]; then
        check_lines ${stem}.tsv \$(( n_pairs * ${reps} + 1 )) "\${n_pairs} pairs x ${reps} simulations"
        gzip -f ${stem}.tsv
    fi
    if [ -f ${stem}.partial.tsv ]; then
        check_lines ${stem}.partial.tsv \$(( n_pairs + 1 )) "\${n_pairs} pairs"
        # Every pair's counts must cover all of its simulations.
        awk -F '\\t' 'NR == 1 { for (i = 1; i <= NF; i++) col[\$i] = i; next }
            \$col["n_simulations"] != ${reps} { bad++ }
            END { if (bad) { print "ERROR: " bad " pair(s) do not have ${reps} simulations" > "/dev/stderr"; exit 1 } }' \\
            ${stem}.partial.tsv
        gzip -f ${stem}.partial.tsv
    fi
    """

    stub:
    def stem = "${split.baseName}_es${effect_size}_rep${rep_offset}"
    out_name = "${stem}.tsv.gz"
    counts_name = "${stem}.partial.tsv.gz"
    def whole = (params.reps_per_chunk as int) == (params.num_replicates as int)
    def write_rows = params.keep_per_simulation || !whole
    """
    ${write_rows ? "printf '' | gzip > ${out_name}" : ''}
    ${whole ? "printf '' | gzip > ${counts_name}" : ''}
    """
}
