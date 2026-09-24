// Step 3 -- the simulation. One task per (split, effect size, replicate chunk).
//
// This is where essentially all of the compute goes.
//
// WHY THE FAN-OUT IS OVER (split x effect size) AND NOT A LOOP OVER EFFECT SIZES
//
// Looping inside the task would share the per-task setup across effect sizes. That setup is a
// second or so against a task measured in minutes, while the loop would multiply each task's
// duration by the number of effect sizes and divide the parallelism by the same factor.
//
// THERE IS NO --null-precomputations, AND THAT IS THE POINT. R hoisted the per-gene null model out
// of the simulation because refitting it inside every call cost 4.3x. The Python path fits each
// gene's null on that replicate's own simulated counts as a matter of course, which is the faithful
// configuration R's hoist was built to approximate -- so FIT_NULL_MODELS and MERGE_NULL_MODELS are
// gone, along with the seed-matching guard between a bundle and the run that consumes it.

process POWER_SIMULATION {
    tag "${meta.id} ${split.baseName} es${effect_size} reps ${rep_offset + 1}-${rep_offset + reps}"

    // NOT published. These per-split files are what makes the run resumable and
    // preemption-tolerant -- many tasks write concurrently and no single-file format supports that
    // -- but as a published artefact they cost 30-90 s of parsing per read. CONSOLIDATE_REPLICATES
    // turns each effect size into one Parquet file and publishes that instead. They survive in the
    // work directory until Nextflow cleans it, so a failed consolidation cannot lose data.

    input:
    // analysis_mode.tsv is read by the CLI from --prepared (here, the task directory). It was not
    // staged until 2026-09-24, so every real run of this process died with FileNotFoundError and
    // every Python number so far came from direct CLI calls; the stub never reads it.
    tuple val(meta), path(sim_input), path(grna_targets), path(analysis_mode), path(pairs_with_info),
          path(split), val(effect_size), val(rep_offset), val(reps)

    output:
    tuple val(meta), val(effect_size), path(out_name), emit: sim

    script:
    // The replicate range is in the filename so chunks of one split cannot collide and a stray file
    // is attributable.
    out_name = "${split.baseName}_es${effect_size}_rep${rep_offset}.tsv.gz"
    """
    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        watteg-run-power-simulation \\
            --prepared . \\
            --pairs ${split} \\
            --effect-size ${effect_size} \\
            --reps ${reps} \\
            --rep-offset ${rep_offset} \\
            --guide-sd ${params.guide_sd} \\
            --seed ${params.seed} \\
            --n-jobs ${task.cpus} \\
            --expression-model ${params.expression_model} \\
            --out ${out_name.replace('.gz', '')}
    gzip -f ${out_name.replace('.gz', '')}

    # A short file means replicate rows were lost, which compute_power would otherwise absorb as a
    # smaller denominator for the affected pairs.
    n_pairs=\$(( \$(wc -l < ${split}) - 1 ))
    expected=\$(( n_pairs * ${reps} + 1 ))
    actual=\$(gzip -cd ${out_name} | wc -l)
    if [ "\${actual}" -ne "\${expected}" ]; then
        echo "ERROR: wrote \${actual} lines, expected \${expected} (\${n_pairs} pairs x ${reps} replicates + header)." >&2
        echo "A target that perturbs no cell is skipped and writes no rows; look for 'skipped' above." >&2
        exit 1
    fi
    """

    stub:
    out_name = "${split.baseName}_es${effect_size}_rep${rep_offset}.tsv.gz"
    """
    printf '' | gzip > ${out_name}
    """
}
