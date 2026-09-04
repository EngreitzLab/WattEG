// Step 4 -- the simulation. One task per (split, effect size, replicate chunk).
//
// This is where essentially all of the compute goes: 635 CPU-hours per effect size on
// day0_grna20_no_shuffle at 100 replicates, against a few CPU-hours for everything else combined.
//
// WHY THE FAN-OUT IS OVER (split x effect size) AND NOT A LOOP OVER EFFECT SIZES
//
// Looping inside the task would share the per-task setup -- the 16 MB sim_input.rds and 47 MB
// sceptre_template.rds reads, plus container and R startup -- across effect sizes. That setup is
// ~1 second against a task measured in tens of minutes, while the loop would multiply each task's
// duration by the number of effect sizes and divide the parallelism by the same factor. With
// thousands of cores available on `owners`, repeating a little I/O is much the better trade.
//
// --null-precomputations is NOT optional here. Without it run_power_simulation.R either refits the
// null model inside every call (4.3x the cost) or -- if the template still carried the inherited
// real-data cache -- would silently test simulated counts against real-data coefficients and
// understate power. The template has that slot cleared by slim_sceptre_object(), and the script
// errors rather than falling back, so the failure mode that went unnoticed through the whole
// refactor cannot recur here.

process POWER_SIMULATION {
    tag "${meta.id} ${split.baseName} es${effect_size} reps ${rep_offset + 1}-${rep_offset + reps}"

    // Published for the same reason the sbatch runner keeps them: they are the only record of the
    // per-replicate p-values, they let a run be compared split by split against the other runner,
    // and re-deriving one costs a full task. This is the bulky output -- ~1,000 files per effect
    // size, a few hundred MB in total.
    // per_replicate, not sim: a row here is one (pair, replicate) test result, and "sim" said nothing
// about that. This is the only output from which power can be RE-derived -- subsampling replicates
// to study reduced designs, or bootstrapping -- so it is worth keeping and worth naming clearly.
    publishDir "${params.outdir}/${meta.id}/per_replicate/es${effect_size}", mode: params.publish_mode

    input:
    tuple val(meta), path(sim_input), path(sceptre_template), path(grna_targets),
          path(null_models), path(split), val(effect_size), val(rep_offset), val(reps)

    output:
    tuple val(meta), val(effect_size), path(out_name), emit: sim

    script:
    // The replicate range is part of the filename so chunks of one split cannot collide, and so a
    // stray file is attributable. With no chunking (reps_per_chunk == num_replicates) there is
    // exactly one per split, which is what every measured run has done.
    out_name = "${split.baseName}_es${effect_size}_rep${rep_offset}.tsv.gz"
    """
    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        Rscript ${projectDir}/src/run_power_simulation.R \\
            --sim-input ${sim_input} \\
            --sceptre-template ${sceptre_template} \\
            --pairs ${split} \\
            --grna-targets ${grna_targets} \\
            --null-precomputations ${null_models} \\
            --effect-size ${effect_size} \\
            --reps ${reps} \\
            --rep-offset ${rep_offset} \\
            --guide-sd ${params.guide_sd} \\
            --seed ${params.seed} \\
            --out ${out_name}

    # A short file means replicate rows were lost, which compute_power.R would otherwise absorb as
    # a smaller denominator for the affected pairs.
    n_pairs=\$(( \$(wc -l < ${split}) - 1 ))
    expected=\$(( n_pairs * ${reps} + 1 ))
    actual=\$(gzip -cd ${out_name} | wc -l)
    if [ "\${actual}" -ne "\${expected}" ]; then
        echo "ERROR: wrote \${actual} lines, expected \${expected} (\${n_pairs} pairs x ${reps} reps + header)." >&2
        exit 1
    fi
    """

    stub:
    out_name = "${split.baseName}_es${effect_size}_rep${rep_offset}.tsv.gz"
    """
    printf '' | gzip > ${out_name}
    """
}
