#!/usr/bin/env nextflow

// Element-gene power analysis.
//
//   nextflow run . -profile sherlock -params-file config/config.yml
//
// The DAG, and why it has the shape it has:
//
//   samplesheet -> PREPARE_SIM_INPUT -> SPLIT_PAIRS
//                        |                   -> POWER_SIMULATION (split x effect size x rep chunk)
//                        +-> FIT_NULL_MODELS ---^   -> per-pair counts -> COMPUTE_POWER
//                                                   -> [CONSOLIDATE_REPLICATES]   -> SUMMARIZE_POWER
//
// Power is computed inside each simulation task when the task holds all simulations of its pairs
// (reps_per_chunk == num_replicates, the default): it writes per-pair counts and COMPUTE_POWER adds
// them up, which gives the table the per-simulation rows give, byte for byte. The rows, and
// CONSOLIDATE_REPLICATES that publishes them as Parquet, run only under keep_per_simulation or when
// simulations are chunked across tasks, where power is computed from the rows as before.
//
// FIT_NULL_MODELS runs under --driver fast --null-fits reuse: each gene's null model is fitted once
// per simulation, on an independent draw with no knockdown, and every simulation task reuses that
// fit for every target the gene is tested with, instead of refitting it per target. It is R's old
// FIT_NULL_MODELS approximation, rebuilt in Python and measured before it was adopted
// (docs/pysceptre-backend.md, section 13, item 3c). --null-fits refit skips it and keeps the exact
// per-target refit.

nextflow.enable.dsl = 2

include { PREPARE_SIM_INPUT } from './modules/local/prepare_sim_input'
include { SPLIT_PAIRS       } from './modules/local/split_pairs'
include { FIT_NULL_MODELS   } from './modules/local/fit_null_models'
include { POWER_SIMULATION  } from './modules/local/power_simulation'
include { CONSOLIDATE_REPLICATES } from './modules/local/consolidate_replicates'
include { COMPUTE_POWER     } from './modules/local/compute_power'
include { SUMMARIZE_POWER   } from './modules/local/summarize_power'

// Resolve a samplesheet path against the repository root rather than the launch directory, so a
// run's validity does not depend on where it was started from. A scheme-prefixed URI (gs://, s3://,
// az://) is absolute in the same sense a leading '/' is, and must not be prefixed either.
//
// Test the STRING, not file(p).isAbsolute(): Nextflow's file() resolves a relative path against
// launchDir and hands back an absolute path, so isAbsolute() is always true and the projectDir
// fallback would never fire -- silently making resolution depend on the launch directory, which is
// the thing this is here to prevent.
//
// A function, not a closure assigned with `def`: Nextflow 26.04's strict syntax does not see the
// latter from inside a workflow body.
def resolve(p) {
    p.toString().startsWith('/') || p.toString().matches('(?i)^[a-z][a-z0-9+.-]*://.*')
        ? file(p)
        : file("${projectDir}/${p}")
}

workflow {
    main:

    // ---- parameter checks -----------------------------------------------------------------
    //
    // Cheap to check here, expensive to discover 40 minutes into a 1,000-task run.
    // guide_sd was an absolute spread (0.13). Its replacement, guide_spread_c, is a coefficient on
    // es * (1 - es), so reading an old 0.13 as c would shrink the spread fivefold. Refuse it.
    if (params.guide_sd != null) {
        error "guide_sd was replaced by guide_spread_c (default 0.65) on 2026-09-25; the spread is " +
              "now c * es * (1 - es), not an absolute sd. Remove guide_sd from the params."
    }
    if (!params.effect_sizes || params.effect_sizes.size() == 0) {
        error "effect_sizes is empty -- nothing to simulate."
    }
    if (params.reps_per_chunk != params.num_replicates && !params.keep_per_simulation) {
        log.warn "reps_per_chunk (${params.reps_per_chunk}) < num_replicates " +
                 "(${params.num_replicates}): a pair's simulations span several tasks, so power " +
                 "is computed from the per-simulation rows, which are written and consolidated " +
                 "as if keep_per_simulation were true."
    }
    if (params.num_replicates % params.reps_per_chunk != 0) {
        error "num_replicates (${params.num_replicates}) must be a multiple of reps_per_chunk " +
              "(${params.reps_per_chunk}); otherwise the last chunk is short and the power " +
              "denominators differ between pairs."
    }
    // The fast driver shares one set of permutation matrices across a target's simulations, which
    // is the engine's test only when they share one permutation set; reused fits are read only by
    // the fast driver. Both combinations would fail every simulation task.
    if (params.driver == 'fast' && params.permutations != 'per-target') {
        error "driver 'fast' needs permutations 'per-target' (got '${params.permutations}'); " +
              "set driver 'engine' to run per-replicate permutations."
    }
    if (params.null_fits == 'reuse' && params.driver != 'fast') {
        error "null_fits 'reuse' needs driver 'fast'; the engine refits every gene itself. " +
              "Set null_fits 'refit' with driver 'engine'."
    }
    // The gcb profile has no default container_image -- see conf/gcb.config for why -- so a run
    // that forgets --container_image would otherwise fail 5-30 minutes in, on the first task,
    // with "pixi: command not found" and no indication the image was ever the problem.
    if (workflow.profile.tokenize(',').contains('gcb') && !params.container_image) {
        error "profile 'gcb' requires --container_image (a Wave/Docker image with the pinned " +
              "pixi environment baked in -- see conf/gcb.config)."
    }

    // ---- inputs ---------------------------------------------------------------------------
    //
    // Samplesheet paths are resolved against the repository root, not the launch directory. The
    // committed samplesheets use repo-relative paths, so resolving against the launch directory
    // would make a run's validity depend on where it was started from.
    // Test the *string*, not file(p).isAbsolute(). Nextflow's file() resolves a relative path
    // against launchDir and hands back an absolute path, so isAbsolute() is always true and the
    // projectDir fallback would never fire -- the resolution would silently depend on the launch
    // directory, which is the thing this is here to prevent.
    //
    // A scheme-prefixed URI (gs://, s3://, az://) is absolute in the same sense a leading '/' is --
    // it already names a full location, not one relative to the repo. Without this check,
    // '${projectDir}/gs://bucket/obj' is nonsense and never exists.

    // Emptiness is checked on the file, eagerly, rather than with .ifEmpty on the channel.
    // ifEmpty's closure is invoked while the DAG is being built, not when the channel turns out to
    // be empty, so `.ifEmpty { error ... }` aborts every run and -- worse -- reports the empty-
    // samplesheet message in place of whatever the real error was.
    def sheet = resolve(params.samplesheet)
    if (!sheet.exists()) {
        error "samplesheet not found: ${sheet}"
    }
    def sheet_rows = sheet.readLines().findAll { it.trim() }
    if (sheet_rows.size() < 2) {
        error "samplesheet ${sheet} has a header but no sample rows."
    }

    ch_samples = Channel
        .fromPath(sheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample?.trim() || !row.dataset?.trim()) {
                error "samplesheet ${params.samplesheet} needs non-empty 'sample' and " +
                      "'dataset' columns; got: ${row}"
            }
            // A .h5mu from pysceptre's export, written with --all-genes --all-cells. The sceptre
            // object is no longer an input to this pipeline at all: converting it is a one-off
            // step that happens outside, which is what lets the environment hold no R.
            def dataset = resolve(row.dataset.trim())
            if (!dataset.exists()) {
                error "sample '${row.sample}': dataset not found at ${dataset}"
            }
            [ [id: row.sample.trim()], dataset ]
        }

    // ---- step 1: reduce the sceptre object -------------------------------------------------
    PREPARE_SIM_INPUT(ch_samples)

    // ---- step 2: split pairs into per-task chunks ------------------------------------------
    SPLIT_PAIRS(PREPARE_SIM_INPUT.out.pairs)

    // One item per split, carrying its sample. `flatten` on the path list would lose the meta, so
    // transpose the [meta, [split, split, ...]] tuple instead.
    //
    // `take` gets params.test_max_splits directly. Neither a `def` local nor an `as int` cast
    // survives Nextflow's operator dispatch -- both arrive wrapped in a PojoWrapper and fail with
    // "Missing process or function take(...)", which reads like a missing operator rather than an
    // argument-type problem. A literal or a params value works; nothing else here does.
    ch_splits_all = SPLIT_PAIRS.out.splits.transpose()

    ch_splits = params.test_max_splits
        ? ch_splits_all.take(params.test_max_splits)
        : ch_splits_all

    if (params.test_max_splits) {
        log.warn "test_max_splits=${params.test_max_splits}: simulating only the first " +
                 "${params.test_max_splits} of ${params.n_splits} splits. NOT a complete run."
    }

    // ---- step 3: the simulation ------------------------------------------------------------
    //
    // The per-sample inputs are one item; the fan-out is the cross product of splits, effect sizes
    // and replicate chunks. Joining the sample's inputs first keeps them together, so `combine`
    // only ever multiplies out the things that genuinely vary.
    //
    // pairs_with_info is optional -- an export with no discovery_pairs_with_info produces none --
    // so it is mixed in with a default rather than joined, which would drop the sample entirely.
    //
    // The null-model fits: FIT_NULL_MODELS' file under null_fits 'reuse', an empty placeholder
    // otherwise. The process is always wired, and fed nothing when it is not wanted, so there is
    // one channel shape either way.
    FIT_NULL_MODELS(PREPARE_SIM_INPUT.out.sim_input
                        .join(PREPARE_SIM_INPUT.out.pairs)
                        .join(PREPARE_SIM_INPUT.out.analysis_mode)
                        .filter { params.null_fits == 'reuse' })
    ch_null_fits = params.null_fits == 'reuse'
        ? FIT_NULL_MODELS.out.fits
        : PREPARE_SIM_INPUT.out.sim_input.map { meta, _sim_input -> [meta, []] }

    ch_sim_inputs = PREPARE_SIM_INPUT.out.sim_input
        .join(PREPARE_SIM_INPUT.out.grna_targets)
        .join(PREPARE_SIM_INPUT.out.analysis_mode)
        .join(PREPARE_SIM_INPUT.out.pairs_with_info, remainder: true)
        .map { meta, sim_input, grna_targets, analysis_mode, info ->
            [meta, sim_input, grna_targets, analysis_mode, info ?: []] }
        .join(ch_null_fits)
        .join(PREPARE_SIM_INPUT.out.threshold)

    ch_rep_chunks = Channel
        .of(0..<(params.num_replicates.intdiv(params.reps_per_chunk)))
        .map { i -> [ i * params.reps_per_chunk, params.reps_per_chunk ] }

    // combine on the meta key so a multi-sample run pairs each sample with its own splits rather
    // than with every sample's.
    // No trailing map: `combine` flattens the [offset, reps] pair into two elements rather than
    // keeping it as one, so this already emits the eleven fields POWER_SIMULATION declares --
    // meta, sim_input, grna_targets, analysis_mode, pairs_with_info, null_fits, threshold, split,
    // effect_size, rep_offset, reps.
    ch_sim_tasks = ch_sim_inputs
        .combine(ch_splits, by: 0)
        .combine(Channel.fromList(params.effect_sizes))
        .combine(ch_rep_chunks)

    POWER_SIMULATION(ch_sim_tasks)

    // ---- step 4: the per-simulation rows, when they are written --------------------------------
    //
    // groupTuple over (meta, effect size) collects every split's output for one effect size. With
    // no rows written (the default) the channel is empty and CONSOLIDATE_REPLICATES never runs.
    // Deliberately not collectFile: the power step takes the file list itself and validates it.
    CONSOLIDATE_REPLICATES(POWER_SIMULATION.out.sim.groupTuple(by: [0, 1]))

    // ---- step 5: power, per (sample, effect size) -------------------------------------------
    //
    // From the tasks' per-pair counts when each task held all simulations of its pairs; from the
    // consolidated rows otherwise. `whole` mirrors the test POWER_SIMULATION makes.
    //
    // Each sample's own threshold, joined on the meta key -- not `.first()`, which scored every
    // sample of a multi-sample run against sample 1's.
    ch_power_in = (params.reps_per_chunk as int) == (params.num_replicates as int)
        ? POWER_SIMULATION.out.partials.groupTuple(by: [0, 1])
        : CONSOLIDATE_REPLICATES.out.parquet.map { meta, es, f -> [meta, es, [f]] }
    COMPUTE_POWER(ch_power_in.combine(PREPARE_SIM_INPUT.out.threshold, by: 0))

    // ---- step 6: one row per pair across every effect size ---------------------------------
    ch_all_power = COMPUTE_POWER.out.power
        .map { meta, es, power -> [meta, power] }
        .groupTuple()

    // sim_input rides along so the summary can carry the gene's dispersion and normalised mean.
    // Joined at this step rather than emitted per replicate: it is a per-GENE constant, and adding
    // it to the simulation output would only help sweeps run after the change, whereas summarising
    // again from stored power tables costs seconds and works on sweeps already finished.
    //
    // `.join` on the meta key, NOT two separate input channels. Nextflow pairs multiple input
    // channels positionally, so with more than one sample the power tables of one could be matched
    // to the sim_input of another -- silently, producing a summary whose dispersions belong to a
    // different dataset. Joining on the key makes that impossible.
    ch_summary_in = ch_all_power.join(PREPARE_SIM_INPUT.out.sim_input)

    SUMMARIZE_POWER(ch_summary_in)
}

// The `workflow.onComplete { ... }` summary that used to sit here is gone: Nextflow 26.04's strict
// syntax rejects top-level statements, and everything it printed -- success, duration, outdir and
// the command line -- Nextflow's own completion summary and `-with-report` already carry.
