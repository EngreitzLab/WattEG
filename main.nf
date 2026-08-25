#!/usr/bin/env nextflow

// Element-gene power analysis.
//
//   nextflow run . -profile sherlock -params-file config/config.yml
//
// The DAG, and why it has the shape it has:
//
//   samplesheet -> PREPARE_SIM_INPUT ---+-> SPLIT_PAIRS -----------------+
//                                       |                               |
//                                       +-> FIT_NULL_MODELS (x reps)    |
//                                           -> MERGE_NULL_MODELS -------+
//                                                                       v
//                                       POWER_SIMULATION (split x effect size)
//                                            -> collectFile by (sample, effect size)
//                                            -> COMPUTE_POWER -> SUMMARIZE_POWER
//
// FIT_NULL_MODELS hangs off PREPARE_SIM_INPUT rather than off SPLIT_PAIRS because a gene's null
// model is independent of both the target and the effect size: it is fitted on a null simulation
// with no knockdown. One set of fits therefore serves every split and every effect size in a sweep
// -- 100 replicates, not 100 x however many effect sizes. Making it a sibling of SPLIT_PAIRS rather
// than a descendant is what expresses that.

nextflow.enable.dsl = 2

include { PREPARE_SIM_INPUT } from './modules/local/prepare_sim_input'
include { SPLIT_PAIRS       } from './modules/local/split_pairs'
include { FIT_NULL_MODELS   } from './modules/local/fit_null_models'
include { MERGE_NULL_MODELS } from './modules/local/merge_null_models'
include { POWER_SIMULATION  } from './modules/local/power_simulation'
include { COMPUTE_POWER     } from './modules/local/compute_power'
include { SUMMARIZE_POWER   } from './modules/local/summarize_power'

workflow {
    main:

    // ---- parameter checks -----------------------------------------------------------------
    //
    // Cheap to check here, expensive to discover 40 minutes into a 1,000-task run.
    if (!params.effect_sizes || params.effect_sizes.size() == 0) {
        error "effect_sizes is empty -- nothing to simulate."
    }
    if (params.num_replicates % params.reps_per_chunk != 0) {
        error "num_replicates (${params.num_replicates}) must be a multiple of reps_per_chunk " +
              "(${params.reps_per_chunk}); otherwise the last chunk is short and the power " +
              "denominators differ between pairs."
    }
    if (params.num_replicates % params.reps_per_null_chunk != 0) {
        error "num_replicates (${params.num_replicates}) must be a multiple of " +
              "reps_per_null_chunk (${params.reps_per_null_chunk}), or some replicate will have " +
              "no null model fitted for it."
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
    def resolve = { p ->
        p.toString().startsWith('/') || p.toString().matches('(?i)^[a-z][a-z0-9+.-]*://.*') ? file(p) : file("${projectDir}/${p}")
    }

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
            if (!row.sample?.trim() || !row.sceptre_object?.trim()) {
                error "samplesheet ${params.samplesheet} needs non-empty 'sample' and " +
                      "'sceptre_object' columns; got: ${row}"
            }
            def obj = resolve(row.sceptre_object.trim())
            if (!obj.exists()) {
                error "sample '${row.sample}': sceptre object not found at ${obj}"
            }
            // Optional: only set when the sceptre object's response matrix is odm-backed. Absent
            // for every existing samplesheet, which is why the column and the check are optional.
            def odm = row.response_odm?.trim() ? resolve(row.response_odm.trim()) : []
            if (odm && !odm.exists()) {
                error "sample '${row.sample}': --response-odm file not found at ${odm}"
            }
            [ [id: row.sample.trim()], obj, odm ]
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

    // ---- step 2b: null models, one task per replicate chunk --------------------------------
    //
    // Fanned out over replicate offsets and joined back to the sample's prepared inputs. Divisibility
    // of num_replicates by reps_per_null_chunk is checked above, so every chunk is full width.
    // The chunk count is bounded when the range is built rather than with `take` afterwards, both
    // because it avoids the operator-argument problem above and because it is what it means: there
    // are only this many chunks, not "there are 100 and we ignore most of them".
    def reps_to_fit    = params.test_max_null_reps ?: params.num_replicates
    def n_null_chunks  = reps_to_fit.intdiv(params.reps_per_null_chunk)

    ch_null_offsets = Channel
        .of(0..<n_null_chunks)
        .map { i -> [ i * params.reps_per_null_chunk, params.reps_per_null_chunk ] }

    ch_prepared = PREPARE_SIM_INPUT.out.sim_input
        .join(PREPARE_SIM_INPUT.out.template)
        .join(PREPARE_SIM_INPUT.out.grna_targets)

    FIT_NULL_MODELS(ch_prepared.combine(ch_null_offsets))

    // ---- step 2c: merge the chunks ---------------------------------------------------------
    //
    // groupTuple with an explicit size would deadlock if a chunk failed; the default waits for the
    // channel to close instead, and merge_null_models.R independently checks the replicate count.
    ch_chunks = FIT_NULL_MODELS.out.chunk.groupTuple()

    def merged_reps = params.test_max_null_reps ?: params.num_replicates
    MERGE_NULL_MODELS(ch_chunks, merged_reps)

    // ---- step 4: the simulation ------------------------------------------------------------
    //
    // The per-sample inputs are one item; the fan-out is the cross product of splits, effect sizes
    // and replicate chunks. Joining the null models in first keeps the sample's five inputs
    // together, so `combine` only ever multiplies out the things that genuinely vary.
    ch_sim_inputs = ch_prepared.join(MERGE_NULL_MODELS.out.null_models)

    ch_rep_chunks = Channel
        .of(0..<(params.num_replicates.intdiv(params.reps_per_chunk)))
        .map { i -> [ i * params.reps_per_chunk, params.reps_per_chunk ] }

    // combine on the meta key so a multi-sample run pairs each sample with its own splits rather
    // than with every sample's.
    // No trailing map: `combine` flattens the [offset, reps] pair into two elements rather than
    // keeping it as one, so this already emits the nine fields POWER_SIMULATION declares --
    // meta, sim_input, template, grna_targets, null_models, split, effect_size, rep_offset, reps.
    ch_sim_tasks = ch_sim_inputs
        .combine(ch_splits, by: 0)
        .combine(Channel.fromList(params.effect_sizes))
        .combine(ch_rep_chunks)

    POWER_SIMULATION(ch_sim_tasks)

    // ---- step 5: power, per (sample, effect size) -------------------------------------------
    //
    // groupTuple over (meta, effect size) collects every split's output for one effect size.
    // Deliberately not collectFile: compute_power.R takes the file list itself and validates the
    // replicate count per pair, which a concatenation would hide.
    ch_by_es = POWER_SIMULATION.out.sim.groupTuple(by: [0, 1])

    COMPUTE_POWER(ch_by_es, PREPARE_SIM_INPUT.out.threshold.first())

    // ---- step 6: one row per pair across every effect size ---------------------------------
    ch_all_power = COMPUTE_POWER.out.power
        .map { meta, es, power -> [meta, power] }
        .groupTuple()

    SUMMARIZE_POWER(ch_all_power)
}

workflow.onComplete {
    log.info(
        """
        ${workflow.success ? 'Completed' : 'FAILED'}: ${workflow.runName}
          duration : ${workflow.duration}
          outdir   : ${params.outdir}
          command  : ${workflow.commandLine}
        """.stripIndent()
    )
}
