// Step 2b -- fit the per-gene null model on a null simulation of each replicate.
//
// WHY THIS EXISTS AT ALL
//
// Before testing a pair, sceptre fits a Poisson GLM of the gene's counts on the cell covariates --
// the null model the perturbed cells are compared against. It *skips* that fit whenever
// @response_precomputations already holds an entry for the response_id, and sceptre_template.rds
// inherited 272 such entries from the real discovery analysis, covering all 237 genes that have
// QC-passing pairs. So every simulated count was being tested against coefficients fitted to *real*
// counts, which understates power. Measured: +0.0063 mean power once corrected, concentrated in the
// transition band. See docs/status.md.
//
// WHY IT FANS OUT OVER REPLICATES AND NOTHING ELSE
//
// A gene's null model is fitted on a null simulation -- no knockdown -- so it depends on neither the
// target nor the effect size. One fit per (gene, simulation) therefore serves every target and every
// effect size in a sweep: 100 simulations, not 100 x targets x effect sizes. Fitting inside each
// simulation task would pay for it n_splits times over.
//
// The seed is derived from (seed, rep) only, deliberately not (seed, target, rep, effect_size), for
// exactly that reason.

process FIT_NULL_MODELS {
    tag "${meta.id} reps ${rep_offset + 1}-${rep_offset + reps}"

    input:
    tuple val(meta), path(sim_input), path(sceptre_template), path(grna_targets), val(rep_offset), val(reps)

    output:
    tuple val(meta), path("chunk_${String.format('%04d', rep_offset)}.rds"), emit: chunk

    script:
    def chunk = "chunk_${String.format('%04d', rep_offset)}.rds"
    """
    pixi run --frozen --manifest-path ${projectDir}/pixi.toml \\
        Rscript ${projectDir}/src/fit_null_models.R \\
            --sim-input ${sim_input} \\
            --sceptre-template ${sceptre_template} \\
            --grna-targets ${grna_targets} \\
            --reps ${reps} \\
            --rep-offset ${rep_offset} \\
            --seed ${params.seed} \\
            --out ${chunk}
    """

    stub:
    """
    touch chunk_${String.format('%04d', rep_offset)}.rds
    """
}
