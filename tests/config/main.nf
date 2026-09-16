// Harness for the resource closures in conf/base.config.
//
// WHAT THIS CATCHES that nothing else does. A directive closure is evaluated when Nextflow submits
// a task, which is BEFORE any task record exists. If it throws, the run stops with every task so
// far marked COMPLETED and nothing in the task list to inspect -- see docs/status.md, trap 5, where
// a String-typed memory param aborted a production run in exactly that way, 30 seconds after
// MERGE_NULL_MODELS, having already spent nothing and explained nothing.
//
// So the closures need a test that evaluates them without running the pipeline. The process names
// here must match the `withName:` selectors in conf/base.config, and the input names must match
// what those closures reference -- `split` for POWER_SIMULATION, `simulations` for the other two.
// The scripts do nothing but echo the memory they were given, which is also the assertion: a
// wrong number is as much a failure as an exception.

process POWER_SIMULATION {
    input:  path split
    output: path 'ps.txt'
    script: """echo "POWER_SIMULATION memory=${task.memory}" | tee ps.txt"""
}

process CONSOLIDATE_REPLICATES {
    input:  path simulations, stageAs: 'sim/*'
    output: path 'cr.txt'
    script: """echo "CONSOLIDATE_REPLICATES memory=${task.memory}" | tee cr.txt"""
}

process COMPUTE_POWER {
    input:  path simulations, stageAs: 'sim/*'
    output: path 'cp.txt'
    script: """echo "COMPUTE_POWER memory=${task.memory}" | tee cp.txt"""
}

workflow {
    // Sizes chosen to reproduce the moi5 trans sweep, the largest input the pipeline has run on:
    // a 30 KB split, 1.90 GB of gzipped per-split output going into consolidation, and the 1.40 GiB
    // Parquet it produces going into compute_power. run.sbatch creates them as sparse files.
    POWER_SIMULATION(Channel.fromPath("${projectDir}/inputs/split_0001.tsv"))
    CONSOLIDATE_REPLICATES(Channel.fromPath("${projectDir}/inputs/trans_sim.tsv.gz").collect())
    COMPUTE_POWER(Channel.fromPath("${projectDir}/inputs/replicates.parquet").collect())
}
