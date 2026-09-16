#!/bin/bash
#SBATCH --job-name=watteg_config_test
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=00:25:00
#SBATCH --output=logs/%x-%A.out
#SBATCH --error=logs/%x-%A.err
#
# Evaluate conf/base.config's resource closures against String-typed params -- see main.nf for why
# this needs its own harness and tests/config/README.md for the expected output.
#
# It runs as a job for two reasons, neither optional: a Nextflow launcher is a JVM and does not
# belong on a login node, and the local executor refuses to start a task requesting more memory
# than the node has -- so asserting the real computed requests (8 + 32 + 47 GB, concurrently)
# needs a node that can hold them. Not `bigmem`: that partition enforces a 256 GB minimum request.
#
# Usage:
#   cd <repo root> && sbatch tests/config/run.sbatch

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null \
    | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)")/../.." && pwd)"
TEST_DIR="${REPO_ROOT}/tests/config"

# Sparse files: the closures read only their SIZE, never their contents, so there is no reason to
# write 3.3 GB of bytes. Sizes match the moi5 trans sweep -- see main.nf.
mkdir -p "${TEST_DIR}/inputs"
truncate -s 30000      "${TEST_DIR}/inputs/split_0001.tsv"
truncate -s 1899680000 "${TEST_DIR}/inputs/trans_sim.tsv.gz"
truncate -s 1501940740 "${TEST_DIR}/inputs/replicates.parquet"

export NXF_HOME="${SCRATCH}/nxf_home_config_test"
export NXF_WORK="${SCRATCH}/watteg_config_test_work"
ml load biology nextflow

cd "${TEST_DIR}"
nextflow run . -ansi-log false

echo
echo "=================================================================="
echo "  computed requests"
echo "=================================================================="
cat "${NXF_WORK}"/*/*/ps.txt "${NXF_WORK}"/*/*/cr.txt "${NXF_WORK}"/*/*/cp.txt
