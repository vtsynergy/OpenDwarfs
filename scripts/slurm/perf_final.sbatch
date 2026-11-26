#!/bin/bash
# Final pinned-configuration run on hudson, in 4-hour chunks. Each chunk resumes from DIR/done/<app>, and resubmits itself when the
# runner exits with 75 (time ran out, work remains). Submit once:   sbatch scripts/slurm/perf_final.sbatch
# Results: omp/results/perf-final/{raw.csv,results.csv,FINAL-REPORT.md,environment.txt,LAUNCH-CONFIGS.md}
#SBATCH --job-name=odwarfs-final
#SBATCH --partition=hudson
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gres=gpu:2
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/home/nabayanc/SOURCES/OpenDwarfs2.0/omp/results/slurm-final-%j.out
REPO=${REPO:-/home/nabayanc/SOURCES/OpenDwarfs2.0}
OUT=$REPO/omp/results/perf-final
CHAIN_LEFT=${CHAIN_LEFT:-6}
cd "$REPO" || exit 1
mkdir -p "$OUT"
END=$(date -d "$(scontrol show job "$SLURM_JOB_ID" | grep -o 'EndTime=[^ ]*' | cut -d= -f2)" +%s)
echo "job $SLURM_JOB_ID on $(hostname) $(date) chain_left=$CHAIN_LEFT end=$(date -d @$END +%T)"
scripts/perf_final.sh --outdir "$OUT" --end-epoch "$END" --reserve 900
rc=$?
echo "runner exit $rc at $(date)"
if [ $rc -eq 75 ] && [ "$CHAIN_LEFT" -gt 0 ]; then
  sbatch --export=ALL,CHAIN_LEFT=$((CHAIN_LEFT-1)) --dependency=afterany:$SLURM_JOB_ID "$REPO/scripts/slurm/perf_final.sbatch"
fi
exit 0
