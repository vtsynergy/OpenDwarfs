#!/bin/bash
# CPU thread sweep (PROTOCOL-final.md Addendum 3): points 8, 32, 128, 1 threads; 4 CPU configurations; reduced repetitions. Resumable per app and per thread point;
# resubmits itself when time runs out. Submit once:  sbatch scripts/slurm/perf_sweep.sbatch      Results: omp/results/perf-sweep-t<N>/
#SBATCH --job-name=odwarfs-sweep
#SBATCH --partition=hudson
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gres=gpu:2
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/home/nabayanc/SOURCES/OpenDwarfs2.0/omp/results/slurm-sweep-%j.out
REPO=${REPO:-/home/nabayanc/SOURCES/OpenDwarfs2.0}
CHAIN_LEFT=${CHAIN_LEFT:-8}
cd "$REPO" || exit 1
END=$(date -d "$(scontrol show job "$SLURM_JOB_ID" | grep -o 'EndTime=[^ ]*' | cut -d= -f2)" +%s)
echo "job $SLURM_JOB_ID on $(hostname) $(date) chain_left=$CHAIN_LEFT end=$(date -d @$END +%T)"
more=0
for t in 8 32 128 1; do
  scripts/perf_final.sh --outdir "$REPO/omp/results/perf-sweep-t$t" --threads $t --devices cpu --cfgs "omp-clang omp-nvhpc opencl opencl-intel" --reps 5 --warmup 1 --timeout 300 --end-epoch "$END" --reserve 600
  rc=$?; echo "point t=$t exit $rc at $(date)"; [ $rc -eq 75 ] && { more=1; break; }
done
if [ $more -eq 1 ] && [ "$CHAIN_LEFT" -gt 0 ]; then
  sbatch --export=ALL,CHAIN_LEFT=$((CHAIN_LEFT-1)) --dependency=afterany:$SLURM_JOB_ID "$REPO/scripts/slurm/perf_sweep.sbatch"
fi
exit 0
