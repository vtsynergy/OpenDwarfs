#!/bin/bash
# =============================================================================
# Is the machine quiet enough for performance comparisons?
#   scripts/system_clear.sh [--gpu N] [--max-load L] [--quiet]
# Exit 0 = clear, 1 = busy. A GPU counts as busy if it has any compute process from
# another user or utilization > 5%. The CPU counts as busy if the 1-min load average
# exceeds --max-load (default: 25% of the logical CPUs). Prints a one-line snapshot that
# the run scripts also store next to every result (see scripts/load_snapshot.sh).
# =============================================================================
GPU=${GPU_INDEX:-0}; MAXLOAD=""; QUIET=0
while [ $# -gt 0 ]; do case $1 in --gpu) GPU=$2; shift 2;; --max-load) MAXLOAD=$2; shift 2;; --quiet) QUIET=1; shift;; *) shift;; esac; done
ncpu=$(nproc); [ -z "$MAXLOAD" ] && MAXLOAD=$(( ncpu / 4 ))
load=$(cut -d' ' -f1 /proc/loadavg)
uuid=$(nvidia-smi --query-gpu=index,uuid --format=csv,noheader | awk -F', ' -v g="$GPU" '$1==g{print $2}')
util=$(nvidia-smi --query-gpu=index,utilization.gpu --format=csv,noheader,nounits | awk -F', ' -v g="$GPU" '$1==g{print $2}')
others=0
while IFS=, read -r pid u; do
  u=${u// /}; [ "$u" = "$uuid" ] || continue
  owner=$(ps -o user= -p "$pid" 2>/dev/null | tr -d ' '); [ "$owner" = "$USER" ] || others=$((others+1))
done < <(nvidia-smi --query-compute-apps=pid,gpu_uuid --format=csv,noheader)
busy=0; why=""
[ "${util:-0}" -gt 5 ] && { busy=1; why="$why gpu_util=${util}%"; }
[ "$others" -gt 0 ] && { busy=1; why="$why foreign_gpu_procs=$others"; }
awk -v l="$load" -v m="$MAXLOAD" 'BEGIN{exit !(l>m)}' && { busy=1; why="$why load=$load>$MAXLOAD"; }
[ $QUIET = 1 ] || echo "gpu=$GPU util=${util}% foreign_procs=$others load1=$load ncpu=$ncpu -> $([ $busy = 0 ] && echo CLEAR || echo "BUSY ($why )")"
exit $busy
