#!/bin/bash
# =============================================================================
# Is the machine quiet enough for performance comparisons?
#   scripts/system_clear.sh [--gpu N] [--max-load L | --no-load] [--max-foreign-cores C] [--quiet]
# Exit 0 = clear, 1 = busy. A GPU counts as busy if it has any compute process from
# another user or utilization > 5%. The CPU counts as busy if the 1-min load average
# exceeds --max-load (default: 25% of the logical CPUs). Prints a one-line snapshot that
# the run scripts also store next to every result (see scripts/load_snapshot.sh).
# =============================================================================
GPU=${GPU_INDEX:-0}; MAXLOAD=""; QUIET=0; NOLOAD=0; MAXFC=""
while [ $# -gt 0 ]; do case $1 in --gpu) GPU=$2; shift 2;; --max-load) MAXLOAD=$2; shift 2;; --no-load) NOLOAD=1; shift;; --max-foreign-cores) MAXFC=$2; shift 2;; --quiet) QUIET=1; shift;; *) shift;; esac; done
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
[ $NOLOAD = 0 ] && awk -v l="$load" -v m="$MAXLOAD" 'BEGIN{exit !(l>m)}' && { busy=1; why="$why load=$load>$MAXLOAD"; }
# CPU used by OTHER users right now (two /proc samples 1 s apart), in cores. The load average cannot be used while our own
# multi-threaded benchmark runs: it counts our threads too.
fc=$(python3 -I - <<'PY'
import os, time
me = os.getuid(); hz = os.sysconf('SC_CLK_TCK')
def sample():
    d = {}
    for p in os.listdir('/proc'):
        if not p.isdigit(): continue
        try:
            if os.stat('/proc/' + p).st_uid == me: continue
            f = open('/proc/%s/stat' % p).read().rsplit(')', 1)[1].split()
            d[p] = int(f[11]) + int(f[12])
        except Exception: pass
    return d
a = sample(); time.sleep(1.0); b = sample()
print("%.1f" % (sum(max(0, b[p] - a.get(p, b[p])) for p in b) / hz))
PY
)
[ -n "$MAXFC" ] && awk -v c="$fc" -v m="$MAXFC" 'BEGIN{exit !(c>m)}' && { busy=1; why="$why foreign_cpu=${fc}>${MAXFC}cores"; }
[ $QUIET = 1 ] || echo "gpu=$GPU util=${util}% foreign_procs=$others foreign_cores=$fc load1=$load ncpu=$ncpu -> $([ $busy = 0 ] && echo CLEAR || echo "BUSY ($why )")"
exit $busy
