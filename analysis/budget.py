"""Additive time-budget partition shared by the collectors.

Every instant of the steady-state window is assigned to exactly ONE state, in priority order:
  1. device_kernel : a kernel is executing on the device (device-first rule: time the device is busy is useful work, whatever the host is doing)
  2. device_copy   : a device-side copy is executing (CUDA path only)
  3. one API class : the host is inside a traced runtime/driver call of that class (priority order below) while the device is idle
  4. host_other    : the host is running code outside the traced calls while the device is idle
`init`-class calls are removed from the window altogether. The buckets therefore sum exactly to the steady-state window, and so do their differences
between two toolchains (the old bucket accounting double-counted kernel time that host calls were waiting for).
All intervals are (start_ns, end_ns) in one clock."""
PRIORITY = ["sync", "copy_d2h", "copy_h2d", "copy_other", "alloc", "free", "event", "launch", "other"]

def merge(iv):
    out = []
    for s, e in sorted(iv):
        if e <= s: continue
        if out and s <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], e))
        else: out.append((s, e))
    return out

def clip(iv, lo, hi):
    return [(max(s, lo), min(e, hi)) for s, e in iv if min(e, hi) > max(s, lo)]

def subtract(a, b):
    """a, b: merged, sorted, disjoint interval lists -> a minus b"""
    out = []; j = 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] <= cur: j += 1
        k = j
        while k < len(b) and b[k][0] < e:
            if b[k][0] > cur: out.append((cur, b[k][0]))
            cur = max(cur, b[k][1]); k += 1
        if cur < e: out.append((cur, e))
    return out

def length(iv): return sum(e - s for s, e in iv)

def partition(api, kernel_iv, copy_iv=(), init_class="init", window_iv=None):
    """api: {class: [(s,e),...]} (any classes; unknown ones are treated as 'other'); kernel_iv/copy_iv: device-busy intervals.
    Returns dict(window_ns, init_in_window_ns, buckets_ns{...}, coverage) with sum(buckets_ns) == window_ns (steady state, init removed)."""
    calls = {c: merge(v) for c, v in api.items()}
    real = [iv for c, v in calls.items() if c != init_class for iv in v]
    if not real: return dict(window_ns=0, init_in_window_ns=0, buckets_ns={}, coverage=0.0)
    lo = min(s for s, e in real); hi = max(e for s, e in real)
    window = [(lo, hi)] if not window_iv else merge(window_iv)
    if window_iv: lo, hi = window[0][0], window[-1][1]
    outside = subtract([(lo, hi)], window)
    init = merge(clip(calls.get(init_class, []), lo, hi))
    steady = subtract(window, init)
    taken = []; out = {}
    def take(name, iv):
        nonlocal taken
        iv = subtract(subtract(subtract(merge(clip(iv, lo, hi)), outside), init), taken)
        out[name] = length(iv); taken = merge(taken + iv)
    take("device_kernel", kernel_iv); take("device_copy", list(copy_iv))
    seen = set()
    for c in PRIORITY + sorted(set(calls) - set(PRIORITY) - {init_class}):
        if c in calls and c != init_class and c not in seen: take(c, calls[c]); seen.add(c)
    out["host_other"] = length(steady) - sum(out.values())
    return dict(window_ns=length(steady), init_in_window_ns=length(init), buckets_ns=out, coverage=0.0)
