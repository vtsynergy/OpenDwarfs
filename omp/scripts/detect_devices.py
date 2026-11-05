#!/usr/bin/env python3
"""Auto-detect this machine's cache sizes for the 7-tier sizing config
(devices.json, see devices.example.json), using tools that are already
installed -- no new packages, no source builds.

Sources, in order of trust:

  - CPU:      `lscpu -J -C`            -- exact L1d/L2/L3 sizes from the
                                          kernel's own cacheinfo, aggregated
                                          across all cores (all-size).
  - AMD GPU:  `rocminfo`                -- exact per-CU L1 and total L2 from
                                          the ROCm runtime; L1 is scaled by
                                          the reported Compute Unit count to
                                          match this project's "total across
                                          the device" convention.
  - NVIDIA GPU: no local command exposes real L1/L2 byte sizes (nvidia-smi
                has no such field; OpenCL's reported cache size for the
                CUDA platform does not match published specs -- verified
                against A100's known 40 MB L2 and found wrong). Identified
                via `nvidia-smi --query-gpu=name,compute_cap`, cache sizes
                then come from a small embedded table of public spec-sheet
                numbers keyed by compute capability. Anything not in that
                table is reported as unknown so it doesn't silently poison
                the tier math -- add it to devices.json by hand instead.

Usage:
    ./detect_devices.py --out devices.json
    ./detect_devices.py --out devices.json --merge existing_devices.json
"""

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path

# Public spec-sheet L1 (combined SM L1/shared, total across the chip) and L2
# (total, shared) sizes in bytes, keyed by CUDA compute capability. This is
# a lookup table, not a measurement -- no local tool exposes NVIDIA cache
# sizes without installing the CUDA toolkit (deviceQuery) or similar.
NVIDIA_CACHE_BY_CC = {
    "7.0": {"sm_count": 80,  "l1_per_sm": 128 * 1024, "l2_total": 6 * 1024**2},    # V100
    "8.0": {"sm_count": 108, "l1_per_sm": 192 * 1024, "l2_total": 40 * 1024**2},   # A100
    "8.6": {"sm_count": 84,  "l1_per_sm": 128 * 1024, "l2_total": 6 * 1024**2},    # A10/A40-class (approx)
    "9.0": {"sm_count": 132, "l1_per_sm": 256 * 1024, "l2_total": 50 * 1024**2},   # H100 / H200
}


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def slugify(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")


def parse_size(s: str) -> int:
    """lscpu's human sizes: '4M' -> 4*1024**2, '512K' -> 512*1024, etc."""
    m = re.match(r"^([\d.]+)([KMGT]?)$", s.strip())
    if not m:
        raise ValueError(f"unrecognized size '{s}'")
    val, unit = float(m.group(1)), m.group(2)
    mult = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}[unit]
    return int(val * mult)


def detect_cpu():
    if not shutil.which("lscpu"):
        return {}
    caches = run(["lscpu", "-J", "-C"])
    model = run(["lscpu", "-J"])
    if caches.returncode != 0:
        return {}
    data = json.loads(caches.stdout)["caches"]
    by_name = {c["name"]: c for c in data}
    l1 = by_name.get("L1d")
    l2 = by_name.get("L2")
    l3 = by_name.get("L3")
    if l1 is None:
        return {}

    cpu_model = "cpu"
    try:
        for field in json.loads(model.stdout)["lscpu"]:
            if field.get("field", "").rstrip(":") == "Model name":
                cpu_model = field["data"]
                break
    except Exception:
        pass

    return {
        slugify(cpu_model) or "cpu": {
            "L1": parse_size(l1["all-size"]),
            "L2": parse_size(l2["all-size"]) if l2 else None,
            "LLC": parse_size(l3["all-size"]) if l3 else None,
            "_source": "lscpu -C (measured, aggregate across all cores)",
        }
    }


def detect_amd_gpus():
    if not shutil.which("rocminfo"):
        return {}
    out = run(["rocminfo"])
    if out.returncode != 0:
        return {}

    devices = {}
    agent = {}

    def flush():
        if agent.get("Device Type") == "GPU" and "Marketing Name" in agent:
            name = agent["Marketing Name"]
            cu = int(agent.get("Compute Unit", 0) or 0)
            l1_per_cu = agent.get("L1_KB")
            l2_total = agent.get("L2_KB")
            if l1_per_cu is None or cu == 0:
                return
            devices[slugify(name)] = {
                "L1": int(l1_per_cu) * 1024 * cu,
                "L2": int(l2_total) * 1024 if l2_total is not None else None,
                "LLC": None,
                "_source": f"rocminfo (measured: {l1_per_cu}KB/CU x {cu} CUs; L2 total as reported)",
            }

    for line in out.stdout.splitlines():
        line = line.strip()
        if re.match(r"^Agent \d+", line):
            flush()
            agent = {}
        elif ":" in line:
            key, _, val = line.partition(":")
            key, val = key.strip(), val.strip()
            if key == "L1":
                m = re.match(r"^(\d+)", val)
                if m:
                    agent["L1_KB"] = m.group(1)
            elif key == "L2":
                m = re.match(r"^(\d+)", val)
                if m:
                    agent["L2_KB"] = m.group(1)
            elif key in ("Marketing Name", "Device Type", "Compute Unit"):
                agent[key] = val
    flush()
    return devices


def detect_nvidia_gpus():
    if not shutil.which("nvidia-smi"):
        return {}
    out = run(["nvidia-smi", "--query-gpu=name,compute_cap", "--format=csv,noheader"])
    if out.returncode != 0:
        return {}

    devices = {}
    unknown = set()
    for line in out.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        name, _, cc = line.partition(",")
        name, cc = name.strip(), cc.strip()
        spec = NVIDIA_CACHE_BY_CC.get(cc)
        key = slugify(name)
        if spec is None:
            unknown.add((name, cc))
            continue
        devices[key] = {
            "L1": spec["l1_per_sm"] * spec["sm_count"],
            "L2": spec["l2_total"],
            "LLC": None,
            "_source": f"nvidia-smi identified compute_cap={cc}; sizes from embedded public-spec "
                       f"table, NOT measured -- no local tool exposes real NVIDIA cache sizes",
        }
    for name, cc in unknown:
        print(f"WARNING: no cache-size entry for '{name}' (compute_cap={cc}); "
              f"add it to NVIDIA_CACHE_BY_CC or to devices.json by hand")
    return devices


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--merge", type=Path, default=None,
                     help="existing devices.json to merge into (detected entries win on key collision)")
    args = ap.parse_args()

    devices = {}
    devices.update(detect_cpu())
    devices.update(detect_amd_gpus())
    devices.update(detect_nvidia_gpus())

    if not devices:
        raise SystemExit("no devices detected -- none of lscpu/rocminfo/nvidia-smi produced usable output")

    if args.merge and args.merge.exists():
        with open(args.merge) as f:
            existing = json.load(f)
        merged = existing.get("devices", existing)
        merged.update(devices)
        devices = merged

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump({"devices": devices}, f, indent=2)

    print(f"Written: {args.out}")
    for name, spec in devices.items():
        src = spec.get("_source", "")
        l1 = spec.get("L1")
        l2 = spec.get("L2")
        llc = spec.get("LLC")
        fmt = lambda b: f"{b / 1024**2:.2f}MB" if b else "-"
        print(f"  {name:24s} L1={fmt(l1):>10s} L2={fmt(l2):>10s} LLC={fmt(llc):>10s}  [{src}]")


if __name__ == "__main__":
    main()
