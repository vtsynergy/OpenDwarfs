#!/usr/bin/env python3
"""Compute 7-tier cache-hierarchy-grounded problem sizes for the streaming
(class M) dwarfs, generically over any set of hardware targets.

This is the portable half of the 7-tier sizing methodology: tier byte
targets are derived purely from the min/max L1/L2/LLC sizes across whatever
devices are listed in a devices.json config (see devices.example.json), so
adding, removing, or replacing hardware targets is a config edit, not a
code change. Each streaming dwarf declares an analytic working-set formula
W(N) (and its inverse, N(W)); this script solves N at each of the 7 tier
byte targets.

The other apps (gemm, cg, lud, gem, nqueens, nw, swat, hmm, bfs, tdm,
cfd) use runtime tiering. Supply candidate sizes, measure them with
run_cache_tier_sweep.sh, and adjust them using the measured runtimes.
The runner measures the candidates; it does not search for sizes.
Merge selected sizes with --seed, using the same manifest shape:
{"dwarfs": {name: {tier: {"size":.., "secondary_size":..}}}}.

Usage:
    ./compute_cache_tier_sizes.py --devices devices.json
    ./compute_cache_tier_sizes.py --devices devices.json --dwarfs atax spmv fft
    ./compute_cache_tier_sizes.py --devices devices.json --seed empirical_seed.json --out out.json

With no --out, the manifest is written under data/generated/ (relative to
this repo's root) with a filename that records which dwarfs it covers, so
picking a --dwarfs subset never silently overwrites a full-suite manifest
(or vice versa).
"""

import argparse
import json
import math
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

# ---------------------------------------------------------------------------
# Memory-streaming dwarfs (class M): analytical W(N) -> N inverses.
# "square": True means secondary_size == size; else secondary_size is fixed.
# ---------------------------------------------------------------------------
STREAMING_DWARFS = {
    "atax":     {"square": True,  "sec": None, "solve": lambda W: int(math.sqrt(W / 4))},
    "crc":      {"square": False, "sec": 128,  "solve": lambda W: W - 512},
    "fdtd2d":   {"square": True,  "sec": None, "solve": lambda W: int(math.sqrt(W / 12))},
    # fft: W = 8N + 256*N*8 = N*(8+2048) ~= 2056N (sec=256 bins, N samples, 4 arrays of floats)
    "fft":      {"square": False, "sec": 256,  "solve": lambda W: int((W / 8 - 256) / 257)},
    "kmeans":   {"square": False, "sec": 16,   "solve": lambda W: int((W - 128) / 12)},
    "pagerank": {"square": False, "sec": 2,    "solve": lambda W: int(math.sqrt(W / 8))},
    "spmv":     {"square": False, "sec": 16,   "solve": lambda W: int((W - 4) / 140)},
    "srad":     {"square": True,  "sec": None, "solve": lambda W: int(math.sqrt(W / 8))},
}


def load_devices(path: Path) -> dict:
    with open(path) as f:
        raw = json.load(f)
    devices = raw["devices"] if "devices" in raw else raw
    if not devices:
        raise SystemExit(f"{path}: no devices listed")
    return devices


def tier_targets(devices: dict) -> dict:
    def vals(key):
        out = [v[key] for v in devices.values() if v.get(key) is not None]
        if not out:
            raise SystemExit(
                f"no device in the config defines a '{key}' cache size; "
                f"every tier needs at least one L1, one L2, and one LLC value "
                f"across the device set"
            )
        return out

    l1_min, l1_max = min(vals("L1")), max(vals("L1"))
    l2_min, l2_max = min(vals("L2")), max(vals("L2"))
    llc_min, llc_max = min(vals("LLC")), max(vals("LLC"))

    return {
        "l1lo":   int(l1_min),    # resident in every device's L1
        "l1hi":   int(l1_max),    # resident in the largest L1
        "l2lo":   int(l2_min),    # crosses the smallest L2
        "l2hi":   int(l2_max),    # crosses the largest L2
        "llclo":  int(llc_min),   # crosses the smallest LLC
        "llchi":  int(llc_max),   # crosses the largest LLC
        "stream": 1024 * 1024 * 1024,  # 1 GB - exceeds every listed cache
    }


def solve_streaming(info: dict, targets: dict) -> dict:
    result = {}
    for tier, target in targets.items():
        n = max(1, info["solve"](target))
        sec = n if info["square"] else info["sec"]
        result[tier] = {"size": n, "secondary_size": sec}
    return result


def default_out_path(dwarf_names) -> Path:
    out_dir = REPO_ROOT / "data" / "generated"
    all_names = sorted(STREAMING_DWARFS)
    if sorted(dwarf_names) == all_names:
        stem = "cache_tier_sizes_7tier"
    else:
        stem = "cache_tier_sizes_7tier__" + "_".join(sorted(dwarf_names))
    return out_dir / f"{stem}.json"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--devices", required=True, type=Path,
                     help="JSON file listing {device: {L1, L2, LLC}} cache sizes in bytes (see devices.example.json)")
    ap.add_argument("--dwarfs", nargs="+", metavar="NAME", default=None,
                     choices=sorted(STREAMING_DWARFS),
                     help="only compute sizes for these streaming dwarfs (default: all of them). "
                          f"Choices: {', '.join(sorted(STREAMING_DWARFS))}")
    ap.add_argument("--seed", type=Path, default=None,
                     help="optional JSON with empirically-calibrated tier sizes for non-streaming dwarfs "
                          "(produced by sweeping run_cache_tier_sweep.sh on real hardware); merged in as-is")
    ap.add_argument("--out", type=Path, default=None,
                     help="output manifest path (default: data/generated/cache_tier_sizes_7tier"
                          "[__<selected-dwarfs>].json under the repo root)")
    args = ap.parse_args()

    devices = load_devices(args.devices)
    targets = tier_targets(devices)
    selected = args.dwarfs if args.dwarfs else sorted(STREAMING_DWARFS)
    out_path = args.out if args.out is not None else default_out_path(selected)

    manifest = {
        "metadata": {
            "tier_targets_bytes": targets,
            "device_cache_sizes_bytes": devices,
            "dwarfs_selected": selected,
        },
        "dwarfs": {},
    }

    for name in selected:
        manifest["dwarfs"][name] = solve_streaming(STREAMING_DWARFS[name], targets)

    if args.seed is not None:
        with open(args.seed) as f:
            seed = json.load(f)
        for name, tiers in seed.get("dwarfs", seed).items():
            if name in selected or args.dwarfs is None:
                manifest["dwarfs"][name] = tiers

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"Written: {out_path}")
    print("Devices used:", ", ".join(devices.keys()))
    print("Dwarfs selected:", ", ".join(selected))
    print("Tier targets (bytes):")
    for tier, b in targets.items():
        print(f"  {tier:8s}: {b / 1024**2:.2f} MB")
    print("\nStreaming dwarf sizes:")
    for name in selected:
        sizes = " | ".join(f"{t}:N={params['size']}"
                           for t, params in manifest["dwarfs"][name].items()
                           if params is not None)
        print(f"  {name:10s}: {sizes}")
    if args.seed is None:
        print(
            "\nNote: no --seed given; only streaming (class M) dwarfs have sizes. "
            "Compute-pinned/latency dwarfs (gemm, lud, gem, nqueens, nw, swat, hmm, "
            "bfs, tdm, cfd, ...) need an empirical sweep -- see run_cache_tier_sweep.sh."
        )


if __name__ == "__main__":
    main()
