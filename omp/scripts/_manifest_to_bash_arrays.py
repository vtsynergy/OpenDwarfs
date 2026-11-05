#!/usr/bin/env python3
"""Validate a tier manifest and emit numeric Bash array assignments."""
import json
import re
import sys


def assignments(data):
    for dwarf, tiers in data["dwarfs"].items():
        if not re.fullmatch(r"[a-z][a-z0-9_]*", dwarf):
            raise ValueError(f"invalid benchmark name: {dwarf!r}")
        for tier, params in tiers.items():
            if tier not in {"l1lo", "l1hi", "l2lo", "l2hi", "llclo", "llchi", "stream"}:
                raise ValueError(f"invalid tier: {tier!r}")
            if params is None:
                continue
            for array, value in (("SIZE", params["size"]),
                                 ("SEC", params.get("secondary_size", params["size"])),
                                 ("ITERS", params.get("iterations", 1))):
                if type(value) is not int or value < (0 if array == "SEC" else 1):
                    raise ValueError(f"{dwarf}:{tier}: invalid {array}: {value!r}")
                yield f"{array}[{dwarf}:{tier}]={value}"


if __name__ == "__main__":
    try:
        with open(sys.argv[1]) as source:
            lines = list(assignments(json.load(source)))
        print("\n".join(lines))
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        sys.exit(f"Invalid tier manifest: {exc}")
