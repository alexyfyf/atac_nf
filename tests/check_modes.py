#!/usr/bin/env python3
"""Assert that every assay mode preset is complete and that nextflow_schema.json agrees.

An assay mode is data, not code (see conf/modes.config), so adding cut&run or cut&tag should
be a new block and nothing else. The cost of that is that a typo or a forgotten attribute has
no compiler to catch it: the lookup just returns null, and the failure surfaces much later as
a MACS3 invocation with a stray "null" in it, or as a Tn5 shift that quietly does not run.

This check closes that gap at the source. It asserts that every mode declares the full
attribute contract with the right type, and that the `mode` enum in nextflow_schema.json lists
exactly the modes that exist -- so a new preset cannot ship undocumented, and a removed one
cannot linger in the schema.

`main.nf` re-checks the selected mode at runtime, which covers a mode added by a user's own
config; this covers every mode in the repository, whether or not CI ever runs it.

Reads `nextflow config` rather than the config file, so ${projectDir} and friends are already
resolved. Needs no pipeline run and touches no network.
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

BOOL, STR, CHOICE = "bool", "str", "choice"

# Keep in step with modeContract() in main.nf and the header of conf/modes.config.
CONTRACT = {
    "description": STR,
    "single_end": BOOL,
    "shift": BOOL,
    "trim": BOOL,
    "adapter": STR,
    "macs_input": CHOICE,
    "macs_extra_args": STR,
    "hmmratac": BOOL,
}
MACS_INPUTS = {"bed", "bam"}
BOOLS = {"true", "false"}

SCHEMA = Path(__file__).resolve().parent.parent / "nextflow_schema.json"


def mode_blocks(text):
    """Yield (name, {key: value}) for each mode in `nextflow config` output."""
    m = re.search(r"^   modes \{$(.*?)^   \}$", text, re.M | re.S)
    if not m:
        sys.exit("could not find a `modes` block in the config output")
    for block in re.finditer(r"^      (\S+) \{$(.*?)^      \}$", m.group(1), re.M | re.S):
        # Strings are quoted, booleans are not, so capture the raw right-hand side.
        entries = dict(re.findall(r"^\s+(\w+) = (.*)$", block.group(2), re.M))
        yield block.group(1), {k: v.strip().strip("'") for k, v in entries.items()}


def raw_values(text, mode):
    """The same entries, unstripped, so `true` can be told from `'true'`."""
    m = re.search(rf"^      {re.escape(mode)} \{{$(.*?)^      \}}$", text, re.M | re.S)
    return dict(re.findall(r"^\s+(\w+) = (.*)$", m.group(1), re.M)) if m else {}


def schema_modes():
    with SCHEMA.open() as fh:
        schema = json.load(fh)
    defs = schema.get("$defs") or schema.get("definitions") or {}
    for group in defs.values():
        prop = group.get("properties", {}).get("mode")
        if prop is not None:
            return prop.get("enum")
    return None


def main():
    # CI has nextflow on PATH; $NEXTFLOW lets a developer point at a local build instead.
    nextflow = os.environ.get("NEXTFLOW", "nextflow")
    out = subprocess.run([nextflow, "config"], capture_output=True, text=True, check=True).stdout

    failures, found = [], []
    for name, entry in mode_blocks(out):
        found.append(name)
        raw = raw_values(out, name)

        missing = [k for k in CONTRACT if k not in entry]
        if missing:
            failures.append(f"{name}: missing {', '.join(missing)}")
        extra = [k for k in entry if k not in CONTRACT]
        if extra:
            failures.append(
                f"{name}: declares {', '.join(extra)}, which no attribute contract covers. "
                "Add it to CONTRACT here and to modeContract() in main.nf, or drop it."
            )

        for key, kind in CONTRACT.items():
            if key not in entry:
                continue
            value = entry[key]
            if kind is BOOL and raw.get(key, "").strip() not in BOOLS:
                failures.append(f"{name}: {key} must be a bare true/false, got {raw.get(key)}")
            elif kind is CHOICE and value not in MACS_INPUTS:
                failures.append(
                    f"{name}: {key} must be one of {'/'.join(sorted(MACS_INPUTS))}, got '{value}'"
                )
            elif kind is STR and not raw.get(key, "").strip().startswith("'"):
                # An unquoted value is read as a boolean or a number, and a mode attribute that
                # reaches a script block as `true` renders literally: `--keep-dup all true`.
                failures.append(f"{name}: {key} must be a quoted string, got {raw.get(key)}")
            elif kind is STR and key == "description" and not value:
                failures.append(f"{name}: description is empty")

        if not failures or found[-1] == name:
            print(f"  {name:<8} shift={entry.get('shift', '?'):<5} "
                  f"trim={entry.get('trim', '?'):<5} "
                  f"macs_input={entry.get('macs_input', '?'):<4} "
                  f"single_end={entry.get('single_end', '?'):<5} {entry.get('description', '')}")

    enum = schema_modes()
    if enum is None:
        failures.append("nextflow_schema.json has no `mode` property with an enum")
    elif sorted(enum) != sorted(found):
        failures.append(
            f"nextflow_schema.json lists modes {sorted(enum)} but conf/modes.config "
            f"defines {sorted(found)}"
        )

    if failures:
        print("\nFAILED:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        sys.exit(1)
    if not found:
        sys.exit("no assay modes were checked")
    print(f"\nOK: {len(found)} assay modes are complete and match the schema")


if __name__ == "__main__":
    main()
