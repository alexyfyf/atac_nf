#!/usr/bin/env python3
"""Assert that every genome preset pairs a blacklist and mitochondrial name whose contig
naming matches its own reference.

Mixing the two is the nastiest failure mode in this pipeline because it is silent: a
chr-prefixed blacklist against an Ensembl-named BAM overlaps nothing, and bedtools reports
that as success, so the blacklist fraction reads 0 and nothing is ever filtered.

Reads `nextflow config` rather than the config file, so ${projectDir} and friends are already
resolved. Needs no pipeline run and touches no network.
"""

import os
import re
import subprocess
import sys
from pathlib import Path

# Ensembl distributes bare contig names (1, 2, MT); UCSC and NCBI iGenomes builds are
# chr-prefixed (chr1, chrM).
BARE, PREFIXED = "bare", "chr-prefixed"
EXPECTED_MITO = {BARE: "MT", PREFIXED: "chrM"}


def genome_blocks(text):
    """Yield (name, {key: value}) for each genome in `nextflow config` output."""
    m = re.search(r"^   genomes \{$(.*?)^   \}$", text, re.M | re.S)
    if not m:
        sys.exit("could not find a `genomes` block in the config output")
    for block in re.finditer(r"^      (\S+) \{$(.*?)^      \}$", m.group(1), re.M | re.S):
        entries = dict(re.findall(r"^\s+(\w+) = '([^']*)'$", block.group(2), re.M))
        yield block.group(1), entries


def style_of_reference(fasta):
    if "/Ensembl/" in fasta:
        return BARE
    if "/UCSC/" in fasta or "/NCBI/" in fasta:
        return PREFIXED
    return None


def style_of_contig(contig):
    return PREFIXED if contig.startswith("chr") else BARE


def main():
    # CI has nextflow on PATH; $NEXTFLOW lets a developer point at a local build instead.
    nextflow = os.environ.get("NEXTFLOW", "nextflow")
    out = subprocess.run([nextflow, "config"], capture_output=True, text=True, check=True).stdout

    failures, checked = [], 0
    for name, entry in genome_blocks(out):
        fasta, blacklist, mito = entry.get("fasta"), entry.get("blacklist"), entry.get("mito_name")
        if not (fasta and blacklist and mito):
            failures.append(f"{name}: missing fasta, blacklist or mito_name")
            continue

        expected = style_of_reference(fasta)
        if expected is None:
            failures.append(f"{name}: cannot tell the naming style from {fasta}")
            continue

        path = Path(blacklist)
        if not path.is_file():
            failures.append(f"{name}: blacklist not found at {blacklist}")
            continue

        with path.open() as fh:
            contigs = {line.split("\t", 1)[0] for line in fh if line.strip() and not line.startswith("#")}
        if not contigs:
            failures.append(f"{name}: blacklist {path.name} is empty")
            continue

        got = {style_of_contig(c) for c in contigs}
        if got != {expected}:
            failures.append(
                f"{name}: reference is {expected} but blacklist {path.name} has "
                f"{'/'.join(sorted(got))} contigs (e.g. {sorted(contigs)[0]})"
            )
        if mito != EXPECTED_MITO[expected]:
            failures.append(
                f"{name}: reference is {expected} so mito_name should be "
                f"{EXPECTED_MITO[expected]}, not {mito}"
            )
        checked += 1
        print(f"  {name:<8} {expected:<13} mito={mito:<5} {path.name}")

    if failures:
        print("\nFAILED:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        sys.exit(1)
    if checked == 0:
        sys.exit("no genome presets were checked")
    print(f"\nOK: {checked} genome presets have matching contig naming")


if __name__ == "__main__":
    main()
