#!/usr/bin/env python3
"""Reshape a deeptools plotProfile table into MultiQC custom-content line-plot data.

The exact layout of plotProfile's --outFileNameData output varies between deeptools
releases, so rather than indexing fixed rows this keeps only the rows whose values are
entirely numeric and rebuilds the x-axis from the window and bin size that were requested.

Output is JSON, not the TSV custom-content format, for two reasons learned the hard way:

  * TSV gives MultiQC no way to know the x axis is numeric. It kept the header values as
    strings and sorted them lexicographically -- -10, -100, -1000, ... 990 -- so the line
    was drawn joining points in alphabetical order, producing a zigzag over the real curve.
    JSON carries genuine integers, which sort numerically.
  * The TSV header lines are parsed as YAML, where an apostrophe inside a single-quoted
    scalar ends the string and crashes MultiQC's whole custom_content module. json.dump
    escapes text properly, so no description can break the report.
"""

import argparse
import json
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="plotProfile --outFileNameData table")
    parser.add_argument("--output", required=True, help="MultiQC custom-content TSV to write")
    parser.add_argument("--window", type=int, required=True, help="bp either side of the TSS")
    parser.add_argument("--bin-size", type=int, required=True, help="bp per bin")
    parser.add_argument("--set-label", default="all",
                        help="biotype set the TSS positions were filtered to; keeps the MultiQC "
                             "section ids distinct when more than one profile is produced")
    return parser.parse_args()


# deeptools writes these bookkeeping rows alongside the per-sample ones. The bin-index row is
# entirely numeric, so it has to be excluded by name rather than by content.
METADATA_LABELS = {"bins", "bin labels", "genes"}


def is_number(cell):
    try:
        float(cell)
    except ValueError:
        return False
    return True


def read_series(path):
    """Return [(sample_label, [values])] for every per-sample row.

    plotProfile writes a variable number of *label* columns before the numbers: with one region
    file the row is `<sample>\t<region group>\t<value>...`, i.e. two. Assuming a single label
    column silently discarded every row -- the "genes" group name is not a float -- so this
    step has been emitting a header with no data since it was written. Leading non-numeric
    cells are therefore skipped, and only the sample name is kept as the label, because that is
    what MultiQC matches against the other sections.
    """
    series = []
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            label = row[0].strip()
            if label.lower() in METADATA_LABELS:
                continue
            index = 1
            while index < len(row) and not is_number(row[index].strip()):
                index += 1
            values = []
            for cell in row[index:]:
                cell = cell.strip()
                if not cell:
                    continue
                if not is_number(cell):
                    values = []
                    break
                values.append(float(cell))
            if label and values:
                series.append((label, values))
    return series


def describe(set_label):
    """Human-readable name for the gene set behind this profile."""
    if set_label == "all":
        return "all annotated genes"
    return "genes with biotype " + set_label


def main():
    args = parse_args()
    series = read_series(args.input)
    what = describe(args.set_label)

    if not series:
        print(f"WARNING: no numeric rows found in {args.input}", file=sys.stderr)

    n_bins = len(series[0][1]) if series else 0
    positions = [-args.window + i * args.bin_size for i in range(n_bins)]

    data = {}
    for label, values in series:
        if len(values) != n_bins:
            print(
                f"WARNING: skipping '{label}' ({len(values)} values, expected {n_bins})",
                file=sys.stderr,
            )
            continue
        data[label] = [[x, y] for x, y in zip(positions, values)]

    payload = {
        # Each biotype set needs its own id, or the second profile overwrites the first.
        "id": f"atac_tss_enrichment_{args.set_label}",
        "section_name": f"TSS enrichment ({what})",
        "description": (
            f"Mean coverage around the transcription start sites of {what}. A sharp peak at 0 "
            "indicates good ATAC signal-to-noise. Profiles over different gene sets are not "
            "comparable with each other: a set including pseudogenes and other biotypes with "
            "little ATAC signal gives a flatter curve than protein-coding genes alone."
        ),
        "plot_type": "linegraph",
        "pconfig": {
            "xlab": "Distance from TSS (bp)",
            "ylab": "Mean coverage",
        },
        "data": data,
    }

    with open(args.output, "w") as out:
        json.dump(payload, out, indent=1)
        out.write("\n")


if __name__ == "__main__":
    main()
