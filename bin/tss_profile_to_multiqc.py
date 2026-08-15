#!/usr/bin/env python3
"""Reshape a deeptools plotProfile table into MultiQC custom-content line-plot data.

The exact layout of plotProfile's --outFileNameData output varies between deeptools
releases, so rather than indexing fixed rows this keeps only the rows whose values are
entirely numeric and rebuilds the x-axis from the window and bin size that were requested.
"""

import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="plotProfile --outFileNameData table")
    parser.add_argument("--output", required=True, help="MultiQC custom-content TSV to write")
    parser.add_argument("--window", type=int, required=True, help="bp either side of the TSS")
    parser.add_argument("--bin-size", type=int, required=True, help="bp per bin")
    return parser.parse_args()


# deeptools writes these bookkeeping rows alongside the per-sample ones. The bin-index row is
# entirely numeric, so it has to be excluded by name rather than by content.
METADATA_LABELS = {"bins", "bin labels", "genes"}


def read_series(path):
    """Return [(sample_label, [values])] for every fully numeric per-sample row."""
    series = []
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            label = row[0].strip()
            if label.lower() in METADATA_LABELS:
                continue
            values = []
            for cell in row[1:]:
                cell = cell.strip()
                if not cell:
                    continue
                try:
                    values.append(float(cell))
                except ValueError:
                    values = []
                    break
            if label and values:
                series.append((label, values))
    return series


HEADER = [
    "# id: 'atac_tss_enrichment'",
    "# section_name: 'TSS enrichment'",
    "# description: 'Mean coverage around annotated transcription start sites. "
    "A sharp peak at 0 indicates good ATAC signal-to-noise.'",
    "# format: 'tsv'",
    "# plot_type: 'linegraph'",
    "# pconfig:",
    "#     namespace: 'ATAC'",
    "#     xlab: 'Distance from TSS (bp)'",
    "#     ylab: 'Mean coverage'",
]


def main():
    args = parse_args()
    series = read_series(args.input)

    with open(args.output, "w") as out:
        out.write("\n".join(HEADER) + "\n")
        if not series:
            print(f"WARNING: no numeric rows found in {args.input}", file=sys.stderr)
            return
        n_bins = len(series[0][1])
        positions = [str(-args.window + i * args.bin_size) for i in range(n_bins)]
        out.write("Sample\t" + "\t".join(positions) + "\n")
        for label, values in series:
            if len(values) != n_bins:
                print(
                    f"WARNING: skipping '{label}' ({len(values)} values, expected {n_bins})",
                    file=sys.stderr,
                )
                continue
            out.write(label + "\t" + "\t".join(str(v) for v in values) + "\n")


if __name__ == "__main__":
    main()
