#!/usr/bin/env python3
"""Assemble a UCSC track hub from the pipeline's bigWigs and peak bigBeds.

Built with the `trackhub` package (https://daler.github.io/trackhub/) rather than by writing
hub.txt/genomes.txt/trackDb.txt by hand: it validates track types and parameters, and it knows
the composite/view grammar that makes a hub with two files per sample usable in the browser
instead of a flat list of 2N tracks.

The layout follows the one already in use in the lab's hand-written hubs: one composite track
per condition, coloured per condition, with the samples of that condition as subtracks sharing a
y-axis (`autoScale group`). Peaks are added as a second view inside the same composite, so a
condition's signal and peaks show and hide together.

The result is a self-contained directory:

    <outdir>/
        hub.txt
        genomes.txt
        <genome>/
            trackDb.txt
            <sample>.bw
            <sample>_peaks.bb

Copy that directory to any web server and give UCSC the URL of hub.txt.
"""

import argparse
import os
import shutil
import sys

from trackhub import CompositeTrack, Genome, GenomesFile, Hub, Track, TrackDb, ViewTrack
from trackhub.helpers import sanitize
from trackhub.upload import stage_hub

# One colour per condition, so the samples of one group read as one group in the browser. This
# is Okabe-Ito rather than the ggplot hue palette the hand-written hubs used: same purpose, but
# it stays distinguishable for red-green colour blindness.
PALETTE = [
    "0,114,178", "213,94,0", "0,158,115", "204,121,167",
    "230,159,0", "86,180,233", "240,228,66", "0,0,0",
]

BIGWIG_SUFFIXES = (".bw", ".bigWig", ".bigwig")
BIGBED_SUFFIXES = ("_narrow_peaks.bb", "_peaks.bb", ".bb", ".bigBed")


class NamedHub(Hub):
    """A Hub whose hub.txt actually carries the hub's name.

    trackhub 1.0's Hub.__str__ hardcodes the first stanza line as the literal `hub hub`
    (`("hub", "hub")` in trackhub/hub.py, where the second element should be `self.hub`), so
    every hub it renders is identified to UCSC as "hub" regardless of what it was called.
    Overriding __str__ is the least invasive fix; drop this class if the upstream bug is fixed.
    """

    def __str__(self):
        lines = Hub.__str__(self).split("\n")
        lines[0] = "hub {}".format(self.hub)
        return "\n".join(lines) + "\n"


def sample_name(path, suffixes):
    """Strip the pipeline's naming conventions off a file name to recover the sample ID."""
    name = os.path.basename(path)
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


def read_metadata(path):
    """sample -> condition, from the samplesheet-derived CSV. Missing file is not an error."""
    conditions = {}
    if not path or not os.path.exists(path):
        return conditions
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split(",")
        try:
            sample_i, condition_i = header.index("sample"), header.index("condition")
        except ValueError:
            return conditions
        for line in handle:
            fields = line.rstrip("\n").split(",")
            if len(fields) > max(sample_i, condition_i) and fields[condition_i]:
                conditions[fields[sample_i]] = fields[condition_i]
    return conditions


def dereference(root):
    """Replace the symlinks stage_hub() leaves behind with the files themselves.

    stage_hub symlinks the bigWigs and bigBeds into the staging directory, which points back
    into the Nextflow work directory. A published hub has to survive `-resume`, work-directory
    cleanup and being copied to a web server, so the data files are materialised here.
    """
    for dirpath, _dirnames, filenames in os.walk(root):
        for filename in filenames:
            path = os.path.join(dirpath, filename)
            if not os.path.islink(path):
                continue
            target = os.path.realpath(path)
            os.remove(path)
            shutil.copyfile(target, path)


def check_name_clashes(parser, names, kind):
    """Fail loudly when two names sanitize to the same track name."""
    clashes = {}
    for name in names:
        clashes.setdefault(sanitize(name), []).append(name)
    for sanitized, originals in sorted(clashes.items()):
        if len(originals) > 1:
            parser.error("{} names {} all reduce to the track name '{}'; rename them in the "
                         "samplesheet so they differ by more than punctuation"
                         .format(kind, ", ".join(originals), sanitized))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genome", required=True,
                        help="UCSC assembly name the data is aligned to, e.g. hg38 or mm10")
    parser.add_argument("--hub-name", default="atac_nf",
                        help="Hub name; whitespace becomes underscores")
    parser.add_argument("--email", required=True, help="Contact address; UCSC requires one")
    parser.add_argument("--bigwig", nargs="*", default=[], help="Per-sample bigWig files")
    parser.add_argument("--bigbed", nargs="*", default=[], help="Per-sample peak bigBed files")
    parser.add_argument("--metadata", default=None,
                        help="CSV with sample,condition,replicate; groups and colours the tracks")
    parser.add_argument("--short-label", default=None)
    parser.add_argument("--long-label", default=None)
    parser.add_argument("--url", default=None,
                        help="Base URL the hub will be served from; only used to print the "
                             "hub.txt URL to paste into UCSC")
    parser.add_argument("--outdir", default="trackhub")
    args = parser.parse_args()

    if not args.bigwig and not args.bigbed:
        parser.error("no bigWig or bigBed files given; there would be nothing in the hub")

    # UCSC hub names are identifiers, so whitespace has to go, but there is no reason to flatten
    # a readable name like "Bcell_ATACseq" to "BcellATACseq": strict=False keeps - . _ as well.
    hub_name = sanitize(args.hub_name, strict=False)
    conditions = read_metadata(args.metadata)

    signal = {sample_name(f, BIGWIG_SUFFIXES): f for f in args.bigwig}
    peaks = {sample_name(f, BIGBED_SUFFIXES): f for f in args.bigbed}
    samples = sorted(set(signal) | set(peaks))
    print("Samples in the hub: {}".format(", ".join(samples)), file=sys.stderr)

    group_of = {s: conditions.get(s, "unknown") for s in samples}
    groups = sorted(set(group_of.values()))

    # Track names have to be alphanumeric, so two samples -- or two conditions -- differing only
    # in punctuation collide into one track, and the hub silently loses one of them.
    check_name_clashes(parser, samples, "sample")
    check_name_clashes(parser, groups, "condition")

    colour_of = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(groups)}

    hub = NamedHub(
        hub=hub_name,
        short_label=args.short_label or hub_name,
        long_label=args.long_label or "{} ATAC-seq ({})".format(hub_name, args.genome),
        email=args.email,
    )
    genomes_file = GenomesFile()
    hub.add_genomes_file(genomes_file)
    genome = Genome(args.genome)
    genomes_file.add_genome(genome)
    trackdb = TrackDb()
    genome.add_trackdb(trackdb)

    # One composite per condition. Each carries two views, because a composite's subtracks have
    # to agree on type and the signal is bigWig while the peaks are bigBed.
    for group in groups:
        group_samples = [s for s in samples if group_of[s] == group]
        colour = colour_of[group]

        composite = CompositeTrack(
            name=sanitize(group),
            short_label=group,
            long_label="{} ATAC-seq".format(group),
            tracktype="bigWig",
            visibility="full",
            maxHeightPixels="100:30:10",
        )

        # autoScale group: every sample in the condition shares a y-axis, which is the only way
        # the tracks can be compared by eye.
        signal_view = ViewTrack(
            name=sanitize(group + "_signalview"),
            view="signal",
            visibility="full",
            tracktype="bigWig",
            short_label="{} signal".format(group),
            autoScale="group",
            maxHeightPixels="100:30:10",
        )
        composite.add_view(signal_view)

        peaks_view = ViewTrack(
            name=sanitize(group + "_peaksview"),
            view="peaks",
            visibility="dense",
            tracktype="bigBed 6 +",
            short_label="{} peaks".format(group),
        )
        composite.add_view(peaks_view)

        for sample in group_samples:
            if sample in signal:
                signal_view.add_tracks(Track(
                    name=sanitize(sample + "_signal"),
                    source=signal[sample],
                    # Keep the sample's own name on the staged file: the hub directory is
                    # something people browse and copy around, not just something UCSC reads.
                    filename=os.path.join(args.genome,
                                          "{}.bw".format(sanitize(sample, strict=False))),
                    short_label=sample,
                    long_label="{} ATAC signal ({})".format(sample, group),
                    tracktype="bigWig",
                    visibility="full",
                    color=colour,
                    smoothingWindow="7",
                ))
            if sample in peaks:
                peaks_view.add_tracks(Track(
                    name=sanitize(sample + "_peaks"),
                    source=peaks[sample],
                    filename=os.path.join(args.genome,
                                          "{}_peaks.bb".format(sanitize(sample, strict=False))),
                    short_label="{} peaks".format(sample),
                    long_label="{} MACS3 narrow peaks ({})".format(sample, group),
                    tracktype="bigBed 6 +",
                    color=colour,
                ))

        trackdb.add_tracks(composite)

    # render() validates the hub and writes hub.txt/genomes.txt/trackDb.txt; stage_hub() then
    # gathers those plus the data files into one directory. render() is given the output
    # directory explicitly because with no argument it renders into a temporary directory that
    # nothing ever reads.
    os.makedirs(args.outdir, exist_ok=True)
    hub.render(staging=args.outdir)
    stage_hub(hub, staging=args.outdir)
    dereference(args.outdir)

    hub_txt = os.path.join(args.outdir, "{}.hub.txt".format(hub_name))
    if not os.path.exists(hub_txt):  # older trackhub releases name it plain hub.txt
        hub_txt = os.path.join(args.outdir, "hub.txt")
    if args.url:
        print("Hub URL for UCSC: {}/{}".format(args.url.rstrip("/"),
                                               os.path.basename(hub_txt)), file=sys.stderr)
    print("Staged hub in {}".format(args.outdir), file=sys.stderr)


if __name__ == "__main__":
    main()
