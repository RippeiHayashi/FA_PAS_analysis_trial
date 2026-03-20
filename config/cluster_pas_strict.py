#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict, Counter


def read_end_records(paths):
    for path in paths:
        with open(path) as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                row["pos0"] = int(row["pos0"])
                yield row


def coarse_cluster_positions(sorted_positions, max_dist=24):
    """
    Single-linkage clustering on already-filtered exact positions.
    """
    if not sorted_positions:
        return

    cluster = [sorted_positions[0]]
    last = sorted_positions[0]

    for p in sorted_positions[1:]:
        if p - last <= max_dist:
            cluster.append(p)
        else:
            yield cluster
            cluster = [p]
        last = p

    yield cluster


def find_local_peak_candidates(cluster_positions, pos_counts, half_window, min_peak_reads, min_peak_fraction):
    """
    Peak candidate = local maximum within +/- half_window.
    Keep only peaks with enough reads and enough fraction of cluster reads.
    """
    total_reads = sum(pos_counts[p] for p in cluster_positions)
    candidates = []

    for p in cluster_positions:
        c = pos_counts[p]
        if c < min_peak_reads:
            continue
        if total_reads == 0 or (c / total_reads) < min_peak_fraction:
            continue

        left = p - half_window
        right = p + half_window
        neighborhood = [q for q in cluster_positions if left <= q <= right]
        max_in_neighborhood = max(pos_counts[q] for q in neighborhood)

        # local maximum; tie broken by keeping all tied for now
        if c == max_in_neighborhood:
            candidates.append(p)

    # collapse near-duplicate candidates: if two candidates are very close,
    # keep the stronger one; if tied, keep the leftmost.
    collapsed = []
    for p in sorted(candidates):
        if not collapsed:
            collapsed.append(p)
            continue

        prev = collapsed[-1]
        if p - prev <= half_window:
            prev_count = pos_counts[prev]
            curr_count = pos_counts[p]
            if curr_count > prev_count:
                collapsed[-1] = p
            # if tie, keep prev (leftmost)
        else:
            collapsed.append(p)

    return collapsed


def assign_positions_to_peaks(cluster_positions, peak_positions):
    """
    Assign each exact position in the coarse cluster to nearest retained peak.
    """
    assignments = defaultdict(list)
    if not peak_positions:
        return assignments

    for p in cluster_positions:
        best_peak = None
        best_dist = None
        for peak in peak_positions:
            d = abs(p - peak)
            if best_dist is None or d < best_dist or (d == best_dist and peak < best_peak):
                best_peak = peak
                best_dist = d
        assignments[best_peak].append(p)

    return assignments


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True,
                    help="Input read_ends.tsv files")
    ap.add_argument("--max-dist", type=int, default=24,
                    help="Max distance for coarse clustering")
    ap.add_argument("--min-pos-reads", type=int, default=2,
                    help="Minimum pooled reads at exact position to keep it before clustering")
    ap.add_argument("--min-pos-samples", type=int, default=2,
                    help="Minimum number of samples supporting exact position to keep it before clustering")
    ap.add_argument("--peak-half-window", type=int, default=12,
                    help="+/- window for local peak calling inside coarse cluster")
    ap.add_argument("--min-peak-reads", type=int, default=3,
                    help="Minimum reads at summit position for local peak candidate")
    ap.add_argument("--min-peak-fraction", type=float, default=0.10,
                    help="Minimum fraction of coarse-cluster reads at summit for local peak candidate")
    ap.add_argument("--min-final-reads", type=int, default=10,
                    help="Minimum total assigned reads for final PAS")
    ap.add_argument("--min-final-samples", type=int, default=2,
                    help="Minimum number of samples supporting final PAS")
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    # Count exact positions
    # key = (chrom, strand, pos0)
    pos_total_counts = Counter()
    pos_sample_counts = defaultdict(Counter)

    for row in read_end_records(args.inputs):
        chrom = row["chrom"]
        strand = row["strand"]
        pos0 = row["pos0"]
        sample_id = row["sample_id"]

        key = (chrom, strand, pos0)
        pos_total_counts[key] += 1
        pos_sample_counts[key][sample_id] += 1

    # Build filtered position lists by chrom/strand
    grouped_positions = defaultdict(list)
    for (chrom, strand, pos0), total_count in pos_total_counts.items():
        n_samples = len(pos_sample_counts[(chrom, strand, pos0)])

        if total_count >= args.min_pos_reads or n_samples >= args.min_pos_samples:
            grouped_positions[(chrom, strand)].append(pos0)

    with open(args.out_tsv, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow([
            "pas_id",
            "chrom",
            "start",
            "end",
            "strand",
            "summit",
            "cluster_start",
            "cluster_end",
            "cluster_width",
            "n_reads_total",
            "n_unique_positions",
            "summit_reads",
            "n_samples_supporting_cluster",
            "n_samples_supporting_summit",
            "pct_reads_at_summit",
            "coarse_cluster_start",
            "coarse_cluster_end",
            "coarse_cluster_width",
            "n_reads_coarse_cluster",
            "n_unique_positions_coarse_cluster",
            "n_peak_candidates_coarse_cluster"
        ])

        pas_idx = 1

        for (chrom, strand), pos_list in sorted(grouped_positions.items()):
            pos_list = sorted(set(pos_list))

            # coarse clusters on filtered positions
            for coarse_positions in coarse_cluster_positions(pos_list, max_dist=args.max_dist):
                coarse_start = min(coarse_positions)
                coarse_end = max(coarse_positions) + 1
                coarse_width = coarse_end - coarse_start

                coarse_total_reads = sum(
                    pos_total_counts[(chrom, strand, p)] for p in coarse_positions
                )

                # find local peaks within coarse cluster
                pos_counts_local = {p: pos_total_counts[(chrom, strand, p)] for p in coarse_positions}
                peak_candidates = find_local_peak_candidates(
                    coarse_positions,
                    pos_counts_local,
                    half_window=args.peak_half_window,
                    min_peak_reads=args.min_peak_reads,
                    min_peak_fraction=args.min_peak_fraction
                )

                # fallback: if nothing passes, use strongest summit of coarse cluster
                if not peak_candidates:
                    peak_candidates = [
                        sorted(coarse_positions, key=lambda p: (-pos_counts_local[p], p))[0]
                    ]

                # split coarse cluster by nearest peak
                peak_to_positions = assign_positions_to_peaks(coarse_positions, peak_candidates)

                for summit, final_positions in sorted(peak_to_positions.items()):
                    final_positions = sorted(final_positions)

                    n_reads_total = sum(
                        pos_total_counts[(chrom, strand, p)] for p in final_positions
                    )
                    if n_reads_total < args.min_final_reads:
                        continue

                    # sample support for final PAS
                    sample_set = set()
                    for p in final_positions:
                        sample_set.update(pos_sample_counts[(chrom, strand, p)].keys())

                    n_samples_supporting_cluster = len(sample_set)
                    if n_samples_supporting_cluster < args.min_final_samples:
                        continue

                    cluster_start = min(final_positions)
                    cluster_end = max(final_positions) + 1
                    cluster_width = cluster_end - cluster_start
                    summit_reads = pos_total_counts[(chrom, strand, summit)]
                    n_samples_supporting_summit = len(pos_sample_counts[(chrom, strand, summit)])
                    pct_reads_at_summit = summit_reads / n_reads_total if n_reads_total > 0 else 0.0

                    pas_id = f"PAS{pas_idx:08d}"
                    pas_idx += 1

                    w.writerow([
                        pas_id,
                        chrom,
                        summit,
                        summit + 1,
                        strand,
                        summit,
                        cluster_start,
                        cluster_end,
                        cluster_width,
                        n_reads_total,
                        len(final_positions),
                        summit_reads,
                        n_samples_supporting_cluster,
                        n_samples_supporting_summit,
                        f"{pct_reads_at_summit:.4f}",
                        coarse_start,
                        coarse_end,
                        coarse_width,
                        coarse_total_reads,
                        len(coarse_positions),
                        len(peak_candidates)
                    ])


if __name__ == "__main__":
    main()
