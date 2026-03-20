#!/usr/bin/env python3

import argparse
import csv
from bisect import bisect_left
from collections import defaultdict


def load_pas_catalog(path):
    pas_by_chr_strand = defaultdict(list)
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            row["summit"] = int(row["summit"])
            pas_by_chr_strand[(row["chrom"], row["strand"])].append(row)

    for key in pas_by_chr_strand:
        pas_by_chr_strand[key].sort(key=lambda x: x["summit"])
    return pas_by_chr_strand


def load_read_ends(path):
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            row["pos0"] = int(row["pos0"])
            yield row


def nearest_pas(pos, pas_list, max_dist):
    coords = [p["summit"] for p in pas_list]
    i = bisect_left(coords, pos)

    best = None
    best_dist = None
    for j in [i - 1, i]:
        if 0 <= j < len(pas_list):
            d = abs(pos - pas_list[j]["summit"])
            if best_dist is None or d < best_dist:
                best = pas_list[j]
                best_dist = d

    if best is not None and best_dist <= max_dist:
        return best, best_dist
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--read-ends", required=True)
    ap.add_argument("--pas-catalog", required=True)
    ap.add_argument("--max-assign-dist", type=int, default=24)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    pas_by_chr_strand = load_pas_catalog(args.pas_catalog)

    with open(args.out_tsv, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow([
            "read_id", "sample_id", "condition", "chrom", "strand", "read_end",
            "pas_id", "pas_summit", "assign_dist", "pt", "mapq"
        ])

        for row in load_read_ends(args.read_ends):
            key = (row["chrom"], row["strand"])
            pas_list = pas_by_chr_strand.get(key, [])
            if not pas_list:
                continue

            pas, dist = nearest_pas(row["pos0"], pas_list, args.max_assign_dist)
            if pas is None:
                continue

            w.writerow([
                row["read_id"],
                row["sample_id"],
                row["condition"],
                row["chrom"],
                row["strand"],
                row["pos0"],
                pas["pas_id"],
                pas["summit"],
                dist,
                row["pt"],
                row["mapq"],
            ])


if __name__ == "__main__":
    main()
