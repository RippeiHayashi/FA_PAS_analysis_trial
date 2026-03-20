#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assignments", nargs="+", required=True)
    ap.add_argument("--out-counts", required=True)
    ap.add_argument("--out-long", required=True)
    args = ap.parse_args()

    counts = defaultdict(int)
    samples = set()
    pas_ids = set()

    for path in args.assignments:
        with open(path) as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                sample = row["sample_id"]
                pas = row["pas_id"]
                counts[(pas, sample)] += 1
                samples.add(sample)
                pas_ids.add(pas)

    samples = sorted(samples)
    pas_ids = sorted(pas_ids)

    with open(args.out_counts, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["pas_id"] + samples)
        for pas in pas_ids:
            w.writerow([pas] + [counts[(pas, s)] for s in samples])

    with open(args.out_long, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["pas_id", "sample_id", "count"])
        for pas in pas_ids:
            for s in samples:
                w.writerow([pas, s, counts[(pas, s)]])


if __name__ == "__main__":
    main()
