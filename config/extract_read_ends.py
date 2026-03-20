#!/usr/bin/env python3

import argparse
import csv
import pysam


def get_three_prime_end(aln: pysam.AlignedSegment):
    """
    Return BED-style 1bp genomic interval representing the transcript 3' end.
    For a read aligned to + strand transcript/genome orientation:
      3' end = reference_end - 1
    For a read aligned to - strand:
      3' end = reference_start
    """
    chrom = aln.reference_name
    strand = "-" if aln.is_reverse else "+"
    if strand == "+":
        pos0 = aln.reference_end - 1
    else:
        pos0 = aln.reference_start
    pos1 = pos0 + 1
    return chrom, pos0, pos1, strand


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--condition", required=True)
    ap.add_argument("--min-mapq", type=int, default=1)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")

    with open(args.out_tsv, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow([
            "read_id", "sample_id", "condition", "chrom", "pos0", "pos1", "strand",
            "mapq", "read_length", "aligned_ref_len", "pt", "has_pt",
            "is_supplementary", "is_secondary", "is_unmapped"
        ])

        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if aln.mapping_quality < args.min_mapq:
                continue
            if aln.reference_name is None or aln.reference_end is None:
                continue

            chrom, pos0, pos1, strand = get_three_prime_end(aln)

            pt = None
            has_pt = 0
            try:
                pt = aln.get_tag("pt")
                has_pt = 1
            except KeyError:
                pt = ""

            read_length = aln.infer_read_length() or 0
            aligned_ref_len = aln.reference_length or 0

            w.writerow([
                aln.query_name,
                args.sample_id,
                args.condition,
                chrom,
                pos0,
                pos1,
                strand,
                aln.mapping_quality,
                read_length,
                aligned_ref_len,
                pt,
                has_pt,
                int(aln.is_supplementary),
                int(aln.is_secondary),
                int(aln.is_unmapped),
            ])

    bam.close()


if __name__ == "__main__":
    main()
