#!/usr/bin/env python3
"""Inspect the tiny synthetic paired-end FASTQ supplied for Question 2."""

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path

ADAPTER = "AGATCGGAAGAGC"


def read_fastq(path):
    records = []
    with gzip.open(path, "rt") as handle:
        while name := handle.readline().rstrip():
            seq = handle.readline().rstrip()
            plus = handle.readline().rstrip()
            qual = handle.readline().rstrip()
            if not name.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError(f"Malformed FASTQ record in {path}: {name}")
            records.append((name[1:], seq, qual))
    return records


def pair_key(name):
    return name.rsplit(":", 1)[0]


def summarize(records):
    lengths = [len(seq) for _, seq, _ in records]
    adapter_reads = sum(ADAPTER in seq for _, seq, _ in records)
    high_gc_reads = sum((seq.count("G") + seq.count("C")) / len(seq) >= 0.70 for _, seq, _ in records)
    sequence_counts = Counter(seq for _, seq, _ in records)
    exact_duplicates = len(records) - len(sequence_counts)
    largest_duplicate_group = max(sequence_counts.values())
    low_tail_reads = sum(
        sum(ord(char) - 33 for char in qual[50:]) / len(qual[50:]) <= 15
        for _, _, qual in records
    )
    return {
        "reads": len(records),
        "length": f"{min(lengths)}-{max(lengths)} bp",
        "adapter": f"{adapter_reads} ({adapter_reads / len(records):.1%})",
        "high_gc": f"{high_gc_reads} ({high_gc_reads / len(records):.1%})",
        "exact_duplicates": f"{exact_duplicates} ({exact_duplicates / len(records):.1%})",
        "largest_duplicate_group": f"{largest_duplicate_group} ({largest_duplicate_group / len(records):.1%})",
        "low_tail_quality": f"{low_tail_reads} ({low_tail_reads / len(records):.1%})",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("r1", type=Path)
    parser.add_argument("r2", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    r1, r2 = read_fastq(args.r1), read_fastq(args.r2)
    if len(r1) != len(r2):
        raise ValueError(f"Mate counts differ: R1={len(r1)}, R2={len(r2)}")
    matched = sum(pair_key(a[0]) == pair_key(b[0]) for a, b in zip(r1, r2))
    s1, s2 = summarize(r1), summarize(r2)

    rows = [
        ("Read count", s1["reads"], s2["reads"], "Direct FASTQ count"),
        ("Sequence length range", s1["length"], s2["length"], "Minimum-maximum observed length"),
        ("Matched pair IDs", f"{matched}/{len(r1)}", f"{matched}/{len(r2)}", "Mate names match after removing :1/:2"),
        ("Reads containing adapter", s1["adapter"], s2["adapter"], f"Exact {ADAPTER} match"),
        ("High-GC reads", s1["high_gc"], s2["high_gc"], "GC fraction >= 70%"),
        ("Exact duplicate observations", s1["exact_duplicates"], s2["exact_duplicates"], "Read count minus unique sequence count"),
        ("Largest identical-sequence group", s1["largest_duplicate_group"], s2["largest_duplicate_group"], "Exact sequence identity within each mate file"),
        ("Low late-cycle quality", s1["low_tail_quality"], s2["low_tail_quality"], "Mean Phred <= 15 from cycle 51 onward"),
    ]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "R1", "R2", "definition"])
        writer.writerows(rows)

    assert len(rows) == 8 and matched == len(r1), "FASTQ inspection self-check failed"
    print(f"Wrote {args.output} from {len(r1)} read pairs")


if __name__ == "__main__":
    main()
