#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
"""
# File name  : split_fq_bc2.py
# Create date: 2024/06/26
# Author     : xzh
# Description: 根据barcode将fastq文件拆分
"""

import gzip
import argparse
from Bio import Align
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import threading


def count_mismatches(s1, s2):
    # return sum(1 for a, b in zip(s1, s2) if a != b and (a != '-' and b != '-'))
    return sum(1 for a, b in zip(s1, s2) if a != b)


def split_fq_by_bc(input_fq, out_fq, barcode, mismatch=2):
    bc_len = len(barcode)
    aligner = Align.PairwiseAligner(
        match_score=1.0, mismatch_score=-1.0, gap_score=-1.0, extend_gap_score=-1.0
    )
    aligner.mode = "global"
    total_reads_num = 0
    index2_reads_num = 0
    lock = threading.Lock()  # 创建一个锁

    def process_record(record):
        nonlocal total_reads_num, index2_reads_num
        total_reads_num += 1
        temp_seq = str(record.seq)
        seq_bc = temp_seq[:bc_len]
        alignments = aligner.align(barcode, seq_bc)
        for alignment in alignments:
            mismatches = count_mismatches(alignment.query, alignment.target)
            if mismatches <= mismatch:
                sub_rec = record[bc_len:]  # 截取非barcode序列
                with lock:  # 使用锁来同步文件写入
                    SeqIO.write(sub_rec, f2, "fastq")
                index2_reads_num += 1
                break

    with gzip.open(input_fq, "rt") as f, gzip.open(
        out_fq, "wt"
    ) as f2, ThreadPoolExecutor(max_workers=2) as executor:
        for record in SeqIO.parse(f, "fastq"):
            executor.submit(process_record, record)

    rate = index2_reads_num / total_reads_num * 100
    return total_reads_num, index2_reads_num, rate


def main():
    parser = argparse.ArgumentParser(
        description="根据barcode将fastq文件拆分",
    )
    parser.add_argument("-i", "--input", help="输入fastq文件(gzip格式)", required=True)
    parser.add_argument("-o", "--output", help="输出fastq文件(gzip格式)", required=True)
    parser.add_argument("-b", "--barcode", help="barcode序列", required=True)
    parser.add_argument("-m", "--mismatch", help="允许的错配数", default=1)
    args = parser.parse_args()
    t, i, r = split_fq_by_bc(args.input, args.output, args.barcode, args.mismatch)
    print("#Total\tindex2_reads\tRate")
    print(f"{t}\t{i}\t{r}")


if __name__ == "__main__":
    main()
