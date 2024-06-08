import os
import re
import sys
from subprocess import run, PIPE
import mappy as mp
import pandas as pd
import numpy as np
import logging
import argparse
from restrion_site import restrion_site


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


def get_primer(seq, cutter1, cutter2):
    regex_ptr = re.compile(f"([ACGT]+{cutter1}|[ACGT]+{cutter2})")
    primer = re.findall(regex_ptr, seq)
    if len(primer) == 0:
        return ""
    else:
        return primer[0]


def check_primers(fq1, fq2, prim_len=15):
    logger.info("Check the primers of First 10000 reads")
    prim_df = pd.DataFrame(columns=["p1", "p2"], index=range(10000))
    for i, rd1, rd2 in zip(range(10000), mp.fastx_read(fq1), mp.fastx_read(fq2)):
        prim_df.loc[i] = sorted([rd1[1][:prim_len], rd2[1][:prim_len]])
    prim_counts = prim_df.groupby(["p1", "p2"]).size().reset_index(name="counts")
    prim_counts["prop"] = prim_counts["counts"] / prim_counts["counts"].sum()
    prim_counts = prim_counts[prim_counts["prop"] > 0.01]
    prim_counts = prim_counts.sort_values("counts", ascending=False)
    if len(prim_counts) == 0:
        logger.warning("No overpresented primer pairs found!")
    else:
        logger.info("Overpresented primer pairs found!")
        for i, row in prim_counts.iterrows():
            logger.info(f"{row['p1']} - {row['p2']}: {row['prop']:.2%}")


def extract_primer_reads(fq1, fq2, primer1, primer2, output_prefix, gzip=False):
    if len(primer1) > 20:
        trim_len1 = len(primer1) - 20
        primer1 = primer1[-20:]
    else:
        trim_len1 = 0
    if len(primer2) > 20:
        trim_len2 = len(primer2) - 20
        primer2 = primer2[-20:]
    else:
        trim_len2 = 0
    logger.info("Extract reads with primers")
    total, right_p, wrong_p = 0, 0, 0
    with open(f"{output_prefix}_1.fastq", "w") as f1, open(f"{output_prefix}_2.fastq", "w") as f2:
        for rd1, rd2 in zip(mp.fastx_read(fq1), mp.fastx_read(fq2)):
            seq1 = rd1[1][trim_len1:]
            seq2 = rd2[1][trim_len2:]
            qual1 = rd1[2][trim_len1:]
            qual2 = rd2[2][trim_len2:]
            if seq1.startswith(primer1) and seq2.startswith(primer2):
                f1.write(f"@{rd1[0]}\n{seq1}\n+\n{qual1}\n")
                f2.write(f"@{rd2[0]}\n{seq2}\n+\n{qual2}\n")
                right_p += 1
            elif seq1.startswith(primer2) and seq2.startswith(primer1):
                f1.write(f"@{rd2[0]}\n{seq2}\n+\n{qual2}\n")
                f2.write(f"@{rd1[0]}\n{seq1}\n+\n{qual1}\n")
                right_p += 1
            else:
                wrong_p += 1
            total += 1
            if total % 1000000 == 0:
                logger.info(
                    f"Processed {total:,} reads, {right_p:,} ({right_p/total:.2%}) right primers, {wrong_p:,} wrong primers"
                )
    logger.info(
        f"Processed {total:,} reads, {right_p:,} ({right_p/total:.2%}) right primers, {wrong_p:,} wrong primers"
    )
    if gzip:
        logger.info("Gzip output files")
        run(f"gzip -f {output_prefix}_1.fastq", shell=True)
        run(f"gzip -f {output_prefix}_2.fastq", shell=True)


def parse_args():
    parser = argparse.ArgumentParser(description="4C-seq analysis pipeline")
    parser.add_argument("--fastq1", type=str, required=True, help="fastq file 1")
    parser.add_argument("--fastq2", type=str, required=True, help="fastq file 2")
    parser.add_argument(
        "--primer1",
        type=str,
        required=True,
        help="PCR primer 1, should end with the restriction site",
    )
    parser.add_argument(
        "--primer2",
        type=str,
        required=False,
        help="PCR primer 2, should end with the restriction site",
    )
    parser.add_argument("--output_prefix", type=str, required=True, help="output prefix")
    parser.add_argument("--gzip", action="store_true", help="gzip output or not")
    return parser.parse_args()


def main():
    args = parse_args()
    fastq1 = args.fastq1
    fastq2 = args.fastq2
    primer1 = args.primer1
    primer2 = args.primer2
    output_prefix = args.output_prefix
    gzip = args.gzip

    check_primers(fastq1, fastq2)
    extract_primer_reads(fastq1, fastq2, primer1, primer2, output_prefix, gzip)


if __name__ == "__main__":
    main()
