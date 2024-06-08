import os
import re
import sys
from subprocess import run, PIPE
from textwrap import dedent
import mappy as mp
import pandas as pd
import numpy as np
import logging
import argparse
from restrion_site import restrion_site
from pathlib import Path


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


class UMI4C:
    def __init__(self, fq1, fq2, outdir, cutter1):
        self.fq1 = fq1
        self.fq2 = fq2
        self.outdir = outdir
        self.cutter1 = cutter1

    def get_primer(self, seq):
        regex_ptr = re.compile(f"([ACGT]+{self.cutter1}|[ACGT]+{self.cutter2})")
        primer = re.findall(regex_ptr, seq)
        if len(primer) == 0:
            return ""
        else:
            return primer[0]


def parse_args():
    parser = argparse.ArgumentParser(description="UMI4C analysis pipeline")
    parser.add_argument("--fastq1", type=str, required=True, help="fastq file 1")
    parser.add_argument("--fastq2", type=str, required=False, help="fastq file 2")
    parser.add_argument(
        "--genome", type=str, required=True, help="genome file", choices=["hg19", "mm10"]
    )
    parser.add_argument("--enzyme", type=str, required=True, help="restriction enzyme")
    parser.add_argument("--output", type=str, required=True, help="output directory")
    parser.add_argument("--sample", type=str, required=True, help="sample name")
    parser.add_argument("--vp_chr", type=str, required=True, help="chromosome of viewpoint")
    parser.add_argument("--vp_pos", type=int, required=True, help="position of viewpoint")
    return parser.parse_args()
