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


def run_shell_cmd(cmd):
    logger.info(f"Run shell command: {cmd}")
    res = run(cmd, shell=True, stdout=PIPE, stderr=PIPE, check=True)
    if res.returncode != 0:
        logger.error(res.stderr.decode())
        sys.exit(1)


class PeakC:
    def __init__(
        self, fastq1, fastq2, genome, bwaidx, cutter1, cutter2, output, sample, vp_chr, vp_pos
    ):
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.genome = genome
        self.bwaidx = bwaidx
        self.cutter1 = cutter1
        self.cutter2 = cutter2
        self.output = output
        self.sample = sample
        self.vp_chr = vp_chr
        self.vp_pos = vp_pos
        self.current_script_path = os.path.abspath(__file__)
        self.script_dir = Path(self.current_script_path).parent
        self.interim_dir = self.script_dir.parent / "data" / "interim" / sample
        if not self.interim_dir.exists():
            self.interim_dir.mkdir(parents=True)
        self.processed_dir = self.script_dir.parent / "data" / "processed" / sample
        if not self.processed_dir.exists():
            self.processed_dir.mkdir(parents=True)

    def digest_read(self, name, seq, qual):
        cutter1, cutter2 = self.cutter1, self.cutter2
        if cutter1 in seq or cutter2 in seq:
            digest_seq = []
            frg_be = 0
            ith = 0
            for res_enz in re.finditer("|".join([cutter1, cutter2]), seq):
                if ith == 0:
                    ith += 1
                    frg_en = res_enz.end()
                    frg_be = res_enz.start()
                    continue
                frg_en = res_enz.end()
                if len(seq[frg_be:frg_en]) >= 20:
                    digest_seq.append(f"@{name}\n{seq[frg_be:frg_en]}\n+\n{qual[frg_be:frg_en]}\n")
                frg_be = res_enz.start()
                ith += 1
            if len(seq[frg_be:]) >= 20:
                digest_seq.append(f"@{name}\n{seq[frg_be:]}\n+\n{qual[frg_be:]}\n")
            return "\n".join(digest_seq)
        else:
            return ""

    def dist_fq(self):
        logger.info(f"digest fastq of {self.sample}")
        dist_fq = open(self.interim_dir / f"{self.sample}.dist.fq", "w")
        for fq in [self.fastq1, self.fastq2]:
            for read in mp.fastx_read(fq):
                bwa_read = self.digest_read(read[0], read[1], read[2])
                dist_fq.write(dedent(bwa_read))
        dist_fq.close()

    def align(self):
        logger.info(f"align fastq of {self.sample} using bwa samse")
        cmd = f'bwa aln -t 10 {self.bwaidx} {self.interim_dir / f"{self.sample}.dist.fq"} > {self.interim_dir / f"{self.sample}.sai"}'
        run_shell_cmd(cmd)
        cmd = f'bwa samse {self.bwaidx} {self.interim_dir / f"{self.sample}.sai"} {self.interim_dir / f"{self.sample}.dist.fq"} > {self.interim_dir / f"{self.sample}.sam"}'
        run_shell_cmd(cmd)
        cmd = f'samtools view -bS {self.interim_dir / f"{self.sample}.sam"} > {self.interim_dir / f"{self.sample}.bam"}'
        run_shell_cmd(cmd)
        cmd = f'samtools sort {self.interim_dir / f"{self.sample}.bam"} -o {self.interim_dir / f"{self.sample}.sorted.bam"}'
        run_shell_cmd(cmd)
        cmd = f'samtools index {self.interim_dir / f"{self.sample}.sorted.bam"}'
        run_shell_cmd(cmd)
        cmd = f'samtools view -b {self.interim_dir / f"{self.sample}.sorted.bam"} {self.vp_chr} > {self.interim_dir / f"{self.sample}.sorted_{self.vp_chr}.bam"}'
        run_shell_cmd(cmd)
        cmd = f'samtools index {self.interim_dir / f"{self.sample}.sorted_{self.vp_chr}.bam"}'
        run_shell_cmd(cmd)

    def bamcoverage(self):
        logger.info(f"bamcoverage of {self.sample}")
        cmd = f'bamCoverage -b {self.interim_dir / f"{self.sample}.sorted_{self.vp_chr}.bam"} -o {self.interim_dir / f"{self.sample}.bw"} -r {self.vp_chr}'
        run_shell_cmd(cmd)
        cmd = f'bigWigToWig {self.interim_dir / f"{self.sample}.bw"} {self.interim_dir / f"{self.sample}.wig"}'
        run_shell_cmd(cmd)
        df = pd.read_csv(
            self.interim_dir / f"{self.sample}.wig", sep="\t", names=range(4), skiprows=1
        )
        df[[1, 3]].to_csv(
            self.processed_dir / f"{self.sample}.4C.txt", sep="\t", index=False, header=False
        )

    def plot(self, window=1000000):
        logger.info(f"PeakC plot of {self.sample}")
        cmd = f"Rscript {self.script_dir}/peakC.R {self.sample} {self.processed_dir} {self.processed_dir}/{self.sample}.4C.txt {self.vp_pos} {window}"
        run_shell_cmd(cmd)
        logger.info(f"results are saved in {self.processed_dir}")

    def peakc_pipe(self):
        self.dist_fq()
        self.align()
        self.bamcoverage()
        self.plot()


class FourCSeq:
    def __init__(
        self,
        fastq1,
        genome,
        bwaidx,
        enzyme1,
        enzyme2,
        cutter1,
        cutter2,
        output,
        sample,
        vp_chr,
        vp_pos,
    ):
        self.fastq1 = fastq1
        self.genome = genome
        self.bwaidx = bwaidx
        self.enzyme1 = enzyme1
        self.enzyme2 = enzyme2
        self.cutter1 = cutter1
        self.cutter2 = cutter2
        self.output = output
        self.sample = sample
        self.vp_chr = vp_chr
        self.vp_pos = vp_pos
        self.current_script_path = os.path.abspath(__file__)
        self.script_dir = Path(self.current_script_path).parent
        self.interim_dir = self.script_dir.parent / "data" / "interim" / sample
        self.cwd = os.getcwd()
        if not self.interim_dir.exists():
            self.interim_dir.mkdir(parents=True)
        self.processed_dir = self.script_dir.parent / "data" / "processed" / sample
        if not self.processed_dir.exists():
            self.processed_dir.mkdir(parents=True)
        if self.genome == "hg19":
            self.species = "Homo_sapiens"
            self.track_root = "trackdb/hg19_trackdb"
        elif self.genome == "mm9":
            self.species = "Mus_musculus"
            self.track_root = "trackdb/mm9_trackdb"
        else:
            logger.error("Invalid genome build, use hg19 or mm9")
            exit(1)

    def trim_fq(self, keep_len=36):
        logger.info(f"trim fastq of {self.sample}")
        outfq = open(self.interim_dir / f"{self.sample}.trim.fq", "w")
        for read in mp.fastx_read(self.fastq1):
            if len(read[1]) >= keep_len:
                outfq.write(f"@{read[0]}\n{read[1][:keep_len]}\n+\n{read[2][:keep_len]}\n")
        outfq.close()

    def write_4Cseqconf(self):
        conf = f"""\
                index=index.txt
                trackdb_root={self.track_root}
                trackset_name=4C
                binsize=16
                rawdir=rawdata
                Rscript_path=Rscript"""
        with open(f"{self.script_dir}/4Cseqpipe/4cseqpipe.conf", "w") as f:
            f.write(dedent(conf))

    def write_4Cseqindex(self):
        all_index = pd.read_csv(f"{self.script_dir}/4Cseqpipe/index.txt", sep="\t")
        if self.sample in all_index["exp"].values:
            i = all_index[all_index["exp"] == self.sample].iloc[0, 0]
            logger.info(f"{self.sample} already in index.txt, id={i}")
        else:
            i = all_index["id"].max() + 1
            for read in mp.fastx_read(self.fastq1):
                seq = read[1]
                break
            used_prim = seq.split(self.cutter1)[0] + self.cutter1
            index = {
                "id": i,
                "run": self.sample,
                "lane_no": 0,
                "exp": self.sample,
                "primer_seq": used_prim,
                "species_name": self.species,
                "first_cutter_name": self.enzyme1,
                "first_cutter_seq": self.cutter1,
                "sec_cutter_name": self.enzyme2,
                "sec_cutter_seq": self.cutter2,
                "linearization_name": None,
                "linearization_seq": None,
                "bait_chromo": self.vp_chr[3:],
                "bait_coord": self.vp_pos,
                "seq_len": 36,
                "raw_fname": f"{self.sample}.txt",
            }
            all_index.loc[i] = index
            all_index.fillna("NA").to_csv(
                f"{self.script_dir}/4Cseqpipe/index.txt", sep="\t", index=False
            )

    def run_4Cseq(self):
        if not os.path.exists(
            f"{self.script_dir}/4Cseqpipe/{self.track_root}/tracks/re/{self.cutter1}_dist2{self.cutter2}_3prime"
        ):
            os.chdir(f"{self.script_dir}/4Cseqpipe")
            logger.info("build restriction enzyme database")
            cmd = f"""perl 4cseqpipe.pl \
            -build_re_db -first_cutter {self.cutter1} \
            -second_cutters {self.cutter2} \
            -trackdb_root {self.track_root}"""
            run_shell_cmd(cmd)
            os.chdir(self.cwd)

        os.chdir(f"{self.script_dir}/4Cseqpipe")
        all_index = pd.read_csv(f"{self.script_dir}/4Cseqpipe/index.txt", sep="\t")
        sample_id = all_index[all_index["exp"] == self.sample].iloc[0, 0]
        cmd = f"""perl 4cseqpipe.pl \
        -fastq2raw \
        -ids {sample_id} \
        -fastq_fn {self.interim_dir / f"{self.sample}.trim.fq"} \
        -convert_qual 1"""
        run_shell_cmd(cmd)
        cmd = f"""perl 4cseqpipe.pl \
        -map \
        -ids {sample_id}"""
        run_shell_cmd(cmd)
        os.chdir(self.cwd)

    def plot(self, window=1000000):
        logger.info(f"plot 4C-seq of {self.sample}")
        os.chdir(f"{self.script_dir}/4Cseqpipe")
        all_index = pd.read_csv(f"{self.script_dir}/4Cseqpipe/index.txt", sep="\t")
        sample_id = all_index[all_index["exp"] == self.sample].iloc[0, 0]
        cmd = f"""perl 4cseqpipe.pl \
        -nearcis \
        -calc_from {self.vp_pos-window} \
        -calc_to {self.vp_pos+window} \
        -stat_type median \
        -trend_resolution 5000 \
        -ids {sample_id} \
        -figure_fn {self.processed_dir}/{self.sample}_{window}.4cseq.pdf"""
        run_shell_cmd(cmd)
        os.chdir(self.cwd)

    def fourcseqpipe(self):
        self.trim_fq()
        self.write_4Cseqconf()
        self.write_4Cseqindex()
        self.run_4Cseq()
        self.plot()
        logger.info(f"results are saved in {self.processed_dir}")


class Basic4C:
    def __init__(
        self,
        fastq1,
        fastq2,
        genome,
        bwaidx,
        enzyme1,
        enzyme2,
        cutter1,
        cutter2,
        output,
        sample,
        vp_chr,
        vp_pos,
    ):
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.genome = genome
        self.bwaidx = bwaidx
        self.enzyme1 = enzyme1
        self.enzyme2 = enzyme2
        self.cutter1 = cutter1
        self.cutter2 = cutter2
        self.output = output
        self.sample = sample
        self.vp_chr = vp_chr
        self.vp_pos = vp_pos
        self.current_script_path = os.path.abspath(__file__)
        self.script_dir = Path(self.current_script_path).parent
        self.interim_dir = self.script_dir.parent / "data" / "interim" / sample
        self.cwd = os.getcwd()
        if not self.interim_dir.exists():
            self.interim_dir.mkdir(parents=True)
        self.processed_dir = self.script_dir.parent / "data" / "processed" / sample
        if not self.processed_dir.exists():
            self.processed_dir.mkdir(parents=True)
        self.ref_dir = self.script_dir.parent / "data" / "ref"
        if not self.ref_dir.exists():
            self.ref_dir.mkdir(parents=True)
        if self.genome == "hg19":
            self.species = "Homo_sapiens"
        elif self.genome == "mm9":
            self.species = "Mus_musculus"
        else:
            logger.error("Invalid genome build, use hg19 or mm9")
            exit(1)

    def make_bait_bed(self):
        if not os.path.exists(f"{self.ref_dir}/{self.cutter1}_{self.cutter2}.txt"):
            cmd = f"""Rscript \
            {self.script_dir}/Basic4Cseq_digest.R {self.cutter1} {self.cutter2} \
            {self.ref_dir}/{self.cutter1}_{self.cutter2}.txt {self.species}"""
            run_shell_cmd(cmd)

    def align(self):
        if not os.path.exists(f"{self.interim_dir}/{self.sample}.sorted_{self.vp_chr}.bam"):
            peakc = PeakC(
                self.fastq1,
                self.fastq2,
                self.genome,
                self.bwaidx,
                self.cutter1,
                self.cutter2,
                self.output,
                self.sample,
                self.vp_chr,
                self.vp_pos,
            )
            peakc.peakc_pipe()
        else:
            logger.info(f"{self.interim_dir}/{self.sample}.sorted_{self.vp_chr}.bam already exists")

    def plot(self, window=1000000):
        logger.info(f"basic4C plot of {self.sample}")
        cmd = f'echo "{self.vp_chr}	{self.vp_pos}	{self.vp_pos}	VP	black" > {self.interim_dir}/vp.bed'
        run_shell_cmd(cmd)
        cmd = f"""Rscript {self.script_dir}/Basic4Cseq.R \
        {self.ref_dir}/{self.cutter1}_{self.cutter2}.txt \
        {self.interim_dir}/vp.bed {self.interim_dir}/{self.sample}.sorted_{self.vp_chr}.bam \
        {self.vp_chr} {self.vp_pos} {window} {self.sample} {self.processed_dir}"""
        run_shell_cmd(cmd)

    def basic4c_pipe(self):
        self.make_bait_bed()
        self.align()
        self.plot()
        logger.info(f"results are saved in {self.processed_dir}")


def parse_args():
    parser = argparse.ArgumentParser(description="4C-seq analysis pipeline")
    parser.add_argument("--fastq1", type=str, required=True, help="fastq file 1")
    parser.add_argument("--fastq2", type=str, required=False, help="fastq file 2")
    parser.add_argument(
        "--genome", type=str, required=True, help="genome file", choices=["hg19", "mm10"]
    )
    parser.add_argument("--bwaidx", type=str, required=False, help="bwa index file")
    parser.add_argument("--enzyme1", type=str, required=True, help="restriction enzyme 1")
    parser.add_argument("--enzyme2", type=str, required=True, help="restriction enzyme 2")
    parser.add_argument("--output", type=str, required=True, help="output directory")
    parser.add_argument("--sample", type=str, required=True, help="sample name")
    parser.add_argument("--vp_chr", type=str, required=True, help="chromosome of viewpoint")
    parser.add_argument("--vp_pos", type=int, required=True, help="position of viewpoint")
    return parser.parse_args()


def main():
    args = parse_args()
    fastq1 = args.fastq1
    fastq2 = args.fastq2
    genome = args.genome
    bwaidx = args.bwaidx
    enzyme1 = args.enzyme1
    enzyme2 = args.enzyme2
    output = args.output
    sample = args.sample
    vp_chr = args.vp_chr
    vp_pos = args.vp_pos

    if enzyme1 in restrion_site and enzyme2 in restrion_site:
        cutter1 = restrion_site[enzyme1]
        cutter2 = restrion_site[enzyme2]
    else:
        logger.error("Invalid enzyme")
        sys.exit(1)

    peakc = PeakC(fastq1, fastq2, genome, bwaidx, cutter1, cutter2, output, sample, vp_chr, vp_pos)
    peakc.peakc_pipe()
    fourcseq = FourCSeq(
        fastq1, genome, bwaidx, enzyme1, enzyme2, cutter1, cutter2, output, sample, vp_chr, vp_pos
    )
    fourcseq.fourcseqpipe()
    basic4c = Basic4C(
        fastq1, fastq2, genome, bwaidx, enzyme1, enzyme2, cutter1, cutter2, output, sample, vp_chr, vp_pos
    )
    basic4c.basic4c_pipe()


if __name__ == "__main__":
    main()
