# python pcr_qc.py \
# --fastq1 ../data/raw/test_1.fastq.gz \
# --fastq2 ../data/raw/test_2.fastq.gz \
# --primer1 GGAAATAATGAGCCACATTCATG \
# --primer2 ACTCCCTTTTTATCCCAAACGTTCGTAAATTTTGTATCTGATAAAGAGCATACTTCCATCTAATACAAATATGTTCCCCCCTTCAGATCTTCTCAGCATTCGAGAGATCTGTAC \
# --output_prefix ../data/interim/test \
# --gzip

python ./regular_4C.py \
--fastq1 ../data/interim/test_1.fastq.gz \
--fastq2 ../data/interim/test_2.fastq.gz \
--genome hg19 \
--bwaidx ~/REF/hg19/Sequence/BWAIndex/genome.fa \
--enzyme1 NlaIII \
--enzyme2 CviQI \
--output ../data/processed/test \
--sample test \
--vp_chr chr9 \
--vp_pos 21973773