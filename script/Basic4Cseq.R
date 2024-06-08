library(Basic4Cseq)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicAlignments)
library(caTools)

args<-commandArgs(T)

libraryFile <- args[1]
liverPoints <- readPointsOfInterestFile(args[2])
bam <- args[3]

chrom = args[4]
vp = as.numeric(args[5])
window = as.numeric(args[6])
from = vp-window
to = vp+window
sample <- args[7]
outdir <- args[8]

vp_s = vp-5000
vp_e = vp+5000

zkReads <- readGAlignments(bam)
liverData = Data4Cseq(viewpointChromosome = chrom, viewpointInterval = c(vp_s,vp_e),readLength = 50, pointsOfInterest = liverPoints, rawReads = zkReads)
rawFragments(liverData)<-readsToFragments(liverData, libraryFile)
nearCisFragments(liverData)<-chooseNearCisFragments(liverData,regionCoordinates = c(from, to))
nearCisFragments(liverData)<-normalizeFragmentData(liverData)

output <- paste0(outdir,"/",sample,"_",args[6],".Basic4Cseq.pdf")
title <- sample
pdf(file=output, width=10, height=3)
visualizeViewpoint(liverData,
	mainColour = "green",
	plotTitle = title,
	loessSpan = 0.1,
	)