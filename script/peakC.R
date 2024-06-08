library(peakC)
args<-commandArgs(T)

sample <- args[1]
outdir <- args[2]
file <- args[3]
vp <- as.numeric(args[4])
window <- as.numeric(args[5])

id <- sample
output <- paste0(sample,"_",args[5],".peakC.pdf")
data <- readqWig(file.path(file), vp.pos = vp,window=window)
res <- single.analysis(data=data[[1]],vp.pos= vp)
pdf(file.path(outdir,output),height = 3,width = 10)
plot_C(res,num.exp = 0, y.min = 0, y.max = 2000)
dev.off()