#!/usr/bin/Rscript
x = read.table("output_GC_UTR.txt", h=T)
y = read.table("output_GC_ORF.txt", h=T)
z = read.table("output_ENCprime_JoNo.txt", h=T)

res = merge(y, x, by.x=1, by.y=1)

modL = loess(as.numeric(res$nReads)~as.numeric(res$GC12))

GCcorrected_expression = residuals(modL)

res = cbind(res, GCcorrected_expression)

res = merge(res, z, by.x=1, by.y=1)

write.table(res, col.names=T, row.names=F, quote=F, file="output_summarized.txt", sep="\t", na = "nan")

