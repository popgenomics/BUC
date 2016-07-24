#!/usr/bin/Rscript
x = read.table("output_GC_UTR.txt",h=T)
y = read.table("output_GC_ORF.txt",h=T)

z = merge(y,x,by.x=1,by.y=1)

modL = loess(as.numeric(z$nReads)~as.numeric(z$GC12))

GCcorrected_expression = residuals(modL)

z = cbind(z, GCcorrected_expression)

write.table(z, col.names=T, row.names=F, quote=F, file="output_summarized.txt", sep="\t", na = "nan")

