#!/usr/bin/Rscript
# ./statistics.R dataset=Ciona_intestinalis_A NCBI=Ciona_intestinalis_A 
for(i in commandArgs()){
	tmp = strsplit(i, "=")
	if(tmp[[1]][1] == "dataset"){dataset = tmp[[1]][2]}
	if(tmp[[1]][1] == "NCBI"){NCBI = gsub("_", " ", tmp[[1]][2])}
}

x = read.table("output_summarized.txt", h=T)

nLoci = nrow(x)
GCtot_avg = mean(x$GCtot, na.rm = T); GCtot_std = sd(x$GCtot, na.rm = T)
GC12_avg = mean(x$GC12, na.rm = T); GC12_std = sd(x$GC12, na.rm = T)
GC3_avg = mean(x$GC3, na.rm = T); GC3_std = sd(x$GC3, na.rm = T)
GCutr_avg = mean(x$GCutr, na.rm = T); GCutr_std = sd(x$GCutr, na.rm = T)

z = t(t(x))
var1 = which(colnames(z) == "GC3")
var2 = which(colnames(z) == "GCutr")
r_expression_GC3_GCutr = cor.test(as.numeric(z[, var1]), as.numeric(z[, var2]))

var1 = which(colnames(z) == "GCcorrected_expression")
var2 = which(colnames(z) == "Nc_JN")
r_expression_ENc_JN = cor.test(as.numeric(z[, var1]), as.numeric(z[, var2]))

var1 = which(colnames(z) == "GCcorrected_expression")
var2 = which(colnames(z) == "Ncp_JN")
r_expression_ENcP_JN = cor.test(as.numeric(z[, var1]), as.numeric(z[, var2]))

res = c(dataset, NCBI)
res = c(res, nLoci)
res = c(res, round(c(GCtot_avg, GCtot_std, GC12_avg, GC12_std, GC3_avg, GC3_std, GCutr_avg, GCutr_std), 5))
res = c(res, round(c(r_expression_GC3_GCutr$estimate, r_expression_ENc_JN$estimate, r_expression_ENcP_JN$estimate), 5))

names(res) = c("dataset", "NCBIname", "nLoci", "GCtot_avg", "GCtot_std", "GC12_avg", "GC12_std", "GC3_avg", "GC3_std", "GCutr_avg", "GCutr_std", "r_GC3_GCutr", "r_expression_ENc_JN", "r_expression_ENcP_JN")
res = as.matrix(res, nrow=1)

write.table(t(res), file = "output_final_2.txt", quote = F, sep = "\t", row.names = F, col.names = T)

