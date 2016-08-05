#!/usr/bin/Rscript
args = commandArgs(TRUE)

buc = read.table("../final_output_2.txt", sep="\t", header=T)
bio = read.csv("../bioWithPerry.csv")
bio = read.csv("../bio.csv")
biobuc = merge(bio,buc, by.x="focal", by.y="dataset")

attach(biobuc)
names(biobuc)

plot(GC3_avg, GC3_std, sub=cor.test(GC3_avg, GC3_std)$p.val, main=cor.test(GC3_avg, GC3_std)$estimate)

plot(log(juvSize), r_expression_ENcP, sub=cor.test(log(juvSize), r_expression_ENcP)$p.val, main=(cor.test(log(juvSize), r_expression_ENcP)$estimate)^2)
text(log(juvSize), r_expression_ENcP, sub=cor.test(log(juvSize), r_expression_ENcP)$p.val, main=(cor.test(log(juvSize), r_expression_ENcP)$estimate)^2, focal)

plot(piS, r_expression_ENcP, sub=cor.test(piS, r_expression_ENcP)$p.val, main=(cor.test(piS, r_expression_ENcP)$estimate)^2)
text(log(piS), r_expression_ENcP, sub=cor.test(log(piS), r_expression_ENcP)$p.val, main=cor.test(log(piS), r_expression_ENcP)$estimate, focal)

plot(log(longevity), r_expression_ENcP, sub=cor.test(log(longevity), r_expression_ENcP)$p.val, main=(cor.test(log(longevity), r_expression_ENcP, method="pearson")$estimate)^2, bg='grey60', pch=21, cex=1.5)
abline(lm(r_expression_ENcP~log(longevity)), col='red')
text(log(longevity), r_expression_ENcP, sub=cor.test(log(longevity), r_expression_ENcP)$p.val, main=(cor.test(log(longevity), r_expression_ENcP)$estimate)^2, focal)
dev.print(pdf,'../figures/longevityBuc.pdf')

plot(log(sexmat), r_expression_ENcP, sub=cor.test(log(sexmat), r_expression_ENcP)$p.val, main=(cor.test(log(sexmat), r_expression_ENcP, method="pearson")$estimate)^2, bg='grey60', pch=21, cex=1.5)
