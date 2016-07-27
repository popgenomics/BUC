#!/usr/bin/Rscript
args = commandArgs(TRUE)

buc = read.table("../final_output_2.txt", sep="\t", header=T)
bio = read.csv("../life_history_traits.csv")
biobuc = merge(bio,buc, by.x="focal", by.y="dataset")

attach(biobuc)
names(biobuc)

# plot 1: %GC3_avg, %GC3_std
plot(GC3_avg, GC3_std, sub=round(cor.test(GC3_avg, GC3_std)$p.val, 5), main=round(cor.test(GC3_avg, GC3_std)$estimate, 5))

# plot 2: log(juvenile's size), correlation(ENc' ~ expression)
plot(log(juvSize), r_expression_ENcP, sub=round(cor.test(log(juvSize), r_expression_ENcP)$p.val, 5), main=(round(cor.test(log(juvSize), r_expression_ENcP)$estimate^2, 5)))
text(log(juvSize), r_expression_ENcP, sub=cor.test(log(juvSize), r_expression_ENcP)$p.val, main=(cor.test(log(juvSize), r_expression_ENcP)$estimate)^2, focal)

# plot 3: piS, correlation(ENc' ~ expression)
plot(piS, r_expression_ENcP, sub=round(cor.test(piS, r_expression_ENcP)$p.val, 5), main=(round(cor.test(piS, r_expression_ENcP)$estimate^2, 5)))
text(log(piS), r_expression_ENcP, sub=cor.test(log(piS), r_expression_ENcP)$p.val, main=cor.test(log(piS), r_expression_ENcP)$estimate, focal)

# plot 4: longevity, correlation(ENc' ~ expression)
plot(log(longevity), r_expression_ENcP, sub=round(cor.test(log(longevity), r_expression_ENcP)$p.val, 5), main=(round(cor.test(log(longevity), r_expression_ENcP, method="kendall")$estimate^2, 5)), bg='grey60', pch=21, cex=1.5)
abline(lm(r_expression_ENcP~log(longevity)), col='red')

dev.print(pdf,'../figures/longevityBuc.pdf')

