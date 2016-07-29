#!/usr/bin/Rscript
correlation = function(x, y){
	# function that performs cor.test with the 3 methods
	# and returns the $estimate and the $p.value
	res = NULL
	for(i in c("pearson", "spearman", "kendall")){
		tmp = cor.test(x, y, method=i)
		res = rbind(res, c(tmp$estimate, tmp$p.value))
	}
	rownames(res) = c("pearson", "spearman", "kendall")
	colnames(res) = c("coefficient", "pvalue")
	return(round(res, 5))
}

args = commandArgs(TRUE)

buc = read.table("../final_output_2.txt", sep="\t", header=T)
bio = read.csv("../life_history_traits.csv")
biobuc = merge(bio,buc, by.x="focal", by.y="dataset")

popphyl = which(biobuc$source == "This study")
biobuc = biobuc[popphyl, ]

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

# plot 4-a: longevity, correlation(ENc' ~ expression)
plot(log(longevity), r_expression_ENcP, sub=round(cor.test(log(longevity), r_expression_ENcP)$p.val, 5), main=(round(cor.test(log(longevity), r_expression_ENcP, method="kendall")$estimate^2, 5)), bg='grey60', pch=21, cex=1.5)
abline(lm(r_expression_ENcP~log(longevity)), col='red')

# plot 4-b:
inver = which(biobuc$taxonomy1=="invertebrate")
ver = which(biobuc$taxonomy1=="vertebrate")

dark = grey(0.2)
green = rgb(0, 1, 0, 0.75)
red = rgb(1, 0, 0, 0.75)
diam = 1.25 # diametre of the points

plot(log(longevity), r_expression_ENcP, col="white", xlab = "log(longevity)", ylab = "cor (gene expression, ENc')", main="", cex.lab = 1.25, cex.axis = 1.2)
points(log(longevity)[ver], r_expression_ENcP[ver], pch = 21, col = dark, bg = green, cex = diam)
points(log(longevity)[inver], r_expression_ENcP[inver], pch = 21, col = dark, bg = red, cex = diam)

legend("bottomright", c("vertebrate", "invertebrate"), pch=21, col = dark, cex = diam, pt.bg = c(green, red))

dev.print(pdf,'../figures/longevityBuc.pdf')

