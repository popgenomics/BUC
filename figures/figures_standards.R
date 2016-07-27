library(ggplot2)
library(ggrepel)

x=read.csv("final_output.csv", sep="\t", h=T)
x = as.data.frame(x)

# add a stat
#stdGC3_over_avgGC3 = x$GC3_std/x$GC3_avg
#x = cbind(x, stdGC3_over_avgGC3)

# remove dataset with too much or too few loci
x = x[x$nLoci > 2000 & x$nLoci < 30000, ]

ggplot(x, aes(r_GC3_GCutr, r_expression_ENcP)) + geom_point(color = "red") +
geom_text_repel(aes(label = x$dataset), size = 3) + geom_label_repel(aes(fill = x$taxonomy1), colour = "white", fontface = "bold") +
theme_classic(base_size = 16)

dev.print(pdf, "figure_isochores_BUC_phylogenie_2D.pdf",bg="white")
dev.off()

ggplot(x, aes(r_GC3_GCutr, r_expression_ENcP)) + geom_point(color = "red") +
geom_text_repel(aes(label = x$dataset), size = 3) + geom_label_repel(aes(fill = x$taxonomy1), colour = "white", fontface = "bold") +
theme_classic(base_size = 16)


p = ggplot(x, aes(r_GC3_GCutr, r_expression_ENcP, label = x$dataset))
p + geom_point(aes(colour = factor(taxonomy1))) + geom_text_repel(aes(colour = factor(taxonomy1))) + theme_classic(base_size = 16)
p + geom_point(aes(colour = factor(taxonomy2)), size = 2) + theme_classic(base_size = 16)

p = ggplot(x, aes(GC3_std, r_expression_ENcP, label = x$dataset))
p + geom_point(aes(colour = factor(taxonomy2)), size = 2) + theme_classic(base_size = 16)


