library(ade4)
library(ape)

# return a color for one variable in one individual relatively to values measured in other individuals 
couleurStat = function(individual_value, all_values){
	# all_values = rbeta(1000, 1, 2); individual_value = all_values[123]
	values = seq(range(all_values, na.rm=T)[1], range(all_values, na.rm=T)[2], length.out=100)
	tmp = which(abs(individual_value-values)==min(abs(individual_value-values)))
	if(individual_value == "NaN"){return("black")}
	else{
		return(colorRampPalette(c("blue", "white", "orange"))(100)[tmp])
#		return(rev(heat.colors(100))[tmp])
	}
}

setwd("/home/croux/Documents/BUC_data_OLD_cleaned/figures")

# tree made from NCBI
phyl = read.tree("phyliptree.phy") 

# resBUC.txt -> final.txt modified with an aditionnal column [,2] corresponding to names in [,1] that are recognized by NCBI
x = read.table("/home/croux/Documents/BUC_data_OLD_cleaned/final.txt", h=T, sep="\t")

# transform names in [,2] in order to exactly fit with names in the tree.
# the read.tree function removes " " and add some "'" around label tips
newNames = NULL
for(i in x[,2]){
	tmp = gsub("'", "", as.character(i))
	newNames = c(newNames, paste("'",gsub(" ", "", tmp),"'", sep=""))
}

x[,2] = newNames

# rename the tip labels with x[,1]
for(i in 1:length(phyl$tip.label)){
	tmp = which(x[,2] == phyl$tip.label[i])
	phyl$tip.label[i] = as.character(x[tmp, 1])

}


x11()
# plot the tree
par(mfrow=c(1,2), mar=c(0,4,0,0))
plot.phylo(phyl, use.edge.length=F, cex=0.5, align.tip.label=T)

# plot the statistics
plot.new()
par(mar=c(0,0,0,0))
step = 1/length(phyl$tip.label)

nStats = 4 # number of statistics to show
stats = 1/nStats
for(i in 1:length(phyl$tip.label)){ # loop over the 1:nIndividuals
	tmp = which(x[,1]==phyl$tip.label[i])
#	rect(0, (i-1)*step, 1, i*step, col=couleurStat(x$nLoci[tmp], x$nLoci))
	rect(0*stats, (i-1)*step, 1*stats, i*step, col=couleurStat(x$GC3_avg[tmp], x$GCutr_avg), border=F) # GC3_mean
	rect(1*stats, (i-1)*step, 2*stats, i*step, col=couleurStat(x$GC3_std[tmp], x$GC3_std), border=F) # GC3_std
	rect(2*stats, (i-1)*step, 3*stats, i*step, col=couleurStat(x$tau.nReads_ENcPFDS[tmp], x$tau.expression_ENcPFDS), border=F) # one stat of BUC
	rect(3*stats, (i-1)*step, 4*stats, i*step, col=couleurStat(x$rho.nReads_ENcPFDS[tmp], x$r.expression_ENcPFDS), border=F) # one stat of BUC
}

