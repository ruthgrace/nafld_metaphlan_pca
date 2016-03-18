# load the data and the colours
d.pro.0 <- read.table("bbv_probiotic_samples.txt", header=T, row.names=1)
# remove awkward values from the names
rn <- gsub("_",".", rownames(d.pro.0))
rownames(d.pro.0) <- rn
# the first two rows and three columns of the data looks like this:
d.pro.0[1:2,1:3]
## B208_bv A208_n B210_bv
## Actinobacteria:Actinomyces 1 11 8
## Actinobacteria:Arcanobacterium 1 0 2
# a correspondence table of taxa and colours
col.tax <- read.table("bbv_colours.txt", header=T, row.names=1, comment.char="")
# again, change awkward characters in the row names
rownames(col.tax) <- gsub("_",".", rownames(col.tax))
# replace 0 values with the count zero multiplicative method and output counts
#
# this function expects the samples to be in rows and OTUs to be in columns
# so the dataset is turned sideways on input, and then back again on output
# you need to know which orientation your data needs to be in for each tool
d.pro <- t(cmultRepl(t(d.pro.0), method="CZM", output="counts"))
## No. corrected values: 42
# convert to proportions by sample (columns) using the apply function
d.pro.prop <- apply(d.pro, 2, function(x){x/sum(x)})
#####
# Make a dataset where the taxon is more abundant than 0.1% in all samples
# remove all taxa that are less than 0.1\% abundant in any sample
d.pro.abund.unordered <- d.pro[apply(d.pro.prop, 1, min) > 0.001,]
# add in the names again and sort by abundance
d.names <- rownames(d.pro.abund.unordered)[
order(apply(d.pro.abund.unordered, 1, sum), decreasing=T) ]
# make a standard list of colours for plotting
colours <- as.character(col.tax[d.names,])
6
# get the taxa in the reduced dataset by name
d.pro.abund_unordered <- d.pro.abund.unordered[d.names,]
# order the taxa by their diagnosis bv, n or i
d.pro.abund <- data.frame(d.pro.abund_unordered[,grep("_bv", colnames(d.pro.abund_unordered))],
d.pro.abund_unordered[,grep("_n", colnames(d.pro.abund_unordered))],
d.pro.abund_unordered[,grep("_i", colnames(d.pro.abund_unordered))])
# make our compositional dataset
d.clr.abund <- t(apply(d.pro.abund, 2, function(x){log(x) - mean(log(x))}))
# more name plumbing!
colnames(d.clr.abund) <- gsub("\\w+:", "", colnames(d.clr.abund))

# Singlular value decompositon method of making a PCA (base R)
pcx.abund <- prcomp(d.clr.abund)
# getting info to color the samples
conds <- data.frame(c(rep(1,length(grep("_bv", rownames(d.clr.abund)))),
rep(2, length(grep("_n", rownames(d.clr.abund)))),
rep(3, length(grep("_i", rownames(d.clr.abund)))) ))
colnames(conds) <- "cond"
palette=palette(c(rgb(1,0,0,0.6), rgb(0,0,1,0.6), rgb(0,1,1,0.6)))

layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
par(mgp=c(2,0.5,0))
# make a covariance biplot of the data with compositions function
coloredBiplot(pcx.abund, col="black", cex=c(0.6, 0.6), xlabs.col=conds$cond,
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(pcx.abund$sdev[1]^2)/mvar(d.clr.abund),3), sep=""),
ylab=paste("PC2 ", round (sum(pcx.abund$sdev[2]^2)/mvar(d.clr.abund),3), sep=""),
expand=0.8,var.axes=T, scale=1, main="Biplot")
barplot(pcx.abund$sdev^2/mvar(d.clr.abund),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot