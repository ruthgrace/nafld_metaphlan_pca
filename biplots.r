library(zCompositions)
library(randomcoloR)

d <- read.table("data/summary_all_count.txt", header=T, row.names=1, sep="\t",quote="",comment.char="",stringsAsFactors=FALSE)

d <- d[,(grepl("CL*", colnames(d)) | grepl("HLD*", colnames(d)))]

# for some reason the unclassified species has a bunch of dashes instead of numbers
d <- d[which(rownames(d)!="s__unclassified"),]
d.rownames <- rownames(d)
d <- apply(d,2,function(x) as.numeric(x))
rownames(d) <- d.rownames

# remove all features with zero counts for all samples
d.sum <- apply(d,1,sum)

d <- d[which(d.sum>0),]

original.data <- d

sample.sum <- apply(d,2,sum)
one.percent <- sample.sum*0.01

d.adj.zero <- t(cmultRepl(t(d),method="CZM"))

filter <- apply(d,1,function(x) length(which(x > one.percent)))
d.filter <- d.adj.zero[which(filter > 0),]

taxa.col <- data.frame(as.character(rownames(d)),rownames(d))
colnames(taxa.col) <- c("taxon","color")
taxa.col[,2] <- distinctColorPalette(length(taxa.col[,2]))

d.prop <- apply(d.adj.zero,2,function(x){x/sum(x)})
d.filter.prop <- apply(d.filter,2,function(x) {x/sum(x)})

d.clr <- t(apply(d.prop,2,function(x){log(x) - mean(log(x))}))
d.filter.clr <- t(apply(d.filter.prop,2,function(x){log(x) - mean(log(x))}))

d.pcx <- prcomp(d.clr)
d.filter.pcx <- prcomp(d.filter.clr)

conds <- data.frame(c(rep("NASH",10),rep("Healthy",10)))
colnames(conds) <- "cond"

palette=palette(c(rgb(1,0,0,0.6), rgb(0,0,1,0.6), rgb(0,1,1,0.6)))

pdf("biplots.pdf")

layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,2), heights=c(8,3))
par(mgp=c(2,0.5,0))
# make a covariance biplot of the data with compositions function
coloredBiplot(d.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.pcx$sdev[1]^2)/mvar(d.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.pcx$sdev[2]^2)/mvar(d.clr),3), sep=""),
expand=0.8,var.axes=T, scale=1, main="Biplot")
barplot(d.pcx$sdev^2/mvar(d.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

coloredBiplot(d.filter.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.filter.pcx$sdev[1]^2)/mvar(d.filter.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.filter.pcx$sdev[2]^2)/mvar(d.filter.clr),3), sep=""),
expand=0.8,var.axes=T, scale=1, main="Biplot")
barplot(d.filter.pcx$sdev^2/mvar(d.filter.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

dev.off()
