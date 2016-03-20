library(zCompositions)
library(randomcoloR)
library(compositions)
library(ALDEx2)
library(stringr)

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

rownames(d) <- gsub("s__","",rownames(d))
rownames(d) <- sub("_"," ",rownames(d))

sample.sum <- apply(d,2,sum)
one.percent <- sample.sum*0.01

# get genus level for metaphlan
species <- rownames(d)
genus <- str_extract(species,"^[A-Za-z]*")

d.genus <- aggregate(d,list(genus),sum)
rownames(d.genus) <- d.genus$Group.1
d.genus <- d.genus[,c(2:ncol(d.genus))]
colnames(d.genus) <- str_extract(colnames(d.genus),"^[A-Z]*_[0-9]*")

# get genus level for 16S
otu.tab <- read.table("data/td_OTU_tag_mapped_lineage_working.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

colnames(otu.tab) <- gsub("[.]","_",colnames(otu.tab))

taxonomy <- otu.tab$taxonomy

otu.tab <- otu.tab[,c(1:(ncol(otu.tab)-1))]
otu.genus <- c(as.character(taxonomy))

for (i in c(1:length(taxonomy))) {
  otu.genus[i] <- strsplit(otu.genus[i],c(";"))[[1]][6]
}

otu.tab.genus <- aggregate(otu.tab,list(otu.genus),sum)
rownames(otu.tab.genus) <- otu.tab.genus$Group.1
otu.tab.genus <- otu.tab.genus[,c(2:ncol(otu.tab.genus))]

# adjust zeros
d.adj.zero <- t(cmultRepl(t(d),method="CZM"))
d.genus.adj.zero <- t(cmultRepl(t(d.genus),method="CZM"))
otu.tab.genus.adj.zero <- t(cmultRepl(t(otu.tab.genus),method="CZM"))

filter <- apply(d,1,function(x) length(which(x > one.percent)))
d.filter <- d.adj.zero[which(filter > 0),]

d.filter.counts <- d[which(filter>0),]

d.adj.zero <- d.adj.zero[order(apply(d.adj.zero,1,sum),decreasing=TRUE),]
d.filter <- d.filter[order(apply(d.filter,1,sum),decreasing=TRUE),]
d.genus.adj.zero <- d.genus.adj.zero[order(apply(d.genus.adj.zero,1,sum),decreasing=TRUE),]
otu.tab.genus.adj.zero <- otu.tab.genus.adj.zero[order(apply(otu.tab.genus.adj.zero,1,sum),decreasing=TRUE),]

d.names <- rownames(d.adj.zero)
d.filter.names <- rownames(d.filter)
d.genus.names <- rownames(d.genus.adj.zero)
otu.tab.genus.names <- rownames(otu.tab.genus.adj.zero)

taxa.col <- data.frame(as.character(rownames(d)),rownames(d))
colnames(taxa.col) <- c("taxon","color")
taxa.col[,2] <- distinctColorPalette(length(taxa.col[,2]))

taxa.filter.col <- data.frame(as.character(rownames(d.filter)),rownames(d.filter))
colnames(taxa.filter.col) <- c("taxon","color")
taxa.filter.col[,2] <- taxa.col[match(taxa.filter.col[,1],taxa.col[,1]),2]

all.genus <- unique(c(rownames(d.genus),rownames(otu.tab.genus)))
all.genus.colors <- distinctColorPalette(length(all.genus))

taxa.d.genus.col <- data.frame(rownames(d.genus),rownames(d.genus))
colnames(taxa.d.genus.col) <- c("taxon","color")
taxa.d.genus.col[,2] <- all.genus.colors[match(rownames(d.genus),all.genus)]

taxa.otu.tab.genus.col <- data.frame(rownames(otu.tab.genus),rownames(otu.tab.genus))
colnames(taxa.otu.tab.genus.col) <- c("taxon","color")
taxa.otu.tab.genus.col[,2] <- all.genus.colors[match(rownames(otu.tab.genus),all.genus)]


d.prop <- apply(d.adj.zero,2,function(x){x/sum(x)})
d.filter.prop <- apply(d.filter,2,function(x) {x/sum(x)})
d.genus.prop <- apply(d.genus.adj.zero, 2,function(x) {x/sum(x)})
otu.tab.genus.prop <- apply(otu.tab.genus.adj.zero, 2,function(x) {x/sum(x)})


d.clr <- t(apply(d.prop,2,function(x){log(x) - mean(log(x))}))
d.filter.clr <- t(apply(d.filter.prop,2,function(x){log(x) - mean(log(x))}))
d.genus.clr <- t(apply(d.genus.prop,2,function(x){log(x) - mean(log(x))}))
otu.tab.genus.clr <- t(apply(otu.tab.genus.prop,2,function(x){log(x) - mean(log(x))}))


d.pcx <- prcomp(d.clr)
d.filter.pcx <- prcomp(d.filter.clr)
d.genus.pcx <- prcomp(d.genus.clr)
otu.tab.genus.pcx <- prcomp(otu.tab.genus.clr)

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
xlabs.col=c(rep("red",10),rep("black",10)),
expand=0.8,var.axes=T, scale=1, main="Biplot")
barplot(d.pcx$sdev^2/mvar(d.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

coloredBiplot(d.filter.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.filter.pcx$sdev[1]^2)/mvar(d.filter.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.filter.pcx$sdev[2]^2)/mvar(d.filter.clr),3), sep=""),
xlabs.col=c(rep("red",10),rep("black",10)),
expand=0.8,var.axes=T, scale=1, main="Biplot")
barplot(d.filter.pcx$sdev^2/mvar(d.filter.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

coloredBiplot(d.genus.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(d.genus.pcx$sdev[1]^2)/mvar(d.genus.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(d.genus.pcx$sdev[2]^2)/mvar(d.genus.clr),3), sep=""),
xlabs.col=c(rep("red",10),rep("black",10)),
expand=0.8,var.axes=T, scale=1, main="Biplot")
barplot(d.genus.pcx$sdev^2/mvar(d.genus.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

coloredBiplot(otu.tab.genus.pcx, cex=c(0.6, 0.6),
arrow.len=0.05,
xlab=paste("PC1 ", round (sum(otu.tab.genus.pcx$sdev[1]^2)/mvar(otu.tab.genus.clr),3), sep=""),
ylab=paste("PC2 ", round (sum(otu.tab.genus.pcx$sdev[2]^2)/mvar(otu.tab.genus.clr),3), sep=""),
xlabs.col=c(rep("red",10),rep("black",10)),
expand=0.8,var.axes=T, scale=1, main="Biplot")
barplot(otu.tab.genus.pcx$sdev^2/mvar(otu.tab.genus.clr),  ylab="variance explained", xlab="Component", main="Scree plot") # scree plot

dev.off()


# generate the distance matrix
d.dist <- dist(d.clr, method="euclidian")
d.filter.dist <- dist(d.filter.clr, method="euclidian")
d.genus.dist <- dist(d.genus.clr, method="euclidian")
otu.tab.genus.dist <- dist(otu.tab.genus.clr, method="euclidian")

# cluster the data
d.hc <- hclust(d.dist, method="ward.D2")
d.filter.hc <- hclust(d.filter.dist, method="ward.D2")
d.genus.hc <- hclust(d.genus.dist, method="ward.D2")
otu.tab.genus.hc <- hclust(otu.tab.genus.dist, method="ward.D2")

# now re-order the data to plot the barplot in the same order
d.order <- d.adj.zero[,d.hc$order]
d.filter.order <- d.filter[,d.filter.hc$order]
d.genus.order <- d.genus.adj.zero[,d.genus.hc$order]
otu.tab.genus.order <- otu.tab.genus.adj.zero[,otu.tab.genus.hc$order]

d.acomp <- acomp(t(d.order))
d.filter.acomp <- acomp(t(d.filter.order))
d.genus.acomp <- acomp(t(d.genus.order))
otu.tab.genus.acomp <- acomp(t(otu.tab.genus.order))

pdf("dendogram_barplot.pdf")

layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(8,10), height=c(4,4))
par(mar=c(2,1,1,1)+0.1)
# plot the dendrogram
plot(d.hc, cex=0.6)
# plot the barplot below
barplot(d.acomp, legend.text=F, col=as.character(taxa.col[,2]), axisnames=F, border=NA, xpd=T)
par(mar=c(0,1,1,1)+0.1)
# and the legend
plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE)
legend(x="center", legend=d.names, col=as.character(taxa.col[,2]), lwd=5, cex=.3, border=NULL,ncol=3)

layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(6,4), height=c(4,4))
# plot the dendrogram
plot(d.filter.hc, cex=0.6)
# plot the barplot below
barplot(d.filter.acomp, legend.text=F, col=as.character(taxa.filter.col[,2]), axisnames=F, border=NA, xpd=T)
par(mar=c(0,1,1,1)+0.1)
# and the legend
plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE)
legend(x="center", legend=d.filter.names, col=as.character(taxa.filter.col[,2]), lwd=5, cex=.5, border=NULL)

layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(8,6), height=c(4,4))
# plot the dendrogram
plot(d.genus.hc, cex=0.6)
# plot the barplot below
barplot(d.genus.acomp, legend.text=F, col=as.character(taxa.d.genus.col[,2]), axisnames=F, border=NA, xpd=T)
par(mar=c(0,1,1,1)+0.1)
# and the legend
plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE)
legend(x="center", legend=d.genus.names, col=as.character(taxa.d.genus.col[,2]), lwd=5, cex=.5, border=NULL,ncol=2)

layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(8,6), height=c(4,4))
# plot the dendrogram
plot(otu.tab.genus.hc, cex=0.6)
# plot the barplot below
barplot(otu.tab.genus.acomp, legend.text=F, col=as.character(taxa.otu.tab.genus.col[,2]), axisnames=F, border=NA, xpd=T)
par(mar=c(0,1,1,1)+0.1)
# and the legend
plot(1,2, pch = 1, lty = 1, ylim=c(-20,20), type = "n", axes = FALSE, ann = FALSE)
legend(x="center", legend=otu.tab.genus.names, col=as.character(taxa.otu.tab.genus.col[,2]), lwd=5, cex=.5, border=NULL,ncol=2)

dev.off()

# generate the dataset by making a data frame of
d.h <- colnames(d)[grep("HLD*", colnames(d))] # Before samples
d.n <- colnames(d)[grep("CL*", colnames(d))] # After samples
d.aldex <- data.frame(d[,d.h], d[,d.n]) # make a data frame
# make the vector of set membership in the same order as
conds.aldex <- c(rep("Healthy", 10), rep("NASH", 10))
# generate 128 Dirichlet Monte-Carlo replicates
x <- aldex.clr(d.aldex, mc.samples=128, verbose=FALSE)
## [1] "operating in serial mode"
# calculate p values for each replicate and report the mean
x.t <- aldex.ttest(x, conds.aldex)
# calculate mean effect sizes
x.e <- aldex.effect(x, conds.aldex, verbose=FALSE)
## [1] "operating in serial mode"
# save it all in a data frame
x.all <- data.frame(x.e,x.t)

# generate the dataset by making a data frame of
d.filter.h <- colnames(d.filter.counts)[grep("HLD*", colnames(d.filter.counts))] # Before samples
d.filter.n <- colnames(d.filter.counts)[grep("CL*", colnames(d.filter.counts))] # After samples
d.filter.aldex <- data.frame(d.filter.counts[,d.filter.h], d.filter.counts[,d.filter.n]) # make a data frame
# make the vector of set membership in the same order as
conds.aldex <- c(rep("Healthy", 10), rep("NASH", 10))
# generate 128 Dirichlet Monte-Carlo replicates
x.filter <- aldex.clr(d.filter.aldex, mc.samples=128, verbose=FALSE)
## [1] "operating in serial mode"
# calculate p values for each replicate and report the mean
x.filter.t <- aldex.ttest(x.filter, conds.aldex)
# calculate mean effect sizes
x.filter.e <- aldex.effect(x.filter, conds.aldex, verbose=FALSE)
## [1] "operating in serial mode"
# save it all in a data frame
x.filter.all <- data.frame(x.filter.e,x.filter.t)

# generate the dataset by making a data frame of
d.genus.h <- colnames(d.genus)[grep("HLD*", colnames(d.genus))] # Before samples
d.genus.n <- colnames(d.genus)[grep("CL*", colnames(d.genus))] # After samples
d.genus.aldex <- data.frame(d.genus[,d.genus.h], d.genus[,d.genus.n]) # make a data frame
# make the vector of set membership in the same order as
conds.aldex <- c(rep("Healthy", 10), rep("NASH", 10))
# generate 128 Dirichlet Monte-Carlo replicates
x.genus <- aldex.clr(d.genus.aldex, mc.samples=128, verbose=FALSE)
## [1] "operating in serial mode"
# calculate p values for each replicate and report the mean
x.genus.t <- aldex.ttest(x.genus, conds.aldex)
# calculate mean effect sizes
x.genus.e <- aldex.effect(x.genus, conds.aldex, verbose=FALSE)
## [1] "operating in serial mode"
# save it all in a data frame
x.genus.all <- data.frame(x.genus.e,x.genus.t)

# read metadata
MyMeta<- read.table("../exponentUnifrac/data/nash_data/metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
metadata <- MyMeta[grepl("a$",rownames(MyMeta)),]
rownames(metadata) <- gsub("a$","",rownames(metadata))
samples <- str_extract(rownames(metadata), "^[A-Z]*-[0-9]*")
samples <- gsub("-","_",samples)
unique.samples <- unique(samples)
metadata <- metadata[match(unique.samples,samples),]
rownames(metadata) <- unique.samples

# conditions: Originally 0 meant steatohepatosis, and 1 meant NASH
groups <- metadata$SSvsNASH[match(colnames(otu.tab.genus),rownames(metadata))]
originalgroups <- groups

# Make healthy represented by 0, SS by 1, NASH by 2
groups <- groups + 1;
groups[which(is.na(groups))] <- 0

# make healthy 1, ss 2, nash 3 (healthy metagenomic will be 0 and nash metagenomic will be 4)
groups <- groups + 1

# mark healthy samples selected for metagenomic study
groups[which(colnames(otu.tab.genus) %in% metagenomic_samples & groups == 1)] <- 0

# mark nash samples selected for metagenomic study
groups[which(colnames(otu.tab.genus) %in% metagenomic_samples & groups == 3)] <- 4

groups[which(groups == 0)] <- "Healthy Metagenomic"
groups[which(groups == 1)] <- "Healthy"
groups[which(groups == 2)] <- "SS"
groups[which(groups == 3)] <- "NASH"
groups[which(groups == 4)] <- "NASH Metagenomic"

groups <- as.factor(groups)

# generate the dataset by making a data frame of
otu.tab.genus.h <- colnames(otu.tab.genus)[grepl("^Healthy*",as.character(groups))] # Before samples
otu.tab.genus.n <- colnames(otu.tab.genus)[grepl("^NASH*",as.character(groups))] # After samples
otu.tab.genus.aldex <- data.frame(otu.tab.genus[,otu.tab.genus.h], otu.tab.genus[,otu.tab.genus.n]) # make a data frame
# make the vector of set membership in the same order as
otu.tab.genus.conds.aldex <- c(rep("Healthy", length(otu.tab.genus.h)), rep("NASH", length(otu.tab.genus.n)))
# generate 128 Dirichlet Monte-Carlo replicates
x.otu.tab.genus <- aldex.clr(otu.tab.genus.aldex, mc.samples=128, verbose=FALSE)
## [1] "operating in serial mode"
# calculate p values for each replicate and report the mean
x.otu.tab.genus.t <- aldex.ttest(x.otu.tab.genus, otu.tab.genus.conds.aldex)
# calculate mean effect sizes
x.otu.tab.genus.e <- aldex.effect(x.otu.tab.genus, otu.tab.genus.conds.aldex, verbose=FALSE)
## [1] "operating in serial mode"
# save it all in a data frame
x.otu.tab.genus.all <- data.frame(x.otu.tab.genus.e,x.otu.tab.genus.t)


pdf("aldex_plots.pdf")

layout(matrix(c(1,2,3,1,2,3),2,3, byrow=T), widths=c(5,2,2), height=c(4,4))
par(mar=c(5,4,4,1)+0.1)
aldex.plot(x.all, test="wilcox", cutoff=0.05, all.cex=0.8, called.cex=1)
plot(x.all$effect, x.all$wi.eBH, log="y", pch=19, main="Effect",
cex=0.5, xlab="Effect size", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)
plot(x.all$diff.btw, x.all$wi.eBH, log="y", pch=19, main="Volcano",
cex=0.5, xlab="Difference", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)

aldex.plot(x.filter.all, test="wilcox", cutoff=0.05, all.cex=0.8, called.cex=1)
plot(x.filter.all$effect, x.filter.all$wi.eBH, log="y", pch=19, main="Effect",
cex=0.5, xlab="Effect size", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)
plot(x.filter.all$diff.btw, x.filter.all$wi.eBH, log="y", pch=19, main="Volcano",
cex=0.5, xlab="Difference", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)

aldex.plot(x.genus.all, test="wilcox", cutoff=0.05, all.cex=0.8, called.cex=1)
plot(x.genus.all$effect, x.genus.all$wi.eBH, log="y", pch=19, main="Effect",
cex=0.5, xlab="Effect size", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)
plot(x.genus.all$diff.btw, x.genus.all$wi.eBH, log="y", pch=19, main="Volcano",
cex=0.5, xlab="Difference", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)

aldex.plot(x.otu.tab.genus.all, test="wilcox", cutoff=0.05, all.cex=0.8, called.cex=1)
plot(x.otu.tab.genus.all$effect, x.otu.tab.genus.all$wi.eBH, log="y", pch=19, main="Effect",
cex=0.5, xlab="Effect size", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)
plot(x.otu.tab.genus.all$diff.btw, x.otu.tab.genus.all$wi.eBH, log="y", pch=19, main="Volcano",
cex=0.5, xlab="Difference", ylab="Expected Benjamini-Hochberg P")
abline(h=0.05, lty=2)

dev.off()

# COMPARE EFFECT SIZE WITH 16S

metagenomic_samples <- c("CL_119", "CL_139", "CL_141", "CL_144", "CL_160", "CL_165", "CL_166", "CL_169", "CL_173", "CL_177", "HLD_100", "HLD_102", "HLD_111", "HLD_112", "HLD_23", "HLD_28", "HLD_47", "HLD_72", "HLD_80", "HLD_85")

## sanity check to make sure all your counts have metadata
# which(!(colnames(otu.tab) %in% rownames(metadata)))

h.metnash <- otu.tab.genus
h.metnash.cond <- groups
h.metnash <- h.metnash[,which(h.metnash.cond == "Healthy Metagenomic" | h.metnash.cond == "NASH Metagenomic")]
h.metnash.cond <- h.metnash.cond[which(h.metnash.cond == "Healthy Metagenomic" | h.metnash.cond == "NASH Metagenomic")]

h.metnash.aldex <- aldex(data.frame(h.metnash),as.character(h.metnash.cond))


d.groups <- metadata$SSvsNASH[match(colnames(d.genus),rownames(metadata))]
d.originalgroups <- d.groups

d.cond <- d.groups
d.cond[which(is.na(d.cond))] <- "Healthy Metagenomic"
d.cond[which(d.cond==1)] <- "NASH Metagenomic"
d.cond <- as.factor(d.cond)

h.metnash.d <- d.genus

h.metnash.d.aldex <- aldex(data.frame(h.metnash.d),as.character(d.cond))

d.select <- d.genus
suspect.samples <- c("HLD_80","HLD_85","CL_165")
d.select <- d.genus[,which(!(colnames(d) %in% suspect.samples))]
d.select <- d.select[which(apply(d.select,1,sum)>0),]
d.select.cond <- colnames(d.select)
d.select.cond <- sub("CL.*$","NASH Metagenomic",d.select.cond)
d.select.cond <- sub("HLD.*$","Healthy Metagenomic",d.select.cond)

h.metnash.d.select.aldex <- aldex(data.frame(d.select),d.select.cond)

mycolor <- c(col2rgb("turquoise4"))
red <- mycolor[1]
green <- mycolor[2]
blue <- mycolor[3]
mycolor <- rgb(red/255, green/255, blue/255, 0.3)

pdf("metaphlan_vs_16S_effect_sizes.pdf")
common.genus <- rownames(h.metnash.d.aldex)[which(rownames(h.metnash.d.aldex) %in% rownames(h.metnash.aldex))]
otu.tab.common.effect <- h.metnash.aldex$effect[match(common.genus,rownames(h.metnash.aldex))]
d.common.effect <- h.metnash.d.aldex$effect[match(common.genus,rownames(h.metnash.d.aldex))]
plot(otu.tab.common.effect, d.common.effect, pch=19,col=mycolor, main="Effect sizes of healthy vs extreme NASH\nfor MetaPhlAn results vs. 16S sequencing",xlab="16S rRNA gene tag sequencing",ylab="MetaPhlAn results from metagenomic sequencing")
cor(otu.tab.common.effect, y = d.common.effect, use = "everything", method = "spearman")
# [1] 0.4456304

common.genus <- rownames(h.metnash.d.select.aldex)[which(rownames(h.metnash.d.select.aldex) %in% rownames(h.metnash.aldex))]
otu.tab.common.effect <- h.metnash.aldex$effect[match(common.genus,rownames(h.metnash.aldex))]
d.select.common.effect <- h.metnash.d.select.aldex$effect[match(common.genus,rownames(h.metnash.d.select.aldex))]
plot(otu.tab.common.effect, d.select.common.effect, pch=19,col=mycolor, main="Effect sizes of healthy vs extreme NASH\nfor MetaPhlAn results vs. 16S sequencing",xlab="16S rRNA gene tag sequencing",ylab="MetaPhlAn results from metagenomic sequencing")
cor(otu.tab.common.effect, y = d.select.common.effect, use = "everything", method = "spearman")
# [1] 0.5193331
dev.off()


