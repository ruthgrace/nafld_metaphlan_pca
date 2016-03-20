library(zCompositions)
library(randomcoloR)
library(compositions)
library(ALDEx2)
library(stringr)

metagenomic_samples <- c("CL_119", "CL_139", "CL_141", "CL_144", "CL_160", "CL_165", "CL_166", "CL_169", "CL_173", "CL_177", "HLD_100", "HLD_102", "HLD_111", "HLD_112", "HLD_23", "HLD_28", "HLD_47", "HLD_72", "HLD_80", "HLD_85")

# read metadata for 16S samples
MyMeta<- read.table("../exponentUnifrac/data/nash_data/metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
metadata <- MyMeta[grepl("a$",rownames(MyMeta)),]
rownames(metadata) <- gsub("a$","",rownames(metadata))
samples <- str_extract(rownames(metadata), "^[A-Z]*-[0-9]*")
samples <- gsub("-","_",samples)
unique.samples <- unique(samples)
metadata <- metadata[match(unique.samples,samples),]
rownames(metadata) <- unique.samples

# get metaphlan data
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


sample.sum <- apply(d,2,sum)
one.percent <- sample.sum*0.01

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

# add condition onto labels of 16S hclust data
otu.tab.conditions <- as.character(groups)
otu.tab.conditions <- gsub(" ","_",otu.tab.conditions)
attributes(otu.tab.genus.dist)$Labels <- paste(otu.tab.conditions, attributes(otu.tab.genus.dist)$Labels, sep="_")

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

layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(10,6), height=c(4,4))
# plot the dendrogram
plot(otu.tab.genus.hc, cex=0.4, hang=-1)
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

# see if 16S effect sizes correspond with qPCR
x.otu.tab.genus.all <- x.otu.tab.genus.all[order(abs(x.otu.tab.genus.all$effect),decreasing=TRUE),]
head(x.otu.tab.genus.all)

# rab.all rab.win.Healthy rab.win.NASH   diff.btw diff.win
# Adlercreutzia         3.175086      3.99640882     2.340650 -1.2016842 2.523590
# Odoribacter          -3.277856     -0.07782502    -4.971920 -3.0712957 6.373980
# Escherichia-Shigella  2.756946      1.81004950     3.722981  1.7656502 3.741580
# Subdoligranulum       7.237639      7.51347381     6.921751 -0.6670286 1.681645
# Faecalibacterium      8.628496      8.77647376     8.449255 -0.4848740 1.366769
# Howardella           -4.864411     -3.57357161    -5.685705 -2.3342699 6.163914
#                     diff.btw.025 diff.btw.975 diff.win.025 diff.win.975
# Adlercreutzia           -5.256415     3.712586    0.4053187     5.447547
# Odoribacter            -15.379227     9.001643    0.8151467    16.409252
# Escherichia-Shigella    -5.311607    13.235725    0.5209411    12.426127
# Subdoligranulum         -3.925335     2.493043    0.3513319     3.955797
# Faecalibacterium        -3.114079     1.929234    0.1946130     3.316228
# Howardella             -13.326391     9.818731    0.9930043    15.988690
#                         effect effect.025 effect.975   overlap       we.ep
# Adlercreutzia        -0.4835327  -5.458831   2.521393 0.3077446 0.020343489
# Odoribacter          -0.4649340  -7.072444   3.021696 0.2764946 0.021992210
# Escherichia-Shigella  0.3956649  -2.254910   8.163925 0.2683424 0.006706071
# Subdoligranulum      -0.3719839  -4.507363   2.402492 0.3387644 0.055698109
# Faecalibacterium     -0.3610408  -5.838529   2.404305 0.3179348 0.032446750
# Howardella           -0.3532042  -4.687023   2.905916 0.3462322 0.109705221
#                        we.eBH       wi.ep    wi.eBH
# Adlercreutzia        0.3534599 0.021369055 0.3306036
# Odoribacter          0.2994716 0.009201213 0.2211472
# Escherichia-Shigella 0.2918079 0.004772191 0.2176577
# Subdoligranulum      0.3971124 0.067241972 0.3990891
# Faecalibacterium     0.3692215 0.046623776 0.3620974
# Howardella           0.4603037 0.108507134 0.4397481

h.metnash.aldex <- h.metnash.aldex[order(abs(h.metnash.aldex$effect),decreasing=TRUE),]
head(h.metnash.aldex)

# rab.all rab.win.Healthy.Metagenomic rab.win.NASH.Metagenomic
# Ruminococcus     7.6816501                    8.254224               6.45055079
# Adlercreutzia    3.0391414                    4.578077               2.32654167
# Alistipes        2.9798894                    4.037848               2.51044936
# Subdoligranulum  7.2174157                    7.419332               6.89971491
# Paraprevotella  -3.3494022                   -4.737476               0.45134788
# Sporobacter      0.5569811                    1.317570               0.08132112
# diff.btw diff.win diff.btw.025 diff.btw.975 diff.win.025
# Ruminococcus    -1.3305060 1.436077    -3.753465    1.5297000    0.1742363
# Adlercreutzia   -1.7579496 2.280052    -5.265028    3.1525008    0.2578323
# Alistipes       -1.5810377 2.630745    -5.203863    3.2883952    0.2746771
# Subdoligranulum -0.7711102 1.174769    -3.655815    0.8264022    0.2182949
# Paraprevotella   4.5098653 7.139519    -8.950650   16.1562824    0.9054702
# Sporobacter     -1.7950275 3.583661   -14.269794    2.4021013    0.4350850
# diff.win.975     effect effect.025 effect.975   overlap
# Ruminococcus        3.077371 -0.9082133 -10.360034  1.2750943 0.1578126
# Adlercreutzia       5.574447 -0.7285971  -7.970118  2.7028970 0.2375001
# Alistipes           5.785567 -0.6776068  -8.399277  2.5259105 0.2296876
# Subdoligranulum     3.257063 -0.6216454  -7.206900  0.8898567 0.2230890
# Paraprevotella     16.621770  0.5720271  -2.058568  7.9548203 0.2750001
# Sporobacter        13.317951 -0.5357797 -12.280273  1.2931115 0.2500001
#   we.ep    we.eBH       wi.ep    wi.eBH
# Ruminococcus    0.00832152 0.4587131 0.008848072 0.4150885
# Adlercreutzia   0.05264239 0.5405332 0.060235865 0.5371605
# Alistipes       0.06207430 0.5427525 0.050206252 0.5280711
# Subdoligranulum 0.03060380 0.5237721 0.052333382 0.5096242
# Paraprevotella  0.11230424 0.5885137 0.119365626 0.5934511
# Sporobacter     0.03873904 0.5324252 0.071479326 0.5453548