library(dada2)
library(ShortRead)
library(Biostrings)
library(rlang)
install.packages("rlang", version = "1.1.0")

path <- "C:/Users/vredko/Desktop/NGS/16S/All_in_one/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
# print(sample.names)

plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))

# Dapters used by Genomed
# CCTACGGGAGGCAGCAG -- F_adapt
# CAGTGCCTTTGCAACTGT -- R_adapt

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = 2,minQ = 20,# trimLeft = 40,
                     truncQ = 2, minLen = 50, rm.phix = TRUE, 
                     compress = F, multithread = FALSE) # simillar results were gained with minQ - 15!!!
head(out)


errF <- learnErrors(filtFs, multithread=T)
errR <- learnErrors(filtRs, multithread=T)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


dadaFs <- dada(filtFs, err=errF, multithread=T, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=T, pool="pseudo")

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=F) 

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(mergers, rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
# write.table(track,"C:/Users/vredko/Desktop/NGS/Drugi_poziom_fastq/trimmed/Drugi_poziom.txt",sep="\t",col.names = NA)

taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/vredko/Desktop/NGS/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE) # Use it for 18S

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# write.csv(taxa.print, "C:/Users/vredko/Desktop/NGS/Drugi_poziom_fastq/Drugi_poziom_fastq_taxa.csv", row.names=FALSE)

# write.csv(taxa, "C:/Users/vredko/Desktop/NGS/Drugi_poziom_fastq/Drugi_poziom_fastq_taxa_with_sqe.csv", row.names=TRUE)

#-----------------------#
# https://yulab-smu.top/treedata-book/chapter4.html
# https://joey711.github.io/phyloseq/plot_tree-examples.html
library(phyloseq); 
packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library("gridExtra")
library(ggtree)
library("knitr")
library("BiocStyle")
library(DECIPHER)
library(phangorn)
library(ape)
library(scales)
library(phyloseqGraphTest)
library(igraph)
library(gridExtra)
library(microbiome)
library(vegan)

theme_set(theme_bw())
sequences<-getSequences(seqtab.nochim)
names(sequences)<-sequences
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
# map <- import_qiime_sample_data("~/DADA2_Tutorial/Tutorial_Map.txt")
# Abundance normalization#: Total sum scaling (TSS) - was performed/Relative log expression (RLE)/Cumulative sum scaling (CSS)
samples.out <- rownames(seqtab.nochim)
samdf <- data.frame(Sample = samples.out)
rownames(samdf) <- samples.out
samdf

otu_table <- apply(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 2, function(x) x/sum(x))
ps2<- phyloseq(transform_sample_counts(otu_table(seqtab.nochim,taxa_are_rows=FALSE), function(x) x/sum(x)),
               tax_table(taxa),
               phy_tree(fitGTR$tree))
ps2

# The ps(phyloseq file) works for sure
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa),phy_tree(fitGTR$tree),
               sample_data(samdf))
ps

ps %>% 
  psmelt %>%
  group_by(Phylum) %>% 
  summarize(OTUn = (unique(OTU) %>% length)) %>% 
  arrange(desc(OTUn))

ntaxa(ps)

ordinate(ps, "PCoA", "bray") %>% 
  plot_ordination(ps, ., color = "Sample", title = "Bray-Curtis")



plot_bar(ps2, fill = "Phylum")
plot_bar(ps, fill = "Phylum")
plot_bar(ps, fill = "Genus")
plot_bar(ps, x="Sample", fill="Phylum" , facet_grid= ~ Phylum)
plot_bar(ps2, x="Sample", fill="Phylum" , facet_grid= ~ Phylum)


# plot_tree(ps, nodelabf=nodeplotblank, ladderize="left", color="Sample")
# plot_tree(ps, color="Phylum", ladderize="left") + coord_polar(theta="y")
# plot_tree(ps, nodelabf=nodeplotblank, color="Genus", ladderize="left") + coord_polar(theta="y")

ps                      
rank_names(ps)

#---------------------#
# Taxonomy tree - try to modify this
ggtree(ps2) + 
  aes(color= Sample) + geom_tippoint(aes(size= Abundance), alpha=.6) + 
  geom_tiplab(aes(label= Genus), offset=.05, size=1.75) +
  xlim(0,3) 

ggtree(ps) + 
  aes(color= Sample) + geom_tippoint(aes(size= Abundance), alpha=.6) + 
  geom_tiplab(aes(label= Genus), offset=.05, size=1.75) +
  xlim(0,3)
ggsave("C:/Users/vredko/Desktop/NGS/ITS/All_in_one/First_tour/tree_large2_phylum_normalized.png", width = 50, height = 50, units = "cm", limitsize = FALSE)

ggtree(ps2) + 
  aes(color= Sample) + geom_tippoint(aes(size= Abundance), alpha=.6) + 
  geom_tiplab(aes(label= Genus), offset=.05, size=1.75) + coord_polar(theta="y") +
  xlim(0,3) 

physeq = prune_taxa(taxa_names(ps)[1:500], ps)

GPUF <- UniFrac(ps)
ps_pcoa = ordinate(ps, method="PCoA", distance=GPUF)
plot_scree(ps_pcoa, "Scree plot for Global Patterns, UniFrac/PCoA")



# Nest one app.
windows(height=8,width=8)
ggtree(ps) + 
  geom_nodelab(aes(label=Genus), hjust=-.05, size=3) +
  geom_point(aes(x=x+hjust, color=Sample, shape= Sample, 
                 size=Abundance), na.rm=TRUE) +
  geom_tiplab(aes(label=Genus), hjust=-.85, size = 1) +                
  scale_size_continuous(trans=log_trans(10)) +
  theme(legend.position="right") + hexpand(.4) + xlim(0,3)
ggsave("C:/Users/vredko/Desktop/NGS/ITS/All_in_one/est.jpg", width = 50, height = 50, units = "cm", limitsize = FALSE)

#--------------analysis with R and phyloseq---------------#

library("phyloseq")
library("ggplot2")
library("vegan")

remove.packages("rlang")
remove.packages("phyloseq")

install.packages("rlang")
install.packages("phyloseq")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rlang")
library(rlang)
library(phyloseq)

rarecurve(t(otu_table(ps.rarefied)), step=50, cex=10.5)
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=0.9*min(sample_sums(ps)), replace=F)
ps.rarefied
plot_bar(ps, fill= "Phylum") + facet_wrap(~Phylum, scales = "free_x", nrow = 1)
#---------------------Alpha deversity------------------#
# Disclanmer: To interpret alpha and beta diversity in DNA sequence analyses, 
# you should consider the specific indices used to measure them, as well as 
# the context in which the analysis was performed. High alpha diversity may indicate 
# a rich and diverse microbial community within a single location or sample, 
# while high beta diversity may indicate differences in microbial community 
# composition between different locations or samples.

# It is important to note that alpha and beta diversity 
# are not always independent measures, and changes in one 
# may affect the other. For example, an increase in alpha diversity 
# within a location may result in a decrease in beta diversity 
# between different locations, if the new species added to the community 
# are already present in other locations. Therefore, it is important 
# to carefully consider the biological and ecological factors that may 
# be driving changes in alpha and beta diversity in DNA sequence analyses.

plot_heatmap(physeq = ps, method = "MDS", distance = "wUniFrac", 
             title = "weighted-UniFrac", taxa.label = FALSE)


plot_net(ps, "WUniFrac", color = "Sample")
plot_net(ps, "WUniFrac", color = "Sample", laymeth = "circle")


alpha_div <- estimate_richness(ps, measures=NULL)
summary(alpha_div)
plot_richness(ps, measures=NULL)

ps.ord <- ordinate(ps, "NMDS", "bray")
p2 = plot_ordination(ps, ps.ord, type="samples", color="Sampleі") 
p2 + geom_polygon(aes(fill=Sample)) + geom_point(size=5) + ggtitle("samples")

p2
# PCoA plot using the unweighted UniFrac as distance
# wunifrac_dist = distance(ps.rarefied, method="unifrac", weighted=F)
ordination = ordinate(ps, "PCoA", "bray", weighted=TRUE)#, distance= wunifrac_dist)
plot_ordination(ps, ordination, type="taxa", color="Phylum", title="taxa") + theme(aspect.ratio=1)

test_ps <- make_network(ps, max.dist= 1000,distance='bray')
plot_network(test_ps,ps,label= 'Sample')

library(metacoder)
install.packages("microbiomeSeq")
ord.res <- ordinate(ps,distance="bray",method="NMDS",grouping_column="Sample")
p <- plot_ordination(ps, ord.res)
print(p)


env.taxa.cor <- env_taxa_correlation(ps, grouping_column="Sample", method="pearson", pvalue.threshold=0.05,
                                     padjust.method="BH", adjustment=5, num.taxa=50, select.variables=NULL)
heat_tree(ps, node_size = n_obs, node_label = name, node_color = n_obs)
