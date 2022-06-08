# Libraries used for microbiome analysis- vegan, phyloseq, microbiome, microbiomeSeq, microbiomeutilities

# Read in sample files made with QIIME

sample_info<- read.table("meta.all.tsv", header=T, row.names=1,
                         check.names=F, sep="\t")

tax_tab <- as.matrix(read.table("QIIME_all_tax.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

count_tab <- read.table("QIIME_all_OTU_tab.tsv", header=T, 
                        row.names=1, check.names=F, sep="\t")

ps_allQ <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), sample_data(sample_info), 
                         tax_table(tax_tab))

# Check
sample_sums(ps_allQ)
tax_table(ps_allQ)

# Remove taxa not assigned at Domain level (Kingdom)

rank_names(ps_allQ)
ps_taxa_filt <- subset_taxa(ps_allQ, !is.na(domain) & !domain %in% c("", "uncharacterized"))
table(tax_table(ps_allQ)[, "domain"], exclude = NULL)

prevdf = apply(X = otu_table(ps_taxa_filt),
               MARGIN = ifelse(taxa_are_rows(ps_taxa_filt), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps_taxa_filt),
                    tax_table(ps_taxa_filt))

plyr::ddply(prevdf, "domain", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
filterDomain = c("unidentified")
ps_filtDom = subset_taxa(ps_taxa_filt, !domain %in% filterDomain)
prevdf1 = subset(prevdf, domain %in% get_taxa_unique(ps_filtDom, "domain"))

keepTaxa = rownames(prevdf1)
ps_filt_2 = prune_taxa(keepTaxa, ps_taxa_filt)

# filter taxa (if any) with 0 counts
ps_filt_2 = prune_taxa(taxa_sums(ps_filt_2) > 0, ps_filt_2)
## filter out single OTUs
ps_filt_2 = filter_taxa(ps_filt_2, function (x) {sum(x > 0) > 1}, prune=TRUE)

# Check
sample_data(ps_filt_2)
tax_table(ps_filt_2)

sample_data(ps_filt_2)$sample_type <- as.factor(sample_data(ps_filt_2)$sample_type)

# Make seperate phyloseq objects for gut and human samples

ps_subG = subset_samples(ps_filt_2, sample_type %in% c("Norm. Gut", "PDAC Gut"))
ps_subT = subset_samples(ps_filt_2, sample_type %in% c("Norm. Panc.", "PDAC Panc."))

# Aggregate to top 10 Genera

ps.gen10T <- aggregate_top_taxa2(ps_subT, "genus", top = 10) 
ps.gen10G <- aggregate_top_taxa2(ps_subG, "genus", top = 10)

# Plot read counts

plot.composition.COuntAbun <- plot_composition(ps.gen10T, group_by = "sample_type", otu.sort = "abundance",
                                               sample.sort = "Malassezia", 
                                             x.label = "") + theme(legend.position = "bottom") +
  scale_fill_lancet_adaptive() + theme_bw(base_size = 17) +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.text.x = element_text(angle = 90)) + #ylim(0,40) +
  #ggtitle("Total abundance; Malassezia OTU's in Pancreas Tissue") +  
  theme(legend.title = element_text(size = 15)) +
  facet_grid(~Group, scales = "free", space = "free_x") +
  guides(fill=guide_legend(title="Malassezia OTUs"))
#plot.composition.COuntAbun$data$Group <- factor(plot.composition.COuntAbun$data$Group, levels = desired_order)
plot.composition.COuntAbun

# Subset gut and pancreas samples to Domain and Malassezia only and plot

ps.dom <- aggregate_taxa(ps_subG, "domain", verbose = FALSE) # ps_subT
ps.mal <- subset_taxa(ps.gen10G, genus == "Malassezia") # ps_sub10T
sample_sums(ps.mal)
meltsub <- psmelt(ps.dom) # psmelt(ps.mal)
meltsub

psmelt = ggplot(meltsub, aes(x = sample_type, y = Abundance)) +  ylim(0,2000) +
  geom_boxplot(outlier.shape  = NA, size = 1.5) +
  stat_compare_means(label.y = Inf, vjust = 1.5, size = 7, label.x.npc = 0.135) +
  geom_jitter(aes(color=sample_type, fill=sample_type), shape=21, size = 7,
              color = "black", position = position_jitterdodge(1))+
  #scale_fill_manual(values=c( "#FFCD00FF", "#FFCD00FF")) + # For Malassezia
  scale_fill_manual(values=c("green3", "red2")) + # blue and orange for pancreas
  theme_bw(base_size = 27) +
  labs(x = "", y = "Read Counts\n") +
  theme(legend.position = "none") +
  facet_wrap(~ OTU, scales = "free") 
psmelt

# Transform to relative abundance and plot

ps.gen.relT <- microbiome::transform(ps.gen10T, "compositional")
ps.gen.relG <- microbiome::transform(ps.gen10G, "compositional")

plot.composition.relAbunT <- plot_composition(ps.gen.relT, group_by = "sample_type", #otu.sort = "abundance", 
                                             otu.sort = c("Candida", "Malassezia", "Cystobasidium", "Aureobasidium",
                                                          "Aspergillus", "Tilletia", "Ustilago", "Exophiala", "Sclerotinia","Unknown","Other"),
                                             sample.sort = "Malassezia", 
                                             x.label = "") + theme(legend.position = "right") +
  scale_fill_manual(values = c( "skyblue3", "#FFCD00FF","#7CB22BFF", "#885682FF" ,"#CC0C00FF",
                                "plum3" ,"#00AF66FF" ,"steelblue1","tomato","grey37", "grey87")) +
  theme_bw(base_size = 17) + 
  geom_bar(colour = "black", stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "", y = "Read Counts\n") +
  theme(legend.text = element_text(colour="black", size = 15, face = "italic")) +
  guides(fill=guide_legend(title="Genera"))
plot.composition.relAbunT


plot.composition.relAbunG <- plot_composition(ps.gen.relG, group_by = "sample_type", #otu.sort = "abundance",
                                              otu.sort = c("Candida", "Saccharomyces", "Malassezia", "Debaryomyces", "Piskurozyma", "Ramularia",
                                                           "Filobasidium","Vishniacozyma","Hanseniaspora", "Unknown","Other"),
                                             sample.sort = "Malassezia", 
                                             x.label = "") + theme(legend.position = "right") +
  scale_fill_manual(values = c( "skyblue3", "green3","#FFCD00FF", "#885682FF" ,"red2",
                                "cyan3" ,"plum" ,"royalblue3","tomato","grey37", "grey87")) +
  theme_bw(base_size = 17) + 
  geom_bar(colour = "black", stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "", y = "Read Counts\n") +
  theme(legend.text = element_text(colour="black", size = 15, face = "italic")) +
  guides(fill=guide_legend(title="Genera"))
plot.composition.relAbunG

# Make Alpha diveristy plots with wilcoxon derived p-values
# Remove samples with 0 ASVs

ps_sub0T = prune_samples(sample_sums(ps_subT) > 0, ps_subT)
ps_sub0G = prune_samples(sample_sums(ps_subG) > 0, ps_subG)

p1 = plot_richness(ps_sub0G, x="sample_type", measures=c("Observed", "ACE", "Chao1", "Shannon", "Simpson" )) + #ps-sub0T
  geom_boxplot(size = 1, color="black") +
  theme(text = element_text(size = 21, color = "black")) + 
  theme_bw(base_size = 19) +
  theme(axis.text.x = element_text(angle =45, vjust = 1, hjust=1)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  stat_compare_means(method = "wilcox", label.y = Inf, vjust = 1.5, label = "p", size = 5) 
p1 + labs(x = "", y = "Alpha Diversity Measure") 

# Make ordination plot (PCoA) using Bray-Curtis distance matrix and PERMANOVA derived p-values

otu_r <- abundances(ps_sub0G) # ps-sub0T
meta_r <- meta(ps_sub0G) # ps-sub0T
permanova <- adonis(t(otu_r) ~ sample_type,
                    data = meta_r, permutations=999, method = "bray") 
permanova

ps_nmdsB <- ordinate(physeq = ps_sub0G, method = "PCoA", distance = "bray") # ps-sub0T

plot_ordination(physeq = ps_sub0G, ordination = ps_nmdsB, # ps-sub0T
                color = "sample_type",title = "p = 0.47") + theme_bw()+   # Use p-value derived from adonis for title
  scale_colour_manual(values=c("Norm. Gut"="green3", "PDAC Gut"="red2")) + # blue and orange for pancreas tissue
  geom_point(aes(color = sample_type), size = 5) + 
  theme_classic(base_size = 23) + 
  theme(legend.title = element_blank()) +
  theme(plot.title = element_text(vjust =  -7, hjust = 1.0, size = 23)) +
  stat_ellipse(aes(group = sample_type), size = 1.3, type="t")

