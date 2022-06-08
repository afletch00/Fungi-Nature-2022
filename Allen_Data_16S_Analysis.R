# read in files to make a phyloseq object

sample_info<- read.table("pilot_meta.tsv", header=T, row.names=1,
                         check.names=F, sep="\t")
sample_info$Kit[sample_info$Kit == "Power Soil"] <- "PowerSoil"

tax_tab <- as.matrix(read.table("coreV4V5_tax.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

count_tab <- read.table("coreV4V5_counts.tsv", header=T, 
                        row.names=1, check.names=F, sep="\t")

ps_all <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), sample_data(sample_info), 
                         tax_table(tax_tab))

sample_sums(ps_all)
tax_table(ps_all)

# Remove taxa not assigned at Kingdom level

rank_names(ps_all)
ps_taxa_filt <- subset_taxa(ps_all, !is.na(Kingdom) & !Kingdom %in% c("", "uncharacterized"))
table(tax_table(ps_all)[, "Kingdom"], exclude = NULL)

prevdf = apply(X = otu_table(ps_taxa_filt),
               MARGIN = ifelse(taxa_are_rows(ps_taxa_filt), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps_taxa_filt),
                    tax_table(ps_taxa_filt))

plyr::ddply(prevdf, "Kingdom", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
filterDomain = c("Eukaryota")
ps_filtDom = subset_taxa(ps_taxa_filt, !Kingdom %in% filterDomain)
prevdf1 = subset(prevdf, Kingdom %in% get_taxa_unique(ps_filtDom, "Kingdom"))

keepTaxa = rownames(prevdf)
ps_filt_2 = prune_taxa(keepTaxa, ps_taxa_filt)

# double check

sample_data(ps_filt_2)
tax_table(ps_filt_2)

ps_16S = ps_filt_2
sample_data(ps_16S)$Group <- as.factor(sample_data(ps_sub)$Group)

# Aggregate taxa to kingdom level and top Genera

ps.dom <- aggregate_taxa(ps_16S, "Kingdom")
ps.gen10 <- aggregate_top_taxa2(ps_16S, "Genus", top = 10) 
ps.gen10 <- subset_taxa(ps.gen10, !(Genus == "hgcI clade"))

# Transform to relative abundance and make composition plots.

ps.gen.rel <- microbiome::transform(ps.gen10, "compositional")

plot.composition.COuntAbun <- plot_composition(ps.gen.rel, group_by = "Group1", 
                                               otu.sort = c("Sediminibacterium", "Paraburkholderia","Vibrio","Flavobacterium", 
                                                           "Erwinia", "Uruburuella", "Paracoccus","Halomonas","Streptococcus",  
                                                            "Other"), 
                                               sample.sort = "Kit", 
                                             x.label = "Kit") + theme(legend.position = "bottom") +
  scale_fill_manual(values = c("lightskyblue3", "green4","tomato2", "#885682FF", "#6492AEFF", "yellowgreen",
                              "#FFCD00FF","#B0A255FF", "slateblue","#00B2B0FF", "grey79")) + 
  theme_bw(base_size = 17) + 
  guides(fill=guide_legend(ncol=1)) +
  geom_bar(colour = "black", stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x = "", y = "Relative Abundance\n") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + #ylim(0,40) +
  theme(legend.title = element_text(size = 15)) +
  facet_grid(~Group, scales = "free", space = "free_x") +
  theme(legend.text = element_text(colour="black", size = 15, face = "italic")) +
  guides(fill=guide_legend(title="Genera"))
plot.composition.COuntAbun

