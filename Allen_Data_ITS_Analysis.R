
sample_info<- read.table("pilot_meta.tsv", header=T, row.names=1,
                         check.names=F, sep="\t")
sample_info$Kit[sample_info$Kit == "Power Soil"] <- "PowerSoil"

tax_tab <- as.matrix(read.table("coreITS_tax.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

count_tab <- read.table("coreITS_counts.tsv", header=T, 
                        row.names=1, check.names=F, sep="\t")

ps_all <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE), sample_data(sample_info), 
                         tax_table(tax_tab))

sample_sums(ps_all)
tax_table(ps_all)
sample_data(ps_all)

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
filterDomain = c("Metazoa", "Viridiplantae")
ps_filtDom = subset_taxa(ps_taxa_filt, !Kingdom %in% filterDomain)
prevdf1 = subset(prevdf, Kingdom %in% get_taxa_unique(ps_filtDom, "Kingdom"))

plyr::ddply(prevdf1, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
filterPhyla = c("unidentified")
ps_filtPhy = subset_taxa(ps_filtDom, !Phylum %in% filterPhyla)
prevdf2 = subset(prevdf1, Phylum %in% get_taxa_unique(ps_filtPhy, "Phylum"))

keepTaxa = rownames(prevdf2)
ps_filt_2 = prune_taxa(keepTaxa, ps_taxa_filt)

sample_data(ps_filt_2)
tax_table(ps_filt_2)
sample_sums(ps_filt_2)
tax_table(ps_filt_2)[tax_table(ps_filt_2) == "unidentified"] <- "Unidentified"
##################### Subset samples ###########################

ps_sub = ps_filt_2
sample_data(ps_sub)
sample_sums(ps_sub)
sample_data(ps_sub)$Group <- as.factor(sample_data(ps_sub)$Group)

tax_table(ps_sub)
ps.dom <- aggregate_taxa(ps_sub, "Kingdom")
ps.gen10 <- aggregate_top_taxa2(ps_sub, "Genus", top = 10) 
sample_sums(ps.gen10)

mypal = pal_startrek("uniform", alpha = 1)(7)
mypal
#CC0C00FF" "#5C88DAFF" "#84BD00FF" "#FFCD00FF" "#7C878EFF" "#00B5E2FF" "#00AF66FF"

plot.composition.COuntAbun <- plot_composition(ps.gen10, group_by = "Kit", 
                                               otu.sort = c("Cladosporium", "Alternaria", "Naganishia", "Rhodosporidiobolus",
                                                           "Penicillium", "Aspergillus", "Unidentified"), 
                                               #sample.sort = "Kit", 
                                             x.label = "ID") + theme(legend.position = "bottom") +
  scale_fill_manual(values = c("#CC0C00FF", "#5C88DAFF","#00AF66FF", "#00B5E2FF",
    "#7C878EFF", "orange","grey89")) + 
  theme_bw(base_size = 17) +
  guides(fill=guide_legend(ncol=1)) +
  geom_bar(colour = "black", stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  #ggtitle("DNA Extraction Kit") + 
  labs(x = "", y = "Read Counts\n") +
  theme(legend.title = element_text(size = 15)) +
  facet_grid(~Group, scales = "free", space = "free_x") +
  guides(fill=guide_legend(title="Genera")) +
  theme(legend.text = element_text(colour="black", size = 15, face = "italic")) 
#plot.composition.COuntAbun$data$Group <- factor(plot.composition.COuntAbun$data$Group, levels = desired_order)
plot.composition.COuntAbun

