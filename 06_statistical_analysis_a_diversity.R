# load the packges ----
pacman::p_load(
    rio, tidyverse, magrittr,
    doParallel, foreach,
    phyloseq, microbiome, vegan, pairwiseAdonis,
    RColorBrewer, pheatmap,
    ggsci, ggalluvial
)

# alpha statistics ----

phyloseq_obj_rar <- rio::import(snakemake@input$phyloseq_obj_rar)
phyloseq_obj_rar_a_div <- alpha(phyloseq_obj_rar, index = "all") %>%
    rownames_to_column("Sample_name")
# get the metadata out as seprate object
phyloseq_obj_rar_meta <- meta(phyloseq_obj_rar) %>%
    rownames_to_column("Sample_name")
# merge these two data frames into one

# phyloseq_obj_rar_meta

phyloseq_obj_rar_a_div_meta <- left_join(phyloseq_obj_rar_a_div, phyloseq_obj_rar_meta, by = "Sample_name") %>%
    as_tibble()


# for (alpha in snakemake@params$alpha_dis) {
anova_result <- aov(as.formula(paste0(snakemake@params$alpha_dis, " ~ ", snakemake@params$meta_fct)),
    data = phyloseq_obj_rar_a_div_meta
) %T>%
    broom::tidy() %T>%
    rio::export(snakemake@output$alpha_dis_stat)
TukeyHSD(anova_result) %>%
    broom::tidy() %T>%
    rio::export(snakemake@output$alpha_dis_stat_tukey)
# }
