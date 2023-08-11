# load the packges ----
pacman::p_load(
    rio, tidyverse, magrittr,
    doParallel, foreach,
    phyloseq, microbiome, vegan, pairwiseAdonis,
    RColorBrewer, pheatmap,
    ggsci, ggalluvial
)

# load the Phyloseq object and define the params -----
phyloseq_obj <- rio::import(snakemake@input$phyloseq_obj)
phyloseq_obj

ab <- tax_glom(phyloseq_obj, taxrank = l)
rel_ab <- transform_sample_counts(ab, function(x) x / sum(x))
rel_ab_df <- psmelt(rel_ab) %>%
    as_tibble() %>%
    mutate(OTU = as_factor(OTU)) %>%
    arrange(OTU)

# l_pos <- match(l, lev_tax)

rio::export(rel_ab_df, str_glue("01_Taxa/CRC_rel_abundance_level_{l}.xlsx"))
abs_ab_df <- psmelt(ab)
rio::export(abs_ab_df, str_glue("01_Taxa/CRC_abs_abundance_level_{l}.xlsx"))
