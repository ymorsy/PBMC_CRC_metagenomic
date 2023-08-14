# load the packges ----
pacman::p_load(
    rio, tidyverse, magrittr,
    doParallel, foreach,
    phyloseq, microbiome, vegan, pairwiseAdonis,
    RColorBrewer, pheatmap,
    ggsci, ggalluvial
)

phyloseq_obj_rar <- rio::import(snakemake@input$phyloseq_obj_rar)


# beta diversity statistics ----

snakemake@params$betas_dis
as.symbol(snakemake@params$meta_fct)
# for (dist in snakemake@params$beta_dis) {
dist_matrix <- distance(phyloseq_obj_rar, method = snakemake@params$betas_dis)
result <- adonis2(
    dist_matrix ~ Status,
    as(sample_data(phyloseq_obj_rar), "data.frame")
) %>%
    broom::tidy()

# export(result, str_glue("./03_Diversity/{dist}_permanova.xlsx"))

pairwise_results <- pairwise.adonis2(
    dist_matrix ~ Status,
    as(sample_data(phyloseq_obj_rar), "data.frame")
)

zl <- list()
for (i in seq_along(pairwise_results)[-1]) {
    tbl <- tibble(contrast = names(pairwise_results[i]), broom::tidy(pairwise_results[[i]]))
    zl[[i]] <- tbl
}

pairwise_results_df <- reduce(zl, bind_rows)

rio::export(pairwise_results_df, snakemake@output$betas_dis_stat)
# }
