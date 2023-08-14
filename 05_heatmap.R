# load the packges ----
pacman::p_load(
    rio, tidyverse, magrittr,
    doParallel, foreach,
    phyloseq, microbiome, vegan, pairwiseAdonis,
    RColorBrewer, pheatmap,
    ggsci, ggalluvial
)

# load the Phyloseq object and define the params -----

df_meta <- rio::import(snakemake@input$metadata) %>%
    select(
        Sample_name = 1, snakemake@params$meta_fct
    ) %>%
    column_to_rownames("Sample_name")

head(df_meta)
df_meta$`snakemake@params$meta_fct` %>% levels()
palette <- rio::import(snakemake@input$clrs)
palette
names(palette) <- pull(df_meta, snakemake@params$meta_fct) %>% levels()
palette

colrs_list <- list(palette)
names(colrs_list)[1] <- snakemake@params$meta_fct
colrs_list

rio::import(snakemake@input$metadata) %>%
    as_tibble()
if (!is.null(snakemake@input$df_rel)) {
    l <- str_replace(snakemake@input$df_rel, "(.+_level_)(.+)(\\.xlsx)", "\\2")
} else {
    l <- str_replace(snakemake@input$df_abs, "(.+_level_)(.+)(\\.xlsx)", "\\2")
}

rel_ab_df <- rio::import(snakemake@input$df_rel) %>%
    as_tibble() %>%
    mutate(OTU = factor(!!as.symbol(l))) %>%
    select(c(Sample, OTU, Abundance)) %>%
    pivot_wider(names_from = Sample, values_from = Abundance) %>%
    column_to_rownames("OTU") %>%
    Filter(function(x) !all(is.na(x)), .) %>%
    dplyr::mutate(sum = rowSums(.)) %>%
    arrange(desc(sum))

newnames <- rownames(rel_ab_df %>% select(-sum) %>% slice(1:30)) %>%
    # str_replace(., "(.+__)(.+)", "\\2") %>%
    lapply(., function(x) bquote(italic(.(x))))

meta <- df_meta

meta$Status
pdf(snakemake@output$hm_rel, width = 12)
print(
    rel_ab_df %>% select(-sum) %>% slice(1:30) %>%
        pheatmap(
            color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256), annotation_col = meta,
            show_colnames = F, cluster_rows = F, cluster_cols = F, scale = "row", annotation_colors = colrs_list,
            labels_row = as.expression(newnames)
        )
)

print(
    rel_ab_df %>% select(-sum) %>% slice(1:30) %>%
        pheatmap(
            color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256), annotation_col = meta,
            show_colnames = F, cluster_rows = T, cluster_cols = T, scale = "row", annotation_colors = colrs_list,
            labels_row = as.expression(newnames)
        )
)
dev.off()

abs_ab_df <- rio::import(snakemake@input$df_abs) %>%
    as_tibble() %>%
    mutate(OTU = factor(!!as.symbol(l))) %>%
    select(c(Sample, OTU, Abundance)) %>%
    pivot_wider(names_from = Sample, values_from = Abundance) %>%
    column_to_rownames("OTU") %>%
    Filter(function(x) !all(is.na(x)), .) %>%
    dplyr::mutate(sum = rowSums(.)) %>%
    arrange(desc(sum))

newnames <- rownames(abs_ab_df %>% select(-sum) %>% slice(1:30)) %>%
    # str_replace(., "(.+__)(.+)", "\\2") %>%
    lapply(., function(x) bquote(italic(.(x))))

meta <- df_meta


pdf(snakemake@output$hm_abs, width = 12)
print(
    abs_ab_df %>% select(-sum) %>% slice(1:30) %>%
        pheatmap(
            color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256), annotation_col = meta,
            show_colnames = F, cluster_rows = F, cluster_cols = F, annotation_colors = colrs_list,
            labels_row = as.expression(newnames)
        )
)

print(
    abs_ab_df %>% select(-sum) %>% slice(1:30) %>%
        pheatmap(
            color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256), annotation_col = meta,
            show_colnames = F, cluster_rows = T, cluster_cols = T, annotation_colors = colrs_list,
            labels_row = as.expression(newnames)
        )
)
dev.off()
