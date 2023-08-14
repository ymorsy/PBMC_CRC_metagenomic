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

palette <- rio::import(snakemake@input$clrs)

if (!is.null(snakemake@input$df_rel)) {
    l <- str_replace(snakemake@input$df_rel, "(.+_level_)(.+)(\\.xlsx)", "\\2")
} else {
    l <- str_replace(snakemake@input$df_abs, "(.+_level_)(.+)(\\.xlsx)", "\\2")
}

if (length(palette) == 2) {
    df_stat_t_test <- rio::import(snakemake@input$df_rel) %>%
        as_tibble() %>%
        group_by(!!as.symbol(l)) %>%
        nest() %>%
        mutate(
            Exp_test = map(
                .x = data,
                ~ wilcox.test(x = .x$Abundance, g = .x$Status) %>%
                    broom::tidy()
            )
        ) %>%
        unnest(Exp_test) %>%
        mutate(
            adj.p.value = p.adjust(p.value, method = "fdr", n = length(.))
        ) %>%
        select(-data)

    rio::export(df_stat_t_test, snakemake@output$stat_rel)
    rio::export(df_stat_t_test %>% filter(adj.p.value < 0.05), snakemake@output$stat_rel_sig)
} else {
    df_stat_anova <- rio::import(snakemake@input$df_rel) %>%
        as_tibble() %>%
        group_by(!!as.symbol(l)) %>%
        nest() %>%
        mutate(Exp_test = map(
            .x = data,
            ~ aov(Abundance ~ !!as.symbol(snakemake@params$meta_fct), data = .x) %>%
                broom::tidy()
        )) %>%
        unnest(Exp_test) %>%
        select(-data)
    rio::export(df_stat_anova, snakemake@output$stat_rel)
    rio::export(df_stat_anova %>% filter(adj.p.value < 0.05), snakemake@output$stat_rel_sig)
}
