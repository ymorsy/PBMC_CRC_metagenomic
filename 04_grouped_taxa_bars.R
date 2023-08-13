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

source("./workflow/stacked_bar_prep.R")

pdf(snakemake@output$grouped_taxa_bars_plot,
    width = 15, height = 9
)

for (l in snakemake@params$taxrank) {
    rel_ab_df <- stacked_bar_prep(phyloseq_obj,
        tax_level = l,
        cutoff = snakemake@params$cutoff,
        meta_fct = snakemake@params$meta_fct,
    )

    no_l <- length(levels(rel_ab_df$ASV))
    ifelse(no_l <= 10,
        palette <- c(pal_aaas()(no_l - 1), "gray20"),
        palette <- c(pal_d3("category20")(no_l - 1), "gray20")
    )

    print(
        rel_ab_df %>%
            group_by(!!as.symbol(snakemake@params$meta_fct), ASV) %>%
            summarise(mean_rel_abund = sum(mean_rel_abund)) %T>%
            export(str_glue("{str_remove(snakemake@output$grouped_taxa_bars_plot, '.pdf')}_{l}.xlsx")) %>%
            ggplot(aes_string(y = "mean_rel_abund", fill = "ASV", x = "Status")) +
            # facet_wrap(~ Status, ncol = 3, scales = "free") +
            geom_col(width = 0.6, color = "white") +
            scale_fill_manual(values = palette) +
            labs(
                title = l,
                y = "Mean Relative Abundance (%)"
            ) +
            theme_classic() +
            geom_text(aes(label = paste(round(mean_rel_abund, 1), "%")),
                position = position_stack(vjust = 0.5), size = 2, color = "white"
            ) +
            theme(
                legend.text = element_text(face = "italic"),
                legend.key.size = unit(10, "pt"),
                axis.text.y = element_text(size = rel(2)),
                axis.title.y = element_text(size = rel(1.5)),
                axis.text.x = element_text(size = rel(1.5), angle = 45, hjust = 1),
            )
    )
}
dev.off()
