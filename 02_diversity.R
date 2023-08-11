# load the packges ----
pacman::p_load(
    rio, tidyverse, magrittr,
    doParallel, foreach,
    phyloseq, microbiome, vegan, pairwiseAdonis,
    RColorBrewer, pheatmap,
    ggsci, ggalluvial,
    plotly, htmlwidgets
)

# load the Phyloseq object and define the params -----
alpha_dis <- snakemake@params$alphas

alpha_dis

phyloseq_obj <- rio::import(snakemake@input$phyloseq_obj)
phyloseq_obj
sort(sample_sums(phyloseq_obj))

set.seed(5234)
rare_depth <- sort(sample_sums(phyloseq_obj))[Position(function(x) x > snakemake@params$depth, sort(sample_sums(phyloseq_obj)))]
phyloseq_obj_rar <- rarefy_even_depth(phyloseq_obj, sample.size = rare_depth)

# alpha
phyloseq_obj_rar_a_div <- alpha(phyloseq_obj_rar, index = "all") %>%
    rownames_to_column("Sample_name")
# get the metadata out as seprate object
phyloseq_obj_rar_meta <- meta(phyloseq_obj_rar) %>%
    rownames_to_column("Sample_name")
# merge these two data frames into one

phyloseq_obj_rar_meta

phyloseq_obj_rar_a_div_meta <- left_join(phyloseq_obj_rar_a_div, phyloseq_obj_rar_meta, by = "Sample_name") %>%
    as_tibble()
# x <- "response"
# facet <- "visit"
# alpha <- "observed"
str(phyloseq_obj_rar_a_div_meta)

# no_l <- length(levels(phyloseq_obj_rar_a_div_meta[[snakemake@params$x]]))

# ifelse(no_l <= 10,
#     palette <- c(pal_aaas()(no_l)),
#     palette <- c(pal_d3("category20")(no_l))
# )
palette <- rio::import(snakemake@input$clrs)

palette

# alpha diversity plotting ----
pdf(snakemake@output$alpha_diversity, height = 7, width = 10)
if (is.null(snakemake@params$facet)) {
    for (alpha in alpha_dis) {
        p <- ggplot(phyloseq_obj_rar_a_div_meta, aes_string(
            x = snakemake@params$x,
            y = alpha,
            color = snakemake@params$x
        )) +
            stat_boxplot(geom = "errorbar", width = 0.15) +
            geom_boxplot(aes_string(fill = snakemake@params$x), width = 0.3, alpha = 0.5) +
            geom_jitter(width = 0.1, size = 2.5) +
            scale_color_manual(values = palette) +
            scale_fill_manual(values = palette) +
            theme_bw() +
            theme(
                legend.key.size = unit(9, "pt"),
                axis.text.x = element_text(size = rel(0.8)),
                axis.text.y = element_text(size = rel(2)),
                axis.title.y = element_text(size = rel(1.5)),
            )

        print(p)
        saveWidget(
            ggplotly(p),
            file = str_glue("{snakemake@params$alpha_diversity_html}_{alpha}.html")
        )
    }
} else {
    for (alpha in alpha_dis) {
        p <- ggplot(phyloseq_obj_rar_a_div_meta, aes_string(
            x = snakemake@params$x,
            y = alpha,
            color = snakemake@params$x
        )) +
            facet_wrap(as.formula(paste("~", snakemake@params$facet))) +
            stat_boxplot(geom = "errorbar", width = 0.15) +
            geom_boxplot(aes_string(fill = snakemake@params$x), width = 0.3, alpha = 0.5) +
            geom_jitter(width = 0.1, size = 2.5) +
            scale_color_manual(values = palette) +
            scale_fill_manual(values = palette) +
            theme_bw() +
            theme(
                legend.key.size = unit(9, "pt"),
                axis.text.x = element_text(size = rel(0.8)),
                axis.text.y = element_text(size = rel(2)),
                axis.title.y = element_text(size = rel(1.5)),
            )

        print(p)
        saveWidget(ggplotly(p),
            file = str_glue("{snakemake@params$alpha_diversity_html}_{alpha}.html")
        )
    }
}
dev.off()
# beta diversity plotting ----

betas_distances <- c("bray", "jaccard")

betas <- snakemake@params$betas %>%
    purrr::set_names() %>%
    purrr::map(~ ordinate(phyloseq_obj_rar, "PCoA", .))


pdf(snakemake@output$beta_diversity, height = 7, width = 10)
if (is.null(snakemake@params$facet)) {
    for (i in seq_along(betas)) {
        p2 <- plot_ordination(phyloseq_obj_rar, betas[[i]], color = snakemake@params$x)
        p2 <- p2 + ggtitle(names(betas[i])) + geom_point(size = 2)
        p2 <- p2 + theme_classic() +
            scale_color_manual(values = palette) +
            theme(
                legend.key.size = unit(9, "pt"),
                axis.text.x = element_text(size = rel(0.8)),
                axis.text.y = element_text(size = rel(2)),
                axis.title.y = element_text(size = rel(1.5)),
            )
        print(p2)
        beta_name <- names(betas[i])
        saveWidget(
            ggplotly(p2),
            file = str_glue("{snakemake@params$beta_diversity_html}_{beta_name}.html")
        )
    }
    for (i in seq_along(betas)) {
        p2 <- plot_ordination(phyloseq_obj_rar, betas[[i]], color = snakemake@params$x)
        p2 <- p2 + ggtitle(names(betas[i])) + geom_point(size = 2)
        p2 <- p2 + theme_classic() + stat_ellipse() +
            scale_color_manual(values = palette) +
            theme(
                legend.key.size = unit(9, "pt"),
                axis.text.x = element_text(size = rel(0.8)),
                axis.text.y = element_text(size = rel(2)),
                axis.title.y = element_text(size = rel(1.5)),
            )
        print(p2)
        beta_name <- names(betas[i])
        saveWidget(
            ggplotly(p2),
            file = str_glue("{snakemake@params$beta_diversity_html}_{beta_name}.html")
        )
    }
} else {
    for (i in seq_along(betas)) {
        p2 <- plot_ordination(phyloseq_obj_rar, betas[[i]], color = snakemake@params$x, shape = snakemake@params$facet)
        p2 <- p2 + ggtitle(names(betas[i])) + geom_point(size = 2)
        p2 <- p2 + theme_classic() +
            scale_color_manual(values = palette) +
            theme(
                legend.key.size = unit(9, "pt"),
                axis.text.x = element_text(size = rel(0.8)),
                axis.text.y = element_text(size = rel(2)),
                axis.title.y = element_text(size = rel(1.5)),
            )
        print(p2)
        beta_name <- names(betas[i])
        saveWidget(
            ggplotly(p2),
            file = str_glue("{snakemake@params$beta_diversity_html}_{beta_name}.html")
        )
    }
    for (i in seq_along(betas)) {
        p2 <- plot_ordination(phyloseq_obj_rar, betas[[i]], color = snakemake@params$x, shape = snakemake@params$x) +
            facet_wrap(as.formula(paste("~", snakemake@params$facet)))
        p2 <- p2 + ggtitle(names(betas[i])) + geom_point(size = 2)
        p2 <- p2 + theme_classic() + stat_ellipse() +
            scale_color_manual(values = palette) +
            theme(
                legend.key.size = unit(9, "pt"),
                axis.text.x = element_text(size = rel(0.8)),
                axis.text.y = element_text(size = rel(2)),
                axis.title.y = element_text(size = rel(1.5)),
            )
        print(p2)
        beta_name <- names(betas[i])
        saveWidget(
            ggplotly(p2),
            file = str_glue("{snakemake@params$beta_diversity_html}_{beta_name}.html")
        )
    }
}
dev.off()
