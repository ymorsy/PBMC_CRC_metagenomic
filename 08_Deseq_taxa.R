# Purpose: Differential abundance analysis of taxa
pacman::p_load(
    rio, tidyverse, magrittr, ggsci, optparse, GenomicFeatures,
    Rsamtools, DESeq2, EnhancedVolcano, magrittr, PCAtools, pheatmap, RColorBrewer,
    foreach, doParallel
)

n_cores <- 8


my_cluster <- parallel::makeCluster(
    n_cores,
    type = "FORK"
)

doParallel::registerDoParallel(cl = my_cluster)
# prepare the count file

counts <- rio::import(snakemake@input$taxa_count_abs_in) %>%
    as_tibble() %>%
    dplyr::select(!!as.symbol(snakemake@params$taxrank), Sample, Abundance) %>%
    pivot_wider(
        id_cols = !!as.symbol(snakemake@params$taxrank),
        names_from = Sample, values_from = Abundance,
        values_fn = {
            mean
        }
    ) %>%
    dplyr::filter((rowSums(across(where(is.numeric))) != 0)) %>%
    column_to_rownames(snakemake@params$taxrank)

# dplyr::rename(Taxon = !!as.symbol(snakemake@params$taxrank)) %>%


counts %>%
    dplyr::slice(5) %>%
    export(snakemake@output$taxa_count_abs_out)


metadata <- rio::import(snakemake@input$metadata) %>%
    arrange(!!as.symbol(snakemake@params$fct))



comps <- combn(levels(metadata[[snakemake@params$fct]]), 2)
comps
number_of_comp <- length(comps)
number_of_comp <- number_of_comp / 2
print(paste("The number of comparisons is :", number_of_comp))


palette <- rio::import(snakemake@input$clrs)
palette
names(palette) <- pull(metadata, snakemake@params$fct) %>% levels()
palette


foreach(i = 1:number_of_comp) %dopar% {
    comp <- comps[, i]
    # name the comparison
    comp2 <- str_glue("{comp[1]}_vs_{comp[2]}")
    str_glue("{snakemake@params$exp}/08_Read_norm_Deseq/{comp2}")

    out_path <- str_glue("{snakemake@params$exp}/08_Read_norm_Deseq/00{i}_{comp2}")
    dir.create(out_path)
    dir.create(str_glue("{out_path}/01_Counts"))
    dir.create(str_glue("{out_path}/02_PCA"))
    dir.create(str_glue("{out_path}/03_Differentially abundant taxa"))
    dir.create(str_glue("{out_path}/05_Heatmap"))
    dir.create(str_glue("{out_path}/04_Volcano plots"))

    meta_f <- metadata %>%
        dplyr::rename(Sample_name = 1) %>%
        dplyr::mutate(Sample_name = as.character(Sample_name)) %>%
        dplyr::filter(!!as.symbol(snakemake@params$fct) %in% comp)

    colrs <- palette[comp]
    print(colrs)
    export(meta_f, str_glue("{out_path}/{comp2}_metadata.xlsx"))

    counts_f <- counts %>%
        dplyr::select(meta_f$Sample_name)


    dds <- DESeqDataSetFromMatrix(
        countData = counts_f,
        colData = meta_f,
        design = as.formula(str_glue("~ {snakemake@params$fct}"))
    )
    DAT_deseq <- DESeq(dds)

    nor_count <- counts(DAT_deseq, normalized = T) %>%
        as_tibble(rownames = "Taxa_name")


    export(
        counts_f
        %<>% rownames_to_column("Taxa_name"),
        str_glue("{out_path}/01_Counts/Counts.xlsx")
    )
    export(
        nor_count,
        str_glue("{out_path}/01_Counts/Normalized_counts.xlsx")
    )
    # PCA

    p <- nor_count %>%
        as_tibble() %>%
        filter(rowSums(across(where(is.numeric))) != 0) %>%
        column_to_rownames("Taxa_name") %>%
        pca(scale = T, metadata = meta_f %>% column_to_rownames("Sample_name"))

    p$label <- p[["yvars"]]
    ### PCA ellipse sample names
    pca_scaled_ellipse <- PCAtools::biplot(p,
        colby = snakemake@params$fct, colkey = colrs,
        pointSize = 3,
        legendPosition = "right",
        boxedLabels = TRUE,
        ellipse = TRUE, ellipseType = "t",
        ellipseLevel = 0.85, ellipseFill = TRUE,
        ellipseAlpha = 0.25, ellipseLineSize = 1.0
    )

    ggsave(str_glue("{out_path}/02_PCA/PCA_ellipse.pdf", plot = pca_scaled_ellipse, height = 7, width = 9, units = "in"))

    # PCA ellipse without sample names
    pca_scaled_ellipse_noNames <- PCAtools::biplot(p,
        colby = snakemake@params$fct, colkey = colrs,
        pointSize = 3,
        legendPosition = "right",
        lab = NULL,
        ellipse = TRUE, ellipseType = "t",
        ellipseLevel = 0.85, ellipseFill = TRUE,
        ellipseAlpha = 0.25, ellipseLineSize = 1.0
    )

    ggsave(str_glue("{out_path}/02_PCA/PCA_ellipse_noNames.pdf", plot = pca_scaled_ellipse_noNames, height = 7, width = 9, units = "in"))



    # PCA without ellipse


    pca_scaled <- PCAtools::biplot(p,
        colby = snakemake@params$fct, colkey = colrs,
        pointSize = 3,
        legendPosition = "right",
        boxedLabels = TRUE
    )

    ggsave(str_glue("{out_path}/02_PCA/PCA.pdf", plot = pca_scaled, height = 7, width = 9, units = "in"))



    # PCA without ellipse witout names


    pca_scaled_noNames <- PCAtools::biplot(p,
        colby = snakemake@params$fct, colkey = colrs,
        pointSize = 3,
        legendPosition = "right",
        lab = NULL
    )


    ggsave(str_glue("{out_path}/02_PCA/PCA_noNames.pdf", plot = pca_scaled_noNames, height = 7, width = 9, units = "in"))


    # DAT

    res <- results(DAT_deseq, contrast = c(snakemake@params$fct, comp[2], comp[1]))

    # res <- results(DAT_deseq)
    # res <- results(DAT_deseq, contrast=c("Condition",comp,ref))


    DAT <- res %>%
        as_tibble(rownames = "Taxa_name") %>%
        filter(padj != is.na(padj)) %>%
        arrange(padj) %>%
        left_join(nor_count, by = "Taxa_name")

    colnames(DAT) <- c(
        colnames(DAT)[1:7],
        str_c(
            meta_f[[snakemake@params$fct]],
            "_",
            meta_f[["Sample_name"]]
        )
    )

    DAT1 <- DAT %>%
        filter(padj < 0.05)

    DAT2 <- DAT %>%
        filter(padj < 0.05) %>%
        filter(log2FoldChange > 2 | log2FoldChange < -2)

    DAT3 <- DAT %>%
        filter(padj < 0.001)

    DAT4 <- DAT %>%
        filter(padj < 0.001) %>%
        filter(log2FoldChange > 2 | log2FoldChange < -2)

    dir.create("03_Differentially abundant taxa")

    export(DAT, str_glue("{out_path}/03_Differentially abundant taxa/DAT.xlsx"))
    export(DAT1, str_glue("{out_path}/03_Differentially abundant taxa/DAT1.xlsx"))
    export(DAT2, str_glue("{out_path}/03_Differentially abundant taxa/DAT2.xlsx"))
    export(DAT3, str_glue("{out_path}/03_Differentially abundant taxa/DAT3.xlsx"))
    export(DAT4, str_glue("{out_path}/03_Differentially abundant taxa/DAT4.xlsx"))



    label <- c("DAT", "DAT1", "DAT2", "DAT3", "DAT4")

    description <- c(
        " all taxa ",
        " taxa with P value less than 0.05 ",
        " taxa with P value less than 0.05 and fold change less than  -2 or more than 2 ",
        " taxa with P value less than 0.001 ",
        " taxa with P value less than 0.001 and fold change less than  -2 or more than 2 "
    )

    taxaNo <- c(nrow(DAT), nrow(DAT1), nrow(DAT2), nrow(DAT3), nrow(DAT4))

    DATStat <- as.data.frame(cbind(label, description, taxaNo))

    colnames(DATStat) <- c("Label", "Description", "taxa numbers")


    export(DATStat, str_glue("{out_path}/03_Differentially abundant taxa/DAT_stat.xlsx"))

    tryCatch(
        {
            jpeg(str_glue("{out_path}/04_Volcano plots/volcano1.jpeg"), width = 15, height = 15, units = "in", res = 1200)
            print(
                EnhancedVolcano(DAT, lab = DAT$Taxa_name, x = "log2FoldChange", y = "padj", pCutoff = 0.05, FCcutoff = 0)
            )
            dev.off()

            jpeg(str_glue("{out_path}/04_Volcano plots/volcano2.jpeg"), width = 25, height = 25, units = "in", res = 1200)
            print(
                EnhancedVolcano(DAT, lab = DAT$Taxa_name, x = "log2FoldChange", y = "padj", pCutoff = 0.05, FCcutoff = 2)
            )
            dev.off()

            jpeg(str_glue("{out_path}/04_Volcano plots/volcano3.jpeg"), width = 25, height = 25, units = "in", res = 1200)
            print(
                EnhancedVolcano(DAT, lab = DAT$Taxa_name, x = "log2FoldChange", y = "padj", pCutoff = 0.001, FCcutoff = 0)
            )
            dev.off()

            jpeg(str_glue("{out_path}/04_Volcano plots/volcano4.jpeg"), width = 25, height = 25, units = "in", res = 1200)
            print(
                EnhancedVolcano(DAT, lab = DAT$Taxa_name, x = "log2FoldChange", y = "padj", pCutoff = 0.001, FCcutoff = 2)
            )
            dev.off()
        },
        error = function(e) {
            cat("ERROR :", conditionMessage(e), "\n")
        }
    )


    # Heatmap

    # meta_m <- as.data.frame(meta_f)
    # colrs_list <- list(colrs)
    # names(colrs_list)[1] <- snakemake@params$meta_fct
    # rownames(meta_m) <- str_c(meta_m$Condition, "_", meta_m$Sample_name)
    meta <- meta_f
    colrs_list <- list(colrs)
    names(colrs_list)[1] <- snakemake@params$fct

    rownames(meta) <- str_c(
        meta_f[[snakemake@params$fct]],
        "_",
        meta_f[["Sample_name"]]
    )
    colrs_list
    meta
    tryCatch(
        {
            pdf(str_glue("{out_path}/05_Heatmap/Hm1.pdf"), width = 20)
            print(
                DAT1 %>%
                    dplyr::select(c(1, 8:length(DAT1))) %>%
                    column_to_rownames("Taxa_name") %>%
                    pheatmap(
                        color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256),
                        scale = "row"
                    )
            )
            dev.off()
            pdf(str_glue("{out_path}/05_Heatmap/Hm1_top50.pdf"), height = 20, width = 20)
            print(
                DAT1 %>%
                    dplyr::select(c(1, 8:length(DAT1))) %>%
                    dplyr::slice(1:50) %>%
                    column_to_rownames("Taxa_name") %>%
                    pheatmap(
                        color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256),
                        scale = "row", cellheight = 9.5
                    )
            )
            dev.off()
            pdf(str_glue("{out_path}/05_Heatmap/Hm2.pdf"), width = 20)
            print(
                DAT2 %>%
                    dplyr::select(c(1, 8:length(DAT2))) %>%
                    column_to_rownames("Taxa_name") %>%
                    pheatmap(
                        color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256),
                        scale = "row"
                    )
            )
            dev.off()
            pdf(str_glue("{out_path}/05_Heatmap/Hm2_top50.pdf"), height = 20, width = 20)
            print(
                DAT2 %>%
                    dplyr::select(c(1, 8:length(DAT2))) %>%
                    dplyr::slice(1:50) %>%
                    column_to_rownames("Taxa_name") %>%
                    pheatmap(
                        color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256),
                        scale = "row", cellheight = 9.5
                    )
            )
            dev.off()
            pdf(str_glue("{out_path}/05_Heatmap/Hm3.pdf"), width = 20)
            print(
                DAT3 %>%
                    dplyr::select(c(1, 8:length(DAT3))) %>%
                    column_to_rownames("Taxa_name") %>%
                    pheatmap(
                        color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256),
                        scale = "row"
                    )
            )
            dev.off()
            pdf(str_glue("{out_path}/05_Heatmap/Hm3_top50.pdf"), height = 20, width = 20)
            print(
                DAT3 %>%
                    dplyr::select(c(1, 8:length(DAT3))) %>%
                    dplyr::slice(1:50) %>%
                    column_to_rownames("Taxa_name") %>%
                    pheatmap(
                        color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256),
                        scale = "row", cellheight = 9.5
                    )
            )
            dev.off()
            pdf(str_glue("{out_path}/05_Heatmap/Hm4.pdf"), width = 20)
            print(
                DAT4 %>%
                    dplyr::select(c(1, 8:length(DAT4))) %>%
                    column_to_rownames("Taxa_name") %>%
                    pheatmap(
                        color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256),
                        scale = "row"
                    )
            )
            dev.off()
            pdf(str_glue("{out_path}/05_Heatmap/Hm4_top50.pdf"), height = 20, width = 20)
            print(
                DAT4 %>%
                    dplyr::select(c(1, 8:length(DAT4))) %>%
                    dplyr::slice(1:50) %>%
                    column_to_rownames("Taxa_name") %>%
                    pheatmap(
                        color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(256),
                        scale = "row", cellheight = 9.5
                    )
            )
            dev.off()
        },
        error = function(e) {
            cat("ERROR :", conditionMessage(e), "\n")
        }
    )
}
# stop the cluster server
parallel::stopCluster(cl = my_cluster)
