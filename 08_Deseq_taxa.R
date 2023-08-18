# Purpose: Differential abundance analysis of taxa
pacman::p_load(
    rio, tidyverse, magrittr, ggsci, optparse, GenomicFeatures,
    Rsamtools, DESeq2, EnhancedVolcano, magrittr, PCAtools, pheatmap, RColorBrewer
)

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



for (i in 1:number_of_comp) {
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

    export(meta_f, str_glue("{out_path}/{comp2}_metadata.xlsx"))

    counts_f <- counts %>%
        dplyr::select(meta_f$Sample_name)


    # dds <- DESeqDataSetFromMatrix(
    #     countData = counts_f,
    #     colData = meta_f,
    #     design = as.formula(str_glue("~ {snakemake@params$fct}"))
    # )
    # DAT_deseq <- DESeq(dds)

    # nor_count <- counts(DAT_deseq, normalized = T) %>%
    #     as_tibble(rownames = "Taxa_name")
}
