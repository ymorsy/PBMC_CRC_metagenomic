# load the packges ----
pacman::p_load(rio, tidyverse, magrittr, phyloseq)

# Phyloseq object -----
# function to convert the metaphylan to phyloseq object
source("./workflow/metaphlanToPhyloseq.R")

# snakemake@input$df
# snakemake@input$metadata
# snakemake@params$fct
# snakemake@params$exp
# snakemake@output$phyloseq_obj
# snakemake@params$col_selection



select_vector <- c(
    snakemake@params$col_selection[5:length(snakemake@params$col_selection)],
    str_glue("{snakemake@params$exp}_{snakemake@params$col_selection[1]}"),
    str_glue("{snakemake@params$exp}_{snakemake@params$col_selection[2]}"),
    str_glue("{snakemake@params$exp}_{snakemake@params$col_selection[3]}"),
    str_glue("{snakemake@params$exp}_{snakemake@params$col_selection[4]}")
)

exp_name <- str_glue("{snakemake@params$exp}_exp")
exp_name
# fct <- str_glue("{snakemake@params$exp}_{snakemake@params$col_selection[1]}")
# fct



# Build the phyloseq object
df_metadata <- rio::import(snakemake@input$metadata) %>%
    rename(Sample_name = 1) %>%
    select(all_of(select_vector)) %>%
    filter(!!as.symbol(exp_name) == "Yes")

sample_names <- df_metadata %>%
    pull(Sample_name)

fct_level <- df_metadata %>%
    arrange(!!as.symbol(str_glue("{snakemake@params$exp}_{snakemake@params$col_selection[4]}"))) %>%
    distinct(!!as.symbol(str_glue("{snakemake@params$exp}_{snakemake@params$col_selection[1]}"))) %>%
    pull()

fct_level
clrs_level <- df_metadata %>%
    arrange(!!as.symbol(str_glue("{snakemake@params$exp}_{snakemake@params$col_selection[4]}"))) %>%
    distinct(!!as.symbol(str_glue("{snakemake@params$exp}_{snakemake@params$col_selection[3]}"))) %>%
    pull()

clrs_level

export(clrs_level, file = snakemake@output$clrs)

df_meta <- df_metadata %>%
    mutate(Status = factor(!!as.symbol(str_glue("{snakemake@params$exp}_{snakemake@params$col_selection[1]}")),
        levels = fct_level
    )) %>%
    column_to_rownames("Sample_name")

df_meta$Status
df_meta
full_df <- rio::import(snakemake@input$df) %>%
    rename(Taxa = 1) %>%
    select("Taxa", all_of(sample_names)) %>%
    column_to_rownames("Taxa")

# check the data match
rownames(full_df)[duplicated(rownames(full_df))]
table(duplicated(rownames(full_df)))
all(colnames(full_df) == rownames(df_meta))

# build and export the phyl_obj object
phyl_obj <- metaphlanToPhyloseq(full_df, df_meta, simplenames = F)
phyl_obj

export(phyl_obj, file = snakemake@output$phyloseq_obj)
export(df_meta %>% rownames_to_column("Sample_name"), file = snakemake@output$metadata)
export(df_meta %>% rownames_to_column("Sample_name"), file = snakemake@output$metadata_rds)
