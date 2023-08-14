pacman::p_load(tidyverse, magrittr, ggsci, rio)

#load the reference map krakeb database
mpa_map <- import(snakemake@input$k2_mpa_map) %>%
    as_tibble() %>%
    select(Taxa = 1) %>%
    mutate(
        Taxa = str_replace_all(Taxa, "(?<!_)_(?!_)", " "),
        match := str_extract(Taxa, "(?<=__)[^__]+$") %>%
            str_replace_all(" sp. ", "_") %>%
            str_replace_all(" ", "_")
    )


# prepare the data ----
import(snakemake@input$taxa_count_rel) %>%
    as_tibble() %>%
    # filter(!!as.symbol(filter_factor) %in% filter_criteria) %>%
    mutate(!!as.symbol(snakemake@params$taxrank) := str_replace_all(!!as.symbol(snakemake@params$taxrank), " sp. ", "_") %>%
        str_replace_all(" ", "_")) %>%
    left_join(., mpa_map, by = setNames("match", snakemake@params$taxrank), copy = TRUE, suffix = c(".x", ".y"), relationship = "many-to-many") %>%
    drop_na() %T>%
    { # notice the T pipe(%T>%) in the previous line
        # this block is used to define the names of ASV to select them later in the loop # nolint
        pull(., Taxa) ->> x
    } %>%
    pivot_wider(id_cols = c(Sample, Status), names_from = Taxa, values_from = Abundance, values_fn = {
        mean
    }) %>%
    select(Sample, Status, unique(x)) %>%
    group_by(Sample, Status) %>%
    summarise(across(where(is.numeric), na.omit), .groups = "drop") %>%
    export(snakemake@output$lefse_df)
