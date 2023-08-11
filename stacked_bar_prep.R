# function to summerize the level of GP for stacked bar plot ----
stacked_bar_prep <- function(phyloseq_obj,
                             tax_level = String,
                             cutoff = numeric,
                             meta_fct = String,
                             meta_fct2 = String) {
    if (missing(meta_fct2)) {
        df <- tax_glom(phyloseq_obj, taxrank = tax_level) %>%
            transform_sample_counts(function(x) x / sum(x)) %>%
            psmelt() %>%
            as_tibble() %>%
            mutate(
                OTU = as.factor(!!as.symbol(tax_level))
            ) %>%
            arrange(OTU) %>%
            group_by(!!as.symbol(meta_fct), OTU) %>%
            summarise(mean_rel_abund = mean(Abundance, na.rm = TRUE) * 100, .groups = "drop") %>%
            group_by(OTU) %>%
            mutate(pool = max(mean_rel_abund) < cutoff) %>%
            mutate(OTU = if_else(pool, "Other", as.character(OTU))) %>%
            ungroup(OTU) %>%
            arrange(desc(mean_rel_abund)) %>%
            mutate(OTU = factor(OTU) %>%
                fct_reorder(mean_rel_abund, max, .desc = TRUE) %>%
                fct_relevel("Other", after = Inf)) %>%
            select(ASV = OTU, !!as.symbol(meta_fct), mean_rel_abund)

        return(df)
    } else {
        df2 <- tax_glom(phyloseq_obj, taxrank = tax_level) %>%
            transform_sample_counts(function(x) x / sum(x)) %>%
            psmelt() %>%
            as_tibble() %>%
            mutate(
                OTU = as.factor(!!as.symbol(tax_level))
            ) %>%
            arrange(OTU) %>%
            group_by(!!as.symbol(meta_fct), !!as.symbol(meta_fct2), , OTU) %>%
            summarise(mean_rel_abund = mean(Abundance, na.rm = TRUE) * 100, .groups = "drop") %>%
            group_by(OTU) %>%
            mutate(pool = max(mean_rel_abund) < cutoff) %>%
            mutate(OTU = if_else(pool, "Other", as.character(OTU))) %>%
            ungroup(OTU) %>%
            arrange(desc(mean_rel_abund)) %>%
            mutate(OTU = factor(OTU) %>%
                fct_reorder(mean_rel_abund, max, .desc = TRUE) %>%
                fct_relevel("Other", after = Inf)) %>%
            select(ASV = OTU, !!as.symbol(meta_fct), !!as.symbol(meta_fct2), , mean_rel_abund)

        return(df2)
    }
}
