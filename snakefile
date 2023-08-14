configfile: "config.yaml"

rule targets:
    input: 
        expand("{exp}/{exp}_phyloseq_obj.rds", exp=config["Exp_filters"][1]),
        expand("{exp}/{exp}_phyloseq_obj_rar.rds", exp=config["Exp_filters"][1]),
        expand("{exp}/{exp}_metadata.xlsx", exp=config["Exp_filters"][1]),
        expand("{exp}/{exp}_metadata.rds", exp=config["Exp_filters"][1]),
        expand("{exp}/{exp}_clrs.rds", exp=config["Exp_filters"][1]),
        expand("{exp}/02_Diversity/01_Alpha_diversity_{exp}.pdf", exp=config["Exp_filters"][1]),
        expand("{exp}/02_Diversity/02_Beta_diversity_{exp}.pdf", exp=config["Exp_filters"][1]),
        expand("{exp}/01_Taxa/{exp}_rel_abundance_level_{l}.xlsx", exp=config["Exp_filters"][1], l=config["lev_tax"][1]),
        expand("{exp}/01_Taxa/{exp}_abs_abundance_level_{l}.xlsx", exp=config["Exp_filters"][1], l=config["lev_tax"][1]),
        expand("{exp}/03_Taxa_bar_grouped/{exp}_grouped_taxa_bars.pdf", exp=config["Exp_filters"][1]),
        expand("{exp}/05_Heatmaps/{exp}_Top_30_abund_rel_{l}.pdf", exp=config["Exp_filters"][1], l=config["lev_tax"][1]),
        expand("{exp}/05_Heatmaps/{exp}_Top_30_abund_abs_{l}.pdf", exp=config["Exp_filters"][1], l=config["lev_tax"][1]),
        expand("{exp}/04_Statistical_analysis/02_Taxa/{exp}_stat_rel_{l}.xlsx", exp=config["Exp_filters"][1], l=config["lev_tax"][1]),
        expand("{exp}/04_Statistical_analysis/02_Taxa/{exp}_stat_rel_sig_{l}.xlsx", exp=config["Exp_filters"][1], l=config["lev_tax"][1]),
        expand("{exp}/04_Statistical_analysis/01_Diversity/{exp}_{alpha}_dis_stat.xlsx", exp=config["Exp_filters"][1], alpha=config["diversity_params"][0]["alpha_distances"]),
        expand("{exp}/04_Statistical_analysis/01_Diversity/{exp}_{alpha}_dis_stat_tukey.xlsx", exp=config["Exp_filters"][1], alpha=config["diversity_params"][0]["alpha_distances"]),
        expand("{exp}/04_Statistical_analysis/01_Diversity/{exp}_{beta}_dis_stat_perm.xlsx", exp=config["Exp_filters"][1], beta=config["diversity_params"][1]["betas_distances"]),

rule create_phyloseq_obj:
    input:
        df=config["df"],
        metadata=config["metadata"],
    params:
        exp="{exp}",
        fct=config["grouped_taxa_bars_params"][1]["meta_fct"],
        col_selection=config["Column_selection"],
    output:
        phyloseq_obj="{exp}/{exp}_phyloseq_obj.rds",
        metadata="{exp}/{exp}_metadata.xlsx",
        metadata_rds="{exp}/{exp}_metadata.rds",
        clrs="{exp}/{exp}_clrs.rds",
    script:
        "01_create_phyloseq_obj.R"

rule diversity:
    input:
        phyloseq_obj="{exp}/{exp}_phyloseq_obj.rds",
        clrs="{exp}/{exp}_clrs.rds",
    output:
        alpha_diversity="{exp}/02_Diversity/01_Alpha_diversity_{exp}.pdf",
        beta_diversity="{exp}/02_Diversity/02_Beta_diversity_{exp}.pdf",
        phyloseq_obj_rar="{exp}/{exp}_phyloseq_obj_rar.rds",
    params:
        alpha_diversity_html="{exp}/02_Diversity/01_Alpha_diversity_",
        beta_diversity_html="{exp}/02_Diversity/02_Beta_diversity_",
        depth=config["analysis_fctrs"][0]["depth"],
        x=config["analysis_fctrs"][1]["x"],
        facet=config["analysis_fctrs"][2]["facet"],
        alphas=config["diversity_params"][0]["alpha_distances"],
        betas=config["diversity_params"][1]["betas_distances"],
    script:
        "02_diversity.R"

rule taxa_count:
    input:
        phyloseq_obj="{exp}/{exp}_phyloseq_obj.rds",
        clrs="{exp}/{exp}_clrs.rds",
    output:
        taxa_count_abs="{exp}/01_Taxa/{exp}_abs_abundance_level_{l}.xlsx",
        taxa_count_rel="{exp}/01_Taxa/{exp}_rel_abundance_level_{l}.xlsx",
    params:
        taxrank=config["lev_tax"][1],
    script:
        "03_taxa_count_export.R"

rule grouped_taxa_bars:
    input:
        phyloseq_obj="{exp}/{exp}_phyloseq_obj.rds",
        clrs="{exp}/{exp}_clrs.rds",
    output:
        grouped_taxa_bars_plot="{exp}/03_Taxa_bar_grouped/{exp}_grouped_taxa_bars.pdf",
    params:
        taxrank=config["lev_tax"][1],
        cutoff=config["grouped_taxa_bars_params"][0]["cutoff"],
        meta_fct=config["grouped_taxa_bars_params"][1]["meta_fct"],
    script:
        "04_grouped_taxa_bars.R"

rule heatmpas:
    input:
        metadata="{exp}/{exp}_metadata.rds",
        clrs="{exp}/{exp}_clrs.rds",
        df_rel="{exp}/01_Taxa/{exp}_rel_abundance_level_{l}.xlsx",
        df_abs="{exp}/01_Taxa/{exp}_abs_abundance_level_{l}.xlsx",
    output:
        hm_rel="{exp}/05_Heatmaps/{exp}_Top_30_abund_rel_{l}.pdf",
        hm_abs="{exp}/05_Heatmaps/{exp}_Top_30_abund_abs_{l}.pdf",
    params:
        taxrank=config["lev_tax"][1],
        exp="{exp}",
        meta_fct=config["grouped_taxa_bars_params"][1]["meta_fct"],
    script:
        "05_heatmap.R"

rule stat:
    input:
        # phyloseq_obj_rar="{exp}/{exp}_phyloseq_obj_rar.rds",
        metadata="{exp}/{exp}_metadata.rds",
        clrs="{exp}/{exp}_clrs.rds",
        df_rel="{exp}/01_Taxa/{exp}_rel_abundance_level_{l}.xlsx",
    output:
        stat_rel="{exp}/04_Statistical_analysis/02_Taxa/{exp}_stat_rel_{l}.xlsx",
        stat_rel_sig="{exp}/04_Statistical_analysis/02_Taxa/{exp}_stat_rel_sig_{l}.xlsx",
        # alpha_dis_stat="{exp}/04_Statistical_analysis/01_Diversity/{exp}_{alpha}_dis_stat.xlsx",
        # alpha_dis_stat_tukey="{exp}/04_Statistical_analysis/01_Diversity/{exp}_{alpha}_dis_stat_tukey.xlsx",
        # betas_dis_stat="{exp}/04_Statistical_analysis/01_Diversity/{exp}_{beta}_dis_stat.xlsx",
    params:
        taxrank=config["lev_tax"][1],
        exp="{exp}",
        meta_fct=config["grouped_taxa_bars_params"][1]["meta_fct"],
        alpha_dis=config["diversity_params"][0]["alpha_distances"],
        betas_dis=config["diversity_params"][1]["betas_distances"],
    script:
        "06_statistical_analysis.R"

rule stat_diversity_alpha:
    input:
        phyloseq_obj_rar="{exp}/{exp}_phyloseq_obj_rar.rds",
        # metadata="{exp}/{exp}_metadata.rds",
        # clrs="{exp}/{exp}_clrs.rds",
    output:
        alpha_dis_stat="{exp}/04_Statistical_analysis/01_Diversity/{exp}_{alpha}_dis_stat.xlsx",
        alpha_dis_stat_tukey="{exp}/04_Statistical_analysis/01_Diversity/{exp}_{alpha}_dis_stat_tukey.xlsx",
        # betas_dis_stat="{exp}/04_Statistical_analysis/01_Diversity/{exp}_{beta}_dis_stat.xlsx",
    params:
        taxrank=config["lev_tax"][1],
        exp="{exp}",
        meta_fct=config["grouped_taxa_bars_params"][1]["meta_fct"],
        # alpha_dis=config["diversity_params"][0]["alpha_distances"],
        alpha_dis="{alpha}",
        # betas_dis=config["diversity_params"][1]["betas_distances"],
    script:
        "06_statistical_analysis_a_diversity.R"

rule stat_diversity_beta:
    input:
        phyloseq_obj_rar="{exp}/{exp}_phyloseq_obj_rar.rds",
        # metadata="{exp}/{exp}_metadata.rds",
        # clrs="{exp}/{exp}_clrs.rds",
    output:
        betas_dis_stat="{exp}/04_Statistical_analysis/01_Diversity/{exp}_{beta}_dis_stat_perm.xlsx",
    params:
        # taxrank=config["lev_tax"][1],
        exp="{exp}",
        meta_fct=config["grouped_taxa_bars_params"][1]["meta_fct"],
        betas_dis="{beta}",
    script:
        "06_statistical_analysis_b_diversity.R"

