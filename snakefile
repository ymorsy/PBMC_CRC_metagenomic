configfile: "config.yaml"

rule targets:
    input: 
        expand("{exp}/{exp}_phyloseq_obj.rds", exp=config["Exp_filters"][1]),
        expand("{exp}/{exp}_metadata.xlsx", exp=config["Exp_filters"][1]),
        expand("{exp}/{exp}_clrs.rds", exp=config["Exp_filters"][1]),
        expand("{exp}/02_Diversity/01_Alpha_diversity_{exp}.pdf", exp=config["Exp_filters"][1]),
        expand("{exp}/02_Diversity/02_Beta_diversity_{exp}.pdf", exp=config["Exp_filters"][1]),
        expand("{exp}/01_Taxa/{exp}_rel_abundance_level_{l}.xlsx", exp=config["Exp_filters"][1], l=config["lev_tax"][1]),
        expand("{exp}/01_Taxa/{exp}_abs_abundance_level_{l}.xlsx", exp=config["Exp_filters"][1], l=config["lev_tax"][1]),

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

