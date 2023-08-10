configfile: "config.yaml"

rule targets:
    input: 
        expand("{exp}/{exp}_phyloseq_obj.rds", exp=config["Exp_filters"]),
        expand("{exp}/{exp}_metadata.xlsx", exp=config["Exp_filters"]),
        expand("{exp}/{exp}_clrs.rds", exp=config["Exp_filters"]),
        expand("{exp}/{exp}_01_Alpha_diversity.pdf", exp=config["Exp_filters"]),
        expand("{exp}/{exp}_02_Beta_diversity.pdf", exp=config["Exp_filters"]),

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
        alpha_diversity="{exp}/{exp}_01_Alpha_diversity.pdf",
        beta_diversity="{exp}/{exp}_02_Beta_diversity.pdf",
    params:
        depth=config["analysis_fctrs"][0]["depth"],
        x=config["analysis_fctrs"][1]["x"],
        facet=config["analysis_fctrs"][2]["facet"],
        alphas=config["diversity_params"][0]["alpha_distances"],
        betas=config["diversity_params"][1]["betas_distances"],
    script:
        "02_diversity.R"