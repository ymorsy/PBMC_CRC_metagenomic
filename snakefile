configfile: "config.yaml"

rule targets:
    input: 
        expand("{exp}/{exp}_phyloseq_obj.rds", exp=config["Exp_filters"]),
        expand("{exp}/{exp}_metadata.xlsx", exp=config["Exp_filters"]),
        expand("{exp}/{exp}_clrs.rds", exp=config["Exp_filters"]),


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
