import glob

configfile: "config/config.yaml"

singularity: config['singularity_image']

wildcard_object=glob_wildcards(os.path.join('/medstore/logs/pipeline_logfiles/parca', "{run, \w+}", "main_page_stats_{case_control}.tsv"))

print(wildcard_object)

rule all:
    input: expand("{webinterface}/main_page/main_page_stats_all.tsv", webinterface='/medstore/logs/pipeline_logfiles/parca')

rule concat_main_page_stats:
    input: 
        expand("{{webinterface}}/{run}/main_page_stats_{case_control}.tsv", 
            zip,
            run=wildcard_object.run,
            case_control=wildcard_object.case_control)
    output: "{webinterface}/main_page/main_page_stats_all.tsv"
    shell: 
        """
        cat {input} > {output}
        """