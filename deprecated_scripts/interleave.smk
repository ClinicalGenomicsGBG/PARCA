configfile: "config.yaml"
import glob

rule all:
    input:
        expand("{results_path}/interleaved_samples/{sample_ids}_interleaved.fq", results_path=config['results_path'], sample_ids=config['sample_ids'])

rule interleave_paired_reads:
    input: 
        fwd= lambda wildcards: glob.glob("{paired_file_path}/{sample_ids}_R1_val_1.fq".format(paired_file_path=config['paired_file_path'], sample_ids=wildcards.sample_ids)),
        rev= lambda wildcards: glob.glob("{paired_file_path}/{sample_ids}_R2_val_2.fq".format(paired_file_path=config['paired_file_path'], sample_ids=wildcards.sample_ids))
    output: 
        interleaved="{results_path}/interleaved_samples/{sample_ids}_interleaved.fq"
    log: "{results_path}/logs/interleave_paired_reads/{sample_ids}_interleave.log"
    #benchmark: "{results_path}/benchmark/interleave_paired_reads/{sample_ids}_interleave.benchmark.txt"
    shell:
        """
        /apps/bio/software/bbmap/37.68/bbmap/reformat.sh in={input.fwd} in2={input.rev} out={output.interleaved} > {log} 2> {log};
        """
    
