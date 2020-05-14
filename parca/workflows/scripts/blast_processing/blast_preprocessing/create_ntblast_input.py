# Script for extracting reads from a multifasta file where classed reads are removed.
# Author: Pernilla Ericsson (pernilla.ericsson@gu.se)
# Date: 2020-05-06

"""
Input: 
    classed_path = File with classifies read ids by kraken, kaiju and blast.
    fasta_path = Multifasta file to be filtered.
Params:
    chunk_size = The maximum number of fastas in a file.
Output:
    blast_infiles = The name of the output directory where the files should be placed.
    count = File containing stats about how many files that were printed and the total length printed.
"""

import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from create_blast_input import ExtractRead

classed_path = snakemake.input['kmer_blast']
fasta_path = snakemake.input['kmer_input']

chunk_size = snakemake.params['chunk_size']

out_dir= snakemake.output['blast_infiles']
stats_file=snakemake.output['count']

os.system(f"if [ ! -d {out_dir} ]; then mkdir {out_dir};fi;") 

def main():
    ER=ExtractRead(classed_path, fasta_path)

    classed_list=ER.read_classed_file()
    above_dict=ER.create_classed_dict(classed_list, key_column=0, value_column=0)
    fastas=ER.create_filtered_fasta_dict(above_dict)

    fastas_key_list=list(fastas.keys())
    assembledreadlength=0
    printedfiles=0
    iterator=0
    for chunk in range(0,len(fastas),chunk_size):
            fasta_list=[]
            for i in range(chunk, chunk+chunk_size):
                if i == len(fastas):
                    break
                seq_id=fastas_key_list[i]
                sequence_record = ER.create_seq_record(fastas[seq_id],seq_id)
                fasta_list.append(sequence_record)
                assembledreadlength+=len(fastas[seq_id])
            out_name=os.path.join(out_dir, "ntblast"+"__"+str(iterator)+".fasta")
            ER.write_to_fasta(fasta_list, out_name)
            printedfiles+=1
            iterator+=1

    with open(stats_file,"w") as stats:
        stats.write(f"type\tcount\nPrinted files\t{printedfiles}\nAssembled Read Length\t{assembledreadlength}\n")
        stats.close()

if __name__=="__main__":
    main()
