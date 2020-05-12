# Script for extracting reads with the same taxonomic id from a fasta file and write to files called after the taxonomic id the reads are classed to in chunks of a specified size, e.g. "taxid__chunk".
# Author: Pernilla Ericsson (pernilla.ericsson@gu.se)
# Date: 2020-05-06

"""
Input: 
    classed_path = File with read ids classified higher than species by kraken and kaiju.
    fasta_path = Multifasta file to be filtered for higher classed reads.
Params:
    chunk_size = The maximum number of fastas in a file.
Output:
    blast_infiles = The name of the output directory where the files should be placed.
    count = File containing stats about how many files that were printed and the total length printed.
"""


import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

classed_path = snakemake.input['higher']
fasta_path = snakemake.input['kmer_input']

chunk_size = snakemake.params['chunk_size']

out_dir= snakemake.output['blast_infiles']
stats_file=snakemake.output['count']

os.system(f"if [ ! -d {out_dir} ]; then mkdir {out_dir};fi;") 

class ExtractRead:
    def __init__(self, classed_path, fasta_path):
        self.classed_path = classed_path
        self.fasta_path = fasta_path

    def read_classed_file(self):
        with open(self.classed_path, 'r') as open_file:
            classed = [line.strip().split() for line in open_file.readlines() if not line.strip().startswith("seq_id")]
        return(classed)

    def create_classed_dict(self, classed_list):
        above_dict={}
        for classed_read in classed_list:
            try:
                    tmp_list=above_dict[classed_read[1]]
            except:
                    tmp_list=[]

            tmp_list.append(classed_read[0])
            above_dict[classed_read[1]]=tmp_list
        return(above_dict)

    def create_fasta_dict(self):
        fastas = {}
        for fasta in SeqIO.parse(self.fasta_path, 'fasta'):
            fastas[fasta.id] = fasta.seq
        return fastas

    def create_seq_record(self, seq_obj,seq_id):
        sequence_record = SeqRecord(seq_obj, id=seq_id,description= "")
        return(sequence_record)

    def write_to_fasta(self, fasta_list, out_name):
        with open(out_name, "w") as output_handle:
            SeqIO.write(fasta_list, output_handle, "fasta")
    
def main():
    ER=ExtractRead(classed_path, fasta_path)

    classed_list=ER.read_classed_file()
    above_dict=ER.create_classed_dict(classed_list)
    fastas=ER.create_fasta_dict()
    print(list(above_dict.keys()))

    assembledreadlength=0
    printedfiles=0
    for tax_id in above_dict:
        seq_id_list=above_dict[tax_id]
        
        fasta_list=[]
        #chunk_size=2
        iterator=0
        for chunk in range(0,len(seq_id_list),chunk_size):
            for i in range(chunk, chunk+chunk_size):
                if i == len(seq_id_list):
                    break
                else:
                    seq_id=seq_id_list[i]
                    sequence_record = ER.create_seq_record(fastas[seq_id],seq_id)
                fasta_list.append(sequence_record)
                assembledreadlength+=len(fastas[seq_id])
            
            out_name=os.path.join(out_dir, str(tax_id)+"__"+str(iterator))
            ER.write_to_fasta(fasta_list, out_name)
            printedfiles+=1
            iterator+=1

    with open(stats_file,"w") as stats:
        stats.write(f"type\tcount\nPrinted files\t{printedfiles}\nAssembled Read Length\t{assembledreadlength}\n")
        stats.close()

if __name__=="__main__":
    main()
