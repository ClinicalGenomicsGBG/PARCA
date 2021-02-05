#
# Author: Pernilla Ericsson (pernilla.ericsson@gu.se)
# Date: 2020-05-12

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class ExtractRead:
    """
    Class for extracting reads from a fasta file and write to files in chunks of a specified size, e.g. "taxid__chunk".
    """
    def __init__(self, classed_path, fasta_path):
        self.classed_path = classed_path
        self.fasta_path = fasta_path

    def read_classed_file(self):
        with open(self.classed_path, 'r') as open_file:
            classed = [line.strip().split() for line in open_file.readlines() if not line.strip().startswith("seq_id")]
        return(classed)

    def create_classed_dict(self, classed_list,key_column=1,value_column=0):
        above_dict={}
        for classed_read in classed_list:
            try:
                    tmp_list=above_dict[classed_read[key_column]]
            except:
                    tmp_list=[]

            tmp_list.append(classed_read[value_column])
            above_dict[classed_read[key_column]]=tmp_list
        return(above_dict)

    def create_fasta_dict(self):
        fastas = {}
        for fasta in SeqIO.parse(self.fasta_path, 'fasta'):
            fastas[fasta.id] = fasta.seq
        return fastas

    def create_filtered_fasta_dict(self, classed_dict):
        fastas = {}
        classed_keys = list(classed_dict.keys())
        for fasta in SeqIO.parse(self.fasta_path, 'fasta'):
            if fasta.id not in classed_keys:
                fastas[fasta.id] = fasta.seq
        return fastas

    def create_seq_record(self, seq_obj,seq_id):
        sequence_record = SeqRecord(seq_obj, id=seq_id,description= "")
        return(sequence_record)

    def write_to_fasta(self, fasta_list, out_name):
        with open(out_name, "w") as output_handle:
            SeqIO.write(fasta_list, output_handle, "fasta")