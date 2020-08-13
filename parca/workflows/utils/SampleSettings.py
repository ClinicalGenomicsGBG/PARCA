# -*- coding: utf-8 -*-

# Written by Pernilla Ericsson <pernilla.ericsson@gu.se> <clinicalgenomics@gu.se>

import glob
import re
import sys
import os
from ErrorMessages import InputError

"""
Class for determining the run setup parameters.
"""

class SampleSettings:
    def __init__(self, sample_path_list=[], RNA="", get_sample_id=True, suffix_regex='(.fq.gz$)|(.fastq.gz$)|(.fastq$)|(.fq$)'):
        """
        suffix regex is only used if get_sample_id is true
        """
        self.sample_path_list=sample_path_list
        self.RNA=RNA
        self.get_sample_id=get_sample_id
        self.suffix_regex=suffix_regex
    
        self.input_count=len(self.sample_path_list)

        if self.input_count == 2:
            self.fwd_path=self.sample_path_list[0]
            self.rev_path=self.sample_path_list[1]
            
            if isinstance(self.fwd_path, str) and isinstance(self.rev_path, str):
                pass
            else:
                raise InputError(f'Input is not of type character.')

        elif self.input_count == 1:
            self.fwd_path=self.sample_path_list[0]

            if isinstance(self.fwd_path, str):
                pass
            else:
                raise InputError(f'Input is not of type character.')
        else:
            raise InputError(f'Input can only have a length of 1 or 2, not {self.input_count}.')

    def findSampleIdPE(self):
        fwd_basename=list(os.path.basename(self.fwd_path))
        rev_basename=list(os.path.basename(self.rev_path))

        fwd_len=len(fwd_basename)
        sample_id_list=[]
        for i in range(fwd_len):
            fwd_letter=fwd_basename[i]
            rev_letter=rev_basename[i]

            if fwd_letter==rev_letter:
                sample_id_list.append(fwd_letter)
            else:
                break
        sample_id_raw="".join(sample_id_list)
        sample_id=re.sub('_$', '', sample_id_raw)

        return sample_id
    
    def findSampleIdSE(self):
        fwd_basename=os.path.basename(self.fwd_path)
        sample_id=re.sub(self.suffix_regex, '', fwd_basename)

        return sample_id

    def getNucleotide(self):
        try:
            self.RNA = self.RNA.capitalize() 
        except: 
            pass
        
        if self.RNA == "" or self.RNA == None or self.RNA == False or self.RNA != True or self.RNA=="False" or self.RNA=="No":
            nucleotide="DNA"
        else:
            nucleotide="RNA"

        return nucleotide
      
    def settings(self):
        
        nucleotide = self.getNucleotide()

        sample_id=""
        if self.input_count == 2:
            sample_type="PE"
            
            if self.get_sample_id == True:
                sample_id=self.findSampleIdPE()
        else:
            sample_type="SE"

            if self.get_sample_id == True:
                sample_id=self.findSampleIdSE()

        return sample_id, sample_type, nucleotide
