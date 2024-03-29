# -*- coding: utf-8 -*-

# Written by Pernilla Ericsson <pernilla.ericsson@gu.se> <clinicalgenomics@gu.se>

import glob
import re
import sys
import os
from workflows.utils.ErrorMessages import InputError

"""
Class for determining the run setup parameters for one sample.
"""

class SampleSettings:
    def __init__(self, sample_path_list=[], get_sample_id=True, suffix_regex='(.fq.gz$)|(.fastq.gz$)|(.fastq$)|(.fq$)'):
        """
        Purpose: Checks validity of input files. The sample_path_list has to consist of 1 or 2 items in a list and be of type str.
        Parameters: 
            sample_path_list - A list with either 2 or 3 items with the full sample path and RNA or DNA stated last.
            get_sample_id - Set to True for searcing for a sample id that does not contain an extension and is similar between two sequences if sample_path_list has 3 items.
            suffix_regex - suffix_regex is only used if get_sample_id is True and consists of a regex that will match the file type extension suffix.
        Returns:
        Comments:
        """
        self.sample_path_list=sample_path_list
        #self.RNA=RNA
        self.get_sample_id=get_sample_id
        self.suffix_regex=suffix_regex
    
        self.input_count=len(self.sample_path_list)

        if self.input_count == 3:
            self.fwd_path=self.sample_path_list[0]
            self.rev_path=self.sample_path_list[1]
            self.nt=self.sample_path_list[2]
            
            if os.path.isfile(self.fwd_path) and os.path.isfile(self.rev_path) and isinstance(self.nt, str):
                pass
            else:
                raise InputError('Input is not of type character.')

        elif self.input_count == 2:
            self.fwd_path=self.sample_path_list[0]
            self.nt=self.sample_path_list[1]

            if os.path.isfile(self.fwd_path) and isinstance(self.nt, str):
                pass
            else:
                raise InputError('Input is not of type character.')

        else:
            raise InputError(f'Input can only have a length of 2 or 3, not {self.input_count}. Input example [<read>, <DNA or RNA>] or [<read fwd>, <read rev>, <DNA or RNA>].')

    def findSampleIdPE(self):
        """
        Purpose: Generates a sample id that is the same for two sample paths by searching for identical letters and stopping when there is a difference.
        Parameters: 
            self.fwd_path - The path to the forward read.
            self.rev_path - The path to the reverse read.
        Returns: A sample id that is the same between two sample paths.
        Comments:
        """
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
        """
        Purpose: Generates a sample id that does not contain the file type extension
        Parameters: 
            self.fwd_path - The path to the forward read.
            self.suffix_regex - A regex that will match the file type extension suffix.
        Returns: A sample id that does not contain the file type extension.
        Comments:
        """
        fwd_basename=os.path.basename(self.fwd_path)
        sample_id=re.sub(self.suffix_regex, '', fwd_basename)

        return sample_id

    def getNucleotide(self):
        """
        Purpose: Determine whether the samples should be interpreted as DNA or RNA.
        Parameters: 
            self.nt - True or False depending on whether the samples are of DNA or RNA.
        Returns: "DNA" if RNA variable is empty, None, False, not True, "False" or "No", otherwise "RNA" is returned.
        Comments:
        """
        try:
            self.nt = self.nt.upper() 
        except: 
            raise InputError('Nucleotide should be specified as RNA or DNA.')

        if self.nt == "RNA":
            nucleotide="RNA"
        elif self.nt == "DNA":
            nucleotide="DNA"
        else:
            raise InputError('Nucleotide should be specified as RNA or DNA.')

        return nucleotide
      
    def settings(self):
        """
        Purpose: Call other functions to generate settings.
        Parameters: 
            self.getNucleotide() - Function explained above.
            self.findSampleIdPE() - Function explained above.
            self.findSampleIdSE() - Function explained above. 
            self.input_count - The number of items in self.sample_path_list.
            self.get_sample_id -  True if the script should look for a sample id, otherwise False.
        Returns:
            sample_id - If self.get_sample_id is False, an empty string is returned, otherwise an inferred sample id.
            sample_type - "PE" or "SE"
            nucleotide - "DNA" or "RNA"
        Comments:
        """
        nucleotide = self.getNucleotide()

        sample_id=""
        if self.input_count == 3:
            sample_type="PE"
            
            if self.get_sample_id == True:
                sample_id=self.findSampleIdPE()

        elif self.input_count == 2:
            sample_type="SE"

            if self.get_sample_id == True:
                sample_id=self.findSampleIdSE()

        return sample_id, sample_type, nucleotide
