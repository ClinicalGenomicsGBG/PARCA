# -*- coding: utf-8 -*-

# Written by Pernilla Ericsson <pernilla.ericsson@gu.se> <clinicalgenomics@gu.se>

import glob
import re
import sys
import os
from workflows.utils.SampleSettings import SampleSettings

"""
Class for determining the run setup parameters for all samples in a dictionary.
"""

class Setup:
    def __init__(self, sample_dictionary={}, RNA="", get_sample_id=True, suffix_regex='(.fq.gz$)|(.fastq.gz$)|(.fastq$)|(.fq$)'):
        """
        Parameters: 
            self.sample_dictionary - A dictionary with identifiers (could be used as sample id if get_sample_id is set to False) as keys and sample paths in a list as values.
            self.RNA - Set to True if RNA and False if DNA.
            self.get_sample_id - Set to True if the sample id should be inferred from sample paths and False if the dictionary key should be used as sample id.
            self.suffix_regex - suffix_regex is only used if get_sample_id is True and consists of a regex that will match the file type extension suffix.
        Returns:
        Comments:
        """
        self.sample_dictionary=sample_dictionary
        self.RNA=RNA
        self.get_sample_id=get_sample_id
        self.suffix_regex=suffix_regex
    
    def keysList(self, dictionary):
        """
        Purpose: Generate a list of the key names.
        Parameters: 
            dictionary - a dictionary with identifiers as values and sample path lists as values.
        Returns:
            A list of key names from the dictionary.
        Comments:
        """
        keys_list=list(dictionary.keys())

        return keys_list
    
    def generateSettingsLists(self):
        """
        Purpose: Generate run settings for all samples in a dictionary.
        Parameters: 
            self.keysList - A list of keys from the dictionary.
            self.sample_dictionary - A dictionary with identifiers as values and sample path lists as values.
            self.RNA - Set to True if RNA and False if DNA.
            self.get_sample_id - Set to True if the sample id should be inferred from sample paths and False if the dictionary key should be used as sample id.
            self.suffix_regex - self.suffix_regex - suffix_regex is only used if get_sample_id is True and consists of a regex that will match the file type extension suffix.
        Returns:
            sample_id_list - A list of sample ids.
            sample_type_list - A list of sample types, either "PE" or "SE".
            nucleotide_list - A list of "DNA" or "RNA"
        Comments:
        """
        keys_list=self.keysList(self.sample_dictionary)

        sample_id_list=[]
        sample_type_list=[]
        nucleotide_list=[]
        for key in keys_list:
            samples=self.sample_dictionary[key]
            #print(samples)
            
            SS=SampleSettings(samples, self.RNA, self.get_sample_id, self.suffix_regex)
            (sample_id, sample_type, nucleotide) = SS.settings()

            if sample_id=="" and self.get_sample_id!=True:
                sample_id=key
            
            sample_id_list.append(sample_id)
            sample_type_list.append(sample_type)
            nucleotide_list.append(nucleotide)

        return sample_id_list, sample_type_list, nucleotide_list