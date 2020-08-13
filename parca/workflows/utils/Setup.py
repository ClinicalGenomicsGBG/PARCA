# -*- coding: utf-8 -*-

# Written by Pernilla Ericsson <pernilla.ericsson@gu.se> <clinicalgenomics@gu.se>

import glob
import re
import sys
import os

"""
Class for determining the run setup parameters.
"""

class Setup:
    def __init__(self,RNA="", raw_sample_dir="", samples="", suffix_fwd=".fastq", suffix_rev=""):
        """
        Parameters:
            self.RNA - The parameter for selecting RNA or not.
            self.raw_sample_dir - The directory name containing the raw samples.
            self.samples - The sample names.
            self.suffix_fwd - The suffix of the sample name in the forward read.
            self.suffix_rev - The suffix of the sample name in the reverse read.
        """
        self.RNA=RNA
        self.raw_sample_dir = raw_sample_dir
        self.samples = samples
        self.suffix_fwd=suffix_fwd
        self.suffix_rev=suffix_rev
    
    def findSamples(self):
        """ 
        Purpose: Generate a sample list.
        Parameters: 
            self.raw_sample_dir - The directory name containing the raw samples.
            self.samples - The sample names.
            self.suffix_fwd - The suffix of the sample name in the forward read.
        Returns: 
            A list of samplenames (only basename without extension).
             - If self.samplenames is empty it searches the self.raw_sample_dir with the specified suffix.
        Comments: 
            Only uses suffix_fwd for determining. The suffix_rev is not needed for determining the samplenames.
        """ 
        if self.samples == "" or self.samples == None:
            SAMPLENAMES, = glob_wildcards(self.raw_sample_dir+"/{sample}"+self.suffix_fwd)
        else:
            SAMPLENAMES=[re.sub(f'{self.suffix_fwd}$', '', sample) for sample in self.samples]

        return SAMPLENAMES

    def settings(self):
        """ 
        Purpose: Determining the parameters from a variable input to a unified parameter naming.
        Parameters: 
            self.RNA - The parameter for selecting RNA or not.
            self.suffix_rev - The suffix of the sample name in the reverse read. Set to None or "" if single end.
            self.findSamples() - The function fo determining sample names.
        Returns: 
            nucleotide - either "DNA" or "RNA".
            sample_type - either "SE" or "PE".
            sample_ids - a list of samplenames (only basename without extension).
        Comments: 
        """ 
        try:
            self.RNA = self.RNA.capitalize() 
        except: 
            pass
        
        if self.RNA == "" or self.RNA == None or self.RNA == False or self.RNA != True or self.RNA=="False" or self.RNA=="No":
            nucleotide="DNA"
        else:
            nucleotide="RNA"
        
        if self.suffix_rev == "" or self.suffix_rev == None:
            sample_type="SE"
            sample_ids=self.findSamples()
        else:
            sample_type="PE"
            sample_ids=self.findSamples()

        return nucleotide, sample_type, sample_ids


# class SetUp:

#     def __init__(self, raw_sample_dir = "", filenames=""):
#         self.raw_sample_dir = raw_sample_dir
#         self.filenames = filenames
        
#     def detect_samples(self):

#         sample_list = []

#         if self.filenames == "search":
#             self.filenames = [os.path.basename(x) for x in glob.glob(self.raw_sample_dir+"/*.f*")]

#         for file in self.filenames:
#             pattern=r'(.fq.gz$)|(.fastq.gz$)|(.fastq$)|(.fq$)'
            
#             matches=re.split(pattern,file)
#             cleaned_matches=[match for match in matches if match]

#             if len(cleaned_matches) <= 1:
#                 print("Error! \n Input sample name did not match input restrictions. \n Input should be in the format \n \t \t <samplename>.<fq or fastq> \n \t \t or <samplename>.<fq or fastq>.gz")
#                 raise SystemExit
#             else:
#                 sample_list.append(cleaned_matches[0])

#         return sample_list


## Previously placed in main
# def find_samples(raw_sample_dir, samples, suffix_fwd):
#     if samples == "" or samples == None:
#         SAMPLENAMES, = glob_wildcards(raw_sample_dir+"/{sample}"+suffix_fwd)
#     else:
#         SAMPLENAMES=[re.sub(f'{suffix_fwd}$', '', sample) for sample in samples]

#     return SAMPLENAMES

# def settings(RNA="", raw_sample_dir="", samples="", suffix_fwd=".fastq", suffix_rev=""):
#     try:
#         RNA = RNA.capitalize() 
#     except: 
#         pass
    
#     if RNA == "" or RNA == None or RNA == False or RNA != True or RNA=="False" or RNA=="No":
#         nucleotide="DNA"
#     else:
#         nucleotide="RNA"
    
#     if suffix_rev == "" or suffix_rev == None:
#         sample_type="SE"
#         sample_ids=find_samples(raw_sample_dir, samples, suffix_fwd)
#     else:
#         sample_type="PE"
#         sample_ids=find_samples(raw_sample_dir, samples, suffix_fwd)

#     return nucleotide, sample_type, sample_ids