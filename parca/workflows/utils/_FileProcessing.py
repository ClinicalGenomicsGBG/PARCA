
# -*- coding: utf-8 -*-

# Written by Pernilla Ericsson <pernilla.ericsson@gu.se> <clinicalgenomics@gu.se>

"""
The class ProcessFiles contains functions for reading files into an object.
"""

import yaml

class ProcessFiles:
    def __init__(self, filename):
        """
        Parameters: 
            self.filename - The file of interest.
        """
        self.filename=filename
    
    def readYaml(self):
        """
        Purpose: Create a dictionary from a yaml-file.
        Parameters: 
            self.filename - The file that is read.
        Returns: A dictionary.
        Comments: 
        """
        with open(self.filename, 'r') as yamlFile:
            yamlDict = yaml.safe_load(yamlFile)

        return yamlDict
