import glob
import re
import sys
import os

class SetUp:

    def __init__(self, raw_sample_dir = "", filenames=""):
        self.raw_sample_dir = raw_sample_dir
        self.filenames = filenames
        
    def detect_samples(self):

        sample_list = []

        if self.filenames == "search":
            self.filenames = [os.path.basename(x) for x in glob.glob(self.raw_sample_dir+"/*.f*")]

        for file in self.filenames:
            pattern=r'(.fq.gz$)|(.fastq.gz$)|(.fastq$)|(.fq$)'
            
            matches=re.split(pattern,file)
            cleaned_matches=[match for match in matches if match]

            if len(cleaned_matches) <= 1:
                print("Error! \n Input sample name did not match input restrictions. \n Input should be in the format \n \t \t <samplename>.<fq or fastq> \n \t \t or <samplename>.<fq or fastq>.gz")
                raise SystemExit
            else:
                sample_list.append(cleaned_matches[0])

        return sample_list