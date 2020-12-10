# -*- coding: utf-8 -*-

# Written by Pernilla Ericsson <pernilla.ericsson@gu.se> <clinicalgenomics@gu.se>

class Error(Exception):
    pass

class InputError(Error):
    """
    Purpose: Input error exception
    Parameters:
        message - Error message
    """

    def __init__(self, message):
        self.message = message