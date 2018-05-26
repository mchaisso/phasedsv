# Author: David Alexander
from __future__ import absolute_import, division, print_function

__all__ = ["loadCmpH5", "loadBam"]

import os.path
from pbcore.io import AlignmentSet


def loadCmpH5(filename, referenceFname, disableChunkCache=False):
    """
    Get a CmpH5Reader object, disabling the chunk cache if requested.
    """
    filename = os.path.abspath(os.path.expanduser(filename))
    return AlignmentSet(filename)

def loadBam(filename, referenceFname):
    filename = os.path.abspath(os.path.expanduser(filename))
    aln = AlignmentSet(filename, referenceFastaFname=referenceFname)
    return aln
