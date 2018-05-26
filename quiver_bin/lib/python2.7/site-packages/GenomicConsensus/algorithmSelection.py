# Author: David Alexander
from __future__ import absolute_import, division, print_function

from .utils import die

def bestAlgorithm_(sequencingChemistries):
    """
    Identify the (de novo) consensus algorithm we expect to deliver
    the best results, given the sequencing chemistries represented in
    an alignment file.

    We key off the sequencing chemistries as follows:

    - Just RS chemistry data?  Then use quiver (at least for now, until
      we get arrow > quiver on P6-C4)
    - Else (either all Sequel data, or a mix of Sequel and RS data),
      use arrow.
    - Unknown chemistry found? Return None; we should abort if this is found

    Note that the handling/rejection of chemistry mixtures (including
    mixtures of Sequel and RS data) is left to the algorithm itself.
    """
    if len(sequencingChemistries) == 0:
        raise ValueError("sequencingChemistries must be nonempty list or set")
    chems = set(sequencingChemistries)
    anyUnknown    = "unknown" in chems
    allRS         = all(not(chem.startswith("S/")) for chem in chems) and (not anyUnknown)

    if anyUnknown:
        return None
    elif allRS:
        return "quiver"
    else:
        return "arrow"

def bestAlgorithm(sequencingChemistries):
    ba = bestAlgorithm_(sequencingChemistries)
    if ba is None:
        die("Unidentifiable sequencing chemistry present in dataset.  " +
            "Check if your SMRTanalysis installation is out-of-date.")
    else:
        return ba
