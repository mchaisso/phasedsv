# Authors: David Alexander, Lance Hepler
from __future__ import absolute_import, division, print_function

from GenomicConsensus.arrow.utils import allSingleBaseMutations
from GenomicConsensus.variants import Variant
from GenomicConsensus.quiver.diploid import variantsFromAlignment

import numpy as np
import ConsensusCore2 as cc

# IUPAC reference:
#   http://www.bioinformatics.org/sms/iupac.html

_packIupac =   { ("A", "G") : "R" ,
                 ("G", "A") : "R" ,
                 ("C", "T") : "Y" ,
                 ("T", "C") : "Y" ,
                 ("G", "C") : "S" ,
                 ("C", "G") : "S" ,
                 ("A", "T") : "W" ,
                 ("T", "A") : "W" ,
                 ("G", "T") : "K" ,
                 ("T", "G") : "K" ,
                 ("A", "C") : "M" ,
                 ("C", "A") : "M" }

_unpackIupac = { "R" : ("A", "G") ,
                 "Y" : ("C", "T") ,
                 "S" : ("G", "C") ,
                 "W" : ("A", "T") ,
                 "K" : ("G", "T") ,
                 "M" : ("A", "C") }

def packIUPAC(bases):
    return _packIupac[bases]

def unpackIUPAC(iupacCode):
    return _unpackIupac[iupacCode]

def isHeterozygote(base):
    return (base in _unpackIupac)

def packMuts(cssBase, mut1, mut2):
    # Turn two muts (with same Start, End, LengthDiff) into a single mutation to
    # IUPAC.  The no-op mutation is coded as None.
    #
    # Example1: (_, Subs A, Subs T) -> Subs W
    # Example2: (_, Ins A, Ins T)   -> Ins W
    # Example3: (A, None, Subs T)   -> Subs W
    #
    nonNullMut = mut1 or mut2
    start   = nonNullMut.Start()
    mutType = nonNullMut.Type()
    newBase1 = mut1.Bases() if mut1 else cssBase
    newBase2 = mut2.Bases() if mut2 else cssBase
    newBasePacked = packIUPAC((newBase1, newBase2))
    return cc.Mutation(mutType, start, newBasePacked)


def scoresForPosition(ai, pos):
    muts = allSingleBaseMutations(str(ai), positions=[pos])
    noMutScore = [0] * ai.NumReads()
    mutScores_ = [ ai.ReadLLs(mut)
                   for mut in muts ]
    mutScores = np.column_stack([noMutScore] + mutScores_).astype(np.float32)
    return mutScores


def variantsFromConsensus(refWindow, refSequenceInWindow, cssSequenceInWindow,
                          cssQvInWindow=None, siteCoverage=None, aligner="affine",
                          ai=None):
    """
    Compare the consensus and the reference in this window, returning
    a list of variants.

    Uses the integrator to identify heterozygous variants.
    """
    assert (cssQvInWindow is None) == (siteCoverage is None)  # Both or none

    refId, refStart, refEnd = refWindow

    if ai is not None:
        #
        # Hunting diploid variants:
        # 1. find confident heterozygous sites;
        # 2. build a "diploid consensus" using IUPAC encoding
        #    for het sites; mark cssQv accordingly
        # 3. align diploid consensus to reference
        # 4. extract and decorate variants
        #
        assert str(ai) == cssSequenceInWindow
        iupacMutations = []  # List of (Mutation, confidence)
        for pos in xrange(0, ai.Length()):
            ds = cc.IsSiteHeterozygous(scoresForPosition(ai, pos), 40)
            if ds:
                muts = [None] + list(allSingleBaseMutations(cssSequenceInWindow, positions=[pos]))
                mut0 = muts[ds.Allele0]
                mut1 = muts[ds.Allele1]
                cssBase = cssSequenceInWindow[pos]
                packedMut = packMuts(cssBase, mut0, mut1)
                iupacMutations.append((packedMut, 40))

        # Create diploidCss by applying mutations, meanwhile updating the
        # confidence vector accordingly.
        diploidCss = cc.ApplyMutations([pair[0] for pair in iupacMutations],
                                       cssSequenceInWindow)

        diploidQv  = list(cssQvInWindow) if cssQvInWindow is not None else None

        runningLengthDiff = 0
        for (mut, conf) in iupacMutations:
            start = mut.Start() + runningLengthDiff
            end   = mut.End() + runningLengthDiff
            diploidQv[start:end] = [conf]
        assert len(diploidCss) == len(diploidQv)

        cssSequenceInWindow = diploidCss
        cssQvInWindow = diploidQv

    vars = variantsFromAlignment(refWindow,
                                 refSequenceInWindow, cssSequenceInWindow,
                                 cssQvInWindow, siteCoverage)
    return vars
