# Author: David Alexander
from __future__ import absolute_import, division, print_function

import numpy as np

__all__ = [ "Consensus",
            "QuiverConsensus",
            "ArrowConsensus",
            "totalLength",
            "areContiguous",
            "join" ]

class Consensus(object):
    """
    A multiple sequence consensus corresponding to a
    (reference/scaffold) coordinate region
    """
    def __init__(self, refWindow, sequence, confidence):
        assert (len(sequence) ==
                len(confidence))
        self.refWindow  = refWindow
        self.sequence   = sequence
        self.confidence = confidence

    def __cmp__(self, other):
        return cmp(self.refWindow, other.refWindow)

    #
    # Functions for calling the consensus for regions of inadequate
    # coverage
    #

    @classmethod
    def nAsConsensus(cls, refWin, referenceSequence):
        length = len(referenceSequence)
        seq = np.empty(length, dtype="S1")
        seq.fill("N")
        conf = np.zeros(length, dtype=np.uint8)
        return cls(refWin, seq.tostring(), conf)

    @classmethod
    def referenceAsConsensus(cls, refWin, referenceSequence):
        conf = np.zeros(len(referenceSequence), dtype=np.uint8)
        return cls(refWin, referenceSequence, conf)

    @classmethod
    def lowercaseReferenceAsConsensus(cls, refWin, referenceSequence):
        conf = np.zeros(len(referenceSequence), dtype=np.uint8)
        return cls(refWin, referenceSequence.lower(), conf)

    @classmethod
    def noCallConsensus(cls, noCallStyle, refWin, refSequence):
        d = { "nocall"             : cls.nAsConsensus,
              "reference"          : cls.referenceAsConsensus,
              "lowercasereference" : cls.lowercaseReferenceAsConsensus}
        factory = d[noCallStyle]
        return factory(refWin, refSequence)

    @property
    def hasEvidence(self):
        if isinstance(self, ArrowConsensus)  and self.ai  is None:
            return False
        if isinstance(self, QuiverConsensus) and self.mms is None:
            return False
        return True


class QuiverConsensus(Consensus):
    """
    A QuiverConsensus object carries an additional field, `mms`, which
    is the ConsensusCore MultiReadMutationScorer object, which can be
    used to perform some post-hoc analyses (diploid, sample mixture, etc)
    """
    def __init__(self, refWindow, sequence, confidence, mms=None):
        super(QuiverConsensus, self).__init__(refWindow, sequence, confidence)
        self.mms = mms


class ArrowConsensus(Consensus):
    """
    An ArrowConsensus object carries an additional field, `ai`, which
    is the ConsensusCore2 abstract integrator object, which can be used
    to perform some post-hoc analyses (diploid, sample mixture, etc)
    """
    def __init__(self, refWindow, sequence, confidence, ai=None):
        super(ArrowConsensus, self).__init__(refWindow, sequence, confidence)
        self.ai = ai


def totalLength(consensi):
    """
    Total length of reference/scaffold coordinate windows
    """
    return sum(cssChunk.refWindow[2] - cssChunk.refWindow[1]
               for cssChunk in consensi)

def areContiguous(refWindows):
    """
    Predicate that determines whether the reference/scaffold windows
    are contiguous.
    """
    lastEnd = None
    lastId  = None
    for refWin in refWindows:
        id, start, end = refWin
        if ((lastId is not None and id != lastId) or
            (lastEnd is not None and start != lastEnd)):
            return False
        lastEnd = end
        lastId  = id
    return True

def join(consensi):
    """
    [Consensus] -> Consensus

    String together all the consensus objects into a single consensus.
    Will raise a ValueError if the reference windows are not
    contiguous.
    """
    assert len(consensi) >= 1
    sortedConsensi = sorted(consensi)
    if not areContiguous([cssChunk.refWindow for cssChunk in sortedConsensi]):
        raise ValueError("Consensus chunks must be contiguous")

    joinedRefWindow  = (sortedConsensi[0].refWindow[0],
                        sortedConsensi[0].refWindow[1],
                        sortedConsensi[-1].refWindow[2])
    joinedSeq        = "".join([cssChunk.sequence for cssChunk in sortedConsensi])
    joinedConfidence = np.concatenate([cssChunk.confidence for cssChunk in sortedConsensi])

    return Consensus(joinedRefWindow,
                     joinedSeq,
                     joinedConfidence)


#
# Naming convention for consensus contigs
#
def consensusContigName(referenceName, algorithmName):
    return "%s|%s" % (referenceName, algorithmName)
