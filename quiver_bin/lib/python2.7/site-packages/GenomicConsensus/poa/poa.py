# Author: David Alexander
from __future__ import absolute_import, division, print_function

import itertools, logging, math, random
from collections import Counter

import ConsensusCore as cc, numpy as np

from ..utils import readsInWindow, snd, third
from .. import reference
from ..options import options
from ..consensus import Consensus, join
from ..windows import kSpannedIntervals, holes, subWindow
from ..variants import Variant, filterVariants, annotateVariants
from ..Worker import WorkerProcess, WorkerThread
from ..ResultCollector import ResultCollectorProcess, ResultCollectorThread
from ..arrow.utils import variantsFromAlignment

#
# --------------- Configuration ----------------------
#

class PoaConfig(object):
    """
    Poa configuration options
    """
    def __init__(self,
                 aligner="affine",
                 minMapQV=10,
                 minPoaCoverage=3,
                 maxPoaCoverage=100,
                 noEvidenceConsensus="nocall",
                 readStumpinessThreshold=0.1,
                 minReadScore=0.75,
                 minHqRegionSnr=3.75):
        self.aligner                    = aligner
        self.minMapQV                   = minMapQV
        self.minPoaCoverage             = minPoaCoverage
        self.maxPoaCoverage             = maxPoaCoverage
        self.noEvidenceConsensus        = noEvidenceConsensus
        self.readStumpinessThreshold    = readStumpinessThreshold
        self.minReadScore               = minReadScore
        self.minHqRegionSnr             = minHqRegionSnr


#
# -----------  The actual algorithm code -------------
#

def filterAlns(alns, poaConfig):
    """
    Given alns (already clipped to the window bounds), filter out any
    that are deemed insufficiently high-quality for POA.

    By and large we avoid doing any filtering to avoid potential
    reference bias in variant calling.

    However at the moment the POA (and potentially other components)
    breaks when there is a read of zero length.  So we throw away
    reads that are "stumpy", where the aligner has inserted a large
    gap, such that while the alignment technically spans the window,
    it may not have any read content therein:

          Ref   ATGATCCAGTTACTCCGATAAA
          Read  ATG---------------TA-A
          Win.     [              )
    """
    return [ a for a in alns
             if a.readLength >= (poaConfig.readStumpinessThreshold * a.referenceSpan) and
                min(a.hqRegionSnr) >= poaConfig.minHqRegionSnr and
                a.readScore >= poaConfig.minReadScore ]


def variantsAndConfidence(refWindow, refSequence, cssSequence, aligner="affine"):
    """
    Compute the confidence for each position, and compare
    the consensus and reference in this window, returning a list of variants
    """
    refId, refStart, refEnd = refWindow

    if aligner == "affine":
        align = cc.AlignAffine
    else:
        align = cc.Align

    ga = align(refSequence, cssSequence)
    confidence = np.ones((len(cssSequence),), dtype=np.uint8) * 20
    variants = variantsFromAlignment(ga, refWindow, confidence)
    return (confidence, variants)


def consensusAndVariantsForAlignments(refWindow, refSequence, alns, poaConfig):
    """
    Call consensus on this interval---without subdividing the interval
    further.

    Testable!

    Clipping has already been done!
    """
    _, refStart, refEnd = refWindow

    # Compute the POA consensus, which is our initial guess, and
    # should typically be > 99.5% accurate
    fwdSequences = [ a.read(orientation="genomic", aligned=False)
                     for a in alns
                     if a.spansReferenceRange(refStart, refEnd) ]

    try:
        assert len(fwdSequences) >= poaConfig.minPoaCoverage
        p = cc.PoaConsensus.FindConsensus(fwdSequences[:poaConfig.maxPoaCoverage])
    except:
        logging.info("%s: POA could not be generated" % (refWindow,))
        css = Consensus.noCallConsensus(poaConfig.noEvidenceConsensus,
                                        refWindow, refSequence)
        return (css, [])

    poaCss = p.Sequence
    confidence, variants = \
        variantsAndConfidence(refWindow, refSequence, poaCss, poaConfig.aligner)
    css = Consensus(refWindow, poaCss, confidence)

    return (css, variants)


def poaConsensusAndVariants(alnFile, refWindow, referenceContig,
                            depthLimit, poaConfig):
    """
    High-level routine for calling the consensus for a
    window of the genome given an alignment.

    Identifies the coverage contours of the window in order to
    identify subintervals where a good consensus can be called.
    Creates the desired "no evidence consensus" where there is
    inadequate coverage.
    """
    winId, winStart, winEnd = refWindow
    logging.info("POA operating on %s" %
                 reference.windowToString(refWindow))

    if options.fancyChunking:
        # 1) identify the intervals with adequate coverage for poa
        #    consensus; restrict to intervals of length > 10
        alnHits = readsInWindow(alnFile, refWindow,
                                depthLimit=20000,
                                minMapQV=poaConfig.minMapQV,
                                strategy="longest",
                                stratum=options.readStratum,
                                barcode=options.barcode)
        starts = np.fromiter((hit.tStart for hit in alnHits), np.int)
        ends   = np.fromiter((hit.tEnd   for hit in alnHits), np.int)
        intervals = kSpannedIntervals(refWindow, poaConfig.minPoaCoverage,
                                      starts, ends, minLength=10)
        coverageGaps = holes(refWindow, intervals)
        allIntervals = sorted(intervals + coverageGaps)
        if len(allIntervals) > 1:
            logging.info("Usable coverage in %s: %r" %
                         (reference.windowToString(refWindow), intervals))

    else:
        allIntervals = [ (winStart, winEnd) ]

    # 2) pull out the reads we will use for each interval
    # 3) call consensusForAlignments on the interval
    subConsensi = []
    variants = []

    for interval in allIntervals:
        intStart, intEnd = interval
        intRefSeq = referenceContig[intStart:intEnd]
        subWin = subWindow(refWindow, interval)

        windowRefSeq = referenceContig[intStart:intEnd]
        alns = readsInWindow(alnFile, subWin,
                             depthLimit=depthLimit,
                             minMapQV=poaConfig.minMapQV,
                             strategy="longest",
                             stratum=options.readStratum,
                             barcode=options.barcode)
        clippedAlns_ = [ aln.clippedTo(*interval) for aln in alns ]
        clippedAlns = filterAlns(clippedAlns_, poaConfig)

        if len([ a for a in clippedAlns
                 if a.spansReferenceRange(*interval) ]) >= poaConfig.minPoaCoverage:

            logging.debug("%s: Reads being used: %s" %
                          (reference.windowToString(subWin),
                           " ".join([str(hit.readName) for hit in alns])))

            css, variants_ = \
                    consensusAndVariantsForAlignments(subWin, intRefSeq,
                                                      clippedAlns, poaConfig)

            filteredVars =  filterVariants(options.minCoverage,
                                           options.minConfidence,
                                           variants_)
            # Annotate?
            if options.annotateGFF:
                annotateVariants(filteredVars, clippedAlns)

            variants += filteredVars

            # Dump?
            shouldDumpEvidence = \
                ((options.dumpEvidence == "all") or
                 (options.dumpEvidence == "variants") and (len(variants) > 0))
            if shouldDumpEvidence:
                logging.info("POA does not yet support --dumpEvidence")
#                 dumpEvidence(options.evidenceDirectory,
#                              subWin, windowRefSeq,
#                              clippedAlns, css)
        else:
            css = Consensus.noCallConsensus(poaConfig.noEvidenceConsensus,
                                            subWin, intRefSeq)
        subConsensi.append(css)

    # 4) glue the subwindow consensus objects together to form the
    #    full window consensus
    css = join(subConsensi)

    # 5) Return
    return css, variants


#
# --------------  Poa Worker class --------------------
#

class PoaWorker(object):

    @property
    def poaConfig(self):
        return self._algorithmConfig

    def onStart(self):
        random.seed(42)

    def onChunk(self, workChunk):
        referenceWindow  = workChunk.window
        refId, refStart, refEnd = referenceWindow

        refSeqInWindow = reference.sequenceInWindow(referenceWindow)

        # Quick cutout for no-coverage case
        if not workChunk.hasCoverage:
            noCallCss = Consensus.noCallConsensus(self.poaConfig.noEvidenceConsensus,
                                                  referenceWindow, refSeqInWindow)
            return (referenceWindow, (noCallCss, []))

        # General case
        eWindow = reference.enlargedReferenceWindow(referenceWindow,
                                                    options.referenceChunkOverlap)
        _, eStart, eEnd = eWindow

        # We call consensus on the enlarged window and then map back
        # to the reference and clip the consensus at the implied
        # bounds.  This seems to be more reliable thank cutting the
        # consensus bluntly
        refContig = reference.byName[refId].sequence
        refSequenceInEnlargedWindow = refContig[eStart:eEnd]

        #
        # Get the consensus for the enlarged window.
        #
        css_, variants_ = \
            poaConsensusAndVariants(self._inAlnFile, eWindow, refContig,
                                    options.coverage, self.poaConfig)

        #
        # Restrict the consensus and variants to the reference window.
        #
        ga = cc.Align(refSequenceInEnlargedWindow, css_.sequence)
        targetPositions = cc.TargetToQueryPositions(ga)
        cssStart = targetPositions[refStart-eStart]
        cssEnd   = targetPositions[refEnd-eStart]

        cssSequence    = css_.sequence[cssStart:cssEnd]
        cssQv          = css_.confidence[cssStart:cssEnd]
        variants       = [ v for v in variants_
                           if refStart <= v.refStart < refEnd ]

        consensusObj = Consensus(referenceWindow,
                                 cssSequence,
                                 cssQv)

        return (referenceWindow, (consensusObj, variants))


# define both process and thread-based plurality callers
class PoaWorkerProcess(PoaWorker, WorkerProcess): pass
class PoaWorkerThread(PoaWorker, WorkerThread): pass

#
#  --------------------- Plugin API --------------------------------
#

# Pluggable module API for algorithms:
#  - Algorithm lives in a package
#  - Package must never fail to import, even if some of
#    its dependencies are not installed.
#  - Package must provide a main module exposing these top level
#    variables/methods:
#    - name                            = str
#    - availability                    = (bool, str)
#    - configure                      -> options -> cmph5 -> algorithm specific config object;
#                                        (can raise IncompatibleDataException)
#    - slaveFactories                 -> bool -> (class, class)

__all__ = [ "name",
            "availability",
            "configure",
            "slaveFactories" ]

name = "poa"
availability = (True, "OK")

def slaveFactories(threaded):
    if threaded:
        return (PoaWorkerThread,  ResultCollectorThread)
    else:
        return (PoaWorkerProcess, ResultCollectorProcess)

def configure(options, _):
    poaConfig = PoaConfig(aligner=options.aligner,
                          minMapQV=options.minMapQV,
                          noEvidenceConsensus=options.noEvidenceConsensusCall,
                          minReadScore=options.minReadScore,
                          minHqRegionSnr=options.minHqRegionSnr)
    return poaConfig
