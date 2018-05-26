# Authors: David Alexander, Lance Hepler
from __future__ import absolute_import, division, print_function

import numpy as np, ConfigParser, collections, logging
from glob import glob
from os.path import join
from pkg_resources import resource_filename, Requirement

from GenomicConsensus.utils import die
from GenomicConsensus.arrow.utils import fst, snd
from pbcore.chemistry import ChemistryLookupError
from pbcore.io import CmpH5Alignment
import ConsensusCore2 as cc

__all__ = [ "ArrowConfig" ]


#
#  ArrowConfig: the kitchen sink class of arrow options
#

class ArrowConfig(object):
    """
    Arrow configuration options
    """
    def __init__(self,
                 minMapQV=10,
                 minPoaCoverage=3,
                 maxPoaCoverage=11,
                 mutationSeparation=10,
                 mutationNeighborhood=20,
                 maxIterations=40,
                 noEvidenceConsensus="nocall",
                 computeConfidence=True,
                 readStumpinessThreshold=0.1,
                 minReadScore=0.75,
                 minHqRegionSnr=3.75,
                 minZScore=-3.5,
                 minAccuracy=0.82,
                 maskRadius=0,
                 maskErrorRate=0.5,
                 polishDiploid=False):

        self.minMapQV                   = minMapQV
        self.minPoaCoverage             = minPoaCoverage
        self.maxPoaCoverage             = maxPoaCoverage
        self.mutationSeparation         = mutationSeparation
        self.mutationNeighborhood       = mutationNeighborhood
        self.maxIterations              = maxIterations
        self.noEvidenceConsensus        = noEvidenceConsensus
        self.computeConfidence          = computeConfidence
        self.readStumpinessThreshold    = readStumpinessThreshold
        self.minReadScore               = minReadScore
        self.minHqRegionSnr             = minHqRegionSnr
        self.minZScore                  = minZScore
        self.minAccuracy                = minAccuracy
        self.maskRadius                 = maskRadius
        self.maskErrorRate              = maskErrorRate
        self.polishDiploid              = polishDiploid

    def extractMappedRead(self, aln, windowStart):
        """
        Given a clipped alignment, convert its coordinates into template
        space (starts with 0), bundle it up with its features as a
        MappedRead.
        """
        if isinstance(aln, CmpH5Alignment):
            die("Arrow does not support CmpH5 files!")

        assert aln.referenceSpan > 0

        def baseFeature(featureName):
            if aln.reader.hasBaseFeature(featureName):
                rawFeature = aln.baseFeature(featureName, aligned=False, orientation="native")
                return rawFeature.clip(0,255).astype(np.uint8)
            else:
                return np.zeros((aln.readLength,), dtype=np.uint8)

        name = aln.readName
        chemistry = aln.sequencingChemistry
        strand = cc.StrandType_REVERSE if aln.isReverseStrand else cc.StrandType_FORWARD
        read = cc.Read(name,
                       aln.read(aligned=False, orientation="native"),
                       cc.Uint8Vector(baseFeature("Ipd").tolist()),
                       cc.Uint8Vector(baseFeature("PulseWidth").tolist()),
                       cc.SNR(aln.hqRegionSnr),
                       chemistry)
        return cc.MappedRead(read,
                             strand,
                             int(aln.referenceStart - windowStart),
                             int(aln.referenceEnd   - windowStart))
