# Authors: David Alexander, Lance Hepler
from __future__ import absolute_import, division, print_function

__all__ = [ "ArrowEvidence" ]

import logging, os.path, numpy as np
from collections import namedtuple
from itertools import groupby
from bisect import bisect_left, bisect_right
from pbcore.io import FastaReader, FastaWriter
from .utils import scoreMatrix
from .. import reference

class ArrowEvidence(object):
    """If used at runtime, this will require h5py
    (for storing/loading 'arrow-scores.h5')
    """
    Mutation = namedtuple("Mutation", ("Position", "Type", "FromBase", "ToBase"))

    @staticmethod
    def _parseMutName(mutName):
        fields = mutName.split(" ")
        pos = int(fields[0])
        type, fromBase, _, toBase = fields[1:]
        return ArrowEvidence.Mutation(pos, type, fromBase, toBase)

    def __init__(self, refWindow, consensus, rowNames, colNames, baselineScores, scores):
        assert isinstance(consensus, str)
        self.refWindow      = refWindow # tuple(str, int, int)
        self.consensus      = consensus
        self.rowNames       = rowNames
        self.colNames       = colNames
        self.baselineScores = baselineScores
        self.scores         = scores
        self.muts           = map(ArrowEvidence._parseMutName, self.colNames)

    @staticmethod
    def fromConsensus(css):
        rowNames, colNames, baselineScores, scores = scoreMatrix(css.ai)
        return ArrowEvidence(css.refWindow,
                             css.sequence,
                             rowNames,
                             colNames,
                             baselineScores,
                             scores)
    @property
    def refName(self):
        return self.refWindow[0]

    @property
    def refStart(self):
        return self.refWindow[1]

    @property
    def refEnd(self):
        return self.refWindow[2]

    @property
    def positions(self):
        return  [ mut.Position for mut in self.muts ]

    @property
    def uniquePositions(self):
        return sorted(list(set(self.positions)))

    @property
    def delta(self):
        return self.scores - self.baselineScores[:, np.newaxis]

    @staticmethod
    def load(dir):
        """
        Load an ArrowEvidence from a directory
        """
        if dir.endswith("/"): dir = dir[:-1]
        refStart, refEnd = map(int, dir.split("/")[-1].split("-"))
        refName = dir.split("/")[-2]
        refWindow = (refName, refStart, refEnd)
        with FastaReader(dir + "/consensus.fa") as fr:
            consensus = next(iter(fr)).sequence
        import h5py
        with h5py.File(dir + "/arrow-scores.h5", "r") as f:
            scores   = f["Scores"].value
            baselineScores = f["BaselineScores"].value
            colNames = f["ColumnNames"].value
            rowNames = f["RowNames"].value
            return ArrowEvidence(refWindow, consensus,
                                 rowNames, colNames,
                                 baselineScores, scores)

    def save(self, dir):
        """
        Save this ArrowEvidence to a directory.  The directory will be
        *created* by this method.

        Format of evidence dump:
        evidence_dump/
          ref000001/
            0-1005/
              consensus.fa
              arrow-scores.h5
            995-2005/
            ...
        """
        logging.info("Dumping evidence to %s" % (dir,))
        join = os.path.join
        if os.path.exists(dir):
            raise Exception("Evidence dump does not expect directory %s to exist." % dir)
        os.makedirs(dir)
        #refFasta       = FastaWriter(join(dir, "reference.fa"))
        #readsFasta     = FastaWriter(join(dir, "reads.fa"))
        consensusFasta = FastaWriter(join(dir, "consensus.fa"))
        windowName = self.refName + (":%d-%d" % (self.refStart, self.refEnd))
        #refFasta.writeRecord(windowName, self.refSequence)
        #refFasta.close()

        consensusFasta.writeRecord(windowName + "|arrow", self.consensus)
        consensusFasta.close()

        import h5py
        arrowScoreFile = h5py.File(join(dir, "arrow-scores.h5"))
        arrowScoreFile.create_dataset("Scores", data=self.scores)
        vlen_str = h5py.special_dtype(vlen=str)
        arrowScoreFile.create_dataset("RowNames", data=self.rowNames, dtype=vlen_str)
        arrowScoreFile.create_dataset("ColumnNames", data=self.colNames, dtype=vlen_str)
        arrowScoreFile.create_dataset("BaselineScores", data=self.baselineScores)
        arrowScoreFile.close()

        # for aln in alns:
        #     readsFasta.writeRecord(str(aln.rowNumber),
        #                            aln.read(orientation="genomic", aligned=False))
        # readsFasta.close()


    def forPosition(self, pos):
        posStart = bisect_left(self.positions, pos)
        posEnd   = bisect_right(self.positions, pos)
        return ArrowEvidence(self.refStart,
                             self.consensus,
                             self.rowNames,
                             self.colNames[posStart:posEnd],
                             self.baselineScores,
                             self.scores[:, posStart:posEnd])


    def justSubstitutions(self):
        colMask = np.array(map(lambda s: ("Sub" in s), self.colNames))
        return ArrowEvidence(self.refStart,
                             self.consensus,
                             self.rowNames,
                             self.colNames[colMask],
                             self.baselineScores,
                             self.scores[:, colMask])

    def rowNumbers(self):
        # with FastaReader(self.dir + "/reads.fa") as fr:
        #     return [ int(ctg.name) for ctg in fr ]
        raise NotImplementedError
