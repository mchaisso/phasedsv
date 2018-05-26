# Author: David Alexander
from __future__ import absolute_import, division, print_function

__all__ = [ "dumpEvidence",
            "QuiverEvidence" ]

import logging, os.path, numpy as np
from collections import namedtuple
from itertools import groupby
from bisect import bisect_left, bisect_right
from pbcore.io import FastaReader, FastaWriter
from .utils import scoreMatrix
from .. import reference

def dumpEvidence(evidenceDumpBaseDirectory,
                 refWindow, refSequence, alns,
                 quiverConsensus):
    """This will import h5py at runtime.
    """
    # Format of evidence dump:
    # evidence_dump/
    #   ref000001/
    #     0-1005/
    #       reference.fa
    #       reads.fa
    #       consensus.fa
    #       quiver-scores.h5
    #     995-2005/
    #       ...
    join = os.path.join
    refId, refStart, refEnd = refWindow
    refName = reference.idToName(refId)
    windowDirectory = join(evidenceDumpBaseDirectory,
                           refName,
                           "%d-%d" % (refStart, refEnd))
    logging.info("Dumping evidence to %s" % (windowDirectory,))

    if os.path.exists(windowDirectory):
        raise Exception("Evidence dump does not expect directory %s to exist." % windowDirectory)
    os.makedirs(windowDirectory)
    refFasta       = FastaWriter(join(windowDirectory, "reference.fa"))
    readsFasta     = FastaWriter(join(windowDirectory, "reads.fa"))
    consensusFasta = FastaWriter(join(windowDirectory, "consensus.fa"))

    windowName = refName + (":%d-%d" % (refStart, refEnd))
    refFasta.writeRecord(windowName, refSequence)
    refFasta.close()

    consensusFasta.writeRecord(windowName + "|quiver", quiverConsensus.sequence)
    consensusFasta.close()

    rowNames, columnNames, baselineScores, scores = scoreMatrix(quiverConsensus.mms)
    import h5py
    quiverScoreFile = h5py.File(join(windowDirectory, "quiver-scores.h5"))
    quiverScoreFile.create_dataset("Scores", data=scores)
    vlen_str = h5py.special_dtype(vlen=str)
    quiverScoreFile.create_dataset("RowNames", data=rowNames, dtype=vlen_str)
    quiverScoreFile.create_dataset("ColumnNames", data=columnNames, dtype=vlen_str)
    quiverScoreFile.create_dataset("BaselineScores", data=baselineScores)
    quiverScoreFile.close()
    for aln in alns:
        readsFasta.writeRecord(str(aln.rowNumber),
                               aln.read(orientation="genomic", aligned=False))
    readsFasta.close()


class QuiverEvidence(object):
    """
    An experimental reader class for quiver evidence dumps produced by
    quiver --dumpEvidence

    This will import h5py at runtime.
    """
    Mutation = namedtuple("Mutation", ("Position", "Type", "FromBase", "ToBase"))

    @staticmethod
    def _parseMutName(mutName):
        fields = mutName.split(" ")
        pos = int(fields[0])
        type, fromBase, _, toBase = fields[1:]
        return QuiverEvidence.Mutation(pos, type, fromBase, toBase)

    def __init__(self, path, refStart, consensus, rowNames, colNames, baselineScores, scores):
        self.path           = path
        self.refStart       = refStart
        self.consensus      = consensus
        self.rowNames       = rowNames
        self.colNames       = colNames
        self.baselineScores = baselineScores
        self.scores         = scores
        self.muts           = map(QuiverEvidence._parseMutName, self.colNames)

    @property
    def positions(self):
        return  [ mut.Position for mut in self.muts ]

    @property
    def uniquePositions(self):
        return sorted(list(set(self.positions)))

    @property
    def totalScores(self):
        return self.baselineScores[:, np.newaxis] + self.scores

    @staticmethod
    def load(path):
        if path.endswith("/"): path = path[:-1]

        refWin_ = path.split("/")[-1].split("-")
        refStart = int(refWin_[0])

        with FastaReader(path + "/consensus.fa") as fr:
            consensus = next(iter(fr)).sequence

        import h5py
        with h5py.File(path + "/quiver-scores.h5", "r") as f:
            scores   = f["Scores"].value
            baselineScores = f["BaselineScores"].value
            colNames = f["ColumnNames"].value
            rowNames = f["RowNames"].value
            return QuiverEvidence(path, refStart, consensus,
                                  rowNames, colNames,
                                  baselineScores, scores)

    def forPosition(self, pos):
        posStart = bisect_left(self.positions, pos)
        posEnd   = bisect_right(self.positions, pos)
        return QuiverEvidence(self.path,
                              self.refStart,
                              self.consensus,
                              self.rowNames,
                              self.colNames[posStart:posEnd],
                              self.baselineScores,
                              self.scores[:, posStart:posEnd])


    def justSubstitutions(self):
        colMask = np.array(map(lambda s: ("Sub" in s), self.colNames))
        return QuiverEvidence(self.path,
                              self.refStart,
                              self.consensus,
                              self.rowNames,
                              self.colNames[colMask],
                              self.baselineScores,
                              self.scores[:, colMask])

    def rowNumbers(self):
        with FastaReader(self.path + "/reads.fa") as fr:
            return [ int(ctg.name) for ctg in fr ]
