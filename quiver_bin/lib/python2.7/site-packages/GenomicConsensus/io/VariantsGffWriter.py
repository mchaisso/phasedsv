# Author: David Alexander
from __future__ import absolute_import, division, print_function

import time
from pbcore.io import GffWriter, Gff3Record
from GenomicConsensus import __VERSION__, reference


def gffVariantSeq(var):
    if var.isHeterozygous:
        return "%s/%s" % (var.readSeq1 or ".",
                          var.readSeq2 or ".")
    else:
        return var.readSeq1 or "."

def gffVariantFrequency(var):
    if var.frequency1==None:
        return None
    elif var.isHeterozygous:
        return "{0:.3g}/{1:.3g}".format(var.frequency1, var.frequency2)
    else:
        return "{0:.3g}".format(var.frequency1)

def toGffRecord(var):
    varType  = var.variantType
    gffType  = varType.lower()
    gffStart = (var.refStart + 1) if (var.refSeq != "") else var.refStart
    gffEnd   = var.refEnd         if (var.refSeq != "") else var.refStart
    gffFreq = gffVariantFrequency(var)

    record = Gff3Record(reference.idToFullName(var.refId), gffStart, gffEnd, gffType)
    record.reference  = var.refSeq or "."
    record.variantSeq = gffVariantSeq(var)
    if gffFreq:
        record.frequency  = gffFreq
    record.coverage   = var.coverage
    record.confidence = var.confidence
    if var.annotations:
        for (k, v) in var.annotations:
            record.put(k, v)
    return record

class VariantsGffWriter(object):

    ONTOLOGY_URL = \
        "http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12"

    def __init__(self, f, optionsDict, referenceEntries):
        self._gffWriter = GffWriter(f)
        self._gffWriter.writeHeader("##pacbio-variant-version 2.1")
        self._gffWriter.writeHeader("##date %s" % time.ctime())
        self._gffWriter.writeHeader("##feature-ontology %s" % self.ONTOLOGY_URL)
        self._gffWriter.writeHeader("##source GenomicConsensus %s" % __VERSION__)
        self._gffWriter.writeHeader("##source-commandline %s" % optionsDict["shellCommand"])
        self._gffWriter.writeHeader("##source-alignment-file %s" % optionsDict["inputFilename"])
        self._gffWriter.writeHeader("##source-reference-file %s" % optionsDict["referenceFilename"])
        # Reference groups.
        for entry in referenceEntries:
            self._gffWriter.writeHeader("##sequence-region %s 1 %d" \
                                            % (entry.name, entry.length))

    def writeVariants(self, variants):
        for var in variants:
            self._gffWriter.writeRecord(toGffRecord(var))

    def close(self):
        self._gffWriter.close()
