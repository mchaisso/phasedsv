#################################################################################
# Copyright (c) 2017, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#################################################################################

# Author: Lance Hepler

from __future__ import division, print_function

import time
from textwrap import dedent
from GenomicConsensus import __VERSION__, reference

def vcfVariantFrequency(var, labels):
    if var.frequency1 is None:
        return None
    elif var.isHeterozygous:
        denom = var.frequency1 + var.frequency2
        names = ['frequency{}'.format(label) for label in labels]
        freqs = [getattr(var, name) for name in names if getattr(var, name) is not None]
        return 'AF={}'.format(','.join('{:.3g}'.format(f / denom) for f in freqs))
    else:
        # the frequency is 100%, so no need
        return None

class VariantsVcfWriter(object):

    def __init__(self, f, optionsDict, referenceEntries):
        self._vcfFile = open(f, "w")
        print(dedent('''\
            ##fileformat=VCFv4.3
            ##fileDate={date}
            ##source=GenomicConsensusV{version}
            ##reference={reference}''').format(
                date=time.strftime("%Y%m%d"),
                version=__VERSION__,
                reference="file://" + optionsDict["referenceFilename"],
                ), file=self._vcfFile)
        # reference contigs
        for entry in referenceEntries:
            print("##contig=<ID={name},length={length}>".format(
                name=entry.name,
                length=entry.length
                # TODO(lhepler): evaluate adding md5 hexdigest here on large genomes
                ), file=self._vcfFile)
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=self._vcfFile)

    def writeVariants(self, variants):
        for var in variants:
            pos = var.refStart
            ref = ""
            alt = ""
            labels = (1, 2)
            # insertion or deletion
            if var.refSeq == "" or var.readSeq1 == "" or \
                    (var.isHeterozygous and var.readSeq2 == ""):
                # we're anchored on the previous base so no 0- to 1-indexing
                #   correction required
                ref = var.refPrev + var.refSeq
                if var.isHeterozygous:
                    alt = ",".join(var.readPrev + seq for seq in (var.readSeq1, var.readSeq2))
                else:
                    alt = var.readPrev + var.readSeq1
            # substitution
            else:
                # due to 1-indexing, pos needs to be incremented
                pos += 1
                ref = var.refSeq
                if var.isHeterozygous:
                    alt = ",".join(seq for seq in (var.readSeq1, var.readSeq2))
                    if var.refSeq == var.readSeq1:
                        # first variant is same as wildtype
                        alt = var.readSeq2
                        labels = (2,)
                    elif var.refSeq == var.readSeq2:
                        # second variant is same as wildtype
                        alt = var.readSeq1
                        labels = (1,)
                    else:
                        # both variants differ from wildtype
                        alt = ",".join(seq for seq in (var.readSeq1, var.readSeq2))
                else:
                    alt = var.readSeq1
            freq = vcfVariantFrequency(var=var, labels=labels)
            info = "DP={0}".format(var.coverage)
            if freq:
                info = info + ";" + freq
            print("{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}".format(
                chrom=reference.idToFullName(var.refId),
                pos=pos,
                id=".",
                ref=ref,
                alt=alt,
                qual=var.confidence,
                filter="PASS",
                info=info), file=self._vcfFile)

    def close(self):
        self._vcfFile.close()
