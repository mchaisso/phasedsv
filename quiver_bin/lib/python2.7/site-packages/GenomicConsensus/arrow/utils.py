# Authors: David Alexander, Lance Hepler
from __future__ import absolute_import, division, print_function

import numpy as np, itertools, logging, re, sys
from collections import Counter

from GenomicConsensus.variants import *
from GenomicConsensus.utils import *
from GenomicConsensus.consensus import ArrowConsensus
from pbcore.io.rangeQueries import projectIntoRange
import ConsensusCore2 as cc

def uniqueSingleBaseMutations(templateSequence, positions=None):
    """
    Return an iterator over all single-base mutations of a
    templateSequence that result in unique mutated sequences.
    """
    allBases = "ACGT"
    positions = positions or xrange(0, len(templateSequence))
    for tplStart in positions:
        tplBase     = templateSequence[tplStart]
        prevTplBase = templateSequence[tplStart-1] if (tplStart > 0) else None
        # snvs
        for subsBase in allBases:
            if subsBase != tplBase:
                yield cc.Mutation_Substitution(tplStart, subsBase)
        # Insertions---only allowing insertions that are not cognate
        # with the previous base.
        for insBase in allBases:
            if insBase != prevTplBase:
                yield cc.Mutation_Insertion(tplStart, insBase)
        # Deletion--only allowed if refBase does not match previous tpl base
        if tplBase != prevTplBase:
            yield cc.Mutation_Deletion(tplStart, 1)

def allSingleBaseMutations(templateSequence, positions=None):
    """
    Same as ``uniqueSingleBaseMutations``, but no filtering as to
    whether the mutated sequences are unique.
    """
    allBases = "ACGT"
    positions = positions or xrange(0, len(templateSequence))
    for tplStart in positions:
        tplBase = templateSequence[tplStart]
        # snvs
        for subsBase in allBases:
            if subsBase != tplBase:
                yield cc.Mutation_Substitution(tplStart, subsBase)
        # Insertions
        for insBase in allBases:
            yield cc.Mutation_Insertion(tplStart, insBase)
        # Deletion
        yield cc.Mutation_Deletion(tplStart, 1)

def nearbyMutations(mutations, tpl, neighborhoodSize):
    """
    Return mutations nearby the previously-tried mutations
    """
    mutationPositions = map(cc.Mutation.Start, mutations)
    nearbyPositions = set()
    for mp in mutationPositions:
        nearbyPositions.update(range(max(0, mp - neighborhoodSize),
                                     min(len(tpl), mp + neighborhoodSize)))
    return uniqueSingleBaseMutations(tpl, sorted(nearbyPositions))

def bestSubset(mutationsAndScores, separation):
    """
    Given a list of (mutation, score) tuples, this utility method
    greedily chooses the highest scoring well-separated elements.  We
    use this to avoid applying adjacent high scoring mutations, which
    are the rule, not the exception.  We only apply the best scoring one
    in each neighborhood, and then revisit the neighborhoods after
    applying the mutations.
    """
    input = mutationsAndScores[:]
    output = []

    while input:
        best = max(input, key=snd)
        output.append(best)
        nStart = best[0].Start() - separation
        nEnd   = best[0].Start() + separation
        for t in input[:]:
            if nStart <= t[0].Start() <= nEnd:
                input.remove(t)

    return output

def refineConsensus(ai, arrowConfig, polishDiploid=False):
    """
    Given a MultiReadMutationScorer, identify and apply favorable
    template mutations.  Return (consensus, didConverge) :: (str, bool)
    """
    cfg = cc.PolishConfig(arrowConfig.maxIterations,
                          arrowConfig.mutationSeparation,
                          arrowConfig.mutationNeighborhood,
                          polishDiploid)
    if arrowConfig.maskRadius:
        _ = cc.Polish(ai, cfg)
        ai.MaskIntervals(arrowConfig.maskRadius, arrowConfig.maskErrorRate)
    polishResult = cc.Polish(ai, cfg)
    return str(ai), polishResult.hasConverged

def consensusConfidence(ai, positions=None):
    """
    Returns an array of QV values reflecting the consensus confidence
    at each position specified.  If the `positions` argument is
    omitted, confidence values are returned for all positions in the
    consensus (str(ai)).
    """
    return np.array(np.clip(cc.ConsensusQualities(ai), 0, 93), dtype=np.uint8)

IUPACdict = {
    'N' : ('N'),
    'n' : ('n'),
    'A' : ('A'),
    'a' : ('a'),
    'C' : ('C'),
    'c' : ('c'),
    'G' : ('G'),
    'g' : ('g'),
    'T' : ('T'),
    't' : ('t'),
    'M' : ('A', 'C'),
    'm' : ('a', 'c'),
    'R' : ('A', 'G'),
    'r' : ('a', 'g'),
    'W' : ('A', 'T'),
    'w' : ('a', 't'),
    'S' : ('C', 'G'),
    's' : ('c', 'g'),
    'Y' : ('C', 'T'),
    'y' : ('c', 't'),
    'K' : ('G', 'T'),
    'k' : ('g', 't')}

def splitupIUPAC(css):
    listSeq1 = [IUPACdict[x][0] for x in css]
    listSeq2 = [IUPACdict[x][-1] for x in css]

    if listSeq1 == listSeq2:
        # haploid
        readSeq1 = ''.join(listSeq1)
        readSeq2 = None
        freq1 = None
        freq2 = None
    else:
        # diploid
        readSeq1 = ''.join(listSeq1)
        readSeq2 = ''.join(listSeq2)
        freq1 = 0.5
        freq2 = 0.5

    return readSeq1, readSeq2, freq1, freq2

def variantsFromAlignment(a, refWindow, cssQvInWindow=None, siteCoverage=None, effectiveSiteCoverage=None):
    """
    Extract the variants implied by a pairwise alignment to the
    reference.
    """
    variants = []
    refId, refStart, _ = refWindow
    refPos = refStart
    cssPos = 0
    tbl = zip(a.Transcript(),
              a.Target(),
              a.Query())

    # We don't call variants where either the reference or css is 'N'
    grouper = lambda row: "N" if (row[1]=="N" or row[2]=="N") else row[0]
    runs = itertools.groupby(tbl, grouper)

    # track predecessor "anchor" base for vcf output of indel variants
    refPrev = "N"
    cssPrev = "N"

    for code, run in runs:
        assert code in "RIDMN"
        run = list(run)
        ref = "".join(map(snd, run))
        refLen = len(ref) - Counter(ref)["-"]
        css = "".join(map(third, run))
        cssLen = len(css) - Counter(css)["-"]
        variant = None

        if code == "M" or code == "N":
            pass
        elif code == "R":
            assert len(css)==len(ref)
            css, readSeq2, freq1, freq2 = splitupIUPAC(css)
            variant = Variant(refId=refId, refStart=refPos, refEnd=refPos+len(css),
                              refSeq=ref, readSeq1=css, readSeq2=readSeq2,
                              frequency1=freq1, frequency2=freq2,
                              refPrev=refPrev, readPrev=cssPrev)
        elif code == "I":
            css, readSeq2, freq1, freq2 = splitupIUPAC(css)
            variant = Variant(refId=refId, refStart=refPos, refEnd=refPos,
                              refSeq="", readSeq1=css, readSeq2=readSeq2,
                              frequency1=freq1, frequency2=freq2,
                              refPrev=refPrev, readPrev=cssPrev)
        elif code == "D":
            variant = Variant(refId, refPos, refPos + len(ref), ref, "",
                              refPrev=refPrev, readPrev=cssPrev)

        if variant is not None:
            # HACK ALERT: variants at the first and last position
            # are not handled correctly
            if siteCoverage is not None and np.size(siteCoverage) > 0:
                refPos_ = min(refPos-refStart, len(siteCoverage)-1)
                variant.coverage = siteCoverage[refPos_]
            if effectiveSiteCoverage is not None and np.size(effectiveSiteCoverage) > 0:
                refPos_ = min(refPos-refStart, len(siteCoverage)-1)
                variant.annotate("effectiveCoverage", effectiveSiteCoverage[refPos_])
                #import ipdb
                #ipdb.set_trace()
            if cssQvInWindow is not None and np.size(cssQvInWindow) > 0:
                cssPos_ = min(cssPos, len(cssQvInWindow)-1)
                variant.confidence = cssQvInWindow[cssPos_]
            variants.append(variant)

        refPos += refLen
        cssPos += cssLen
        ref = ref.replace("-", "")
        css = css.replace("-", "")
        refPrev = ref[-1] if ref else refPrev
        cssPrev = css[-1] if css else cssPrev

    return variants

def referenceSpanWithinWindow(referenceWindow, aln):
    """
    Helper function for sorting reads by their reference span
    after restriction to a window.
    """
    _, winStart, winEnd = referenceWindow
    return min(winEnd, aln.referenceEnd) - \
           max(winStart, aln.referenceStart)

def lifted(queryPositions, mappedRead):
    """
    Lift a mappedRead into a new coordinate system by using the
    position translation table `queryPositions`
    """
    newStart = queryPositions[mappedRead.TemplateStart]
    newEnd   = queryPositions[mappedRead.TemplateEnd]
    copy = cc.MappedRead(mappedRead)
    copy.TemplateStart = newStart
    copy.TemplateEnd = newEnd
    return copy


_typeMap = { cc.MutationType_INSERTION    : "Ins",
             cc.MutationType_DELETION     : "Del",
             cc.MutationType_SUBSTITUTION : "Sub" }

def _shortMutationDescription(mut, tpl):
    """
    More compact and uniform mutation description strings
    Examples:

    201 Ins . > G
    201 Sub C > T
    201 Del C > .
    """
    _type = _typeMap[mut.Type()]
    _pos = mut.Start()
    _oldBase = "." if mut.IsInsertion() else tpl[_pos]
    _newBase = "." if mut.IsDeletion()  else mut.Bases()
    return "%d %s %s > %s" % (_pos, _type, _oldBase, _newBase)

def scoreMatrix(ai):
    """
    Returns (rowNames, columnNames, S)

    where:
      - S is a matrix where S_{ij} represents the score delta
        of mutation j against read i
      - rowNames[i] is an identifier name for the the read i---presently
        we use the the row number within the cmp.h5, encoded as a string
      - columnNames[j] is an identifier for mutation j, encoding the
        position, type, and base change
    """
    css = str(ai)
    allMutations = sorted(allSingleBaseMutations(css))
    readNames = list(ai.ReadNames())
    numReads = len(readNames)
    shape = (numReads, len(allMutations))
    scoreMatrix = np.zeros(shape)
    for j, mut in enumerate(allMutations):
        mutScores = ai.LLs(mut)
        scoreMatrix[:, j] = mutScores
    baselineScores =  np.array(ai.LLs())
    columnNames = [ _shortMutationDescription(mut, css)
                    for mut in allMutations ]
    rowNames = readNames
    return (rowNames, columnNames, baselineScores, scoreMatrix)

def constructIUPACfreeConsensus(a):
    targetAln = a.Target()
    queryAln = a.Query()

    assert len(targetAln) == len(queryAln)

    pureCss = []

    for i in xrange(len(queryAln)):
        curBase = queryAln[i]

        if curBase != '-':
            if curBase == 'N' or curBase == 'n':
                newBase = curBase
            elif targetAln[i] in IUPACdict[curBase]:
                # construct new sequence with a
                # minimum of divergence, i.e.
                # target: A
                # query:  R (=A+G)
                # -> new base = A
                newBase = targetAln[i]
            else:
                newBase = IUPACdict[curBase][0]

            pureCss.append(newBase)

    newCss = ''.join(pureCss)

    # Be absolutely sure that *really* all
    # ambiguous bases have been removed.
    for i in ('M', 'm', 'R', 'r', 'W', 'w', 'S', 's', 'Y', 'y', 'K', 'k'):
        assert newCss.find(i) == -1

    return newCss

def variantsFromConsensus(refWindow, refSequenceInWindow, cssSequenceInWindow,
                          cssQvInWindow=None, siteCoverage=None, effectiveSiteCoverage=None,
                          aligner="affine", ai=None, diploid=False):
    """
    Compare the consensus and the reference in this window, returning
    a list of variants.
    """
    refId, refStart, refEnd = refWindow

    if diploid:
        align = cc.AlignAffineIupac
    else:
        newCss = ""
        if aligner == "affine":
            align = cc.AlignAffine
        else:
            align = cc.Align

    ga = align(refSequenceInWindow, cssSequenceInWindow)

    if diploid:
        newCss = constructIUPACfreeConsensus(ga)
        # new de-IUPACed sequence still
        # needs to be of same length
        assert(len(newCss) == len(cssSequenceInWindow))

    return variantsFromAlignment(ga, refWindow, cssQvInWindow, siteCoverage, effectiveSiteCoverage), newCss


def filterAlns(refWindow, alns, arrowConfig):
    """
    Given alns (already clipped to the window bounds), filter out any
    that are incompatible with Arrow.

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
             if a.readLength >= (arrowConfig.readStumpinessThreshold * a.referenceSpan) and
                min(a.hqRegionSnr) >= arrowConfig.minHqRegionSnr and
                a.readScore >= arrowConfig.minReadScore ]


def sufficientlyAccurate(mappedRead, poaCss, minAccuracy):
    if minAccuracy <= 0.0:
        return True
    s, e = mappedRead.TemplateStart, mappedRead.TemplateEnd
    tpl = poaCss[s:e]
    if mappedRead.Strand == cc.StrandType_FORWARD:
        pass
    elif mappedRead.Strand == cc.StrandType_REVERSE:
        tpl = reverseComplement(tpl)
    else:
        return False
    aln = cc.AlignLinear(tpl, mappedRead.Seq)
    nErrors = sum(1 for t in aln.Transcript() if t != 'M')
    tlen = len(tpl)
    acc = 1.0 - 1.0 * min(nErrors, tlen) / tlen
    return acc >= minAccuracy


def poaConsensus(fwdSequences, arrowConfig):
    seqLens = [len(seq) for seq in fwdSequences]
    median = np.median(seqLens)
    ordSeqs = sorted(fwdSequences, key=lambda seq: abs(len(seq) - median))
    ordSeqs = ordSeqs[:arrowConfig.maxPoaCoverage]
    cov = len(ordSeqs)
    minCov = 1 if cov < 5 else ((cov + 1) // 2 - 1)
    poaConfig = cc.DefaultPoaConfig(cc.AlignMode_GLOBAL)
    return cc.PoaConsensus.FindConsensus(ordSeqs, poaConfig, minCov)


def consensusForAlignments(refWindow, refSequence, alns, arrowConfig, draft=None, polish=True,
                           alnsUsed=None):
    """
    Call consensus on this interval---without subdividing the interval
    further.

    Returns an ArrowConsensus object.

    Requires that clipping has already been done.

    If `draft` is provided, it will serve as the starting
    point for polishing.  If not, the POA will be used to generate a
    draft starting point.

    If `polish` is False, the arrow polishing procedure will not be
    used, and the draft consensus will be returned.

    `alnsUsed` is an output parameter; if not None, it should be an
    empty list on entry; on return from this function, the list will
    contain the alns objects that were actually used to compute the
    consensus (those not filtered out).
    """
    _, refStart, refEnd = refWindow

    if alnsUsed is not None:
        assert alnsUsed == []

    if draft is None:
        # Compute the POA consensus, which is our initial guess, and
        # should typically be > 99.5% accurate
        fwdSequences = [ a.read(orientation="genomic", aligned=False)
                         for a in alns
                         if a.spansReferenceRange(refStart, refEnd) ]
        assert len(fwdSequences) >= arrowConfig.minPoaCoverage

        try:
            p = poaConsensus(fwdSequences, arrowConfig)
        except Exception:
            logging.info("%s: POA could not be generated" % (refWindow,))
            return ArrowConsensus.noCallConsensus(arrowConfig.noEvidenceConsensus,
                                                  refWindow, refSequence)
        draft = p.Sequence

    ga = cc.Align(refSequence, draft)

    # Extract reads into ConsensusCore2-compatible objects, and map them into the
    # coordinates relative to the POA consensus
    mappedReads = [ arrowConfig.extractMappedRead(aln, refStart) for aln in alns ]
    queryPositions = cc.TargetToQueryPositions(ga)
    mappedReads = [ lifted(queryPositions, mr) for mr in mappedReads ]

    # Load the mapped reads into the mutation scorer, and iterate
    # until convergence.
    ai = cc.Integrator(draft, cc.IntegratorConfig(arrowConfig.minZScore))
    coverage = 0
    for i, mr in enumerate(mappedReads):
        if (mr.TemplateEnd <= mr.TemplateStart or
            mr.TemplateEnd - mr.TemplateStart < 2 or
            mr.Length() < 2):
            continue
        if not sufficientlyAccurate(mr, draft, arrowConfig.minAccuracy):
            tpl = draft[mr.TemplateStart:mr.TemplateEnd]
            if mr.Strand == cc.StrandType_FORWARD:
                pass
            elif mr.Strand == cc.StrandType_REVERSE:
                tpl = reverseComplement(tpl)
            else:
                tpl = "INACTIVE/UNMAPPED"
            logging.debug("%s: skipping read '%s' due to insufficient accuracy, (poa, read): ('%s', '%s')" % (refWindow, mr.Name, tpl, mr.Seq))
            continue
        if ai.AddRead(mr) == cc.State_VALID:
            coverage += 1
            if alnsUsed is not None:
                alnsUsed.append(alns[i])

    if coverage < arrowConfig.minPoaCoverage:
        logging.info("%s: Inadequate coverage to call consensus" % (refWindow,))
        return ArrowConsensus.noCallConsensus(arrowConfig.noEvidenceConsensus,
                                              refWindow, refSequence)

    if not polish:
        confidence = np.zeros(len(draft), dtype=int)
        return ArrowConsensus(refWindow, draft, confidence, ai)

    # Iterate until covergence
    _, converged = refineConsensus(ai, arrowConfig, polishDiploid=False)
    if converged:
        arrowCss = str(ai)
        if arrowConfig.computeConfidence:
            confidence = consensusConfidence(ai)
        else:
            confidence = np.zeros(shape=len(arrowCss), dtype=int)
    else:
        logging.info("%s: Arrow did not converge to MLE" % (refWindow,))
        return ArrowConsensus.noCallConsensus(arrowConfig.noEvidenceConsensus,
                                              refWindow, refSequence)

    if arrowConfig.polishDiploid:
        # additional rounds of diploid polishing
        _, converged = refineConsensus(ai, arrowConfig, polishDiploid=True)
        if converged:
            arrowCss = str(ai)
            if arrowConfig.computeConfidence:
                confidence = consensusConfidence(ai)
            else:
                confidence = np.zeros(shape=len(arrowCss), dtype=int)
        else:
            logging.info("%s: Arrow (diploid) did not converge to optimal solution" % (refWindow,))

    return ArrowConsensus(refWindow,
                          arrowCss,
                          confidence,
                          ai)


def coverageInWindow(refWin, hits):
    winId, winStart, winEnd = refWin
    a = np.array([(hit.referenceStart, hit.referenceEnd)
                  for hit in hits
                  if hit.referenceName == winId])
    if len(a) == 0:
        return np.zeros(winEnd - winStart, dtype=np.uint)
    else:
        tStart = a[:,0]
        tEnd   = a[:,1]
        cov = projectIntoRange(tStart, tEnd, winStart, winEnd)
        return cov
