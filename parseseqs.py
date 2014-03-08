"""Parse and aligns HA sequences to examine mutations at specific sites.

This script examines mutations at the following sites:

    * 180 (sequential numbering, 166 in H3 numbering)

    * 170 (sequential numbering, 156 in H3 numbering)

    * 171 (sequential numbering, 157 in H3 numbering)

    * 172 (sequential numbering, 158 in H3 numbering)

It performs the following operations.

1) Aligns the sequences

2) Computes the frequencies of mutations at specified sites as a function of date.

3) Aligns a subset of sequences for phylogenetic analysis.

Script by Jesse Bloom, 2014.
"""


import sys
import os
import re
import random
import datetime
import mapmuts.align
import mapmuts.sequtils


def DateToOrdinal(datestring, refyear):
    """Converts a date string to an ordinal date.

    *datestring* is a date given by a string such as '2007/2/13' (for
    Feb-13-2007), or '2007/2//' if no day is specified, or
    '2007//' if no day or month is specified. The '/' characters can
    also be '-'.

    *refdate* is an integer year from the approximate timeframe we are examining
    which is used to anchor the datestring date on the assumption
    that each year has 365.25 days.

    The returned value is a number (decimal) giving the date. If no
    day is specified, the 15th (halfway through the month) is chosen.
    If no month or day is specified, July 1 (halfway through the
    year) is chosen.

    >>> print "%.2f" % DateToOrdinal('2007/4/27', 1968)
    2007.32

    >>> print "%.2f" % DateToOrdinal('2007/4/', 1968)
    2007.29

    >>> print "%.2f" % DateToOrdinal('2007//', 1968)
    2007.50

    >>> print "%.2f" % DateToOrdinal('2007-4-27', 1968)
    2007.32

    """
    if not isinstance(refyear, int):
        raise ValueError('refyear is not an integer')
    refdate = datetime.date(refyear, 1, 1).toordinal()
    try:
        if '/' in datestring:
            (year, month, day) = datestring.split('/')
        else:
            (year, month, day) = datestring.split('-')
    except ValueError:
        raise ValueError("Invalid datestring of: %s" % str(datestring))
    if year and month and day:
        (year, month, day) = (int(year), int(month), int(day))
        date = datetime.date(year, month, day)
    elif year and month:
        (year, month) = (int(year), int(month))
        date = datetime.date(year, month, 15)
    elif year:
        year = int(year)
        date = datetime.date(year, 7, 1)
    else:
        raise ValueError("Invalid datestring of: %s" % str(datestring))
    return (date.toordinal() - refdate) / 365.25 + refyear


def ParseDate(head, fulldateonly):
    """Parses an influenza sequence header for year.

    *head* is a string specifying a sequence header, in a format such as
    these examples::

        cds:CAA24268 A/Puerto Rico/8/1934 1934// NP
        cds:BAN57584 A/muscovy duck/Vietnam/LBM330/2013 2013/01/07 NP

    This function parses the header and returns the date as a floating
    point number, or *None* if no date can be parsed.

    *fulldateonly* specifies we only return a date if it includes month and
    year. This option can be set to *True* or *False*.
    """
    headmatch = re.compile('^\S+ (?P<strain>A(\/[\w\. \-\'\(\)\?\+]*){1,6}) +(?P<date>(\d*|unknown)\/\d*\/\d*) \S+$')
    m = headmatch.search(head)
    if not m:
        raise ValueError("Failed to match header:\n%s" % head)
    date = m.group('date')
    if 'unknown' in date:
        return None
    fulldatematch = re.compile('^\d+\/\d+\/\d+$')
    if fulldateonly and not fulldatematch.search(date):
        return None
    date = DateToOrdinal(date, 2009)
    return date


def NeedleProtAlignments(refseq, seqs, needlecmd, tempfile='_alignments.temp'):
    """Uses EMBOSS needle to align sequences in *seqs* to *refseq*.

    The sequences in *seqs* and *refseq* should be protein sequences.

    *refseq* is a string giving the reference sequence.

    *seqs* is a list of *(header, sequence)* 2-tuples.

    *needlecmd* is the path to the EMBOSS needle program.

    *tempfile* is the name of a file temporarily created to store
    the alignments, and then deleted.

    Returns the list *prot_alignment*.

    *prot_alignment* contains each of the proteins encoded in 
    *seqs* as a *(header, protsequence)* 2-tuple, with the
    protein sequence aligned to that in *refseq* and with all
    gaps relative to *resfeq* stripped away.
    """
    prots = []
    heads = []
    for (head, seq) in seqs:
        heads.append(head)
        prots.append(seq)
    try:
        mapmuts.align.Needle(refseq, prots, needlecmd, 'protein', tempfile)
        alignments = mapmuts.sequtils.ReadFASTA(tempfile)
    finally:
        if os.path.isfile(tempfile):
            os.remove(tempfile)
    assert len(alignments) == 2 * len(prots) == 2 * len(heads) == 2 * len(seqs)
    prot_alignment = []
    for i in range(len(prots)):
        prot = prots[i]
        head = heads[i]
        assert seqs[i][0] == head
        (refa, prota) = (alignments[2 * i][1], alignments[2 * i + 1][1])
        assert len(refa) == len(prota)
        iref = iprot = 0
        alignedprot = []
        for (aa_ref, aa_prot) in zip(refa, prota):
            assert (aa_ref == '-' or aa_ref == refseq[iref]), "iref = %d, aa_ref = %s, refseq[iref] = %s" % (iref, aa_ref, refseq[iref])
            assert (aa_prot == '-' or aa_prot == prot[iprot])
            if aa_ref == '-' and aa_prot != '-':
                iprot += 1
            elif aa_prot == '-' and aa_ref != '-':
                alignedprot.append(aa_prot)
                iref += 1
            elif aa_ref != '-' and aa_prot != '-':
                alignedprot.append(aa_prot)
                iref += 1
                iprot += 1
            else:
                raise ValueError("Both prots in alignment have gap")
        alignedprot = ''.join(alignedprot)
        assert len(alignedprot) == len(refseq)
        prot_alignment.append((head, alignedprot))
    return prot_alignment



def main():
    """Main body of script."""
    random.seed(1)

    # input / output files and script parameters
    sites = [180, 170, 171, 172] # sites in 1, 2, ... numbering
    yearpartition = 4 # partition each year into this many pieces
    nperpartition = 10 # keep this many per partition for selected sequences, plus ref sequence
    lastdate = 2014.0 # only keep sequences with dates <= this exact number
    needlecmd = 'needle' # path to EMBOSS needle
    refseqfile = 'California2009-HA.fasta' # reference sequence to which we align
    infile = 'HAseqs.fasta' # input file with all sequences
    purgeambiguous = True # remove sequences with ambiguous nucleotides?
    requirelengthmatch = True # remove sequences if protein length differs from refseq?
    fulldateonly = True # only keep sequences with full (year, month, day) dates
    anomalies = [line.strip() for line in open('AnomalousSequences.txt').readlines() if not line.isspace()] # anomalous sequences to exclude
    fullprotalignmentfile = 'all_alignedproteins.fasta' # created file
    selectedprotalignmentfile = 'selected_alignedproteins.fasta' # created file
    partitioncounts = 'partitioncounts_site%d.txt' # created file

    # get and clean the sequences
    refseq = mapmuts.sequtils.ReadFASTA(refseqfile)
    if len(refseq) != 1:
        raise IOError("Failed to read exactly one sequence from refseqfile %s" % refseqfile)
    (refhead, refseq) = refseq[0]
    refprot = (refhead, mapmuts.sequtils.Translate([(refhead, refseq)])[0][1])
    print "\nUsing a reference sequence of %s.\nThis sequence is %d amino acids in length." % (refhead, len(refprot))
    seqs = mapmuts.sequtils.ReadFASTA(infile)
    random.shuffle(seqs) # so they are not in order
    print "\nRead %d sequences from %s." % (len(seqs), infile)
    print "\nNow checking sequences..."
    cleanprots = []
    nanomalous = nlengthmismatched = nambiguous = nundated = npastdate = iseq = 0
    unambiguousmatch = re.compile('^[ATCGatcg]+$')
    for seq in seqs:
        if purgeambiguous and not unambiguousmatch.search(seq[1]):
            nambiguous += 1
            continue
        prot = mapmuts.sequtils.Translate([seq])[0]
        if requirelengthmatch and len(refprot[1]) != len(prot[1]):
            nlengthmismatched += 1
            continue
        anomalous = False
        for anomaly in anomalies:
            if anomaly in seq[0]:
                nanomalous += 1
                anomalous = True
                break
        if anomalous:
            continue
        date = ParseDate(seq[0], fulldateonly=fulldateonly)
        if not date:
            nundated += 1
            continue
        if date > lastdate:
            npastdate += 1
        iseq += 1
        cleanprots.append(("%.2f_SEQ%d_%s" % (date, iseq, prot[0]), prot[1])) 
    print "Removed %d anomalous sequences, %d ambiguous sequences, %d without a full date, %d past the last retained date, and %d sequences of the wrong length." % (nanomalous, nambiguous, nundated, npastdate, nlengthmismatched)
    print "Retained %d sequences." % len(cleanprots)
    print "Now aligning..."
    prots = NeedleProtAlignments(refprot[1], cleanprots, needlecmd, tempfile='_alignments.temp')
    print "Aligned all %d sequences." % len(prots)
    print "\nWriting all aligned proteins to %s" % fullprotalignmentfile
    del cleanprots
    mapmuts.sequtils.WriteFASTA(prots, fullprotalignmentfile)

    # get assigned partitions
    #
    def AssignPartition(date):
        """Returns partition of sequence, *(partitionstart, partitionend)*."""
        yearfracs = 1.0 / float(yearpartition)
        year = int(date)
        for ipartition in range(yearpartition):
            start = year + ipartition * yearfracs
            end = start + yearfracs
            if start <= date <= end:
                return (start, end)
        raise ValueError("Failed to assign partition")
    #
    partitionedseqs = {}
    selectedprots = []
    for (head, prot) in prots:
        date = float(head.split('_')[0])
        partition = AssignPartition(date)
        if partition not in partitionedseqs:
            partitionedseqs[partition] = []
        partitionedseqs[partition].append((head, prot))
        if len(partitionedseqs[partition]) <= nperpartition:
            selectedprots.append((head, prot))
    print "\nRetained %d selected proteins with at most %d per partition." % (len(selectedprots), nperpartition)
    print "Writing selected proteins to %s" % (selectedprotalignmentfile)
    mapmuts.sequtils.WriteFASTA(selectedprots, selectedprotalignmentfile)
    aminoacids = mapmuts.sequtils.AminoAcids()
    partitionedseqs = partitionedseqs.items()
    partitionedseqs.sort()
    for site in sites:
        outfile = partitioncounts % site
        print "\nWriting the partition counts for site %d to %s" % (site, outfile)
        f = open(outfile, 'w')
        f.write('#Counts of different amino acids in each date partiton.\n#date\tstart\tend\tnseqs\t%s\n' % '\t'.join(['%s' % aa for aa in aminoacids]))
        for (partition, seqlist) in partitionedseqs:
            f.write('%.2f\t%.2f\t%.2f\t%d' % ((partition[0] + partition[1]) / 2.0, partition[0], partition[1], len(seqlist)))
            counts = dict([(aa, 0) for aa in aminoacids])
            for (head, seq) in seqlist:
                counts[seq[site - 1]] += 1
            f.write('\t%s\n' % '\t'.join(['%d' % counts[aa] for aa in aminoacids]))
        f.close()





main() # run the script
