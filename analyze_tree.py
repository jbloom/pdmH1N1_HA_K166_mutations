"""Analyzes BEAST ``.trees`` file to get ancestral states at specific sites.

This script takes as input ``.trees`` files built by BEAST that reconstructs
ancestral sequences at each node. It first parses out the identity at 
specific sites of interest. It then uses TreeAnnotator to build a maximum
clade credibility tree.

Jesse Bloom, 2014."""


import os
import re


def NameTaxa(treefile, taxa_to_name):
    """Adds names to specified taxa.

    *treefile* is a ``.trees`` file from BEAST. This whole file is
    read into memory, so typically you want to make it a contain just
    a single tree, such as after running TreeAnnotator so it is not
    too huge.

    *taxa_to_name* is a dictionary keyed by the taxa names used
    in the taxlabels block of *treefile*, and with values giving
    the name that should be assigned to these blocks.

    On completion, *treefile* has been modified so that all taxa in
    *taxa_to_name* are annotated with the indicated name, and all
    taxa not in *taxa_to_name* are annotated with a blank name ('').
    """
    lines = open(treefile).readlines()
    newlines = []
    intaxablock = False
    for line in lines:
        if intaxablock:
            if line.strip() == ';':
                intaxablock = False
            else:
                line = line.strip()
                if line[0] == "'":
                    assert line[-1] == "'"
                    taxa = line[1 : -1]
                else:
                    taxa = line
                if taxa in taxa_to_name:
                    line = "\t\t'%s'[&!name=" % taxa + '"%s"]\n' % taxa_to_name[taxa]
                else:
                    line = "\t\t'%s'[&!name=" % taxa + '""]\n'
        elif line.strip().lower() == 'taxlabels':
            intaxablock = True
        newlines.append(line)
    open(treefile, 'w').write(''.join(newlines))



def MakeSiteAnnotatedTrees(intreefile, outtreefile, sites, seqlabel):
    """Makes ``.trees`` file with only node annotations of specified sites.

    * *intreefile* is the name of input ``.trees`` file.

    * *outtreefile* is the name of the created output ``.trees`` file.

    * *sites* is a list of the sites of interest in 1, 2, ... numbering

    * *seqlabel* is the the label for the sequences in the tree.
    """
    ratematch = re.compile('\[\&rate\=\d+(\.\d+(e|E(\-){0,1}\d+){0,1}){0,1}\]')
    seqmatch = re.compile('\[\&(rate\=\d+(\.\d+(e|E(\-){0,1}\d+){0,1}){0,1}\,){0,1}%s\=\"(?P<seq>[A-Z,\-]+)\"]' % seqlabel)
    if not os.path.isfile(intreefile):
        raise IOError("Cannot find intreefile %s" % intreefile)
    f = open(outtreefile, 'w')
    for line in open(intreefile):
        if line[ : 5] != 'tree ':  # not a tree line
            f.write(line)
            continue
        # we have a tree line, parse it
        line = ratematch.sub('', line)
        m = seqmatch.search(line)
        while m:
            seq = m.group('seq')
            sitestring = ['%d="%s"' % (site, seq[site - 1]) for site in sites]
            replacestring = '[&%s]' % ','.join(sitestring)
            line = line.replace(m.group(0), replacestring)
            m = seqmatch.search(line)
        f.write(line)
    f.close()


def AnnotateSites(inmaxtreefile, outmaxtreefile, sites):
    """Annotates sites based on residue identity for FigTree display.

    *inmaxtreefile* is the name of a file containing a maximum clade
    credibility tree with the probability of the different residues
    at *sites* annotated for each node. This would typically be generated
    by running *MakeSiteAnnotatedTrees* and then TreeAnnotator.

    *outmaxtreefile* is the name of a file created by this method.
    It is the same as *inmaxtreefile* except that node annotations
    are remove except for a single annotation for each site at each node,
    as defined by *sites*. This file is overwritten if it already exists.

    *sites* is a dictionary that defines how the annotations are done.
    It is keyed by an integer for each site that we are annotating.
    For each node, looks for the annotation for that site (such as
    *259="L"*). If the annotated residue at that site is in *sites[site]*,
    then sets that annotation to the number specified by *sites[site]*.
    If the annotated residue that that site is not in *sites[site]*,
    then sets the annotation to the number specified by *sites['other']*.
    """
    annotationmatch = re.compile('\[\&.+?\]')
    sitematches = dict([(site, re.compile('%d\=\"([A-Z\-])\"' % site)) for site in sites.keys()])
    if not os.path.isfile(inmaxtreefile):
        raise IOError("Can't find inmaxtreefile of %s" % inmaxtreefile)
    f = open(outmaxtreefile, 'w')
    for line in open(inmaxtreefile):
        if line[ : 5] != 'tree ':
            f.write(line) # not a tree line, so just write it
            continue
        m = annotationmatch.search(line)
        while m:
            annotation = m.group(0)
            i = line.index(annotation)
            if annotation == '[&R]':
                i = i + len(annotation)
                m = annotationmatch.search(line[i : ])
                continue
            assert annotation.count('[') == 1
            newannotation = []
            for (site, site_d) in sites.iteritems():
                sitem = sitematches[site].findall(annotation)
                if len(sitem) != 1:
                    raise ValueError("Failed to find exactly one annotation for site %d in:\n%s" % (site, annotation))
                aa = sitem[0]
                if aa in sites[site]:
                    newannotation.append('%d=%s' % (site, sites[site][aa]))
                else:
                    newannotation.append('%d=%s' % (site, sites[site]['other']))
            newannotation = '[&%s]' % ','.join(newannotation)
            line = line.replace(annotation, newannotation, 1)
            i = i + len(newannotation)
            m = annotationmatch.search(line[i : ])
        f.write(line)
    f.close()



def main():
    """Main body of script."""
    # input / output variables
    sites = { # keyed by residue number, values define labeling. Using 1, 2, ... numbering
            180:{'K':'K', 'Q':'Q', 'other':'other'},
            }
    taxa_to_name = { # assign names to these taxa
        '2013.94_SEQ1002_cds:AHG96349 A/Colorado/3514/2013 2013/12/09 HA':"A/Colorado/3514/2013",
        '2009.25_SEQ2533_cds:ACP41105 A/California/04/2009 2009/04/01 HA':"A/California/04/2009",
        }
    intreefiles = ['selected_alignedproteins.trees'] # the input *.trees files
    combinedtreefile = 'burninremoved.trees' # combined intreefiles with burnin removed
    outtreefile = 'branchannotated.trees' # a created file in which branches are only annotated with the sites of interest
    seqlabel = 'selected_alignedproteins' # label of the sequences in intreefile
    maxcredtree = 'maxcredtree.trees' # created maximum clade credibility tree
    annotatedmaxcredtree = 'annotated_maxcredtree.trees' # created site-annotated maximum clade credibility tree
    treeannotator = '/Users/jbloom/BEASTv1.8.0/bin/treeannotator' # path to treeannotator program
    logcombiner = '/Users/jbloom/BEASTv1.8.0/bin/logcombiner' # path to logcombiner program
    burnin = 1000000 # remove this many initial states as burn-in

    # combine trees with logcombiner
    if os.path.isfile(combinedtreefile):
        print "\nThe combined trees file %s already exists, so using this file. If you want to create a new file, delete the existing one and re-run this script." % combinedtreefile
    else:
        print "\nConstructing the combined trees file %s with LogCombiner..." % combinedtreefile
        os.system("%s -trees -burnin %d %s %s" % (logcombiner, burnin, ' '.join(intreefiles), combinedtreefile))

    # build site annotated trees file
    if os.path.isfile(outtreefile):
        print "\nThe site annotated trees file %s already exists, so using this file. If you want to create a new file, delete the existing one and re-run this script." % outtreefile
    else:
        print "\nConstructing the site annotated trees file %s from the BEAST trees file %s..." % (outtreefile, combinedtreefile)
        MakeSiteAnnotatedTrees(combinedtreefile, outtreefile, sites.keys(), seqlabel)

    # build max clade credibility tree
    if os.path.isfile(maxcredtree):
        print "The maximum clade credibility tree %s already exists, just using this existing file -- delete it if you want to replace it." % maxcredtree
    else:
        print "Building the maximum clade credibility tree %s with TreeAnnotator." % maxcredtree
        # no burn-in here as was done during log combining
        os.system('%s %s %s' % (treeannotator, outtreefile, maxcredtree))

    # annotate the sites in the max clade credibility tree
    if os.path.isfile(annotatedmaxcredtree):
        print "\nThe site-annotated maximum clade credibility tree %s already exists, just using this existing file. If you want to create a new file, delete the existing one and re-run the script." % annotatedmaxcredtree
    else:
        print "\nConstructing the site-annotated maximum clade credibility tree %s." % annotatedmaxcredtree
        AnnotateSites(maxcredtree, annotatedmaxcredtree, sites)
        NameTaxa(annotatedmaxcredtree, taxa_to_name)


main() # run the script
