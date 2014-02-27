"""Makes stacked bar graph of amino-acid frequencies over time.

Jesse Bloom, 2014."""


import math
import os
import pylab
import matplotlib
import mapmuts.sequtils


def main():
    """Main body of script."""
    infile = 'partitioncounts.txt'
    plotfile = 'aafracs.pdf'
    labelfrac = 0.01 # only label residues that rise above this level in legend

    # read input file
    aminoacids = mapmuts.sequtils.AminoAcids()
    colors = ['b', 'r', 'g', 'c', 'm', 'y'] + ['k'] * len(aminoacids) # order corresponds to aa_order
    colors = colors[ : len(aminoacids)]
    colors.reverse()
    aa_order = [  # smallest to largest corresponds to plotting top to bottom
            'K',
            'Q',
            'E',
            'I',
            'N',
            'T',
            'A',
            'C',
            'D',
            'F',
            'G',
            'H',
            'L',
            'M',
            'P',
            'R',
            'S',
            'V',
            'W',
            'Y',
            ]
    aa_order.reverse()
    aa_indices = dict([(aa_order[i], i) for i in range(len(aa_order))])
    lines = [line for line in open(infile).readlines() if not (line.isspace() or line[0] == '#')]
    dates = []
    counts = []
    aafracs = [[] for aa in aminoacids]
    for line in lines:
        entries = line.split('\t')
        dates.append(float(entries[0]))
        counts.append(int(entries[3]))
        i = 0
        for aa in aminoacids:
            frac = float(entries[4 + i]) / counts[-1]
            aafracs[aa_indices[aa]].append(frac)
            i += 1
    labeled_aas = []
    for i in range(len(aa_order)):
        aa = aa_order[i]
        maxfrac = max(aafracs[i])
        if maxfrac >= labelfrac:
            labeled_aas.append(aa)

    # make plot
    matplotlib.rc('font', size=9)
    fig = pylab.figure(figsize=(3, 1.9))
    (lmargin, bmargin, rmargin, tmargin) = (0.15, 0.22, 0.01, 0.18)
    axis = pylab.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - tmargin - bmargin])
    plotareas = pylab.stackplot(dates, aafracs, colors=colors)
    pylab.ylabel('fraction', size=10)
    pylab.xlabel('date', size=10)
    axis.set_ylim([0, 1.0])
    yticker = matplotlib.ticker.FixedLocator([0.0, 0.5, 1.0])
    axis.yaxis.set_major_locator(yticker)
    axis.set_xlim([min(dates), max(dates)])
    dateticks = [year for year in range(int(math.ceil(min(dates))), int(max(dates)) + 1)]
    xticker = matplotlib.ticker.FixedLocator(dateticks)
    axis.xaxis.set_major_locator(xticker)
    yformatter = matplotlib.ticker.FixedFormatter(['%d' % year for year in dateticks])
    axis.xaxis.set_major_formatter(yformatter)
    patches = []
    i = 0
    assert len(aa_order) == len(plotareas) == len(colors)
    for (aa, color) in zip(aa_order, colors):
        if aa in labeled_aas:
            patch = pylab.Rectangle((0, 0), 1, 1, color=color)
            patches.append(patch)
        i += 1
    assert len(labeled_aas) == len(patches)
    patches.reverse()
    labeled_aas.reverse()
    pylab.legend(patches, labeled_aas, ncol=len(labeled_aas), handlelength=0.8, fontsize=10, handletextpad=0.4, columnspacing=1, bbox_to_anchor=(0.5, 0.98), loc='lower center')

    print "Creating the plot file %s" % plotfile
    pylab.savefig(plotfile)

    print "Now making JPG version."
    jpgfile = "%s.jpg" % os.path.splitext(plotfile)[0]
    os.system('convert -density 500 %s %s' % (plotfile, jpgfile))

if __name__ == '__main__':
    main() # run the script
