#!/usr/bin/env python

# DBTW.py - generate DoS from CG Js read in from file
# Down by the water // My lovely daughter // I took her home

# Import our numeric library
import numpy as np
import matplotlib.pyplot as pl

import datetime # current date for log files etc.

from IPython import embed# we do this so we can drop to interactive python for debugging; major Python coolio
 #  # --> embed() <-- just add this to any point in code, and TADA!

### Matplotlib setup
#Pretty colours; via http://blog.olgabotvinnik.com/post/58941062205/prettyplotlib-painlessly-create-beautiful-matplotlib
try:
    import brewer2mpl #'pip install brewer2mpl'
# Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
    colours = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
except ImportError: #If no brewer2mpl library
    #Otherwise, boring built in ones...
    colours='brgcmkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk' # Aaah, we fade to grey (fade to grey)
    print "Hey - no brewer2mpl (try 'pip install brewer2mpl'). Thus a simple palette."

# Matplotlib - initialise figure
fig=pl.figure()
pl.axes().set_aspect('equal') # Square data .'. square figure please

n=1000

# Initialise our Hamiltonian matrix
H = np.zeros ( (n,n) )

# Load off-diagonal elements from text file. Format:-
#  (index_i,int) (index_j,int) (value,double)
filein=np.loadtxt("test.edges",
        dtype=([('f1', '<u4'), ('f2', '<u4'), ('f3',np.float64)]) ) #specify datatypes: 4byte int, 4byte int, float64

for datum in filein:
#    print datum
    H[datum[0],datum[1]]=datum[2] # Populate Hamiltonian with off diagonal elements
    H[datum[1],datum[0]]=datum[2]  # Hermition...

print "Loaded Hamiltonian... "

pl.title("Off-diagonal elements of Hamiltonian")
pl.imshow(H,interpolation='nearest', cmap=pl.cm.PuBuGn) # 2D colourmap of Hamiltonian, nearest interpolation.
pl.colorbar()
pl.show()

# Fill the diagonal elements with site energy; for tight binding
np.fill_diagonal(H, -6.0)

print "Hamiltonian fully setup, time to solve!"
# OK; here we go - let's solve that TB Hamiltonian!
evals,evecs=np.linalg.eigh(H)

#print "Eigenvalues", evals
#print "Eigenvectors", evecs
#print "first Eigenvector..."
#print evecs[0]

pl.title("DoS by TightBinding")
pl.subplot(211) #5 subplots stacked on top of one another

#Plot Eigenvalues with filled-in Eigenvectors^2 / electron probabilities
pl.subplot(211)
for j in range(0,5): #Number of eigenvalues plotted (electron wavefns)
    pl.fill_between(range(n),evals[j],evals[j]+evecs[j]*evecs[j], facecolor=colours[j])
pl.ylabel("Occ %")
#pl.ylim((3.8,5.0))
pl.yticks(fontsize=9)
pl.xticks(visible=False)

#Plot DoS
pl.subplot(212)
pl.hist(evals,100,histtype='stepfilled',color=colours[0])
pl.ylabel("DoS")
pl.yticks(fontsize=9)

pl.tight_layout(pad=0.3)

pl.show() #Displays plots!

print "Lowest Eigenvalue:\n", evals[0]

print "Saving figures...(one moment please)"
now=datetime.datetime.now().strftime("%Y-%m-%d-%Hh%Mm") #String of standardised year-leading time
pl.annotate("%s"%now,xy=(0.75,0.02),xycoords='figure fraction') #Date Stamp in corner

fig.savefig("%s-DBTW.pdf"%now) #Save figures as both PDF and easy viewing PNG (perfect for talks)
fig.savefig("%s-DBTW.png"%now)
#fig.savefig("%s-LongSnakeMoan.ps"%now)    # TODO: check with latest python scripts to see best way to export these for future inclusion in Latex etc.

