#!/usr/bin/env python

# ITIAM.py
# Alone, I emplore ya
# I think I'm a mother
 
print("ITIAM.jl - DoS by TB with Polaron self-interaction")

# Import our numeric library
import numpy as np
import scipy.spatial.distance as distance
# Matplotlib
import matplotlib.pyplot as pl
# Sys for arg passing
import sys

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

# Setup size of system to study...
if len(sys.argv) > 1: n = int(sys.argv[1])
else: n=10

# Initialise our Hamiltonian matrix
H = np.zeros ( (n,n) )

if len(sys.argv) > 2: coordfile = sys.argv[2]
else: coordfile="test.xyz"

#cell=[106.287,106.287,106.287] # Hard coded cell dimensions!!! FIXME
cell=[100,100,100]

# Load C60 locations from coordinate file. Format:-
#  X Y Z
#  Assuming angstroms.
locations=np.loadtxt(coordfile)

locations=locations/cell # transpose to fractional coordinates

distancematrix=locations[:,None,...]-locations[None,...] 
# Calculate distance matrix with Numpy functional programming methods. Probably v. memory heavy.

PBCS=True
if (PBCS==True):
    distancematrix[distancematrix<0.5]+=1.0 #minimum image convention
    distancematrix[distancematrix>0.5]-=1.0 #minimum image convention

distancematrix*=cell # back to real coordinates

H=np.linalg.norm(distancematrix,axis=2) # distances via linalg norm command on suitables axes

J0=10
LAMBDA=0.6
H=J0*np.exp(-LAMBDA*H)

print "Generated Hamiltonian... "

np.fill_diagonal(H, 0.0) #set diagonal elements to zero; so we can always see the off-digaonal elements

pl.title("Off-diagonal elements of Hamiltonian")
pl.imshow(H,interpolation='nearest', cmap=pl.cm.PuBuGn) # 2D colourmap of Hamiltonian, nearest interpolation.
pl.colorbar()
pl.show()

# Fill the diagonal elements with site energy; for tight binding
np.fill_diagonal(H, -6.0)

print "Hamiltonian fully setup, time to solve!"
# OK; here we go - let's solve that TB Hamiltonian!

ALPHA = 0.1 # some kind of effective electron phonon coupling / dielectric of medium

siteEs=[]
polarons=[]
for i in range(40):
    evals,evecs=np.linalg.eigh(H)
    polaron=evecs[:,0]*evecs[:,0] #lowest energy state electron density
    print polaron
    polarons.append(polaron)
    siteEs.append(H.diagonal())
    print H.diagonal()
    np.fill_diagonal(H,-6.0-ALPHA*polaron)

fig=pl.figure()
pl.plot(np.transpose(polarons))
pl.plot(np.transpose(siteEs)+6.0)
pl.show()

evals,evecs=np.linalg.eigh(H)

#print "Eigenvalues", evals
#print "Eigenvectors", evecs
#print "first Eigenvector..."
#print evecs[0]

fig=pl.figure()

pl.title("DoS by TightBinding")
pl.subplot(311) #3 subplots stacked on top of one another

#Plot Eigenvalues with filled-in Eigenvectors^2 / electron probabilities
pl.subplot(311)
for j in [0,n/2]: #range(0,5): #Number of eigenvalues plotted (electron wavefns)
    psi=evecs[:,j]*evecs[:,j]
    pl.fill_between(range(n),0,psi, facecolor=colours[j%8])
pl.ylabel("Occupation")
#pl.ylim((3.8,5.0))
pl.yticks(fontsize=9)
pl.xticks(visible=False)

#Plot cumulative eigenvectors / probability density
pl.subplot(312)
for j in [0,1,2,n/4, n/2]: #range(0,5): #Number of eigenvalues plotted (electron wavefns)
    psi=evecs[:,j]*evecs[:,j]    # expectation value
    pl.fill_between(range(n),0,sorted(psi,reverse=True), facecolor=colours[j%8]) #can't see this anymore on large plots...
    psi=sorted(psi,reverse=True) # expectation value, ranked in order (largest first)
    
    psi_sum=[0.0]
    for i in range(len(psi)): # should be a nicer way to do this with functional programming!
        psi_sum.append(psi_sum[-1]+psi[i])
    
    pl.plot(psi_sum, color=colours[j%8])
    pl.plot(y=0.95) # 2 sigma confidence interval? # TODO: why doesn't this work?
pl.ylabel("Cumulative Density")
#pl.ylim((3.8,5.0))
pl.yticks(fontsize=9)
pl.xticks(visible=False)

#Plot DoS
pl.subplot(313)
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

