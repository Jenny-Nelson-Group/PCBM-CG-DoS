#!/usr/bin/env python

# ITIAM.py
# Alone, I emplore ya
# I think I'm a mother
 
print("ITIAM.py - DoS by TB with Polaron self-interaction")
print("call as: #ITIAM.py (sites) (file.xyz) (CellA) (CellB) (CellC)")

# Import our numeric library
import numpy as np
# Matplotlib
import matplotlib.pyplot as pl
# Sys for arg passing
import sys

import datetime # current date for log files etc.
now=datetime.datetime.now().strftime("%Y-%m-%d-%Hh%Mm") #String of standardised year-leading time

def archivefigure(name="default"):
    fig.savefig("%s-Tris_%s.pdf"%(name,ALPHA)) #Save figures as both PDF and easy viewing PNG (perfect for talks)
#fig.savefig("%s-ITIAM_%s.png"%(now,name))

def savedata(name="default"):
    np.savetxt("%s-ITIAM_%s.dat"%(dx,name),polaron_evals,delimiter=' ',newline='\n')

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

# Setup size of system to study...
# if present, read in number of sites from argument {1}
#  TODO: probably redundant as we can count number of lines in XYZ ?
if len(sys.argv) > 1: n = int(sys.argv[1])
else: n=10

# Initialise our Hamiltonian matrix
H = np.zeros ( (n,n) )

# if present, read in coordinate filename from argument {2}
if len(sys.argv) > 2: coordfile = sys.argv[2]
else: coordfile="test.xyz"

# if present, read in cell coordinates from arguments {3,4,5}
if len(sys.argv) > 3: cell=[float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5])]
else: cell=[100,100,100]
#cell=[106.287,106.287,106.287] # Hard coded cell dimensions!!! FIXME
print("Cell dimensions: ",cell)

# Load C60 locations from coordinate file. Format:-
#  X Y Z
#  Assuming angstroms.
locations=np.loadtxt(coordfile) #this defaults to reading them in as floats, which should be fine

locations=locations/cell # scale to fractional coordinates

distancematrix=locations[:,None,...]-locations[None,...] # rolled over 
# Calculate distance matrix with Numpy functional programming methods. 
#  Probably v. memory heavy.

PBCS=False
if (PBCS==True):
    distancematrix[distancematrix<0.5]+=1.0 #minimum image convention
    distancematrix[distancematrix>0.5]-=1.0 #minimum image convention

distancematrix*=cell # scale back to real coordinates
locations*=cell # scale from fractional coordinates to real distances 

H=np.apply_along_axis(np.linalg.norm,2,distancematrix) # distances via linalg norm command on suitables axes
# elements in H are now euler distances between those sites {i,j}

J0=10
LAMBDA=0.6
H=J0*np.exp(-LAMBDA*H) # calculate transfer integrals with isotropic exponential form

print "Generated Hamiltonian... "

#np.fill_diagonal(H, 0.0) #set diagonal elements to zero; so we can always see the off-digaonal elements

# Matplotlib - initialise figure
fig=pl.figure()
pl.axes().set_aspect('equal') # Square data .'. square figure please

pl.title("Off-diagonal elements of Hamiltonian")
pl.imshow(H,interpolation='nearest', cmap=pl.cm.PuBuGn) # 2D colourmap of Hamiltonian, nearest interpolation.
pl.colorbar()
pl.show()

#archivefigure("H")

# Fill the diagonal elements with site energy; for tight binding

#Generate gaussian noise with variance dx and mean 0

print "Hamiltonian fully setup, time to solve!"
# OK; here we go - let's solve that TB Hamiltonian!

np.fill_diagonal(H,-3.7)

Hp=H+0.0 #no copy

dx=0.0
if dx!=0.0:np.fill_diagonal(H,np.random.normal(loc=-3.7,scale=dx,size=n))

ALPHA = 0.5 # some kind of effective electron phonon coupling / dielectric of medium
SCFSTEPS = 3


siteEs=[]
polarons=[]

for i in range(SCFSTEPS): # Number of SCF steps
    evals,evecs=np.linalg.eigh(Hp)
    polaron=evecs[:,0]*evecs[:,0] #lowest energy state electron density
    polarons.append(polaron)
    np.fill_diagonal(Hp,-3.7-ALPHA*polaron)
    Hp_diagonal = np.diagonal(Hp)
    siteEs.append(Hp_diagonal)

print "Hamiltonian solved"
fig=pl.figure()
pl.plot(np.transpose(polarons)) #transposes appended lists so that data is plotted as fn of site
pl.plot(np.transpose(siteEs)+3.7)
pl.legend(range(len(polarons))+range(len(siteEs)))
pl.show()

#archivefigure("SCF")

evals,evecs=np.linalg.eigh(H) # solve final form of Hamiltonian (always computes here even if no SCF steps)

#evals_range = np.max(evals)-np.min(evals)

#print "Range of eigenvalues is: ", evals_range

# TODO: calculate Js from polaron ensemble orbitals

# FIXME: Probably doesn't calculate anything other than spurious numbers

pvals,pvecs=np.linalg.eigh(Hp) #polaron eigenvectors / values

polarons=[]
overlaps=[]
for state in range(n): #:[0,1,2,3]: #range(n): #[1,2,3,500]:
    psi0=evecs[:,state] #.reshape(1,n)
    psi1=pvecs[:,0].reshape(1,n)

#    print "psi0= ",psi0
#    print "psi1= ",psi1
#    print "H= ", H
    J=np.dot(psi0,np.inner(H,psi1))
#    print state,J
    overlaps.append(J)
    polarons.append(state)

#print overlaps

fig=pl.figure()
pl.title("Js by Polaron Orbital Overlap")
pl.plot(polarons,overlaps)
pl.show()

#archivefigure("POO")

# TODO: calculate Js from polaron ensemble orbitals

#print "Eigenvalues", evals
#print "Eigenvectors", evecs
#print "first Eigenvector..."
#print evecs[0]

#Find effective size of polaron

centre = np.zeros(3)
prob = np.zeros(n)
r = np.zeros(n)
polaron_size=np.zeros(n)
cum_prob=np.zeros(n)
sorted_r=np.zeros(n)
num=np.zeros(n)

for i in range(0,1000):

    for j in range(0,n):
        prob[j]=pvecs[j,i]*pvecs[j,i]
        #        centre[0]+=locations[j,0]*prob[j]
        #        centre[1]+=locations[j,1]*prob[j]
        #        centre[2]+=locations[j,2]*prob[j]
    max_index=np.argmax(prob)
    centre[0]=locations[max_index,0]
    centre[1]=locations[max_index,1]
    centre[2]=locations[max_index,2]

    max_prob=max(prob)

    #Calculating how many molecules polarons localised over

    for j in range(0,n):
        r[j]=np.sqrt((centre[0]-locations[j,0])*(centre[0]-locations[j,0])+(centre[1]-locations[j,1])*(centre[1]-locations[j,1])+(centre[2]-locations[j,2])*(centre[2]-locations[j,2]))

    idx = np.argsort(r)
    sorted_r = r[idx]
    sorted_charge = prob[idx]

    cum_prob = np.cumsum(sorted_charge)

    for j in range (0,n):
        if cum_prob[j]>0.95:break

    polaron_size[i] = sorted_r[j]

    for j in range(0,n):
        if prob[j]/max_prob>0.05:num[i]+=1

    centre = np.zeros(3)            #Reset centre coordinates and num

print num[0:10]

#print polaron_size

polaron_evals=np.zeros((n,2))
for i in range (0,n):
   polaron_evals[i,0]=pvals[i]
   polaron_evals[i,1]=polaron_size[i]

savedata("Bis")

#Print number of molecules polaron localised over for first 10 eigenvalues

fig=pl.figure()
pl.bar(pvals[0:1000],polaron_size[0:1000],0.00001)
pl.title("Size of polaron vs eigenvalue")
pl.xlabel("Eigenvalues")
pl.ylabel("Effective size of polaron")
#pl.xlim(-3.88,-3.47)
pl.show()

#archivefigure("size")

fig=pl.figure()

pl.title("DoS by TightBinding")
pl.subplot(311) #3 subplots stacked on top of one another

#Plot Eigenvalues with filled-in Eigenvectors^2 / electron probabilities
pl.subplot(311)

#for j in range(0,5): #[0,n/2]: #Number of eigenvalues plotted (electron wavefns)


#Plot polaron
psi=pvecs[:,0]*pvecs[:,0]
pl.fill_between(range(n),0,psi, facecolor='k')
pl.plot(range(n),pvecs[:,0],color='k')

for j in [0,n/2]: #range(0,5): #Number of eigenvalues plotted (electron wavefns)
    psi=evecs[:,j]*evecs[:,j]
    pl.fill_between(range(n),0,psi, facecolor=colours[j%8])
    pl.plot(range(n),evecs[:,j],color=colours[j%8])
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
pl.hist(pvals,bins=np.linspace(min(pvals),max(pvals),100),histtype='stepfilled',color='r')
pl.hist(evals,bins= np.linspace(min(pvals),max(pvals),100),histtype='stepfilled',color='b', alpha=0.5)
pl.ylabel("DoS")
pl.yticks(fontsize=9)
#pl.xlim(-6.5,-5)


pl.show() #Displays plots!

print "Lowest Eigenvalue:\n", evals[0]

print "Saving figures...(one moment please)"
pl.annotate("%s"%now,xy=(0.75,0.02),xycoords='figure fraction') #Date Stamp in corner

#archivefigure("3fig")

fp=open('eigenvector_balls_pymol.py','w')
fp.write("from pymol.cgo import *    # get constants \nfrom pymol import cmd \n")

psi=pvecs[:,0]*pvecs[:,0] # scale everything relative to max density on first eigenvector

for ei,colour in zip( [0,5,10,50] , [(0,0,1),(0,1,0),(1,1,0),(1,0,0)]):
    psi=pvecs[:,ei]*pvecs[:,ei]
    maxpsi=max(psi)

    fp.write("obj = [\n")

    psisum=0.0
    for i in reversed(np.argsort(psi)): #magic list of sorted array indices
#       print locations[i]
#       print psi[i]
        weight=float(psi[i])/maxpsi #on interval [0,1]
        fp.write("ALPHA, %f,\n" %(weight))
        weight=1
        fp.write("COLOR, %f, %f, %f,\n" %(colour[0]*weight , colour[1]*weight, colour[2]*weight))
        fp.write("SPHERE, %f, %f, %f, 5.0,\n" %(locations[i][0],locations[i][1],locations[i][2]))
 
 #        psisum+=psi[i]
 #       if (psisum>0.90): #only if this amount of electron density or above
 #           print "Eigvec: %d .95 density at %d" %(ei,i)
 #           break
 
    fp.write(" END \n]\ncmd.load_cgo(obj,'EV_%d') \n" %(ei))

