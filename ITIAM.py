#!/usr/bin/env python

# ITIAM.py
# Alone, I emplore ya
# I think I'm a mother

# Import our numeric library
import numpy as np
# Matplotlib
import matplotlib.pyplot as pl
# Sys for arg passing
import sys

import pickle

import matplotlib.animation as animation

from sys import exit

import datetime # current date for log files etc.
now=datetime.datetime.now().strftime("%Y-%m-%d-%Hh%Mm") #String of standardised year-leading time

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
 
def archivefigure(name="default"):
    fig.savefig("%s_%s_%s_%s.pdf"%(name,dx,ALPHA,state))


def savedata(name="default"):
    np.savetxt("%s_%s_%s_%s.dat"%(name,dx,ALPHA,state),pvals,delimiter=' ',newline='\n')

#Nearest neighbours only
def nearestneighbours(Ham):
    for i in range (0,n):
        for j in range (0,n):
            if Ham[i,j]<0.001: Ham[i,j]=0.0
    return Ham


#Plot off-diagonal elements
def offdiagonal(Ham):
    np.fill_diagonal(Ham, 0.0) #set diagonal elements to zero; so we can always see the off-digaonal elements
    
    # Matplotlib - initialise figure
    fig=pl.figure()
    pl.axes().set_aspect('equal') # Square data .'. square figure please
    pl.title("Off-diagonal elements of Hamiltonian")
    pl.imshow(H,interpolation='nearest', cmap=pl.cm.PuBuGn) # 2D colourmap of Hamiltonian, nearest interpolation.
    pl.colorbar()
    #pl.show()
    
    #archivefigure("H")
    return Ham

# Fill the diagonal elements with site energy; for tight binding
def filldiagonal(Ham):
    np.fill_diagonal(Ham,-3.7)
    Ham_p=H+0.0 #no copy
    if dx!=0.0:np.fill_diagonal(Ham_p,np.random.normal(loc=-3.7,scale=dx,size=n))
    return Ham,Ham_p


#Self-consistent polaron generator to model self-trapping of polaron
def SCpolarongenerator(Ham,Ham_p,):
    SCFSTEPS = 20
    
    siteEs=[]
    polarons=[]
    overlaps=[]
    
    max_overlap_idx=state
    
    for i in range(SCFSTEPS): # Number of SCF steps
        evals,evecs=np.linalg.eigh(Ham_p)
        polaron=evecs[:,max_overlap_idx]*evecs[:,max_overlap_idx] #lowest energy state electron density
        polarons.append(polaron)
        Hp_diagonal = np.diagonal(Ham_p)
        siteEs.append(Hp_diagonal)
        np.fill_diagonal(Ham_p,Hp_diagonal-ALPHA*polaron)
        pvals,pvecs=np.linalg.eigh(Ham_p)
        for j in range(0,n):
            psi0=evecs[:,j]
            psi1=pvecs[:,max_overlap_idx].reshape(1,n)
            J=np.dot(psi0,np.inner(Ham_p,psi1))
            #print J
            overlaps.append(J)
            max_overlap_idx=np.argmax(np.absolute(overlaps))
        #print max_overlap_idx
        overlaps=[]
        for j in range(0,n):
            if pvecs[j,state]*pvecs[j,state]>0.99:break
                #print Ham
    
    return Ham,Ham_p


# solve final form of Hamiltonian (always computes here even if no SCF steps)
def solveHandHp(Ham,Ham_p):
    Evals,Evecs=np.linalg.eigh(Ham)
    Pvals,Pvecs=np.linalg.eigh(Ham_p)
    return Evals,Evecs,Pvals,Pvecs
#print Hp


def plotoverlaps(Evecs,Pvecs):
    polarons=[]
    overlaps=[]
    for state in range(n): #:[0,1,2,3]: #range(n): #[1,2,3,500]:
        psi0=Evecs[:,state] #.reshape(1,n)
        psi1=Pvecs[:,0].reshape(1,n)
        J=np.dot(psi0,np.inner(H,psi1))
        overlaps.append(J)
        polarons.append(state)
    
    fig=pl.figure()
    pl.title("Js by Polaron Orbital Overlap")
    pl.plot(polarons,overlaps)
#pl.show()

#    archivefigure("POO")


#Calculating time evolution of wavefunction
def timeevolution(Evecs,Evals,Pvecs):
    psi0=Evecs[:,0]              #Original wavefunction
    psi1=Pvecs[:,0]              #Localised wavefunction to propagate
    
    timesteps=1000
    
    psi_t_occ=np.zeros((n,timesteps))
    psi0_occ=np.zeros(n)
    psi1_occ=np.zeros(n)
    
    hbar=1
    z=1j
    
    psi_t=np.zeros((n,timesteps))
    
    for t in range(0,timesteps):
        for i in range(0,n):
            a=np.inner(Evecs[:,i],psi1)
            psi_t[:,t]+=a*np.exp(-z*Evals[i]*dt)*Evecs[:,i]

    
        psi_t_norm=np.linalg.norm(psi_t[:,t])      #Normalise
        psi_t[:,t]=psi_t[:,t]/psi_t_norm

    
        for j in range(0,n):
            psi_t_occ[j,t]=psi_t[j,t]*np.conj(psi_t[j,t])
            psi0_occ[j]=psi0[j]*np.conj(psi0[j])
            psi1_occ[j]=psi1[j]*np.conj(psi1[j])
    
    fig=pl.figure()           #Plot original wavefunction, localised wavefunction and time evolved wavefunction
    
    #pl.fill_between(range(n),0,psi0_occ, facecolor='y',alpha=0.5)
    #pl.fill_between(range(n),0,psi1_occ,facecolor='r',alpha=0.5)
    
    for t in range(0,timesteps):
        pl.fill_between(range(n),0,psi_t_occ[:,t],facecolor='b',alpha=0.1)

    fig.savefig("%s_%s_%s_time.pdf"%(dx,ALPHA,state))


def calcocc(T,Evals,Evecs,Pvecs):
    
    psi0=Evecs[:,0]              #Original wavefunction
    psi1=Pvecs[:,0]              #Localised wavefunction to propagate
    
    hbar=1
    z=1j
    
    psi_t=np.zeros(n)
    psi_t_occ=np.zeros(n)
    

    for i in range(0,n):
        a=np.inner(Evecs[:,i],psi1)
        psi_t+=a*np.exp(-z*Evals[i]*T/hbar)*Evecs[:,i]
        
        
    psi_t_norm=np.linalg.norm(psi_t)      #Normalise
    psi_t=psi_t/psi_t_norm
        
        
    for j in range(0,n):
        psi_t_occ[j]=psi_t[j]*np.conj(psi_t[j])

    
    return psi_t_occ


def Animate():

    fig=pl.figure()

    xlim=(0,n)
    ylim=(0,1)

    ax = pl.axes(xlim=xlim,ylim=ylim)
    psi_x_line, = ax.plot([],[],c='r')

    def init():
        psi_x_line.set_data([], [])
        return psi_x_line,

    def animate(i):
        x=np.linspace(0,n,n)
        y=calcocc(dt*i,evals,evecs,pvecs)
        psi_x_line.set_data(x,y)

        return psi_x_line,

    anim=animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=10,blit=False)

    pl.show()

    anim.save('wavefunction_propagation.mp4', writer='ffmpeg', fps=30, extra_args=['-vcodec', 'libx264'])


#Find effective size of polaron and no. of molecules polaron localised over
def sizeofpolaron(Pvals,Pvecs):
    centre = np.zeros(3)
    prob = np.zeros(n)
    r = np.zeros(n)
    Polaron_size=np.zeros(n)
    cum_prob=np.zeros(n)
    sorted_r=np.zeros(n)
    Num=np.zeros(n)
    
    for i in range(0,n):
        
        for j in range(0,n):
            prob[j]=Pvecs[j,i]*Pvecs[j,i]
        max_index=np.argmax(prob)
        centre[0]=locations[max_index,0]
        centre[1]=locations[max_index,1]
        centre[2]=locations[max_index,2]
        
        max_prob=max(prob)
        #print max_prob
        #Calculating how many molecules polarons localised over
        
        for j in range(0,n):
            r[j]=np.sqrt((centre[0]-locations[j,0])*(centre[0]-locations[j,0])+(centre[1]-locations[j,1])*(centre[1]-locations[j,1])+(centre[2]-locations[j,2])*(centre[2]-locations[j,2]))
        
        idx = np.argsort(r)
        sorted_r = r[idx]
        sorted_charge = prob[idx]
        
        cum_prob = np.cumsum(sorted_charge)
        
        for j in range (0,n):
            if cum_prob[j]>0.99:break
        
        Polaron_size[i] = sorted_r[j]
        
        for j in range(0,n):
            if prob[j]/max_prob>0.01:Num[i]+=1
        
        centre = np.zeros(3)            #Reset centre coordinates and num
    
    print Num[0:10]
    
    
    
    polaron_evals=np.zeros((n,2))
    for i in range (0,n):
        polaron_evals[i,0]=Pvals[i]
        polaron_evals[i,1]=Polaron_size[i]
    
    return Polaron_size,Num
    
    
    
    savedata("size")

#Find alpha and disorder that will localise polaron on 1 molecule (99%)
def localisationcriteria(Num):
    if Num[0]==1:print "Localised with alpha= ", ALPHA, "and disorder= ", dx
    else:exit(0)


#Print number of molecules polaron localised over for first 10 eigenvalues
def plotsize(Pvals,Polaron_size):
    fig=pl.figure()
    pl.bar(Pvals[0:10],Polaron_size[0:10],0.00001)
    #pl.title("Size of polaron vs eigenvalue")
    pl.xlabel("Eigenvalues")
    pl.ylabel("Effective size of polaron")
    #pl.xlim(-3.88,-3.47)
    
    #pl.show()
    #fig.savefig("Sizeofpolaron.pdf")
    
    #archivefigure("Size")
    fig.savefig("%s_%s_%s_Size.pdf"%(dx,ALPHA,state))



#Plot Occupation, cumulative probailities and DOS
def plot3fig(Pvals,Pvecs):
    fig=pl.figure()
    
    pl.title("DoS by TightBinding")
    pl.subplot(311) #3 subplots stacked on top of one another
    
    #Plot Eigenvalues with filled-in Eigenvectors^2 / electron probabilities
    pl.subplot(311)
    
    psi=Pvecs[:,state]*Pvecs[:,state]
    pl.fill_between(range(n),0,psi, facecolor='blue')
    #pl.plot(range(n),pvecs[:,j],color=colours[j%8])
    pl.ylabel("Occupation")
    #pl.ylim((3.8,5.0))
    pl.yticks(fontsize=9)
    #pl.xticks(visible=False)
    
    #Plot cumulative eigenvectors / probability density
    pl.subplot(312)
    
    for j in [0,1,2,n/4, n/2]: #range(0,5): #Number of eigenvalues plotted (electron wavefns)
        #psi=Pvecs[:,j]*Pvecs[:,j]    # expectation value
        #pl.fill_between(range(n),0,sorted(psi,reverse=True), facecolor=colours[j%8]) #can't see this anymore on large plots...
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
    
    pl.hist(Pvals,bins=np.linspace(min(Pvals),max(Pvals),100),histtype='stepfilled',color='r')
    pl.hist(evals,bins= np.linspace(min(pvals),max(pvals),100),histtype='stepfilled',color='b', alpha=0.5)
    pl.ylabel("DoS")
    pl.yticks(fontsize=9)
    #pl.xlim(-6.5,-5)
    
    #pl.show()
    
    fig.savefig("%s_%s_%s_3fig.pdf"%(dx,ALPHA,state))

#    archivefigure("3fig")

#Plot DOS
def plotDOS(Evals,Pvals):
    fig=pl.figure()
    pl.hist(Pvals,bins=np.linspace(min(pvals),max(pvals),100),histtype='stepfilled',color='r')
    pl.hist(Evals,bins= np.linspace(min(pvals),max(pvals),100),histtype='stepfilled',color='b', alpha=0.5)
    pl.ylabel("DoS")
    pl.yticks(fontsize=9)
    #pl.xlim(-6.5,-5)
    
    #pl.show()
    fig.savefig("%s_%s_%s_DOS.pdf"%(dx,ALPHA,state))

#    archivefigure("DOS")

#Plot occupation of given state
def plotocc(Pvecs):
    fig=pl.figure()
    psi=Pvecs[:,state]*Pvecs[:,state]
    pl.fill_between(range(n),0,psi, facecolor='b')
    #pl.plot(range(n),pvecs[:,j],color=colours[j%8])
    pl.ylabel("Occupation")
    #pl.ylim((3.8,5.0))
    pl.yticks(fontsize=9)
    #pl.xticks(visible=False)
    
    #pl.show()
    fig.savefig("%s_%s_%s_Occ.pdf"%(dx,ALPHA,state))

def polaronvisualise():
    fp=open('eigenvector_balls_pymol.py','w')
    fp.write("from pymol.cgo import *    # get constants \nfrom pymol import cmd \n")
    
    psi=Pvecs[:,0]*Pvecs[:,0] # scale everything relative to max density on first eigenvector
    
    for ei,colour in zip( [0,1,3,5] , [(0,0,1),(0,1,0),(1,1,0),(1,0,0)]):
        
        psi=Pvecs[:,ei]*Pvecs[:,ei]
        maxpsi=max(psi)
        
        fp.write("obj = [\n")
        
        psisum=0.0
        for i in reversed(np.argsort(psi)): #magic list of sorted array indices
            weight=float(psi[i])/maxpsi #on interval [0,1]
            fp.write("ALPHA, %f,\n" %(weight))
            weight=1
            fp.write("COLOR, %f, %f, %f,\n" %(colour[0]*weight , colour[1]*weight, colour[2]*weight))
            fp.write("SPHERE, %f, %f, %f, 5.0,\n" %(locations[i][0],locations[i][1],locations[i][2]))
        
        fp.write(" END \n]\ncmd.load_cgo(obj,'EV_%d') \n" %(ei))


print("ITIAM.py - DoS by TB with Polaron self-interaction")
print("call as: #ITIAM.py (sites) (file.xyz) (CellA) (CellB) (CellC)")



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

ALPHA = float(sys.argv[6])
dx = float(sys.argv[7])
state = float(sys.argv[8])

timesteps=1000      #no. of timesteps for animation
dt=0.1519           #timestep for animation


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
BETA=0.6
H=J0*np.exp(-BETA*H) # calculate transfer integrals with isotropic exponential form
#H=np.random.normal(loc=H,scale=0.00005) # disorder in couplings

Hp=H

evals=np.zeros(n)
pvals=np.zeros(n)
evecs=np.zeros((n,n))
pvecs=np.zeros((n,n))
polaron_size=np.zeros(n)
fig=pl.figure



#---------------------------------------------------------------------------------#





#H=nearestneighbours(H)

print "Generated Hamiltonian... "

#H=offdiagonal(H)

H,Hp=filldiagonal(H)

print "Hamiltonian fully setup, time to solve!"

H,Hp=SCpolarongenerator(H,Hp)

evals,evecs,pvals,pvecs = solveHandHp(H,Hp)

print "Hamiltonian solved"

#plotoverlaps(evecs,pvecs)

#timeevolution(evecs,evals,pvecs)

Animate()

#polaron_size,num=sizeofpolaron(pvals,pvecs)

#localisationcriteria(num)

#plotsize(pvals,polaron_size)

#plot3fig(pvals,pvecs)

#plotDOS(evals,pvals)

#plotocc(pvecs)

#polaronvisualise()

print "Saving figures...(one moment please)"











