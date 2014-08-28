#!/usr/bin/env python

# ITIAM.py
# Alone, I emplore ya
# I think I'm a mother

# Import our numeric library
import numpy as np
# Matplotlib
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection

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


#Self-consistent polaron generator to model self-trapping of polaron using lowest state electron density

def SCpolarongenerator(Ham,Ham_p,):
    SCFSTEPS = 20
    
    siteEs=[]
    polarons=[]
    overlaps=[]
    pvecs_polaron=np.zeros((n,n))
    pvals_polaron=np.zeros((n,n))
    pvecs_size=np.zeros((n,n))
    Hp_diagonal = np.diagonal(Ham_p)
    

    max_overlap_idx=state
    for i in range(SCFSTEPS): # Number of SCF steps
        evals,evecs=np.linalg.eigh(Ham_p)
        polaron=evecs[:,max_overlap_idx]*evecs[:,max_overlap_idx] #lowest energy state electron density
        np.fill_diagonal(Ham_p,Hp_diagonal-ALPHA*polaron)
        pvals,pvecs=np.linalg.eigh(Ham_p)
        for j in range(0,n):
            psi0=evecs[:,j]
            psi1=pvecs[:,max_overlap_idx]
            J=np.dot(psi0,np.inner(Ham_p,psi1))
            #print J
            overlaps.append(J)
            max_overlap_idx=np.argmax(np.absolute(overlaps))
                #print state, max_overlap_idx
        overlaps=[]
    
    if state==0:pvecs_size=pvecs
    
    pvals_polaron[:,state]=pvals
    pvecs_polaron[:,state]=pvecs[:,0]      #All lowest state polaron wavefuntions for calculation of J


    return Ham,Ham_p,pvals_polaron,pvecs_polaron,pvecs_size


#SC polaron self trapping effect using all states for localising polaron

def SCpolarongenerator_allstates(Ham,Ham_p,):
    SCFSTEPS = 20
    
    siteEs=[]
    polarons=[]
    overlaps=[]
    pvecs_polaron=np.zeros((n,n))
    pvals_polaron=np.zeros((n,n))
    pvecs_size=np.zeros((n,n))
    Hp_diagonal = np.diagonal(Ham_p)
    
    for state in range(0,n):
        max_overlap_idx=state
        for i in range(SCFSTEPS): # Number of SCF steps
            evals,evecs=np.linalg.eigh(Ham_p)
            polaron=evecs[:,max_overlap_idx]*evecs[:,max_overlap_idx] #lowest energy state electron density
            np.fill_diagonal(Ham_p,Hp_diagonal-ALPHA*polaron)
            pvals,pvecs=np.linalg.eigh(Ham_p)
            for j in range(0,n):
                psi0=evecs[:,j]
                psi1=pvecs[:,max_overlap_idx]
                J=np.dot(psi0,np.inner(Ham_p,psi1))
                #print J
                overlaps.append(J)
                max_overlap_idx=np.argmax(np.absolute(overlaps))
            #print state, max_overlap_idx
            overlaps=[]
    
        if state==0:pvecs_size=pvecs
    
        pvals_polaron[:,state]=pvals
        pvecs_polaron[:,state]=pvecs[:,0]      #All lowest state polaron wavefuntions for calculation of J
    
    
    return Ham,Ham_p,pvals_polaron,pvecs_polaron,pvecs_size



# solve final form of Hamiltonian (always computes here even if no SCF steps)
def solveHandHp(Ham):
    Evals,Evecs=np.linalg.eigh(Ham)
    return Evals,Evecs


def localisepolaron(Ham,Ham_p):

    Hp_diagonal = np.diagonal(Ham_p)
    site=np.zeros(n)

    for i in range(0,n):
        site[i]=1
        np.fill_diagonal(Ham_p,Hp_diagonal-ALPHA*site)
    return Ham,Pvecs_polaron


def CalculateJ(Ham,Pvecs_polaron):
    Jif = np.zeros((n,n))
    for i in range(0,n):
        for j in range(0,n):
            if i!=j:
                Jif[i,j]=np.dot(Pvecs_polaron[:,i],np.inner(Ham,Pvecs_polaron[:,j]))

    np.savetxt("J_mono.dat",Jif,delimiter=' ',newline='\n')

    return Jif


#def MarcusRate(Jif,Pvals_polaron):
#    Gamma=np.zeros((n,n))
#    hbar=6.582*10**-16
#    lambda=2*ALPHA
#    kb=8.617*10**-5
#    T=300

#    for i in range(0,n):
#        for j in range(0,n):
#            deltaE=Pvals_polaron[0,j]-Pvals_polaron[0,i]
#            Gamma[i,j]=(2*np.pi*Jif[i,j]*Jif[i,j]/hbar)*(4*np.pi*lambda*kb*T)**-0.5*np.exp(-((deltaE*lambda)**2)/(4*lambda*kb*T))

#    np.savetxt("Marcus_rate",Gamma,delimiter=' ',newline='\n')



def plotoverlaps(Evecs,Pvecs_polaron):
    polarons=[]
    overlaps=[]
    for state in range(n): #:[0,1,2,3]: #range(n): #[1,2,3,500]:
        psi0=Evecs[:,state] #.reshape(1,n)
        psi1=Pvecs_polaron[:,0].reshape(1,n)
        J=np.dot(psi0,np.inner(H,psi1))
        overlaps.append(J)
        polarons.append(state)
    
    fig=pl.figure()
    pl.title("Js by Polaron Orbital Overlap")
    pl.plot(polarons,overlaps)
#pl.show()

#    archivefigure("POO")


def calcocc(T,Evals,Evecs,Pvecs_polaron):
    
    psi0=Evecs[:,0]              #Original wavefunction
    psi1=Pvecs_polaron[:,0]      #Localised wavefunction to propagate
    centre=np.zeros(3)
    psi1_occ=np.zeros(n)
    r=np.zeros(n)
    
    
    for i in range(0,n):
        psi1_occ[i]=psi1[i]*np.conj(psi1[i])
    
    
    max_index=np.argmax(psi1_occ)
    centre[0]=locations[max_index,0]
    centre[1]=locations[max_index,1]
    centre[2]=locations[max_index,2]
    
    
    #Positions of all molecules from max occupied one
    
    for j in range(0,n):
        r[j]=np.sqrt((centre[0]-locations[j,0])*(centre[0]-locations[j,0])+(centre[1]-locations[j,1])*(centre[1]-locations[j,1])+(centre[2]-locations[j,2])*(centre[2]-locations[j,2]))
    
    
    idx=np.argsort(r)      #Sort distances
    sort_r=r[idx]
    
    hbar=6.582*10**-16
    z=1j
    
    dt=5*10**-15           #timestep for plots
    
    psi_t=np.zeros(n)
    psi_t_occ=np.zeros(n)
    
    
    for i in range(0,n):
        a=np.inner(Evecs[:,i],psi1)
        psi_t+=a*np.exp(-z*Evals[i]*dt*T/hbar)*Evecs[:,i]
    
    
    psi_t_norm=np.linalg.norm(psi_t)      #Normalise
    psi_t=psi_t/psi_t_norm
    
    
    for j in range(0,n):
        psi_t_occ[j]=psi_t[j]*np.conj(psi_t[j])

    psi_occ=psi_t_occ[idx]
    
    
    return sort_r, psi_occ



def plotpropagation():
    
    tsteps=10              #No. of plots

    Z=np.zeros((n,tsteps))
    T=np.linspace(0,tsteps,tsteps)
    
    for t in range(0,tsteps):

        X,Z[:,t]=calcocc(t,evals,evecs,pvecs_polaron)
    

    fig=pl.figure()           #Plot time evolved wavefunction for different dts

    ax = fig.add_subplot(111, projection='3d')
    
    verts = []
    for i in xrange(T.shape[0]):
        verts.append(zip(X[0:150], Z[:,i]))
    

    poly = PolyCollection(verts, facecolors=(1,1,1,1), edgecolors=(0,0,1,1))
    ax.add_collection3d(poly, zs=T, zdir='y')
    ax.set_xlim3d(np.min(X), 40)
    ax.set_ylim3d(np.min(T), np.max(T))
    ax.set_zlim3d(np.min(Z[:,0]), np.max(Z[:,0]))

    pl.show()
    

    fig.savefig("%s_%s_time.pdf"%(dx,ALPHA))




def Animate():

    fig=pl.figure()

    xlim=(0,60)
    ylim=(0,1)

    ax = pl.axes(xlim=xlim,ylim=ylim)
    psi_x_line, = ax.plot([],[],c='r')

    def init():
        psi_x_line.set_data([], [])
        return psi_x_line,

    def animate(i):
        sort_r,psi_occ=calcocc(i,evals,evecs,pvecs_polaron)
        x=sort_r
        y=psi_occ
        psi_x_line.set_data(x,y)

        return psi_x_line,

    anim=animation.FuncAnimation(fig,animate,init_func=init,frames=1000,interval=10,blit=False)

    pl.show()

    anim.save('wavefunction_propagation.mp4', writer='ffmpeg', fps=30, extra_args=['-vcodec', 'libx264'])


#Find effective size of polaron and no. of molecules polaron localised over
def sizeofpolaron(Pvals_polaron,Pvecs_size):
    centre = np.zeros(3)
    prob = np.zeros(n)
    r = np.zeros(n)
    Polaron_size=np.zeros(n)
    cum_prob=np.zeros(n)
    sorted_r=np.zeros(n)
    Num=np.zeros(n)
    
    for i in range(0,n):
        
        for j in range(0,n):
            prob[j]=Pvecs_size[j,i]*Pvecs_size[j,i]
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
        polaron_evals[i,0]=Pvals_polaron[i,0]
        polaron_evals[i,1]=Polaron_size[i]
    
    return Polaron_size,Num
    
    savedata("size")

#Find alpha and disorder that will localise polaron on 1 molecule (99%)
def localisationcriteria(Num):
    if Num[0]==1:print "Localised with alpha= ", ALPHA, "and disorder= ", dx
    else:exit(0)


#Print number of molecules polaron localised over for first 10 eigenvalues
def plotsize(Pvals_polaron,Polaron_size):
    fig=pl.figure()
    pl.bar(Pvals_polaron[0:10,0],Polaron_size[0:10],0.00001)
    #pl.title("Size of polaron vs eigenvalue")
    pl.xlabel("Eigenvalues")
    pl.ylabel("Effective size of polaron")
    pl.ylim(0,120)
    #pl.xlim(-3.88,-3.47)
    
    #pl.show()
    #fig.savefig("Sizeofpolaron.pdf")
    
    #archivefigure("Size")
    fig.savefig("%s_%s_Size.pdf"%(dx,ALPHA))



#Plot Occupation, cumulative probailities and DOS
def plot3fig(Pvals_polaron,Pvecs_polaron):
    fig=pl.figure()
    for state in range(0,n):
        psi=Pvecs_polaron[:,state]*Pvecs_polaron[:,state]
        pl.fill_between(range(n),0,psi,facecolor=colours[state%8],alpha=0.5)
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
    
    pl.hist(Pvals_polaron,bins=np.linspace(min(Pvals_polaron[:,0]),max(Pvals_polaron[:,0]),100),histtype='stepfilled',color='r')
    pl.hist(evals,bins= np.linspace(min(Pvals_polaron[:,0]),max(Pvals_polaron[:,0]),100),histtype='stepfilled',color='b', alpha=0.5)
    pl.ylabel("DoS")
    pl.yticks(fontsize=9)
    #pl.xlim(-6.5,-5)
    
    #pl.show()
    
    fig.savefig("%s_%s_%s_3fig.pdf"%(dx,ALPHA,state))

#    archivefigure("3fig")

#Plot DOS
def plotDOS(Evals,Pvals_polaron):
    fig=pl.figure()
    pl.hist(Pvals_polaron[:,0],bins=np.linspace(min(Pvals_polaron[:,0]),max(Pvals_polaron[:,0]),100),histtype='stepfilled',color='r')
    pl.hist(Evals,bins= np.linspace(min(Pvals_polaron[:,0]),max(Pvals_polaron[:,0]),100),histtype='stepfilled',color='b', alpha=0.5)
    pl.xlabel("Energy (eV)")
    pl.ylabel("DoS")
    pl.yticks(fontsize=9)
    #pl.xlim(-6.5,-5)
    
    #pl.show()
    fig.savefig("%s_%s_DOS.pdf"%(dx,ALPHA))

#    archivefigure("DOS")

#Plot occupation of lowest energy for localisation with each state
def plotocc(Pvecs_polaron):
    fig=pl.figure()
    for state in range(0,n):
        psi=Pvecs_polaron[:,state]*Pvecs_polaron[:,state]
        pl.fill_between(range(n),0,psi,facecolor=colours[state%8],alpha=0.5)
    #pl.plot(range(n),pvecs[:,j],color=colours[j%8])
    pl.ylabel("Occupation")
    #pl.ylim((3.8,5.0))
    pl.yticks(fontsize=9)
    #pl.xticks(visible=False)
    
    #pl.show()
    fig.savefig("%s_%s_AllOcc.pdf"%(dx,ALPHA))

def polaronvisualise(Pvecs_size):
    fp=open('eigenvector_balls_pymol.py','w')
    fp.write("from pymol.cgo import *    # get constants \nfrom pymol import cmd \n")
    
    psi=Pvecs_size[:,0]*Pvecs_size[:,0] # scale everything relative to max density on first eigenvector
    
    for ei,colour in zip( [0,1,3,5] , [(0,0,1),(0,1,0),(1,1,0),(1,0,0)]):
        
        psi=Pvecs_size[:,ei]*Pvecs_size[:,ei]
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


#Use this one for localising just with lowest state for calculating size and time evolution

H,Hp,pvals_polaron,pvecs_polaron,pvecs_size=SCpolarongenerator(H,Hp)

#Use this one for localising using all states for calculating J

#H,Hp,pvals_polaron,pvecs_polaron,pvecs_size=SCpolarongenerator_allstates(H,Hp)


evals,evecs = solveHandHp(H)

#print "Calculating J"

#J=CalculateJ(H,pvecs_polaron)

print "Hamiltonian solved"

#plotoverlaps(evecs,pvecs)

#timeevolution(evals,evecs,pvecs_polaron)

plotpropagation()

Animate()

#polaron_size,num=sizeofpolaron(pvals_polaron,pvecs_size)

#localisationcriteria(num)

#plotsize(pvals_polaron,polaron_size)

#plot3fig(pvals_polaron,pvecs_polaron)

plotDOS(evals,pvals_polaron)

plotocc(pvecs_polaron)

#polaronvisualise(pvecs_size)

print "Saving figures...(one moment please)"











