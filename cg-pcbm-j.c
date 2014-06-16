// cg-pcbm-j.c
//  Generate Coarse grain Js - calculate isotropic Js from Coarse Grain PCBM or
//  C60 molecular dynamics simulations
// File originally called 'gen_edges.c'; Jarvist Moore Frost 2013

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// for size of data structure. Nb: check your ulimits for heap size if
// segfault!
#define MAXPOINTS 100000
#define TWO_PI 6.2831853071795864769252866



double uniformRandom()
{
    return ( (double)(rand())+1.)/((double)(RAND_MAX)+1.);
}


double gauss(double x, double dx)
{
    double noise;
    double u1=uniformRandom();
    double u2=uniformRandom();
    
    noise=sqrt(dx*-2.*log(u1))*cos(u2)*TWO_PI;
    
    x=x+noise;
    
    return x;
}


int main()
{
    double coords[MAXPOINTS][3];
    double d,s,s_dis,mind,r[3];
    int n,i,j;
    
    double lambda=0.6;     // DRAGONS!
    double prefactor=10; //check this!
    
    scanf("%d %lf",&n,&mind); // First line of STDIN is 'number of points' 'min distance to consider'
    
    if (n>MAXPOINTS)
        fprintf(stderr,"ERROR: Too many points! Will probably now crash...");
    
    for (i=0;i<n;i++) // read tuple coordinates in
        scanf("%lf %lf %lf",& coords[i][0], & coords[i][1], & coords[i][2]);
    
    for (i=0;i<n;i++) // N^2 search over particle pairs (!!!) - this is why we are using C
        for (j=0;j<i;j++) // Nb: limit to i to just do lower diagonal of n*n
            if (j!=i) // no infinities if self
            {
                r[0]=coords[i][0]-coords[j][0]; r[1]=coords[i][1]-coords[j][1]; r[2]=coords[i][2]-coords[j][2];
                d= sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
                
                s=-6.0;
                s_dis=gauss(s,0);
                
                if (d<=mind)
                    printf("%d\t%d\t%lf\t%lf\n",i,j,prefactor*exp(-d*lambda),s_dis); //increasing factor makes J tail off 'sharper'
            }
    
}
