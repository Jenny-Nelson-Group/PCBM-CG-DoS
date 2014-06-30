# EEK - to hard for me to get my head around! ~Jarv
awk 'BEGIN{for (i=0;i<10;i++) for (j=0;j<10;j++) for(k=0;k<10;k++) print 10*i,5*(i%2)+10*j,5*(i%2)+10*k}' > 10x10x10_fcc.xyz
