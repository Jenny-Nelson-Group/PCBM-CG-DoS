for data 
do

#OK, let's generate our edges (connection data)
# # with our dirty little C code program, and a bit of AWK manipulation
pdbcat -fields "${data}"  | grep " C " | awk 'BEGIN{print 1000, 20}{print $10,$11,$12}'  | ./cg-pcbm-j > ${data%.*}.edges

python DoS_by_TB.py ${data%.*}.edges

#pdbcat -fields bis_5ns.pdb  | grep " C " | awk '{print $10,$11,$12,"-"}' > bis.xyz

#OK, let's generate XYZ coordinates, with definition of the contacts / gate edges
#pdbcat -fields "${data}"  | grep " C " | awk '{if ($12 < 10.0 ) {print $10,$11,$12,"g"} else if ($12 > 90.0) {print $10,$11,$12,"c"} else {print $10,$11,$12,"-"}}' >  ${data%.*}.xyz
#tft tof.sim  ${data%.*}.xyz ${data%.*}.edges | tee  ${data%.*}.tofet
done

#grep "MOBILITY FROM TOTAL DISPLACEMENT AND TOTAL TIME" *.tofet
