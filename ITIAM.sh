for pdb
do
 cell=` head "${pdb}" | grep CRYST1 | awk '{print $2,$3,$4}' ` # extract periodic cell vectors
 grep "^ATOM" "${pdb}" | grep " C " | awk '{print $6,$7,$8}' > tmp.xyz #all the 'C'60 pseudo-atoms
 sites=` cat tmp.xyz | wc -l ` #cat'ing so just have single field of lines
 python ITIAM.py ${sites} tmp.xyz ${cell}
done
