for pdb
do
 cell=` head "${pdb}" | grep CRYST1 | awk '{print $2,$3,$4}' `
 grep "^ATOM" "${pdb}" | grep " C " | awk '{print $6,$7,$8}' > tmp.xyz
 python ITIAM.py 1000 tmp.xyz ${cell}
done
