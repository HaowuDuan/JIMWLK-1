for y in 0 0.1 0.2 0.3 0.4 0.5 1; do
# average <k2>
f=`cat efs/R_*.dat | awk -v y=$y '{if($1==y) {f=f+$3/$2;n=n+1}} END{print f/n}'`
# average X k^2
s=`cat efs/R_*.dat | awk -v y=$y '{if($1==y) {f=f+$3;n=n+1}} END{print f/n}'`
echo $f $s
cat efs/R_*.dat  | awk -v y=$y -v f=$f -v s=$s '{if($1==y) {print $3/s, $3/$2/f}}' > efs/cloud_R_$y.dat
done 
