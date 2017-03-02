for y in 0 0.5 1; do
f=`cat efs/NEW_*.dat | awk -v y=$y '{if($1==y) {f=f+$2;n=n+1}} END{print f/n}'`
s=`cat efs/NEW_*.dat | awk -v y=$y '{if($1==y) {f=f+$3/$2/($4*$4);n=n+1}} END{print f/n}'`
echo $f $s
cat efs/NEW_*.dat  | awk -v y=$y -v f=$f -v s=$s '{if($1==y) {print $2/f, $3/$2/($4*$4)/s}}' > efs/cloud_$y.dat
done 
