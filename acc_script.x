for y in 0 0.1 0.2 0.3 0.4 0.5 1; do
# average <k2>
f=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$3/$2;n=n+1}} END{print f/n}'`
# average X k^2
s=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$3;n=n+1}} END{print f/n}'`
echo $f $s
cat efs/R_28*.dat  | awk -v y=$y -v f=$f -v s=$s '{if($1==y) {print $3/s, $3/$2/f}}' > efs/cloud_R_F_$y.dat
done

for y in 0 0.1 0.2 0.3 0.4 0.5 1; do
# average <k2>
f=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$5/$4;n=n+1}} END{print f/n}'`
# average X k^2
s=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$5;n=n+1}} END{print f/n}'`
echo $f $s
cat efs/R_28*.dat  | awk -v y=$y -v f=$f -v s=$s '{if($1==y) {print $5/s, $5/$4/f}}' > efs/cloud_R_F1_$y.dat
done

for y in 0 0.1 0.2 0.3 0.4 0.5 1; do
# average <k2>
f=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$7/$6;n=n+1}} END{print f/n}'`
# average X k^2
s=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$7;n=n+1}} END{print f/n}'`
echo $f $s
cat efs/R_28*.dat  | awk -v y=$y -v f=$f -v s=$s '{if($1==y) {print $7/s, $7/$6/f}}' > efs/cloud_R_F2_$y.dat
done



for y in 0 0.1 0.2 0.3 0.4 0.5 1; do
# average <X>
f=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$2;n=n+1}} END{print f/n}'`
# average X k^2
s=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$3;n=n+1}} END{print f/n}'`
echo $f $s
cat efs/R_28*.dat  | awk -v y=$y -v f=$f -v s=$s '{if($1==y) {print ($3-s), ($3-s)/($2*$8*$8)}}' > efs/cloud_R_F_$y.dat
done

for y in 0 0.1 0.2 0.3 0.4 0.5 1; do
# average <k2>
f=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$5/$4;n=n+1}} END{print f/n}'`
# average X k^2
s=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$5;n=n+1}} END{print f/n}'`
echo $f $s
cat efs/R_28*.dat  | awk -v y=$y -v f=$f -v s=$s '{if($1==y) {print ($5-s), ($5-s)/($4*$8*$8)}}' > efs/cloud_R_F1_$y.dat
done

for y in 0 0.1 0.2 0.3 0.4 0.5 1; do
# average <k2>
f=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$7/$6;n=n+1}} END{print f/n}'`
# average X k^2
s=`cat efs/R_28*.dat | awk -v y=$y '{if($1==y) {f=f+$7;n=n+1}} END{print f/n}'`
echo $f $s
cat efs/R_28*.dat  | awk -v y=$y -v f=$f -v s=$s '{if($1==y) {print ($7-s)*(25.6*25.6/4.0/3.141592654/2), ($7-s)/($6*$8*$8)}}' > efs/cloud_R_F2_$y.dat
done
