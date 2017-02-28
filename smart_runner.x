j=1
while [ $j -lt 100  ]; do	
	sub=`ps -ax | grep 'jimwlk.x' | wc | awk '{print $1}'`;
	while [ $sub -lt 2 ]; do 
		sub=`ps -ax | grep 'jimwlk.x' | wc | awk '{print $1}'`;
		echo $j| ./jimwlk.x > out_$j.dat &
		echo "Launching instance # " $j 
		let "j = j+1"  ;
		sleep 10;
	done;
	sleep 30;
done;
