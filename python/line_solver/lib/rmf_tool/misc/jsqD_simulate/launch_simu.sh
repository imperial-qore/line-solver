g++ -W -Wall -O3 JSQd_simulate.cc -o simulate_JSQ

number_of_simulation() {
    N=$1
    C=`expr 1000 + $N \* $N / 100`
    echo $C
}

launch_for_rho_d (){
    rho=$1
    d=$2
    finished=()
    for N in 5 10 20 30 50 100; do
	if ! [ ${finished[$N]} ]; then 
	    fileName=results/exp_rho${rho}_d${d}_N${N}
	    if [ -e $fileName ]; then A=`wc -l $fileName`; else A=0; fi
	    simulations_performed=`echo $A | sed 's/ .*//'`
	    simulations_todo=`number_of_simulation $N`
	    if [ $simulations_performed -le $simulations_todo ] ;
	    then
		# We launch 500 additional simulations 
		echo "Simulations for N=$N d=$d rho=0.${rho} not done    ($simulations_performed < $simulations_todo)"
		{ time  ./simulate_JSQ r0.${rho} d${d} N${N} e5 >> ${fileName};} 2>&1 | grep real
	    else
		finished[$N]=1
		echo "N=$N d=$d rho=0.$rho done     ($simulations_performed < $simulations_todo)";
	    fi
	fi
    done
}


for i in `seq 1 100`; do
    for rho in 90 95 70; do
	for d in 4 3 2; do
	    launch_for_rho_d $rho $d
	done
    done
done
