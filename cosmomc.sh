#!/bin/bash

DIR=/sps/snls14/CosmoMC2016/cosmomc_galileon

mpi_env=$2

cd $DIR

if [ "$mpi_env" == "mpich2" ]
then
    MPIEXEC="/usr/local/mpich2/bin/mpiexec"

    export OMP_NUM_THREADS=$((NSLOTS/4))
    NMACHINES=4
    ${MPIEXEC} -rmk sge -iface eth2 -np $NMACHINES $DIR/cosmomc_mpich2 $1
else
    if [ "$mpi_env" == "openmpi" ] || [ "$mpi_env" == "openmpi_4" ] || [ "$mpi_env" == "openmpi_8" ] || [ "$mpi_env" == "openmpi_16" ]
    then 

	while read line
	do
	    if [ "${line:0:1}" = '#' ]
	    then 
		continue
	    else
		var1=`echo $line |awk -F " = " '{print $1}'`
		var2=`echo $line | awk -F " = " '{print $2}'`
		if [ "$var1" == "chain_num" ]
		then
		    NMACHINES=$var2
		fi
	    fi
	done<$1

	if [ $NMACHINES -eq 1 ]
	then
	    $DIR/cosmomc $1
#	    echo "$DIR/cosmomc $1"
	else
#	MPIEXEC="/usr/local/intel/compilers_and_libraries_2017.1.132/linux/mpi/intel64/bin/mpiexec"
#	${MPIEXEC} --mca btl ^udapl,openib --mca btl_tcp_if_include eth2 -np $NMACHINES $DIR/cosmomc $1
	    mpiexec -np $NMACHINES $DIR/cosmomc $1
#	    echo "mpiexec -np $NMACHINES $DIR/cosmomc $1"
	fi

    fi
fi

if [ -z "$mpi_env" ]
then
    $DIR/cosmomc $1
fi
