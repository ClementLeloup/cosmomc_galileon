#!/bin/bash

# this is a script to launch job on the ccin2p3 cluster

DIR=/sps/snls14/CosmoMC2016/cosmomc_galileon

input_file= 
interact_bool=0
type=
project="P_snovae"
ct_loc=
h_rt_loc=
s_rss_loc=
fsize_loc=
name=
path=$DIR/outputs/
nb_core=
mpi_env=
q_name=
cosmomc_input=

#what to do in interactive mode
interactive()
{

read -p "CPU time requested ? " ct_i

# loop to check if valid time in the hh:mm:ss format, in case not empty
while :; do
    if [ ! -z "$ct_i" ]
    then
	if [[ $ct_i =~ ([0-9][0-9]):([0-9][0-9]):([0-9][0-9])$ ]]
	then
	    if [ ${BASH_REMATCH[3]} -lt 60 ] && [ ${BASH_REMATCH[2]} -lt 60 ]
	    then 
		break
	    fi
	fi

    echo >&2 "Error : The format is wrong, please enter a valid time in the hh:mm:ss format !"
    read ct_i

    else
	break
    fi

done


read -p "Elapse time requested ? " h_rt_i

# loop to check if valid time in the hh:mm:ss format
while :; do
    if [ ! -z "$h_rt_i" ]
    then
	if [[ $h_rt_i =~ ([0-9][0-9]):([0-9][0-9]):([0-9][0-9])$ ]]
	then
	    if [ ${BASH_REMATCH[3]} -lt 60 ] && [ ${BASH_REMATCH[2]} -lt 60 ]
	    then 
		break
	    fi
	fi

    echo >&2 "Error : The format is wrong, please enter a valid time in the hh:mm:ss format !"
    read h_rt_i

    else
	break
    fi

done


read -p "RAM used ? " s_rss_i

read -p "Allocatable memory need expected ? " fsize_i

read -p "Name of the job ? " name_i

read -p "Number of cores used ? " nb_core_i

read -p "Which MPI environment ? " mpi_env_i

# loop to check if valid mpi environment
while :; do
    if [ ! -z "$mpi_env_i" ]
    then
	if [ "$mpi_env_i" == "mpich2" ] || [ "$mpi_env_i" == "openmpi" ] || [ "$mpi_env_i" == "openmpi_4" ] || [ "$mpi_env_i" == "openmpi_8" ] || [ "$mpi_env_i" == "openmpi_16" ]
	then
	    break
	else
	    echo >&2 "Error : Invalid MPI environment, please chose a valid MPI environment (mpich2,openmpi,openmpi_4,openmpi_8,openmpi_16)"
	    read mpi_env_i
	fi
    else
	break
    fi

done

read -p "Name of the queue ? " q_name_i

read -p "Input file for cosmomc ? " cosmomc_input_i

# loop to check if valid input file for cosmomc
while :; do
    if [ ! -z "$cosmomc_input_i" ]
    then
	if [ -e "$cosmomc_input_i" ]
	then
	    break
	fi

    echo >&2 "Error : The file does not exist, please enter a valid file !"
    read cosmomc_input_i

    else
	break
    fi

done

}


#flags for execution (input_file & interactive mode)
while getopts "if:" OPTION
do 
    case $OPTION in
	f)
	    input_file=$OPTARG
	    ;;
	i)
	    echo "Interactive mode"
	    interact_bool=1
	    interactive
	    ;;
	?)
	    echo " Error : Parameter unknown" >&2
	    exit
	    ;;
    esac
done

#check if input_file given & exist
if [[ -z "$input_file" ]]
then
    echo "No input file" >&2
    exit
elif [[ ! -e "$input_file" ]]
then
    echo "Error : File doesn't exist" >&2
    exit
fi



#loop that read the Card file

while read line
do
#for i in $(cat $input_file)
#do
#    line=`echo $i`
    if [ "${line:0:1}" = '#' ]
    then 
	continue
    else
	var1=`echo $line |awk -F "=" '{print $1}'`
	var2=`echo $line | awk -F "=" '{print $2}'`
	if [ "$interact_bool" -eq 1 ]
	then
	    case $var1 in
		"Type") 
		    if [ -z "$type_i" ]
		    then
			type=$var2
		    else
			sed -i 's/Type='"$var2"'/Type='"$type_i"'/' $input_file
			type=$type_i
		    fi
		    ;;
		"Project") 
		    if [ -z "$project_i" ]
		    then
			project=$var2
		    else
			sed -i 's/Project='"$var2"'/Project='"$project_i"'/' $input_file
			project=$project_i
		    fi
		    ;;
		"CPU_time")
		    if [ -z "$ct_i" ]
		    then
			ct_loc=$var2
		    else
			sed -i 's/CPU_time='"$var2"'/CPU_time='"$ct_i"'/' $input_file
			ct_loc=$ct_i
		    fi
		    ;;
		"Elapse_time") 
		    if [ -z "$h_rt_i" ]
		    then
			h_rt_loc=$var2
		    else
			sed -i 's/Elapse_time='"$var2"'/Elapse_time='"$h_rt_i"'/' $input_file
			h_rt_loc=$h_rt_i
		    fi
		    ;;
		"Mem_size") 
		    if [ -z "$s_rss_i" ]
		    then
			s_rss_loc=$var2
		    else
			sed -i 's/Mem_size='"$var2"'/Mem_size='"$s_rss_i"'/' $input_file
			s_rss_loc=$s_rss_i
		    fi
		    ;;
		"Alloc_mem") 
		    if [ -z "$fsize_i" ]
		    then
			fsize_loc=$var2
		    else
			sed -i 's/Alloc_mem='"$var2"'/Alloc_mem='"$fsize_i"'/' $input_file
			fsize_loc=$fsize_i
		    fi
		    ;;
		"Job_name") 
		    if [ -z "$name_i" ]
		    then
			name=$var2
		    else
			sed -i 's/Job_name='"$var2"'/Job_name='"$name_i"'/' $input_file
			name=$name_i
		    fi
		    ;;
		"Output_path") 
		    if [ -z "$path_i" ]
		    then
			path=$var2
		    else
			sed -i -r 's#Output_path='"$var2"'#Output_path='"$path_i"'#' $input_file
			path=$path_i
		    fi
		    ;;
		"Num_core")
		    if [ -z "$nb_core_i" ]
		    then
			nb_core=$var2
		    else
			sed -i 's/Num_core='"$var2"'/Num_core='"$nb_core_i"'/' $input_file
			nb_core=$nb_core_i
		    fi
		    ;;
		"MPI_env")
		    if [ -z "$mpi_env_i" ]
		    then
			mpi_env=$var2
		    else
			sed -i 's/MPI_env='"$var2"'/MPI_env='"$mpi_env_i"'/' $input_file
			mpi_env=$mpi_env_i
		    fi
		    ;;
		"Queue_name")
		    if [ -z "$q_name_i" ]
		    then
			q_name=$var2
		    else
			sed -i 's/Queue_name='"$var2"'/Queue_name='"$q_name_i"'/' $input_file
			q_name=$q_name_i
		    fi
		    ;;
		"Cosmomc_input")
		    if [ -z "$cosmomc_input_i" ]
		    then
			cosmomc_input=$var2
		    else
			sed -i 's#Cosmomc_input='"$var2"'#Cosmomc_input='"$cosmomc_input_i"'#' $input_file
			cosmomc_input=$cosmomc_input_i
		    fi
		    ;;
		"")
		    continue
		    ;;
		*) 
		    echo "Warning : unknown variable $var1 !" >&2
		    ;;
	    esac

	else
	    case $var1 in
		"Type") 
		    type=$var2
		    ;;
		"Project") 
	            project=$var2
		    ;;
		"CPU_time") 
	            ct_loc=$var2
		    ;;
		"Elapse_time") 
	            h_rt_loc=$var2
		    ;;
		"Mem_size") 
	            s_rss_loc=$var2
		    ;;
		"Alloc_mem") 
	            fsize_loc=$var2
		    ;;
		"Job_name") 
	            name=$var2
		    ;;
		"Output_path") 
	            path=$var2
		    ;;
		"Num_core")
	            nb_core=$var2
		    ;;
		"MPI_env")
	            mpi_env=$var2
		    ;;
		"Queue_name")
	            q_name=$var2
		    ;;
		"Cosmomc_input")
	            cosmomc_input=$var2
		    ;;
		"")
		    continue
		    ;;
		*) 
		    echo "Warning : unknown variable $var1 !" >&2
		    ;;
	    esac
	
	fi
	
    fi
done<$input_file

#check the format of CPU_time given in input_file
if [[ $ct_loc =~ ([0-9][0-9]):([0-9][0-9]):([0-9][0-9])$ ]]
then
    if [ ${BASH_REMATCH[3]} -ge 60 ] || [ ${BASH_REMATCH[2]} -ge 60 ]
    then 
	echo "Error : not valid CPU_time !" >&2
	exit
    fi
else
    echo "Error : not valid CPU_time format !" >&2
    exit
fi


#check the format of Elapse_time given in input_file
if [[ $h_rt_loc =~ ([0-9][0-9]):([0-9][0-9]):([0-9][0-9])$ ]]
then
    if [ ${BASH_REMATCH[3]} -ge 60 ] || [ ${BASH_REMATCH[2]} -ge 60 ]
    then  
	echo "Error : not valid Elapse_time !" >&2
	exit
    fi
else
    echo "Error : not valid Elapse_time format !" >&2
    exit
fi


#check the format of Mem_size
s_rss_last=`echo ${s_rss_loc:(-1)}`
if [ "$s_rss_last" != "G" ] && [ "$s_rss_last" != "M" ] && [ "$s_rss_last" != "K" ] && [ "$s_rss_last" != "B" ]
then
    echo "Error : invalid Mem_size format !" >&2
    exit
fi

#check the format of Alloc_mem
fsize_last=`echo ${fsize_loc:(-1)}`
if [ "$fsize_last" != "G" ] && [ "$fsize_last" != "M" ] && [ "$fsize_last" != "K" ] && [ "$fsize_last" != "B" ]
then
    echo "Error : invalid Alloc_mem format !" >&2
    exit
fi

#check if input_file given & exist
if [[ ! -e "$cosmomc_input" ]]
then
    echo "Error : Input file for cosmomc doesn't exist" >&2
    exit
fi


echo '###########################################################################################################################################################'
echo "Type of job : $type"
echo "Name of project : $project"
echo "CPU time requested : $ct_loc"
echo "Elapse time requested : $h_rt_loc"
echo "RAM expected : $s_rss_loc"
echo "Allocatable memory needed : $fsize_loc"
echo "Name of the job : $name"
echo "Path to output files : $path"
echo "Number of cores required : $nb_core"
echo "MPI environment : $mpi_env"
echo "Name of the queue : $q_name"
echo "Input card for CosmoMC : $cosmomc_input"
echo '###########################################################################################################################################################'



read -p "Are you sure you want to launch the job ? (yes,no) " last_check

if [ "$last_check" = "yes" ]
then
    echo "Job launched"
    case $type in
	"Basic")
	    if [ -z "$project" ] || [ -z "$ct_loc" ] || [ -z "$h_rt_loc" ] || [ -z "$s_rss_loc" ] || [ -z "$fsize_loc" ] || [ -z "$name" ] || [ -z "$path" ] || [ -z "$cosmomc_input" ]
	    then
			echo "Error : Some parameters are missing !" >&2
	    else
		echo "qsub -b y -P P_$project -l sps=1 -l ct=$ct_loc -l h_rt=$h_rt_loc -l s_rss=$s_rss_loc -l fsize=$fsize_loc -N $name -o $PWD/$path -e $PWD/$path $DIR/cosmomc.sh $DIR/$cosmomc_input"
		qsub -b y -P P_$project -l sps=1 -l ct=$ct_loc -l h_rt=$h_rt_loc -l s_rss=$s_rss_loc -l fsize=$fsize_loc -N $name -o $DIR/$path -e $DIR/$path $DIR/cosmomc.sh $DIR/$cosmomc_input
	    fi
	    ;;
	"Multi-core")
	    if [ -z "$project" ] || [ -z "$ct_loc" ] || [ -z "$h_rt_loc" ] || [ -z "$s_rss_loc" ] || [ -z "$fsize_loc" ] || [ -z "$name" ] || [ -z "$path" ] || [ -z "$nb_core" ] || [ -z "$q_name" ] || [ -z "$cosmomc_input" ]
	    then
			echo "Error : Some parameters are missing !" >&2
	    else
		echo "qsub -cwd -b y -P P_$project -l sps=1 -l ct=$ct_loc -l h_rt=$h_rt_loc -l s_rss=$s_rss_loc -l fsize=$fsize_loc -N $name -o $DIR/$path -e $DIR/$path -pe multicores $nb_core -q $q_name $DIR/cosmomc.sh $DIR/$cosmomc_input"
		qsub -cwd -b y -P P_$project -l sps=1 -l ct=$ct_loc -l h_rt=$h_rt_loc -l s_rss=$s_rss_loc -l fsize=$fsize_loc -N $name -o $DIR/$path -e $DIR/$path -pe multicores $nb_core -q $q_name $DIR/cosmomc.sh $DIR/$cosmomc_input $mpi_env
	    fi
	    ;;
	"Parallel")
            #check the MPI environment
	    if [ "$mpi_env" != "mpich2" ] && [ "$mpi_env" != "openmpi" ] && [ "$mpi_env" != "openmpi_4" ] && [ "$mpi_env" != "openmpi_8" ] && [ "$mpi_env" != "openmpi_16" ]
	    then
		echo "Error : Invalid MPI environment" >&2
		exit
	    fi
	    if [ -z "$project" ] || [ -z "$ct_loc" ] || [ -z "$h_rt_loc" ] || [ -z "$s_rss_loc" ] || [ -z "$fsize_loc" ] || [ -z "$name" ] || [ -z "$path" ] || [ -z "$nb_core" ] || [ -z "$mpi_env" ] || [ -z "$q_name" ] || [ -z "$cosmomc_input" ]
	    then
			echo "Error : Some parameters are missing !" >&2
	    else
		echo "qsub -cwd -b y -P P_$project -l sps=1 -l ct=$ct_loc -l h_rt=$h_rt_loc -l s_rss=$s_rss_loc -l fsize=$fsize_loc -N $name -o $DIR/$path -e $DIR/$path -pe $mpi_env $nb_core -q $q_name $DIR/cosmomc_pa.sh $mpi_env $DIR/$cosmomc_input"
		qsub -cwd -b y -P P_$project -l sps=1 -l ct=$ct_loc -l h_rt=$h_rt_loc -l s_rss=$s_rss_loc -l fsize=$fsize_loc -N $name -o $DIR/$path -e $DIR/$path -pe $mpi_env $nb_core -q $q_name $DIR/cosmomc_pa.sh $mpi_env $DIR/$cosmomc_input
	    fi	
	    ;;
	*)

	    ;;
    esac	

else
    echo "Job launch aborted"
fi











#qsub -P $project -l sps=1 -l ct=XX:XX:XX -l h_rt=XX:XX:XX -l s_rss=XXXM -l fsize=29G -o $output_name -e $erroutput_name job
