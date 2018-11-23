#!/bin/bash

#  Helpers
#-----------------------------------------------------
function last_column_max_row() {
	filename=$1
	awk '!i++{max = $NF;j=NR}{if ($NF > max) {max = $NF; j=NR}}END{printf "%d", j}' $filename
}

function mean_of_col() {
	filename=$1
	col_num=$2
	awk -v col=$col_num 'BEGIN{sum=0}{sum = sum + $col}END{printf "%.4f", sum/NR}' $filename
}

function std_of_col() {
	filename=$1
	col_num=$2
	avg=$3
	if [ -z "$avg" ]; then
		avg=$(mean_of_col $filename $col_num)
	fi
	awk -v col=$col_num -v avg=$avg 'BEGIN{ssum=0}{ssum=ssum+(avg-$col)^2}END{printf "%.4f", sqrt(ssum/(NR-1))}' $filename
}
#-----------------------------------------------------


# 1) Get number of last itteration
#------------------------------------
Ngen=$( awk 'NR==2{print $0} ' runinfo.txt )
Nparam=$(head -1 $(ls curgen_db_0* | tail -1) | awk '{print NF-1}')

echo "Total num. of generations= ${Ngen}" 

# 2) Check if converged (p=1). if yes postprocess the results
#--------------------------------------------------------------
p=$(awk -v rec=$Ngen 'NR==(rec+1005){res=$0}END{printf "%d", res} ' runinfo.txt)
echo "p=${p}"

if [ $p -eq 1 ]; then
	echo " TMCMC converged. Postprocessing results:"
	echo "-----------------------------------------"

	#2.1) Convert samples from posterior from log to real space -> posterior.txt
	#-------------------------------------------------------------------------------
        logPosterior=curgen_db_$(printf "%03d" ${Ngen}).txt
        filename=posterior.txt
	awk -v N=$Nparam '{for(i=1;i<=N;i++) {j=i; $j=exp($j) } print $0 }' ${logPosterior} > ${filename}      
        
	echo "The last generation ${Ngen} where convered from log() in file ${logPosterior} to ${filename}"
        echo "Samples from the posterior are stored in ${filename}"

	# 2.2) Get mean and std of parameters
	#--------------------------------------
	declare -a mean=()
	declare -a std=()
	
	for i in $(seq 1 ${Nparam}); do
	    mean[${i}]=$(mean_of_col $filename ${i})
	    std[${i}]=$(std_of_col $filename ${i} ${mean[$i]})
	done
	
	# 2.3) Get the MAP values
 	#-------------------------
	row=$(last_column_max_row $filename )
	echo "number of row=${row}"
	declare -a map=()
	
	for i in $(seq 1 ${Nparam}); do
		map[${i}]=$(awk -v r=$row -v j=$i 'NR==r{print $j} ' $filename)
	done
	
  	# 2.4) Write statistics to file
	#-------------------------------
        cat > statistics.txt << EOF
MAP  | ${map[@]}
mean | ${mean[@]}
std  | ${std[@]} 
EOF
	
	echo "Inference resutls: (also stored in statistics.txt)"
	echo "---------------------------------------------------------------------------------------"
	echo "MAP  |" ${map[@]}
	echo "mean |" ${mean[@]}
	echo "std  |" ${std[@]}
	echo "---------------------------------------------------------------------------------------"

	# 2.5) Get the MAP tumor
	#-----------------------------
	# 2.5.1 Greate MAP folder,script, input file
	filename=MAP
	mkdir -p ${filename}
        cd $filename
	
	# 2.5.2 get executable + make the run script
	cp ../TMPTMP/brain .
	head -18 ../TMPTMP/runAll.sh > run.sh 	
	sed -i 's/vtk=0/vtk=1/g' run.sh	
	sed -i 's/verbose=0/verbose=1/g' run.sh	
	sed -i 's/adaptive=1/adaptive=0/g' run.sh	
	sed -i 's/dumpfreq=20/dumpfreq=100/g' run.sh	

	# 2.5.3 write simulations input files
	cat > InputParameters.txt << EOF
${map[1]}
${map[2]}
${map[3]}
EOF

	echo "--------------------------------------------------"
	echo "   !!!  Be careful with c and T parameter  !!!    "
	echo "--------------------------------------------------"

	cat > TumorIC.txt << EOF
${map[4]}
${map[5]}
${map[6]}
EOF

	# 2.5.4 Run simulation to get MAP results
        chmod +x run.sh		
	./run.sh

	# 2.6) Convert MAP to nifty


else
        echo "you are screwed"
fi

# 1) apply math to each column and print new output to file -> make file with posteror samples, instead of log()

# awk '{for(i=1;i<4;i++) {j=i; $j=$j*10 } print $0 }' matrix > matrix10

# 2) get statistic - no need for exp() if input is posterior.txt

# 3) write MAP, mean, std to file 
# !!! check if using c or T for the inference

# 4) create MAP file + scripts for getting MAP.dat

# 5) convert MAP.dat to nii

# 6) Visualisations: plot samples, plot MAP overview, input MRI overview -> use vis script

# 7) Collect output to reports.pdf

# 8) Restructure the output files and inform user where to find results
