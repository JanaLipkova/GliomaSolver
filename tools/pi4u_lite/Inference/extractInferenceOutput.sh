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


echo "------------------------------------------------------"
echo "          EXTRACTING THE INFERENCE RESULTS            "
echo "------------------------------------------------------"

# Get number of last itteration, check convergance (p=1) and extract resutls
Ngen=$( awk 'NR==2{print $0} ' runinfo.txt )
Nparam=$(head -1 $(ls curgen_db_0* | tail -1) | awk '{print NF-1}')
p=$(awk -v rec=$Ngen 'NR==(rec+1005){res=$0}END{printf "%d", res} ' runinfo.txt)


if [ $p -eq 1 ]; then
	
	echo " "
	echo ">>> Extracting samples from the posterior PDF <<<"
	echo "-------------------------------------------------"
	echo "TMCMC converged after ${Ngen} generations."
	
	# Convert posterior from log(), T=sqrt(c/rho)
        logPosterior=curgen_db_$(printf "%03d" ${Ngen}).txt
        filename=posterior.txt
	awk -v N=$Nparam '{for(i=1;i<=N;i++) 
                          {j=i; $j=exp($j) } 
                          $3 = sqrt($3/$2);
                          print $0 }' ${logPosterior} > ${filename}      
 	
	echo " "
	echo "---------------------------------------"
	echo ">>> Ploting Posterior PDF manifold  <<<"
	echo "---------------------------------------"
        InputFile=Input.txt
        SolverPath=$(cat ${InputFile} | awk -F '=' '/^SolverPath/ {print $2}')
        MatlabTools="${SolverPath}/tools/UQ_Processing/Matlab/source/"
        MyBase=$(pwd)
        InputPosterior="${MyBase}/posterior.txt"

        cd "${MatlabTools}"
        matlab -nodisplay -nosplash -nodesktop -r "PosteriorPDF('${InputPosterior}'); exit"
        cd "$MyBase"

        echo " "
        echo "---------------------------------------"
        echo ">>> Computing MAP, mean and std      <<<"
        echo "---------------------------------------"
	declare -a mean=()
	declare -a std=()
	declare -a map=()

	for i in $(seq 1 ${Nparam}); do
	    mean[${i}]=$(mean_of_col $filename ${i})
	    std[${i}]=$(std_of_col $filename ${i} ${mean[$i]})
	done
	
	row=$(last_column_max_row $filename )
	for i in $(seq 1 ${Nparam}); do
		map[${i}]=$(awk -v r=$row -v j=$i 'NR==r{print $j} ' $filename)
	done
	
        cat > statistics.txt << EOF
MAP  | ${map[@]}
mean | ${mean[@]}
std  | ${std[@]} 
EOF

	echo " "	
        echo "---------------------------------------------"
        echo ">>> Computing the MAP tumor cell density  <<<"
        echo "---------------------------------------------"
	MAPfolder=MAP
	mkdir -p ${MAPfolder}

        brainHRPath=$( cat ${InputFile} | awk -F '=' '/^DataPath/ {print $2}')
        cp "${brainHRPath}/../brain256" ${MAPfolder}
	cd ${MAPfolder}

	
	cat > InputParameters.txt << EOF
${map[1]}
${map[2]}
${map[3]}
EOF

	cat > TumorIC.txt << EOF
${map[4]}
${map[5]}
${map[6]}
EOF

       head -18 ../TMPTMP/runAll.sh > run.sh
	sed -i 's/program=brain/program=brain256/g' run.sh
        sed -i 's/vtk=0/vtk=1/g' run.sh
        sed -i 's/verbose=0/verbose=1/g' run.sh
        chmod +x run.sh		
	./run.sh
	mv HGG_data.dat MAP.dat	
	mkdir Paraview 
	mv *vtu Paraview/
	
	echo " "
        echo "---------------------------------------------"
        echo ">>> Converting MAP.dat to nii + visualise <<<"
        echo "---------------------------------------------"
	SolverPath=$(cat ../${InputFile} | awk -F '=' '/^SolverPath/ {print $2}')
        MatlabTools="${SolverPath}tools/DataProcessing/source"

	InputNiiPath=$(  cat ../$InputFile | awk -F '=' '/^DataPath/ {print $2}')
	InputNiiFileName=$( ls "${InputNiiPath}" | awk '/\.nii$|\.nii.gz$/ {print $0}' | head -1)
	InputNiiData="${InputNiiPath}/${InputNiiFileName}"        	
	MyBase=$(pwd)
	InputDatData="${MyBase}"/MAP.dat

	cd "${MatlabTools}"
	bRotate=1
	bResize=1
	matlab -nodisplay -nojvm -nosplash -nodesktop -r "dat2nii('${InputDatData}','${InputNiiData}',${bRotate},${bResize}); exit"
	cd "$MyBase"

        echo " "
        echo "---------------------------------------"
        echo ">>> Plotting MAP tumor cell density <<<"
        echo "---------------------------------------"
	VisFolder=Vis
	mkdir "${VisFolder}"
	cp "${InputNiiPath}"/* 	"${VisFolder}"
	cp MAP.nii "${VisFolder}"
	InputNiftyData="${MyBase}/${VisFolder}"
	
	cd "${MatlabTools}"
        matlab -nodisplay -nosplash -nodesktop -r "plotAllNiftyData('${InputNiftyData}'); exit"
        cd "$MyBase"  
	mv Vis/Overview.jpg .
	rm -r Vis      

	echo " "
	echo "---------------------------------------"
	echo ">>>       Writing the report        <<<"
	echo "---------------------------------------"
	cd ..
	echo "Here I am: $(pwd)"
	
        echo " "
        echo "---------------------------------------"
        echo ">>>   Creating the Results folders  <<<"
        echo "---------------------------------------"
	ResultsFolder=Results
	mkdir -p ${ResultsFolder} 
	cp ReportFile.md ${ResultsFolder} 
	cp MAP/MAP.nii ${ResultsFolder}
	
	DetailsFolder=Results/Details		
	mkdir -p ${DetailsFolder}
	cp MAP/Overview.jpg ${DetailsFolder}
	cp PosteriorPDF.jpg ${DetailsFolder}
	cp posterior.txt ${DetailsFolder}
	cp statistics.txt ${DetailsFolder}
	
	# Fill in the ReportFile.pdf
        InputDataDir=$(  cat ../$InputFile | awk -F '=' '/^DataPath/ {print $2}')
        Nsamples=$(  cat ${InputFile} | awk -F '=' '/^Nsamples/ {print $2}')

        cd "${ResultsFolder}"
        sed -i 's|INPUT_DATA_FOLDER|'"${InputDataDir}"'|g' ReportFile.md
        sed -i 's|NUMBER_OF_SAMPLES|'"$Nsamples"'|g' ReportFile.md

	for i in $(seq 1 ${Nparam}); do
	    sed -i 's|_MAP_p'"${i}"'_|'"${map[${i}]}"'|g' ReportFile.md
	    sed -i 's|_mean_p'"${i}"'_|'"${mean[${i}]}"'|g' ReportFile.md
	    sed -i 's|_std_p'"${i}"'_|'"${std[${i}]}"'|g' ReportFile.md
	done

	pandoc ReportFile.md -o ReportFile.pdf
	cd ../
	
	OutputDataDir=$(dirname "${InputDataDir}")
	OutputDataDir="${OutputDataDir}/${ResultsFolder}"	

	if [ ! -d "${OutputDataDir}" ]; then
  	   mv "${ResultsFolder}" "${OutputDataDir}" 
	           echo " The Results are stored in:"
		   echo " ${OutputDataDir}"
	else
           echo "Output folder Results already exist, your results are stored as:"
	   echo "${OutputDataDir}1"
	   mv "${ResultsFolder}" "${OutputDataDir}1"
	fi
	
	echo " "        
	echo "-----------------------------------"
	echo "              DONE                 "
	echo "-----------------------------------"	

else
	echo "============================================================"
        echo ">>>  TMCMC did not yet converged to posterior PDF!!!     <<<" 
        echo "============================================================"
	echo "This can happen because the inference process was terminated or something went wrong. Please just resubmit the inference (run_inference.sh or run_inference_lrz.sh). The inference will continue where it stopped."
	echo "------------------------------------------------------------"

fi

