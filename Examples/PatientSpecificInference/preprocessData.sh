#!/bin/bash

# nifty2dat
#-----------
#InputFolder="$1"
InputFolder='/Volumes/FileStorage/InputData/'
echo ${InputFolder}
echo $InputFolder

echo "Converting nifty to dat"
#matlab -nodisplay -nodesktop -r "convert_folder_content_dat2nii($InputFolder); exit"
matlab -nodisplay -nodesktop -r "run ../../ProcessingData/source/convert_folder_content_dat2nii($InputFolder);exit"
