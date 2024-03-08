#!/bin/bash
text_file="/media/linlab/Data_Lv/XXEYY/BIDS/code/DTIPrep/dwisubj.txt"
rootPath=/media/linlab/Data_Lv/XXEYY/BIDS/derivatives/dwiprep
cd $rootPath
for item in $(cat $text_file); do
	
	ID=$item
	echo -e "\e[33mcurrent subject $ID\e[0m"
	echo Time point of  Start Processing: $(date +%y-%m-%d" "%H:%M:%S)
	cd $ID
	mkdir Connectome
	dwiextract "$ID"_align.mif -bzero - | mrmath -axis 3 - mean Connectome/"$ID"_b0.nii
	
	
	antsRegistrationSyN.sh -d 3 -o ./Connectome/ants -f ./Connectome/"$ID"_b0.nii  -m /media/linlab/Data_Lv/XXEYY/BIDS/code/Reference/MNI152_T1_1mm_brain.nii.gz -n 30
	
	antsApplyTransforms -d 3        \
		-i /media/linlab/Data_Lv/XXEYY/BIDS/code/Reference/relabel_AAL.nii \
		-r ./Connectome/"$ID"_b0.nii\
		-n GenericLabel[Linear] \
		-t ./Connectome/ants1Warp.nii.gz\
		-t ./Connectome/ants0GenericAffine.mat \
		-o ./Connectome/native_AAL.nii.gz
		
	###未加权矩阵
	tck2connectome -symmetric -zero_diagonal -scale_invnodevol ./track/tracks_10M.tck ./Connectome/native_AAL.nii.gz ./Connectome/"$ID"_aal.csv  -nthreads 30 -force 
	
	
	###加权矩阵
	tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in ./track/sift_1M.txt ./track/tracks_10M.tck ./Connectome/native_AAL.nii.gz ./Connectome/"$ID"_aal_weight.csv  -nthreads 30 -force
	cd ..
	
done
	 
