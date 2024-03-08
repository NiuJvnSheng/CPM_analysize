#!/bin/bash

## dim path
rootpath=/media/linlab/Data_Lv/XXEYY/BIDS/
Oppath=/media/linlab/Data_Lv/XXEYY/BIDS/derivatives/dwiprep


text_file="/media/linlab/Data_Lv/XXEYY/BIDS/code/subj.txt"

for item in $(cat $text_file); do
	
	ID=$item

## 添加for循环

dwipath=${rootpath}/${ID}
mkdir ${Oppath}/${ID}

##1	convert nifti to mif


	mrconvert ${rootpath}/${ID}/dwi/"$ID"_dir-AP_dwi.nii.gz ${Oppath}/${ID}/"$ID"_raw.mif -fslgrad ${rootpath}/${ID}/dwi/"$ID"_dir-AP_dwi.bvec ${rootpath}/${ID}/dwi/"$ID"_dir-AP_dwi.bval
	
	
	cd ${Oppath}/${ID}
##2	creat raw mask
	dwi2mask *.mif - | maskfilter - dilate "$ID"_raw_mask.mif -npass 3
	
##3	denoise 
	dwidenoise "$ID"_raw.mif "$ID"_raw_denoise.mif -noise noiselevel.mif -mask "$ID"_raw_mask.mif
	
##4	Gibbs Ring artifact(optional)
	mrdegibbs "$ID"_raw_denoise.mif "$ID"_raw_denoisegibbs.mif
	
##5	Reverse phrase Encoding
	mrconvert ${rootpath}/${ID}/dwi/"$ID"_dir-PA_dwi.nii.gz PA.mif -force
	mrconvert PA.mif -fslgrad ${rootpath}/${ID}/dwi/"$ID"_dir-PA_dwi.bvec ${rootpath}/${ID}/dwi/"$ID"_dir-PA_dwi.bval - | mrmath - mean mean_b0_PA.mif -axis 3 -force
	dwiextract "$ID"_raw_denoisegibbs.mif - -bzero | mrmath - mean mean_b0_AP.mif -axis 3
	mrconvert mean_b0_PA.mif -coord 2 0:69 mean_b0_PA_modified.mif
	mrcat mean_b0_AP.mif mean_b0_PA_modified.mif -axis 3 b0_pair.mif
	
	
##7	headmotion  and distortion
	dwifslpreproc "$ID"_raw_denoisegibbs.mif "$ID"_prep.mif -nocleanup -pe_dir AP -rpe_pair -se_epi b0_pair.mif -eddy_options " --slm=linear --data_is_shelled"  -nthreads 32
	
##8	bias filled correction 可以更换为ANTs
	dwibiascorrect ants "$ID"_prep.mif "$ID"_unbiased.mif -bias bias.mif
	dwi2mask "$ID"_unbiased.mif mask.mif
	
	
	
##9	alignment to T1w("$ID"_align.mif)
	dwiextract "$ID"_unbiased.mif -bzero - | mrmath -axis 3 - mean b0.nii 
	flirt.fsl -dof 6 -cost normmi -ref /media/linlab/Data_Lv/XXEYY/BIDS/derivatives/fmriprep/${ID}/anat/"$ID"_desc-preproc_T1w.nii.gz -in b0.nii -omat T_fsl.txt  
	transformconvert T_fsl.txt b0.nii /media/linlab/Data_Lv/XXEYY/BIDS/derivatives/fmriprep/${ID}/anat/"$ID"_desc-preproc_T1w.nii.gz flirt_import T_DWItoT1.txt
	mrtransform -linear T_DWItoT1.txt "$ID"_unbiased.mif "$ID"_align.mif
	
##10	seg T1
	5ttgen fsl /media/linlab/Data_Lv/XXEYY/BIDS/derivatives/fmriprep/${ID}/anat/"$ID"_desc-preproc_T1w.nii.gz 5tt_coreg.mif -nthreads 32
	5tt2gmwmi 5tt_coreg.mif 5tt_gmwmi.mif	
	
	cd ${rootpath}
done	
