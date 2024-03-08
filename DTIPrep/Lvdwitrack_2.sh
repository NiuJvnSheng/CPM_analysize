#!/bin/bash
text_file="/media/linlab/Data_Lv/XXEYY/BIDS/code/DTIPrep/dwisubj.txt"
rootPath=/media/linlab/Data_Lv/XXEYY/BIDS/derivatives/dwiprep

cd $rootPath
for item in $(cat $text_file); do
	
	ID=$item
	cd $ID
	mkdir track
	echo -e "\e[33mcurrent subject $ID\e[0m"
	echo Time point of  Start Processing: $(date +%y-%m-%d" "%H:%M:%S)

	##	constrained spherical deconvolution(CSD)  反卷积
		
	#	single shell multi tissue
	#	algorithm1	ms_5tt_wm.txt ms_5tt_gm.txt ms_5tt_csf.txt
	dwi2response dhollander "$ID"_align.mif ms_5tt_wm.txt ms_5tt_gm.txt ms_5tt_csf.txt -voxels ms_5tt_voxels.mif -force
		
	#	algorithm2	wm_response.txt
	dwi2response tournier "$ID"_align.mif wm_response.txt -force
	#	estimate FOD
	#	single-shell single-tissue CSD
	dwi2fod csd "$ID"_align.mif wm_response.txt fod.mif -force
	#	mult-shell single-tissue CSD
	dwi2fod msmt_csd "$ID"_align.mif ms_5tt_wm.txt wmfod.mif ms_5tt_gm.txt gmfod.mif ms_5tt_csf.txt csffod.mif -force
	#	whole brain tractography
	tckgen -algo iFOD2 -act 5tt_coreg.mif -backtrack -crop_at_gmwmi -cutoff 0.05 -angle 45 -minlength 25 -maxlength 250 -seed_gmwmi 5tt_gmwmi.mif -select 10000k -nthreads 30 fod.mif ./track/tracks_10M.tck -force
	##	约束
	tcksift2 -act 5tt_coreg.mif -out_mu sift_mu.txt -out_coeffs sift_coeffs.txt -nthreads 30 ./track/tracks_10M.tck fod.mif ./track/sift_1M.txt
    	cd ..
    	echo $PWD
done
