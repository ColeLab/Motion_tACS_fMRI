#!/bin/bash
# Taku Ito
# MotionAdaptation TACS collaboration with Krekelberg lab


##List of AFNI help files: http://afni.nimh.nih.gov/afni/doc/program_help/index.html

##--Analysis parameters--
#Subject numbers: 5 6 7 8 9 10 11 12 14
listOfSubjects="144" 
#List of run numbers, in order
runsPerSubj="1 2 3 4"
#Directory all data is under
basedir=/projects/Collaborations/KrekelbergCollaboration/MotionAdaptation_tACS_FC/
#Location of AFNI installation
AFNI_loc=/usr/local/afni
#Location of Freesurfer installation
FS_loc=/usr/local/freesurfer/bin/freesurfer
#Location of Python installation (must be above version 2.6, under 3.0)
py_loc=/usr/local/anaconda/bin/python
#Set paths
PATH=${py_loc}:${AFNI_loc}:${FS_loc}:${PATH}
#The fMRI TR duration (in seconds)
TR=2s
#The number of TRs per run
numTRs=190
#FWHM smoothing parameter to use
FWHMSmoothing=6
#Name of your current analysis
ANALYSISNAME="tacs_motionadaptation"

# Subject parameters
# T1
seriesDirString="Pos"

# EPIs
seriesDirString1="sess1-0004 sess1-0005 sess1-0007 sess1-0008"


##Running functions on each subject separately
for subjNum in $listOfSubjects
do
	echo "---Subject ${subjNum}---"

	#--Single subject directory setup--
	subjDir=${basedir}/data/${subjNum}/
		if [ "`ls -d $subjDir`" == "" ]; then mkdir $subjDir; fi	#Making directory if it doesn't exist yet
	subjRawDataDIR=${basedir}/data/rawdata/${subjNum}/
	subjfMRIDIR=${subjDir}/fMRI/
		if [ "`ls -d $subjfMRIDIR`" == "" ]; then mkdir $subjfMRIDIR; fi	#Making directory if it doesn't exist yet
	subjMaskDIR=${subjDir}/masks/
		if [ "`ls -d $subjMaskDIR`" == "" ]; then mkdir $subjMaskDIR; fi	#Making directory if it doesn't exist yet
	subjAnalysisDIR=${subjfMRIDIR}/${ANALYSISNAME}Analysis
		if [ "`ls -d $subjAnalysisDIR`" == "" ]; then mkdir $subjAnalysisDIR; fi	#Making directory if it doesn't exist yet
	StimFileDir=${basedir}/data/${subjNum}/sdm/
		#if [ "`ls -d $StimFileDir`" == "" ]; then mkdir $StimFileDir; fi	#Making directory if it doesn't exist yet
	FreesurferDir=${basedir}/data/freesurfer/
		if [ "`ls -d $FreesurferDir`" == "" ]; then mkdir $FreesurferDir; fi	#Making directory if it doesn't exist yet
	scriptDIR=${basedir}/docs/scripts/
	atlasDIR=${basedir}/data/atlases/


	##Preparing MPRAGE file (anatomical image)
	execute=1
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}
		
		echo "Preparing MPRAGE file (anatomical image)"
		
		##MODIFY: Set raw directory structure
		#How to identify MPRAGE series directory
		subDir="DICOM"
		fileFindString="*.dcm"
		
		##Convert MPRAGE raw data to NIFTI format
		#Get MPRAGE directory
		#dirName=`ls -d ${subjRawDataDIR}/${seriesDirString}`
		#Sort DICOM files (to make sure they will be read in the order they were collected in) using Freesurfer
		#dicom-rename ${dirName}/${subDir}/${fileFindString} --o ${subjRawDataDIR}/SortedDICOMs/MPRAGE/MR
		#Convert DICOM files to NIFTI format using Freesurfer
#		mri_convert ${subjRawDataDIR}/dcm/*${seriesDirString}*-00001.dcm --in_type siemens --out_type nii mprage.nii.gz
		#Remove sorted DICOMs
		#rm -rf ${subjRawDataDIR}/SortedDICOMs/MPRAGE

		##Skull strip MPRAGE
		##Use Freesurfer's skullstripping (very slow, but more accurate)
		#recon-all -subject $subjNum -autorecon1 -sd ${FreesurferDir} -force -i mprage.nii.gz
		recon-all -subject $subjNum -all -sd ${FreesurferDir} -i mprage.nii.gz
		mri_convert --in_type mgz --out_type nii ${FreesurferDir}/${subjNum}/mri/brain.mgz mprage_skullstripped.nii
		3dcopy mprage_skullstripped.nii mprage_skullstripped.nii.gz
		rm mprage_skullstripped.nii

		#Running alignment of MPRAGE to Talairach template
		3dcopy mprage_skullstripped.nii.gz anat_mprage_skullstripped
		@auto_tlrc -base ${atlasDIR}/MNI152_1mm_uni+tlrc -input anat_mprage_skullstripped+orig -no_ss -init_xform AUTO_CENTER
		ln -s ${atlasDIR}/MNI152_1mm_uni+tlrc* .

		##Create mask
		3dcalc -overwrite -a anat_mprage_skullstripped+tlrc -expr 'ispositive(a)' -prefix ${subjMaskDIR}/wholebrain_mask+tlrc
		#Link anatomical image to mask directory for checking alignment
		ln -s anat_mprage_skullstripped+tlrc* ${subjMaskDIR}/
		
		popd
	fi



	##Getting raw fMRI data folder names
	#MODIFY: How to identify fMRI series directory in raw data [regular expression used to get series order correct]
	subDir="DICOM"

	##For-loop for functions used across multiple runs (prior to run-concatenation)
	runNumFrom0=0
	runList=" "
	concatString="1D:"
	TRCount=0
        runNum=1
	for epirun in $seriesDirString1
	do
	
		echo "--Run ${runNum}--"
		
		##Convert fMRI data to AFNI format
		execute=1
		if [ $execute -eq 1 ]; then
			pushd ${subjfMRIDIR}
			

			fileFindString="*.dcm"
			#Converting to NIFTI format using Freesurfer
			mri_convert ${subjRawDataDIR}/dcm/*${epirun}*00001.dcm --in_type siemens --out_type nii epi_r${runNum}.nii.gz
			rm epi_r${runNum}+orig*
			3dcopy epi_r${runNum}.nii.gz epi_r${runNum}
			rm epi_r${runNum}.nii.gz

                        rm  epi_short_r${runNum}+orig*

                        # Create task timing files
                        if [ $runNum -eq 1 ]; then stimstring="*_nostim1_*"; fi
                        if [ $runNum -eq 2 ]; then stimstring="*_nostim2_*"; fi
                        if [ $runNum -eq 3 ]; then stimstring="*_stim1_*"; fi
                        if [ $runNum -eq 4 ]; then stimstring="*_stim2_*"; fi

                        # First 10 lines of original stim file is just header information that we don't need
                        tail -n +10 ${StimFileDir}/${stimstring} > ${StimFileDir}/TaskTiming${runNum}.txt

                        # Get number of TRs for each epi
                        tmp=`ls ${subjRawDataDIR}/dcm/*${epirun}* | wc -l`
                        # First 9 TRs are garbage for all EPIs
                        numTRs=$(expr $tmp - 9)
                        echo "Number of TRs for this EPI: ${numTRs}"
		        numTRsPRT=`cat ${StimFileDir}/TaskTiming${runNum}.txt | wc -l`

                        if [ "$numTRs" -le "$numTRsPRT" ]; then
                            numTRs=$numTRs
                        else
                            numTRs=$numTRsPRT
                        fi

                        # Obtain the necessary number of TRs per task
                        # delete if this is the first run (re-running this code block)
                        if [ ${runNum} -eq 1 ]; then rm ${StimFileDir}/allTaskTimings.txt; fi
                        echo "head -n +${numTRs} 'TaskTiming${runNum}.txt' >> ${StimFileDir}/allTaskTimings.txt"

                        head -n +"${numTRs}" "${StimFileDir}/TaskTiming${runNum}.txt" >> ${StimFileDir}/allTaskTimings.txt
                        
                        tmpTRcount=$(expr $numTRs - 1 + 9)
                        # Remove the first 9 TRs of each EPI
                        3dcalc -a epi_r${runNum}+orig[9..${tmpTRcount}] -expr a -prefix epi_short_r${runNum}
			popd

		fi
		nextInputFilename="epi_short"

		##Slice time correction
		execute=1
		if [ $execute -eq 1 ]; then
			pushd ${subjfMRIDIR}

			echo "-SliceTimeCorrection-"
			#Select the slice acquisition order. For Siemens TRIO (standard EPI sequence): alt+z when you have an odd number of slices, alt+z2 when you have an even number of slices
			tpattern="alt+z"
			3dTshift -overwrite -Fourier -TR ${TR} -tpattern ${tpattern} -prefix stc_${nextInputFilename}_r${runNum} ${nextInputFilename}_r${runNum}+orig
			#Remove intermediate analysis file to save disk space
			rm -v ${nextInputFilename}_r${runNum}+????.????.gz

			popd
		fi
		nextInputFilename="stc_"${nextInputFilename}

		#Keep track of run names
		runList="${runList} ${subjfMRIDIR}/${nextInputFilename}_r${runNum}+orig"

		#Keep track of TRs of run onsets for GLM analysis
		concatString="${concatString} ${TRCount}"
		TRCount=$(expr $TRCount + $numTRs)

		#Increment run count
		let runNumFrom0++
                let runNum++
	done

	##Concatenate runs
	execute=1
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}

		echo "-Concatenating runs-"
		echo "Run list: $runList"
		echo "Concatenation string (onset times of each run): $concatString"
		3dTcat -prefix ${nextInputFilename}_allruns ${runList}
		#Remove intermediate analysis file to save disk space
		rm -v ${nextInputFilename}_r*+????.????.gz

		popd
	fi
	nextInputFilename=${nextInputFilename}_allruns


	##Run align_epi_anat.py to align EPIs to MPRAGE, motion correct, and Talairach transform EPIs
	execute=1
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}

		echo "-Run align_epi_anat.py to align EPIs to MPRAGE, motion correct, and Talairach transform EPIs (output in 333 space)-"
		echo "Make sure Python is version 2.6 or greater"
		python -V
		#Correcting for motion, aligning fMRI data to MPRAGE, and aligning fMRI data to Talairach template [applying all transformation at once reduces reslicing artifacts]
		#[You could alternatively analyze all of the data, then Talairach transform the statistics (though this would make extraction of time series based on Talairached ROIs difficult)]
		#Visit for more info: http://afni.nimh.nih.gov/pub/dist/doc/program_help/align_epi_anat.py.html
		align_epi_anat.py -overwrite -anat anat_mprage_skullstripped+orig -epi ${nextInputFilename}+orig -epi_base 10 -epi2anat -anat_has_skull no -AddEdge -epi_strip 3dAutomask -ex_mode quiet -volreg on -deoblique on -tshift off -tlrc_apar anat_mprage_skullstripped+tlrc -master_tlrc ${atlasDIR}/MNI_EPI_333+tlrc

# 		if [ $runNum -eq 1 ]; then
# 			cat ${nextInputFilename}_vr_motion.1D > allruns_motion_params.1D
# 		else
# 			cat ${nextInputFilename}_vr_motion.1D >> allruns_motion_params.1D
# 		fi
		cp ${nextInputFilename}_vr_motion.1D allruns_motion_params.1D
		
		popd
	fi
	nextInputFilename=${nextInputFilename}"_tlrc_al"

	##Check motion parameters
	execute=1
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}
	
		#Plotting motion parameters
		1dplot -sep_scl -plabel ${subjNum}Motion -volreg allruns_motion_params.1D'[0..5]' &

		echo "Mean, standard deviation, and absolute deviation of subject's motion in mm (left to right), by x,y,z direction (top to bottom):" > MotionInfo.txt
		3dTstat -mean -stdev -absmax -prefix stdout: allruns_motion_params.1D'[0..2]'\' >> MotionInfo.txt
		cat MotionInfo.txt

		popd
	fi


	##Create a graymatter mask (based on Freesurfer output), and dilate it
	execute=1
	if [ $execute -eq 1 ]; then
		pushd ${subjMaskDIR}

 		echo "------------------------"
 		echo "Creating graymatter mask based on Freesurfer autosegmentation for subject " $subjNum
 		cd $subjMaskDIR
 		
 		cp ${FreesurferDir}/${subjNum}/mri/aseg.mgz ./${subjNum}_aseg.mgz
 		mri_convert -i ${subjNum}_aseg.mgz -ot nii ${subjNum}_aseg.nii
 		3dcalc -overwrite -a ${subjNum}_aseg.nii -expr 'a' -prefix ${subjNum}_aseg.nii.gz
 		rm ${subjNum}_aseg.nii
 		rm ${subjNum}_aseg.mgz
 		
 		#Using aseg (not aparc+aseg)
 		maskValSet="8 9 10 11 12 13 16 17 18 19 20 26 27 28 47 48 49 50 51 52 53 54 55 56 58 59 60 96 97 3 42"
 		#Add segments to mask
 		maskNum=1
 		for maskval in $maskValSet
 		do
 			if [ ${maskNum} = 1 ]; then
 				3dcalc -a ${subjNum}_aseg.nii.gz -expr "equals(a,${maskval})" -prefix ${subjNum}mask_temp.nii.gz -overwrite
 			else
 				3dcalc -a ${subjNum}_aseg.nii.gz -b ${subjNum}mask_temp.nii.gz -expr "equals(a,${maskval})+b" -prefix ${subjNum}mask_temp.nii.gz -overwrite
 			fi
 			let maskNum++
 		done
 		#Make mask binary
 		3dcalc -a ${subjNum}mask_temp.nii.gz -expr 'ispositive(a)' -prefix ${subjNum}_gmMask+orig -overwrite
		#Transform to TLRC space
		@auto_tlrc -apar ${subjfMRIDIR}/anat_mprage_skullstripped+tlrc -input ${subjNum}_gmMask+orig
 		#Resample to functional space
 		3dresample -overwrite -master ${subjfMRIDIR}/${nextInputFilename}+tlrc -inset ${subjNum}_gmMask+tlrc -prefix ${subjNum}_gmMask_func
		#Dilate mask by 1 functional voxel (just in case the resampled anatomical mask is off by a bit)
		3dLocalstat -overwrite -nbhd 'SPHERE(-1)' -stat 'max' -prefix ${subjNum}_gmMask_func_dil1vox+tlrc ${subjNum}_gmMask_func+tlrc
 		rm ${subjNum}mask_temp.nii.gz
 		
		popd
	fi


	##Create a white matter mask (based on Freesurfer output)
	execute=1
	if [ $execute -eq 1 ]; then
		pushd ${subjMaskDIR}

 		echo "------------------------"
 		echo "Create white matter mask, and erode it for subject $subjNum (MAKE SURE EROSION DOESN'T REMOVE ALL VENTRICLE VOXELS)"
 		cd $subjMaskDIR
 		
 		cp ${FreesurferDir}/${subjNum}/mri/aseg.mgz ./${subjNum}_aseg.mgz
 		mri_convert -i ${subjNum}_aseg.mgz -ot nii ${subjNum}_aseg.nii
 		3dcalc -overwrite -a ${subjNum}_aseg.nii -expr 'a' -prefix ${subjNum}_aseg.nii.gz
 		rm ${subjNum}_aseg.nii
 		rm ${subjNum}_aseg.mgz
		rm ${subjNum}mask_temp.nii.gz
 		
 		#Using aseg (not aparc+aseg)
 		#including all white matter
 		maskValSet="2 7 41 46"
 		#Add segments to mask
 		maskNum=1
 		for maskval in $maskValSet
 		do
 			if [ ${maskNum} = 1 ]; then
 				3dcalc -a ${subjNum}_aseg.nii.gz -expr "equals(a,${maskval})" -prefix ${subjNum}mask_temp.nii.gz -overwrite
 			else
 				3dcalc -a ${subjNum}_aseg.nii.gz -b ${subjNum}mask_temp.nii.gz -expr "equals(a,${maskval})+b" -prefix ${subjNum}mask_temp.nii.gz -overwrite
 			fi
 			let maskNum++
 		done
 		#Make mask binary
 		3dcalc -a ${subjNum}mask_temp.nii.gz -expr 'ispositive(a)' -prefix ${subjNum}_wmMask+orig -overwrite
		#Transform to TLRC space
		@auto_tlrc -apar ${subjfMRIDIR}/anat_mprage_skullstripped+tlrc -input ${subjNum}_wmMask+orig
 		#Resample to functional space
 		3dresample -overwrite -master ${subjfMRIDIR}/${nextInputFilename}+tlrc -inset ${subjNum}_wmMask+tlrc -prefix ${subjNum}_wmMask_func
		#Subtract graymatter mask from white matter mask (avoiding negative #s)
		3dcalc -a ${subjNum}_wmMask_func+tlrc -b ${subjNum}_gmMask_func_dil1vox+tlrc -expr 'step(a-b)' -prefix ${subjNum}_wmMask_func_eroded -overwrite
 		rm ${subjNum}mask_temp.nii.gz
 		
		popd
	fi


	##Create a ventricle mask (based on Freesurfer output)
	execute=1
	if [ $execute -eq 1 ]; then
		pushd ${subjMaskDIR}

 		echo "------------------------"
 		echo "Create ventricle mask, and erode it for subject $subjNum (MAKE SURE EROSION DOESN'T REMOVE ALL VENTRICLE VOXELS)"
 		cd $subjMaskDIR
 		
 		cp ${FreesurferDir}/${subjNum}/mri/aseg.mgz ./${subjNum}_aseg.mgz
 		mri_convert -i ${subjNum}_aseg.mgz -ot nii ${subjNum}_aseg.nii
 		3dcalc -overwrite -a ${subjNum}_aseg.nii -expr 'a' -prefix ${subjNum}_aseg.nii.gz
 		rm ${subjNum}_aseg.nii
 		rm ${subjNum}_aseg.mgz
		rm ${subjNum}mask_temp.nii.gz
 		
 		#Using aseg (not aparc+aseg)
 		#including ventricles
 		maskValSet="4 43 14 15"
 		#Add segments to mask
 		maskNum=1
 		for maskval in $maskValSet
 		do
 			if [ ${maskNum} = 1 ]; then
 				3dcalc -a ${subjNum}_aseg.nii.gz -expr "equals(a,${maskval})" -prefix ${subjNum}mask_temp.nii.gz -overwrite
 			else
 				3dcalc -a ${subjNum}_aseg.nii.gz -b ${subjNum}mask_temp.nii.gz -expr "equals(a,${maskval})+b" -prefix ${subjNum}mask_temp.nii.gz -overwrite
 			fi
 			let maskNum++
 		done
 		#Make mask binary
 		3dcalc -a ${subjNum}mask_temp.nii.gz -expr 'ispositive(a)' -prefix ${subjNum}_ventricles+orig -overwrite
		#Transform to TLRC space
		@auto_tlrc -apar ${subjfMRIDIR}/anat_mprage_skullstripped+tlrc -input ${subjNum}_ventricles+orig
 		#Resample to functional space
 		3dresample -overwrite -master ${subjfMRIDIR}/${nextInputFilename}+tlrc -inset ${subjNum}_ventricles+tlrc -prefix ${subjNum}_ventricles_func
		#Subtract graymatter mask from ventricle mask (avoiding negative #s)
		3dcalc -a ${subjNum}_ventricles_func+tlrc -b ${subjNum}_gmMask_func_dil1vox+tlrc -expr 'step(a-b)' -prefix ${subjNum}_ventricles_func_eroded -overwrite
 		rm ${subjNum}mask_temp.nii.gz
 		
		popd
	fi

	#Extract time series from white matter, ventricle masks, whole brain
	execute=1
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}

 		echo "--Extract time series from white matter, ventricle masks--"
 		3dmaskave -quiet -mask ${subjMaskDIR}/${subjNum}_wmMask_func_eroded+tlrc ${nextInputFilename}+tlrc > ${subjfMRIDIR}/${subjNum}_WM_timeseries_task.1D
 		3dmaskave -quiet -mask ${subjMaskDIR}/${subjNum}_ventricles_func_eroded+tlrc ${nextInputFilename}+tlrc > ${subjfMRIDIR}/${subjNum}_ventricles_timeseries_task.1D
 
 		echo "--Extract whole brain signal--"
 		cd ${subjMaskDIR}
		#Transform aseg to TLRC space
		3dcopy ${subjNum}_aseg.nii.gz aseg+orig
		@auto_tlrc -apar ${subjfMRIDIR}/anat_mprage_skullstripped+tlrc -input aseg+orig
 		3dcalc -overwrite -a aseg+tlrc -expr 'ispositive(a)' -prefix ${subjNum}_wholebrainmask+tlrc
		#Resample to functional space
 		3dresample -overwrite -master ${subjfMRIDIR}/${nextInputFilename}+tlrc -inset ${subjNum}_wholebrainmask+tlrc -prefix ${subjNum}_wholebrainmask_func+tlrc
		#Dilate mask by 1 functional voxel (just in case the resampled anatomical mask is off by a bit)
		3dLocalstat -overwrite -nbhd 'SPHERE(-1)' -stat 'max' -prefix ${subjNum}_wholebrainmask_func_dil1vox+tlrc ${subjNum}_wholebrainmask_func+tlrc
 		cd $subjfMRIDIR
 		3dmaskave -quiet -mask ${subjMaskDIR}/${subjNum}_wholebrainmask_func_dil1vox+tlrc ${nextInputFilename}+tlrc > ${subjfMRIDIR}/${subjNum}_wholebrainsignal_timeseries_task.1D

		popd
	fi


	##Run GLM
	execute=0
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}

		echo "Run GLM"
		3dDeconvolve \
			-input ${nextInputFilename}+tlrc \
			-mask ${subjMaskDIR}/${subjNum}_gmMask_func_dil1vox+tlrc \
			-concat "${concatString}" \
			-polort 3 \
			-num_stimts 360 \
			-stim_times 1 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV1_Trial1.1D.01.1D 'BLOCK(1,1)' -stim_label 1 Trial1 \
			-stim_times 2 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV2_Trial2.1D.01.1D 'BLOCK(1,1)' -stim_label 2 Trial2 \
			-stim_times 3 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV3_Trial3.1D.01.1D 'BLOCK(1,1)' -stim_label 3 Trial3 \
			-stim_times 4 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV4_Trial4.1D.01.1D 'BLOCK(1,1)' -stim_label 4 Trial4 \
			-stim_times 5 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV5_Trial5.1D.01.1D 'BLOCK(1,1)' -stim_label 5 Trial5 \
			-stim_times 6 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV6_Trial6.1D.01.1D 'BLOCK(1,1)' -stim_label 6 Trial6 \
			-stim_times 7 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV7_Trial7.1D.01.1D 'BLOCK(1,1)' -stim_label 7 Trial7 \
			-stim_times 8 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV8_Trial8.1D.01.1D 'BLOCK(1,1)' -stim_label 8 Trial8 \
			-stim_times 9 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV9_Trial9.1D.01.1D 'BLOCK(1,1)' -stim_label 9 Trial9 \
			-stim_times 10 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV10_Trial10.1D.01.1D 'BLOCK(1,1)' -stim_label 10 Trial10 \
			-stim_times 11 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV11_Trial11.1D.01.1D 'BLOCK(1,1)' -stim_label 11 Trial11 \
			-stim_times 12 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV12_Trial12.1D.01.1D 'BLOCK(1,1)' -stim_label 12 Trial12 \
			-stim_times 13 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV13_Trial13.1D.01.1D 'BLOCK(1,1)' -stim_label 13 Trial13 \
			-stim_times 14 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV14_Trial14.1D.01.1D 'BLOCK(1,1)' -stim_label 14 Trial14 \
			-stim_times 15 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV15_Trial15.1D.01.1D 'BLOCK(1,1)' -stim_label 15 Trial15 \
			-stim_times 16 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV16_Trial16.1D.01.1D 'BLOCK(1,1)' -stim_label 16 Trial16 \
			-stim_times 17 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV17_Trial17.1D.01.1D 'BLOCK(1,1)' -stim_label 17 Trial17 \
			-stim_times 18 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV18_Trial18.1D.01.1D 'BLOCK(1,1)' -stim_label 18 Trial18 \
			-stim_times 19 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV19_Trial19.1D.01.1D 'BLOCK(1,1)' -stim_label 19 Trial19 \
			-stim_times 20 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV20_Trial20.1D.01.1D 'BLOCK(1,1)' -stim_label 20 Trial20 \
			-stim_times 21 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV21_Trial21.1D.01.1D 'BLOCK(1,1)' -stim_label 21 Trial21 \
			-stim_times 22 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV22_Trial22.1D.01.1D 'BLOCK(1,1)' -stim_label 22 Trial22 \
			-stim_times 23 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV23_Trial23.1D.01.1D 'BLOCK(1,1)' -stim_label 23 Trial23 \
			-stim_times 24 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV24_Trial24.1D.01.1D 'BLOCK(1,1)' -stim_label 24 Trial24 \
			-stim_times 25 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV25_Trial25.1D.01.1D 'BLOCK(1,1)' -stim_label 25 Trial25 \
			-stim_times 26 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV26_Trial26.1D.01.1D 'BLOCK(1,1)' -stim_label 26 Trial26 \
			-stim_times 27 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV27_Trial27.1D.01.1D 'BLOCK(1,1)' -stim_label 27 Trial27 \
			-stim_times 28 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV28_Trial28.1D.01.1D 'BLOCK(1,1)' -stim_label 28 Trial28 \
			-stim_times 29 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV29_Trial29.1D.01.1D 'BLOCK(1,1)' -stim_label 29 Trial29 \
			-stim_times 30 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV30_Trial30.1D.01.1D 'BLOCK(1,1)' -stim_label 30 Trial30 \
			-stim_times 31 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV31_Trial31.1D.01.1D 'BLOCK(1,1)' -stim_label 31 Trial31 \
			-stim_times 32 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV32_Trial32.1D.01.1D 'BLOCK(1,1)' -stim_label 32 Trial32 \
			-stim_times 33 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV33_Trial33.1D.01.1D 'BLOCK(1,1)' -stim_label 33 Trial33 \
			-stim_times 34 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV34_Trial34.1D.01.1D 'BLOCK(1,1)' -stim_label 34 Trial34 \
			-stim_times 35 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV35_Trial35.1D.01.1D 'BLOCK(1,1)' -stim_label 35 Trial35 \
			-stim_times 36 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV36_Trial36.1D.01.1D 'BLOCK(1,1)' -stim_label 36 Trial36 \
			-stim_times 37 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV37_Trial37.1D.01.1D 'BLOCK(1,1)' -stim_label 37 Trial37 \
			-stim_times 38 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV38_Trial38.1D.01.1D 'BLOCK(1,1)' -stim_label 38 Trial38 \
			-stim_times 39 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV39_Trial39.1D.01.1D 'BLOCK(1,1)' -stim_label 39 Trial39 \
			-stim_times 40 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV40_Trial40.1D.01.1D 'BLOCK(1,1)' -stim_label 40 Trial40 \
			-stim_times 41 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV41_Trial41.1D.01.1D 'BLOCK(1,1)' -stim_label 41 Trial41 \
			-stim_times 42 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV42_Trial42.1D.01.1D 'BLOCK(1,1)' -stim_label 42 Trial42 \
			-stim_times 43 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV43_Trial43.1D.01.1D 'BLOCK(1,1)' -stim_label 43 Trial43 \
			-stim_times 44 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV44_Trial44.1D.01.1D 'BLOCK(1,1)' -stim_label 44 Trial44 \
			-stim_times 45 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV45_Trial45.1D.01.1D 'BLOCK(1,1)' -stim_label 45 Trial45 \
			-stim_times 46 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV46_Trial46.1D.01.1D 'BLOCK(1,1)' -stim_label 46 Trial46 \
			-stim_times 47 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV47_Trial47.1D.01.1D 'BLOCK(1,1)' -stim_label 47 Trial47 \
			-stim_times 48 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV48_Trial48.1D.01.1D 'BLOCK(1,1)' -stim_label 48 Trial48 \
			-stim_times 49 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV49_Trial49.1D.01.1D 'BLOCK(1,1)' -stim_label 49 Trial49 \
			-stim_times 50 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV50_Trial50.1D.01.1D 'BLOCK(1,1)' -stim_label 50 Trial50 \
			-stim_times 51 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV51_Trial51.1D.01.1D 'BLOCK(1,1)' -stim_label 51 Trial51 \
			-stim_times 52 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV52_Trial52.1D.01.1D 'BLOCK(1,1)' -stim_label 52 Trial52 \
			-stim_times 53 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV53_Trial53.1D.01.1D 'BLOCK(1,1)' -stim_label 53 Trial53 \
			-stim_times 54 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV54_Trial54.1D.01.1D 'BLOCK(1,1)' -stim_label 54 Trial54 \
			-stim_times 55 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV55_Trial55.1D.01.1D 'BLOCK(1,1)' -stim_label 55 Trial55 \
			-stim_times 56 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV56_Trial56.1D.01.1D 'BLOCK(1,1)' -stim_label 56 Trial56 \
			-stim_times 57 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV57_Trial57.1D.01.1D 'BLOCK(1,1)' -stim_label 57 Trial57 \
			-stim_times 58 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV58_Trial58.1D.01.1D 'BLOCK(1,1)' -stim_label 58 Trial58 \
			-stim_times 59 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV59_Trial59.1D.01.1D 'BLOCK(1,1)' -stim_label 59 Trial59 \
			-stim_times 60 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV60_Trial60.1D.01.1D 'BLOCK(1,1)' -stim_label 60 Trial60 \
			-stim_times 61 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV61_Trial61.1D.01.1D 'BLOCK(1,1)' -stim_label 61 Trial61 \
			-stim_times 62 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV62_Trial62.1D.01.1D 'BLOCK(1,1)' -stim_label 62 Trial62 \
			-stim_times 63 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV63_Trial63.1D.01.1D 'BLOCK(1,1)' -stim_label 63 Trial63 \
			-stim_times 64 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV64_Trial64.1D.01.1D 'BLOCK(1,1)' -stim_label 64 Trial64 \
			-stim_times 65 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV65_Trial65.1D.01.1D 'BLOCK(1,1)' -stim_label 65 Trial65 \
			-stim_times 66 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV66_Trial66.1D.01.1D 'BLOCK(1,1)' -stim_label 66 Trial66 \
			-stim_times 67 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV67_Trial67.1D.01.1D 'BLOCK(1,1)' -stim_label 67 Trial67 \
			-stim_times 68 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV68_Trial68.1D.01.1D 'BLOCK(1,1)' -stim_label 68 Trial68 \
			-stim_times 69 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV69_Trial69.1D.01.1D 'BLOCK(1,1)' -stim_label 69 Trial69 \
			-stim_times 70 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV70_Trial70.1D.01.1D 'BLOCK(1,1)' -stim_label 70 Trial70 \
			-stim_times 71 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV71_Trial71.1D.01.1D 'BLOCK(1,1)' -stim_label 71 Trial71 \
			-stim_times 72 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV72_Trial72.1D.01.1D 'BLOCK(1,1)' -stim_label 72 Trial72 \
			-stim_times 73 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV73_Trial73.1D.01.1D 'BLOCK(1,1)' -stim_label 73 Trial73 \
			-stim_times 74 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV74_Trial74.1D.01.1D 'BLOCK(1,1)' -stim_label 74 Trial74 \
			-stim_times 75 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV75_Trial75.1D.01.1D 'BLOCK(1,1)' -stim_label 75 Trial75 \
			-stim_times 76 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV76_Trial76.1D.01.1D 'BLOCK(1,1)' -stim_label 76 Trial76 \
			-stim_times 77 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV77_Trial77.1D.01.1D 'BLOCK(1,1)' -stim_label 77 Trial77 \
			-stim_times 78 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV78_Trial78.1D.01.1D 'BLOCK(1,1)' -stim_label 78 Trial78 \
			-stim_times 79 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV79_Trial79.1D.01.1D 'BLOCK(1,1)' -stim_label 79 Trial79 \
			-stim_times 80 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV80_Trial80.1D.01.1D 'BLOCK(1,1)' -stim_label 80 Trial80 \
			-stim_times 81 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV81_Trial81.1D.01.1D 'BLOCK(1,1)' -stim_label 81 Trial81 \
			-stim_times 82 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV82_Trial82.1D.01.1D 'BLOCK(1,1)' -stim_label 82 Trial82 \
			-stim_times 83 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV83_Trial83.1D.01.1D 'BLOCK(1,1)' -stim_label 83 Trial83 \
			-stim_times 84 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV84_Trial84.1D.01.1D 'BLOCK(1,1)' -stim_label 84 Trial84 \
			-stim_times 85 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV85_Trial85.1D.01.1D 'BLOCK(1,1)' -stim_label 85 Trial85 \
			-stim_times 86 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV86_Trial86.1D.01.1D 'BLOCK(1,1)' -stim_label 86 Trial86 \
			-stim_times 87 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV87_Trial87.1D.01.1D 'BLOCK(1,1)' -stim_label 87 Trial87 \
			-stim_times 88 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV88_Trial88.1D.01.1D 'BLOCK(1,1)' -stim_label 88 Trial88 \
			-stim_times 89 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV89_Trial89.1D.01.1D 'BLOCK(1,1)' -stim_label 89 Trial89 \
			-stim_times 90 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV90_Trial90.1D.01.1D 'BLOCK(1,1)' -stim_label 90 Trial90 \
			-stim_times 91 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV91_Trial91.1D.01.1D 'BLOCK(1,1)' -stim_label 91 Trial91 \
			-stim_times 92 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV92_Trial92.1D.01.1D 'BLOCK(1,1)' -stim_label 92 Trial92 \
			-stim_times 93 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV93_Trial93.1D.01.1D 'BLOCK(1,1)' -stim_label 93 Trial93 \
			-stim_times 94 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV94_Trial94.1D.01.1D 'BLOCK(1,1)' -stim_label 94 Trial94 \
			-stim_times 95 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV95_Trial95.1D.01.1D 'BLOCK(1,1)' -stim_label 95 Trial95 \
			-stim_times 96 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV96_Trial96.1D.01.1D 'BLOCK(1,1)' -stim_label 96 Trial96 \
			-stim_times 97 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV97_Trial97.1D.01.1D 'BLOCK(1,1)' -stim_label 97 Trial97 \
			-stim_times 98 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV98_Trial98.1D.01.1D 'BLOCK(1,1)' -stim_label 98 Trial98 \
			-stim_times 99 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV99_Trial99.1D.01.1D 'BLOCK(1,1)' -stim_label 99 Trial99 \
			-stim_times 100 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV100_Trial100.1D.01.1D 'BLOCK(1,1)' -stim_label 100 Trial100 \
			-stim_times 101 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV101_Trial101.1D.01.1D 'BLOCK(1,1)' -stim_label 101 Trial101 \
			-stim_times 102 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV102_Trial102.1D.01.1D 'BLOCK(1,1)' -stim_label 102 Trial102 \
			-stim_times 103 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV103_Trial103.1D.01.1D 'BLOCK(1,1)' -stim_label 103 Trial103 \
			-stim_times 104 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV104_Trial104.1D.01.1D 'BLOCK(1,1)' -stim_label 104 Trial104 \
			-stim_times 105 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV105_Trial105.1D.01.1D 'BLOCK(1,1)' -stim_label 105 Trial105 \
			-stim_times 106 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV106_Trial106.1D.01.1D 'BLOCK(1,1)' -stim_label 106 Trial106 \
			-stim_times 107 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV107_Trial107.1D.01.1D 'BLOCK(1,1)' -stim_label 107 Trial107 \
			-stim_times 108 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV108_Trial108.1D.01.1D 'BLOCK(1,1)' -stim_label 108 Trial108 \
			-stim_times 109 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV109_Trial109.1D.01.1D 'BLOCK(1,1)' -stim_label 109 Trial109 \
			-stim_times 110 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV110_Trial110.1D.01.1D 'BLOCK(1,1)' -stim_label 110 Trial110 \
			-stim_times 111 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV111_Trial111.1D.01.1D 'BLOCK(1,1)' -stim_label 111 Trial111 \
			-stim_times 112 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV112_Trial112.1D.01.1D 'BLOCK(1,1)' -stim_label 112 Trial112 \
			-stim_times 113 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV113_Trial113.1D.01.1D 'BLOCK(1,1)' -stim_label 113 Trial113 \
			-stim_times 114 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV114_Trial114.1D.01.1D 'BLOCK(1,1)' -stim_label 114 Trial114 \
			-stim_times 115 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV115_Trial115.1D.01.1D 'BLOCK(1,1)' -stim_label 115 Trial115 \
			-stim_times 116 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV116_Trial116.1D.01.1D 'BLOCK(1,1)' -stim_label 116 Trial116 \
			-stim_times 117 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV117_Trial117.1D.01.1D 'BLOCK(1,1)' -stim_label 117 Trial117 \
			-stim_times 118 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV118_Trial118.1D.01.1D 'BLOCK(1,1)' -stim_label 118 Trial118 \
			-stim_times 119 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV119_Trial119.1D.01.1D 'BLOCK(1,1)' -stim_label 119 Trial119 \
			-stim_times 120 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV120_Trial120.1D.01.1D 'BLOCK(1,1)' -stim_label 120 Trial120 \
			-stim_times 121 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV121_Trial121.1D.01.1D 'BLOCK(1,1)' -stim_label 121 Trial121 \
			-stim_times 122 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV122_Trial122.1D.01.1D 'BLOCK(1,1)' -stim_label 122 Trial122 \
			-stim_times 123 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV123_Trial123.1D.01.1D 'BLOCK(1,1)' -stim_label 123 Trial123 \
			-stim_times 124 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV124_Trial124.1D.01.1D 'BLOCK(1,1)' -stim_label 124 Trial124 \
			-stim_times 125 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV125_Trial125.1D.01.1D 'BLOCK(1,1)' -stim_label 125 Trial125 \
			-stim_times 126 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV126_Trial126.1D.01.1D 'BLOCK(1,1)' -stim_label 126 Trial126 \
			-stim_times 127 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV127_Trial127.1D.01.1D 'BLOCK(1,1)' -stim_label 127 Trial127 \
			-stim_times 128 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV128_Trial128.1D.01.1D 'BLOCK(1,1)' -stim_label 128 Trial128 \
			-stim_times 129 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV129_Trial129.1D.01.1D 'BLOCK(1,1)' -stim_label 129 Trial129 \
			-stim_times 130 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV130_Trial130.1D.01.1D 'BLOCK(1,1)' -stim_label 130 Trial130 \
			-stim_times 131 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV131_Trial131.1D.01.1D 'BLOCK(1,1)' -stim_label 131 Trial131 \
			-stim_times 132 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV132_Trial132.1D.01.1D 'BLOCK(1,1)' -stim_label 132 Trial132 \
			-stim_times 133 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV133_Trial133.1D.01.1D 'BLOCK(1,1)' -stim_label 133 Trial133 \
			-stim_times 134 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV134_Trial134.1D.01.1D 'BLOCK(1,1)' -stim_label 134 Trial134 \
			-stim_times 135 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV135_Trial135.1D.01.1D 'BLOCK(1,1)' -stim_label 135 Trial135 \
			-stim_times 136 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV136_Trial136.1D.01.1D 'BLOCK(1,1)' -stim_label 136 Trial136 \
			-stim_times 137 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV137_Trial137.1D.01.1D 'BLOCK(1,1)' -stim_label 137 Trial137 \
			-stim_times 138 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV138_Trial138.1D.01.1D 'BLOCK(1,1)' -stim_label 138 Trial138 \
			-stim_times 139 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV139_Trial139.1D.01.1D 'BLOCK(1,1)' -stim_label 139 Trial139 \
			-stim_times 140 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV140_Trial140.1D.01.1D 'BLOCK(1,1)' -stim_label 140 Trial140 \
			-stim_times 141 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV141_Trial141.1D.01.1D 'BLOCK(1,1)' -stim_label 141 Trial141 \
			-stim_times 142 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV142_Trial142.1D.01.1D 'BLOCK(1,1)' -stim_label 142 Trial142 \
			-stim_times 143 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV143_Trial143.1D.01.1D 'BLOCK(1,1)' -stim_label 143 Trial143 \
			-stim_times 144 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV144_Trial144.1D.01.1D 'BLOCK(1,1)' -stim_label 144 Trial144 \
			-stim_times 145 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV145_Trial145.1D.01.1D 'BLOCK(1,1)' -stim_label 145 Trial145 \
			-stim_times 146 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV146_Trial146.1D.01.1D 'BLOCK(1,1)' -stim_label 146 Trial146 \
			-stim_times 147 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV147_Trial147.1D.01.1D 'BLOCK(1,1)' -stim_label 147 Trial147 \
			-stim_times 148 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV148_Trial148.1D.01.1D 'BLOCK(1,1)' -stim_label 148 Trial148 \
			-stim_times 149 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV149_Trial149.1D.01.1D 'BLOCK(1,1)' -stim_label 149 Trial149 \
			-stim_times 150 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV150_Trial150.1D.01.1D 'BLOCK(1,1)' -stim_label 150 Trial150 \
			-stim_times 151 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV151_Trial151.1D.01.1D 'BLOCK(1,1)' -stim_label 151 Trial151 \
			-stim_times 152 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV152_Trial152.1D.01.1D 'BLOCK(1,1)' -stim_label 152 Trial152 \
			-stim_times 153 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV153_Trial153.1D.01.1D 'BLOCK(1,1)' -stim_label 153 Trial153 \
			-stim_times 154 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV154_Trial154.1D.01.1D 'BLOCK(1,1)' -stim_label 154 Trial154 \
			-stim_times 155 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV155_Trial155.1D.01.1D 'BLOCK(1,1)' -stim_label 155 Trial155 \
			-stim_times 156 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV156_Trial156.1D.01.1D 'BLOCK(1,1)' -stim_label 156 Trial156 \
			-stim_times 157 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV157_Trial157.1D.01.1D 'BLOCK(1,1)' -stim_label 157 Trial157 \
			-stim_times 158 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV158_Trial158.1D.01.1D 'BLOCK(1,1)' -stim_label 158 Trial158 \
			-stim_times 159 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV159_Trial159.1D.01.1D 'BLOCK(1,1)' -stim_label 159 Trial159 \
			-stim_times 160 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV160_Trial160.1D.01.1D 'BLOCK(1,1)' -stim_label 160 Trial160 \
			-stim_times 161 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV161_Trial161.1D.01.1D 'BLOCK(1,1)' -stim_label 161 Trial161 \
			-stim_times 162 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV162_Trial162.1D.01.1D 'BLOCK(1,1)' -stim_label 162 Trial162 \
			-stim_times 163 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV163_Trial163.1D.01.1D 'BLOCK(1,1)' -stim_label 163 Trial163 \
			-stim_times 164 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV164_Trial164.1D.01.1D 'BLOCK(1,1)' -stim_label 164 Trial164 \
			-stim_times 165 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV165_Trial165.1D.01.1D 'BLOCK(1,1)' -stim_label 165 Trial165 \
			-stim_times 166 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV166_Trial166.1D.01.1D 'BLOCK(1,1)' -stim_label 166 Trial166 \
			-stim_times 167 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV167_Trial167.1D.01.1D 'BLOCK(1,1)' -stim_label 167 Trial167 \
			-stim_times 168 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV168_Trial168.1D.01.1D 'BLOCK(1,1)' -stim_label 168 Trial168 \
			-stim_times 169 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV169_Trial169.1D.01.1D 'BLOCK(1,1)' -stim_label 169 Trial169 \
			-stim_times 170 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV170_Trial170.1D.01.1D 'BLOCK(1,1)' -stim_label 170 Trial170 \
			-stim_times 171 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV171_Trial171.1D.01.1D 'BLOCK(1,1)' -stim_label 171 Trial171 \
			-stim_times 172 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV172_Trial172.1D.01.1D 'BLOCK(1,1)' -stim_label 172 Trial172 \
			-stim_times 173 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV173_Trial173.1D.01.1D 'BLOCK(1,1)' -stim_label 173 Trial173 \
			-stim_times 174 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV174_Trial174.1D.01.1D 'BLOCK(1,1)' -stim_label 174 Trial174 \
			-stim_times 175 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV175_Trial175.1D.01.1D 'BLOCK(1,1)' -stim_label 175 Trial175 \
			-stim_times 176 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV176_Trial176.1D.01.1D 'BLOCK(1,1)' -stim_label 176 Trial176 \
			-stim_times 177 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV177_Trial177.1D.01.1D 'BLOCK(1,1)' -stim_label 177 Trial177 \
			-stim_times 178 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV178_Trial178.1D.01.1D 'BLOCK(1,1)' -stim_label 178 Trial178 \
			-stim_times 179 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV179_Trial179.1D.01.1D 'BLOCK(1,1)' -stim_label 179 Trial179 \
			-stim_times 180 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV180_Trial180.1D.01.1D 'BLOCK(1,1)' -stim_label 180 Trial180 \
			-stim_times 181 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV181_Trial181.1D.01.1D 'BLOCK(1,1)' -stim_label 181 Trial181 \
			-stim_times 182 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV182_Trial182.1D.01.1D 'BLOCK(1,1)' -stim_label 182 Trial182 \
			-stim_times 183 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV183_Trial183.1D.01.1D 'BLOCK(1,1)' -stim_label 183 Trial183 \
			-stim_times 184 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV184_Trial184.1D.01.1D 'BLOCK(1,1)' -stim_label 184 Trial184 \
			-stim_times 185 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV185_Trial185.1D.01.1D 'BLOCK(1,1)' -stim_label 185 Trial185 \
			-stim_times 186 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV186_Trial186.1D.01.1D 'BLOCK(1,1)' -stim_label 186 Trial186 \
			-stim_times 187 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV187_Trial187.1D.01.1D 'BLOCK(1,1)' -stim_label 187 Trial187 \
			-stim_times 188 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV188_Trial188.1D.01.1D 'BLOCK(1,1)' -stim_label 188 Trial188 \
			-stim_times 189 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV189_Trial189.1D.01.1D 'BLOCK(1,1)' -stim_label 189 Trial189 \
			-stim_times 190 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV190_Trial190.1D.01.1D 'BLOCK(1,1)' -stim_label 190 Trial190 \
			-stim_times 191 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV191_Trial191.1D.01.1D 'BLOCK(1,1)' -stim_label 191 Trial191 \
			-stim_times 192 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV192_Trial192.1D.01.1D 'BLOCK(1,1)' -stim_label 192 Trial192 \
			-stim_times 193 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV193_Trial193.1D.01.1D 'BLOCK(1,1)' -stim_label 193 Trial193 \
			-stim_times 194 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV194_Trial194.1D.01.1D 'BLOCK(1,1)' -stim_label 194 Trial194 \
			-stim_times 195 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV195_Trial195.1D.01.1D 'BLOCK(1,1)' -stim_label 195 Trial195 \
			-stim_times 196 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV196_Trial196.1D.01.1D 'BLOCK(1,1)' -stim_label 196 Trial196 \
			-stim_times 197 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV197_Trial197.1D.01.1D 'BLOCK(1,1)' -stim_label 197 Trial197 \
			-stim_times 198 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV198_Trial198.1D.01.1D 'BLOCK(1,1)' -stim_label 198 Trial198 \
			-stim_times 199 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV199_Trial199.1D.01.1D 'BLOCK(1,1)' -stim_label 199 Trial199 \
			-stim_times 200 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV200_Trial200.1D.01.1D 'BLOCK(1,1)' -stim_label 200 Trial200 \
			-stim_times 201 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV201_Trial201.1D.01.1D 'BLOCK(1,1)' -stim_label 201 Trial201 \
			-stim_times 202 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV202_Trial202.1D.01.1D 'BLOCK(1,1)' -stim_label 202 Trial202 \
			-stim_times 203 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV203_Trial203.1D.01.1D 'BLOCK(1,1)' -stim_label 203 Trial203 \
			-stim_times 204 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV204_Trial204.1D.01.1D 'BLOCK(1,1)' -stim_label 204 Trial204 \
			-stim_times 205 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV205_Trial205.1D.01.1D 'BLOCK(1,1)' -stim_label 205 Trial205 \
			-stim_times 206 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV206_Trial206.1D.01.1D 'BLOCK(1,1)' -stim_label 206 Trial206 \
			-stim_times 207 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV207_Trial207.1D.01.1D 'BLOCK(1,1)' -stim_label 207 Trial207 \
			-stim_times 208 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV208_Trial208.1D.01.1D 'BLOCK(1,1)' -stim_label 208 Trial208 \
			-stim_times 209 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV209_Trial209.1D.01.1D 'BLOCK(1,1)' -stim_label 209 Trial209 \
			-stim_times 210 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV210_Trial210.1D.01.1D 'BLOCK(1,1)' -stim_label 210 Trial210 \
			-stim_times 211 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV211_Trial211.1D.01.1D 'BLOCK(1,1)' -stim_label 211 Trial211 \
			-stim_times 212 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV212_Trial212.1D.01.1D 'BLOCK(1,1)' -stim_label 212 Trial212 \
			-stim_times 213 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV213_Trial213.1D.01.1D 'BLOCK(1,1)' -stim_label 213 Trial213 \
			-stim_times 214 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV214_Trial214.1D.01.1D 'BLOCK(1,1)' -stim_label 214 Trial214 \
			-stim_times 215 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV215_Trial215.1D.01.1D 'BLOCK(1,1)' -stim_label 215 Trial215 \
			-stim_times 216 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV216_Trial216.1D.01.1D 'BLOCK(1,1)' -stim_label 216 Trial216 \
			-stim_times 217 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV217_Trial217.1D.01.1D 'BLOCK(1,1)' -stim_label 217 Trial217 \
			-stim_times 218 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV218_Trial218.1D.01.1D 'BLOCK(1,1)' -stim_label 218 Trial218 \
			-stim_times 219 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV219_Trial219.1D.01.1D 'BLOCK(1,1)' -stim_label 219 Trial219 \
			-stim_times 220 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV220_Trial220.1D.01.1D 'BLOCK(1,1)' -stim_label 220 Trial220 \
			-stim_times 221 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV221_Trial221.1D.01.1D 'BLOCK(1,1)' -stim_label 221 Trial221 \
			-stim_times 222 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV222_Trial222.1D.01.1D 'BLOCK(1,1)' -stim_label 222 Trial222 \
			-stim_times 223 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV223_Trial223.1D.01.1D 'BLOCK(1,1)' -stim_label 223 Trial223 \
			-stim_times 224 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV224_Trial224.1D.01.1D 'BLOCK(1,1)' -stim_label 224 Trial224 \
			-stim_times 225 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV225_Trial225.1D.01.1D 'BLOCK(1,1)' -stim_label 225 Trial225 \
			-stim_times 226 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV226_Trial226.1D.01.1D 'BLOCK(1,1)' -stim_label 226 Trial226 \
			-stim_times 227 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV227_Trial227.1D.01.1D 'BLOCK(1,1)' -stim_label 227 Trial227 \
			-stim_times 228 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV228_Trial228.1D.01.1D 'BLOCK(1,1)' -stim_label 228 Trial228 \
			-stim_times 229 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV229_Trial229.1D.01.1D 'BLOCK(1,1)' -stim_label 229 Trial229 \
			-stim_times 230 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV230_Trial230.1D.01.1D 'BLOCK(1,1)' -stim_label 230 Trial230 \
			-stim_times 231 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV231_Trial231.1D.01.1D 'BLOCK(1,1)' -stim_label 231 Trial231 \
			-stim_times 232 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV232_Trial232.1D.01.1D 'BLOCK(1,1)' -stim_label 232 Trial232 \
			-stim_times 233 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV233_Trial233.1D.01.1D 'BLOCK(1,1)' -stim_label 233 Trial233 \
			-stim_times 234 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV234_Trial234.1D.01.1D 'BLOCK(1,1)' -stim_label 234 Trial234 \
			-stim_times 235 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV235_Trial235.1D.01.1D 'BLOCK(1,1)' -stim_label 235 Trial235 \
			-stim_times 236 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV236_Trial236.1D.01.1D 'BLOCK(1,1)' -stim_label 236 Trial236 \
			-stim_times 237 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV237_Trial237.1D.01.1D 'BLOCK(1,1)' -stim_label 237 Trial237 \
			-stim_times 238 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV238_Trial238.1D.01.1D 'BLOCK(1,1)' -stim_label 238 Trial238 \
			-stim_times 239 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV239_Trial239.1D.01.1D 'BLOCK(1,1)' -stim_label 239 Trial239 \
			-stim_times 240 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV240_Trial240.1D.01.1D 'BLOCK(1,1)' -stim_label 240 Trial240 \
			-stim_times 241 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV241_Trial241.1D.01.1D 'BLOCK(1,1)' -stim_label 241 Trial241 \
			-stim_times 242 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV242_Trial242.1D.01.1D 'BLOCK(1,1)' -stim_label 242 Trial242 \
			-stim_times 243 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV243_Trial243.1D.01.1D 'BLOCK(1,1)' -stim_label 243 Trial243 \
			-stim_times 244 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV244_Trial244.1D.01.1D 'BLOCK(1,1)' -stim_label 244 Trial244 \
			-stim_times 245 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV245_Trial245.1D.01.1D 'BLOCK(1,1)' -stim_label 245 Trial245 \
			-stim_times 246 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV246_Trial246.1D.01.1D 'BLOCK(1,1)' -stim_label 246 Trial246 \
			-stim_times 247 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV247_Trial247.1D.01.1D 'BLOCK(1,1)' -stim_label 247 Trial247 \
			-stim_times 248 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV248_Trial248.1D.01.1D 'BLOCK(1,1)' -stim_label 248 Trial248 \
			-stim_times 249 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV249_Trial249.1D.01.1D 'BLOCK(1,1)' -stim_label 249 Trial249 \
			-stim_times 250 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV250_Trial250.1D.01.1D 'BLOCK(1,1)' -stim_label 250 Trial250 \
			-stim_times 251 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV251_Trial251.1D.01.1D 'BLOCK(1,1)' -stim_label 251 Trial251 \
			-stim_times 252 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV252_Trial252.1D.01.1D 'BLOCK(1,1)' -stim_label 252 Trial252 \
			-stim_times 253 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV253_Trial253.1D.01.1D 'BLOCK(1,1)' -stim_label 253 Trial253 \
			-stim_times 254 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV254_Trial254.1D.01.1D 'BLOCK(1,1)' -stim_label 254 Trial254 \
			-stim_times 255 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV255_Trial255.1D.01.1D 'BLOCK(1,1)' -stim_label 255 Trial255 \
			-stim_times 256 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV256_Trial256.1D.01.1D 'BLOCK(1,1)' -stim_label 256 Trial256 \
			-stim_times 257 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV257_Trial257.1D.01.1D 'BLOCK(1,1)' -stim_label 257 Trial257 \
			-stim_times 258 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV258_Trial258.1D.01.1D 'BLOCK(1,1)' -stim_label 258 Trial258 \
			-stim_times 259 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV259_Trial259.1D.01.1D 'BLOCK(1,1)' -stim_label 259 Trial259 \
			-stim_times 260 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV260_Trial260.1D.01.1D 'BLOCK(1,1)' -stim_label 260 Trial260 \
			-stim_times 261 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV261_Trial261.1D.01.1D 'BLOCK(1,1)' -stim_label 261 Trial261 \
			-stim_times 262 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV262_Trial262.1D.01.1D 'BLOCK(1,1)' -stim_label 262 Trial262 \
			-stim_times 263 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV263_Trial263.1D.01.1D 'BLOCK(1,1)' -stim_label 263 Trial263 \
			-stim_times 264 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV264_Trial264.1D.01.1D 'BLOCK(1,1)' -stim_label 264 Trial264 \
			-stim_times 265 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV265_Trial265.1D.01.1D 'BLOCK(1,1)' -stim_label 265 Trial265 \
			-stim_times 266 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV266_Trial266.1D.01.1D 'BLOCK(1,1)' -stim_label 266 Trial266 \
			-stim_times 267 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV267_Trial267.1D.01.1D 'BLOCK(1,1)' -stim_label 267 Trial267 \
			-stim_times 268 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV268_Trial268.1D.01.1D 'BLOCK(1,1)' -stim_label 268 Trial268 \
			-stim_times 269 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV269_Trial269.1D.01.1D 'BLOCK(1,1)' -stim_label 269 Trial269 \
			-stim_times 270 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV270_Trial270.1D.01.1D 'BLOCK(1,1)' -stim_label 270 Trial270 \
			-stim_times 271 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV271_Trial271.1D.01.1D 'BLOCK(1,1)' -stim_label 271 Trial271 \
			-stim_times 272 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV272_Trial272.1D.01.1D 'BLOCK(1,1)' -stim_label 272 Trial272 \
			-stim_times 273 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV273_Trial273.1D.01.1D 'BLOCK(1,1)' -stim_label 273 Trial273 \
			-stim_times 274 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV274_Trial274.1D.01.1D 'BLOCK(1,1)' -stim_label 274 Trial274 \
			-stim_times 275 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV275_Trial275.1D.01.1D 'BLOCK(1,1)' -stim_label 275 Trial275 \
			-stim_times 276 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV276_Trial276.1D.01.1D 'BLOCK(1,1)' -stim_label 276 Trial276 \
			-stim_times 277 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV277_Trial277.1D.01.1D 'BLOCK(1,1)' -stim_label 277 Trial277 \
			-stim_times 278 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV278_Trial278.1D.01.1D 'BLOCK(1,1)' -stim_label 278 Trial278 \
			-stim_times 279 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV279_Trial279.1D.01.1D 'BLOCK(1,1)' -stim_label 279 Trial279 \
			-stim_times 280 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV280_Trial280.1D.01.1D 'BLOCK(1,1)' -stim_label 280 Trial280 \
			-stim_times 281 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV281_Trial281.1D.01.1D 'BLOCK(1,1)' -stim_label 281 Trial281 \
			-stim_times 282 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV282_Trial282.1D.01.1D 'BLOCK(1,1)' -stim_label 282 Trial282 \
			-stim_times 283 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV283_Trial283.1D.01.1D 'BLOCK(1,1)' -stim_label 283 Trial283 \
			-stim_times 284 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV284_Trial284.1D.01.1D 'BLOCK(1,1)' -stim_label 284 Trial284 \
			-stim_times 285 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV285_Trial285.1D.01.1D 'BLOCK(1,1)' -stim_label 285 Trial285 \
			-stim_times 286 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV286_Trial286.1D.01.1D 'BLOCK(1,1)' -stim_label 286 Trial286 \
			-stim_times 287 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV287_Trial287.1D.01.1D 'BLOCK(1,1)' -stim_label 287 Trial287 \
			-stim_times 288 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV288_Trial288.1D.01.1D 'BLOCK(1,1)' -stim_label 288 Trial288 \
			-stim_times 289 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV289_Trial289.1D.01.1D 'BLOCK(1,1)' -stim_label 289 Trial289 \
			-stim_times 290 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV290_Trial290.1D.01.1D 'BLOCK(1,1)' -stim_label 290 Trial290 \
			-stim_times 291 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV291_Trial291.1D.01.1D 'BLOCK(1,1)' -stim_label 291 Trial291 \
			-stim_times 292 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV292_Trial292.1D.01.1D 'BLOCK(1,1)' -stim_label 292 Trial292 \
			-stim_times 293 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV293_Trial293.1D.01.1D 'BLOCK(1,1)' -stim_label 293 Trial293 \
			-stim_times 294 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV294_Trial294.1D.01.1D 'BLOCK(1,1)' -stim_label 294 Trial294 \
			-stim_times 295 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV295_Trial295.1D.01.1D 'BLOCK(1,1)' -stim_label 295 Trial295 \
			-stim_times 296 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV296_Trial296.1D.01.1D 'BLOCK(1,1)' -stim_label 296 Trial296 \
			-stim_times 297 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV297_Trial297.1D.01.1D 'BLOCK(1,1)' -stim_label 297 Trial297 \
			-stim_times 298 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV298_Trial298.1D.01.1D 'BLOCK(1,1)' -stim_label 298 Trial298 \
			-stim_times 299 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV299_Trial299.1D.01.1D 'BLOCK(1,1)' -stim_label 299 Trial299 \
			-stim_times 300 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV300_Trial300.1D.01.1D 'BLOCK(1,1)' -stim_label 300 Trial300 \
			-stim_times 301 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV301_Trial301.1D.01.1D 'BLOCK(1,1)' -stim_label 301 Trial301 \
			-stim_times 302 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV302_Trial302.1D.01.1D 'BLOCK(1,1)' -stim_label 302 Trial302 \
			-stim_times 303 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV303_Trial303.1D.01.1D 'BLOCK(1,1)' -stim_label 303 Trial303 \
			-stim_times 304 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV304_Trial304.1D.01.1D 'BLOCK(1,1)' -stim_label 304 Trial304 \
			-stim_times 305 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV305_Trial305.1D.01.1D 'BLOCK(1,1)' -stim_label 305 Trial305 \
			-stim_times 306 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV306_Trial306.1D.01.1D 'BLOCK(1,1)' -stim_label 306 Trial306 \
			-stim_times 307 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV307_Trial307.1D.01.1D 'BLOCK(1,1)' -stim_label 307 Trial307 \
			-stim_times 308 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV308_Trial308.1D.01.1D 'BLOCK(1,1)' -stim_label 308 Trial308 \
			-stim_times 309 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV309_Trial309.1D.01.1D 'BLOCK(1,1)' -stim_label 309 Trial309 \
			-stim_times 310 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV310_Trial310.1D.01.1D 'BLOCK(1,1)' -stim_label 310 Trial310 \
			-stim_times 311 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV311_Trial311.1D.01.1D 'BLOCK(1,1)' -stim_label 311 Trial311 \
			-stim_times 312 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV312_Trial312.1D.01.1D 'BLOCK(1,1)' -stim_label 312 Trial312 \
			-stim_times 313 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV313_Trial313.1D.01.1D 'BLOCK(1,1)' -stim_label 313 Trial313 \
			-stim_times 314 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV314_Trial314.1D.01.1D 'BLOCK(1,1)' -stim_label 314 Trial314 \
			-stim_times 315 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV315_Trial315.1D.01.1D 'BLOCK(1,1)' -stim_label 315 Trial315 \
			-stim_times 316 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV316_Trial316.1D.01.1D 'BLOCK(1,1)' -stim_label 316 Trial316 \
			-stim_times 317 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV317_Trial317.1D.01.1D 'BLOCK(1,1)' -stim_label 317 Trial317 \
			-stim_times 318 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV318_Trial318.1D.01.1D 'BLOCK(1,1)' -stim_label 318 Trial318 \
			-stim_times 319 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV319_Trial319.1D.01.1D 'BLOCK(1,1)' -stim_label 319 Trial319 \
			-stim_times 320 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV320_Trial320.1D.01.1D 'BLOCK(1,1)' -stim_label 320 Trial320 \
			-stim_times 321 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV321_Trial321.1D.01.1D 'BLOCK(1,1)' -stim_label 321 Trial321 \
			-stim_times 322 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV322_Trial322.1D.01.1D 'BLOCK(1,1)' -stim_label 322 Trial322 \
			-stim_times 323 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV323_Trial323.1D.01.1D 'BLOCK(1,1)' -stim_label 323 Trial323 \
			-stim_times 324 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV324_Trial324.1D.01.1D 'BLOCK(1,1)' -stim_label 324 Trial324 \
			-stim_times 325 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV325_Trial325.1D.01.1D 'BLOCK(1,1)' -stim_label 325 Trial325 \
			-stim_times 326 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV326_Trial326.1D.01.1D 'BLOCK(1,1)' -stim_label 326 Trial326 \
			-stim_times 327 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV327_Trial327.1D.01.1D 'BLOCK(1,1)' -stim_label 327 Trial327 \
			-stim_times 328 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV328_Trial328.1D.01.1D 'BLOCK(1,1)' -stim_label 328 Trial328 \
			-stim_times 329 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV329_Trial329.1D.01.1D 'BLOCK(1,1)' -stim_label 329 Trial329 \
			-stim_times 330 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV330_Trial330.1D.01.1D 'BLOCK(1,1)' -stim_label 330 Trial330 \
			-stim_times 331 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV331_Trial331.1D.01.1D 'BLOCK(1,1)' -stim_label 331 Trial331 \
			-stim_times 332 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV332_Trial332.1D.01.1D 'BLOCK(1,1)' -stim_label 332 Trial332 \
			-stim_times 333 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV333_Trial333.1D.01.1D 'BLOCK(1,1)' -stim_label 333 Trial333 \
			-stim_times 334 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV334_Trial334.1D.01.1D 'BLOCK(1,1)' -stim_label 334 Trial334 \
			-stim_times 335 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV335_Trial335.1D.01.1D 'BLOCK(1,1)' -stim_label 335 Trial335 \
			-stim_times 336 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV336_Trial336.1D.01.1D 'BLOCK(1,1)' -stim_label 336 Trial336 \
			-stim_times 337 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV337_Trial337.1D.01.1D 'BLOCK(1,1)' -stim_label 337 Trial337 \
			-stim_times 338 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV338_Trial338.1D.01.1D 'BLOCK(1,1)' -stim_label 338 Trial338 \
			-stim_times 339 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV339_Trial339.1D.01.1D 'BLOCK(1,1)' -stim_label 339 Trial339 \
			-stim_times 340 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV340_Trial340.1D.01.1D 'BLOCK(1,1)' -stim_label 340 Trial340 \
			-stim_times 341 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV341_Trial341.1D.01.1D 'BLOCK(1,1)' -stim_label 341 Trial341 \
			-stim_times 342 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV342_Trial342.1D.01.1D 'BLOCK(1,1)' -stim_label 342 Trial342 \
			-stim_times 343 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV343_Trial343.1D.01.1D 'BLOCK(1,1)' -stim_label 343 Trial343 \
			-stim_times 344 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV344_Trial344.1D.01.1D 'BLOCK(1,1)' -stim_label 344 Trial344 \
			-stim_times 345 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV345_Trial345.1D.01.1D 'BLOCK(1,1)' -stim_label 345 Trial345 \
			-stim_times 346 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV346_Trial346.1D.01.1D 'BLOCK(1,1)' -stim_label 346 Trial346 \
			-stim_times 347 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV347_Trial347.1D.01.1D 'BLOCK(1,1)' -stim_label 347 Trial347 \
			-stim_times 348 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV348_Trial348.1D.01.1D 'BLOCK(1,1)' -stim_label 348 Trial348 \
			-stim_times 349 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV349_Trial349.1D.01.1D 'BLOCK(1,1)' -stim_label 349 Trial349 \
			-stim_times 350 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV350_Trial350.1D.01.1D 'BLOCK(1,1)' -stim_label 350 Trial350 \
			-stim_times 351 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV351_Trial351.1D.01.1D 'BLOCK(1,1)' -stim_label 351 Trial351 \
			-stim_times 352 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV352_Trial352.1D.01.1D 'BLOCK(1,1)' -stim_label 352 Trial352 \
			-stim_times 353 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV353_Trial353.1D.01.1D 'BLOCK(1,1)' -stim_label 353 Trial353 \
			-stim_times 354 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV354_Trial354.1D.01.1D 'BLOCK(1,1)' -stim_label 354 Trial354 \
			-stim_times 355 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV355_Trial355.1D.01.1D 'BLOCK(1,1)' -stim_label 355 Trial355 \
			-stim_times 356 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV356_Trial356.1D.01.1D 'BLOCK(1,1)' -stim_label 356 Trial356 \
			-stim_times 357 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV357_Trial357.1D.01.1D 'BLOCK(1,1)' -stim_label 357 Trial357 \
			-stim_times 358 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV358_Trial358.1D.01.1D 'BLOCK(1,1)' -stim_label 358 Trial358 \
			-stim_times 359 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV359_Trial359.1D.01.1D 'BLOCK(1,1)' -stim_label 359 Trial359 \
			-stim_times 360 ${StimFileDir}/stime_${subjNum}_stimfile_RITLBehavCorrTrialCount_EV360_Trial360.1D.01.1D 'BLOCK(1,1)' -stim_label 360 Trial360 \
			-xsave -x1D xmat_rall.x1D -xjpeg xmat_rall.jpg \
			-fout -tout \
			-jobs 4 -float \
			-cbucket GLMTrialCount_coefficients_${nextInputFilename} \
			-bucket GLMTrialCount_outbucket_${nextInputFilename} -overwrite

		
		popd
	fi


	##Spatially smooth the data
	execute=0
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}

		echo "-Spatially smooth data-"
		3dcalc -overwrite -a ${subjMaskDIR}/${subjNum}_gmMask_func_dil1vox+tlrc -b ${nextInputFilename}+tlrc -expr 'a*b' -prefix masked_${nextInputFilename}
		3dmerge -overwrite -doall -quiet -1filter_nzmean $FWHMSmoothing -prefix nzmeanSM_masked_${nextInputFilename} masked_${nextInputFilename}+tlrc
		#3dmerge -overwrite -doall -quiet -1blur_fwhm $FWHMSmoothing -prefix sm_${nextInputFilename} ${nextInputFilename}+tlrc
		#3dBlurInMask -input ${nextInputFilename}+tlrc -FWHM $FWHMSmoothing -mask ${subjMaskDIR}/${subjNum}_gmMask_func_dil1vox+tlrc -prefix smInMask_${nextInputFilename}
		#rm -v ${nextInputFilename}.????

		popd
	fi
	#nextInputFilename="nzmeanSM_masked_"${nextInputFilename}


	##Normalize signal to % signal change (over the entire time series; might be better to do for each run separately, but GLM will account for inter-run differences anyway)
	execute=0
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}

		echo "-Normalize data (% signal change)-"
# 		#Downsample mprage mask to functional space
# 		3dresample -master ${nextInputFilename}+tlrc -prefix ${subjMaskDIR}/wholebrain_mask_func+tlrc -inset ${subjMaskDIR}/wholebrain_mask+tlrc -overwrite -rmode 'NN'
# 		#Dilate mask by 1 functional voxel (just in case the resampled anatomical mask is off by a bit)
# 		3dLocalstat -overwrite -nbhd 'SPHERE(-1)' -stat 'max' -prefix ${subjMaskDIR}/wholebrain_mask_func_dil1vox+tlrc ${subjMaskDIR}/wholebrain_mask_func+tlrc
		#Normalizes to % signal, with mean value of 100 and max value of 200 (like afni_proc.py)
		3dTstat -overwrite -mean -prefix mean_${nextInputFilename} ${nextInputFilename}+tlrc
		3dcalc -overwrite -a ${nextInputFilename}+tlrc -b mean_${nextInputFilename}+tlrc -c ${subjMaskDIR}/${subjNum}_gmMask_func_dil1vox+tlrc -expr 'min(200, (100 * (a/b)))*c' -prefix normPSC_${nextInputFilename}
		#rm -v ${nextInputFilename}.????

		popd
	fi
	nextInputFilename="normPSC_"${nextInputFilename}




	##Convert to NIfTI
	execute=0
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}

		echo "Converting to NIfTI"
		#rm ${nextInputFilename}.nii.gz
		#3dcopy ${nextInputFilename}+tlrc ${nextInputFilename}.nii.gz
		#rm GLMTrialCount_coefficients_${nextInputFilename}.nii.gz
		#3dcopy GLMTrialCount_coefficients_${nextInputFilename}+tlrc ${nextInputFilename}.nii.gz
		
		rm GLMTrialCount_coefficients_${nextInputFilename}.nii.gz
		rm GLMTrialCount_coefficients_normPSC_stc_epi_allruns_tlrc_al_360coefs.nii.gz
		3dcalc -a GLMTrialCount_coefficients_normPSC_stc_epi_allruns_tlrc_al+tlrc.'[40..$]' -expr 'a' -prefix GLMTrialCount_coefficients_normPSC_stc_epi_allruns_tlrc_al_360coefs.nii.gz
		
		3dTcat -prefix GLMTrialCount_coefficients_normPSC_stc_epi_allruns_tlrc_al_360coefs_4d.nii.gz GLMTrialCount_coefficients_normPSC_stc_epi_allruns_tlrc_al_360coefs.nii.gz
		rm GLMTrialCount_coefficients_normPSC_stc_epi_allruns_tlrc_al_360coefs.nii.gz

		popd
	fi

done
