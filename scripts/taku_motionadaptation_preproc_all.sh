#!/bin/bash
# Taku Ito
# MotionAdaptation TACS collaboration with Krekelberg lab


##List of AFNI help files: http://afni.nimh.nih.gov/afni/doc/program_help/index.html

##--Analysis parameters--
#Subject numbers: 5 6 7 8 9 10 11 12 14
listOfSubjects="144"
 #038 069 141 172 173 177 178 083 144 170"
#"" 
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
#FWHM smoothing parameter to use
FWHMSmoothing=6
#Name of your current analysis
ANALYSISNAME="tacs_motionadaptation"


##Running functions on each subject separately
for subjNum in $listOfSubjects
do
        if [ "$subjNum" == "038" ]; then
            # Subject parameters
            # T1
            seriesDirString="sess2-0002"
            # EPIs
            seriesDirString1="sess2-0006 sess2-0007 sess2-0008 sess2-0009"
        elif [ "$subjNum" == "069" ]; then
            # T1
            seriesDirString="sess1-0003"
            # EPIs
            seriesDirString1="sess1-0008 sess1-0009 sess1-0010 sess1-0011"
        elif [ "$subjNum" == "083" ]; then
            # T1
            seriesDirString="Pos"
            # EPIs
            seriesDirString1="sess1-0004 sess1-0005 sess1-0007 sess1-0008"
        elif [ "$subjNum" == "141" ]; then
            # T1
            seriesDirString="sess1-0002"
            # EPIs
            seriesDirString1="sess1-0005 sess1-0006 sess1-0007 sess1-0008"
        elif [ "$subjNum" == "144" ]; then
            # T1
            seriesDirString="Pos"
            # EPIs
            seriesDirString1="sess1-0004 sess1-0005 sess1-0007 sess1-0008"
        elif [ "$subjNum" == "170" ]; then
            # T1
            seriesDirString="Pos"
            # EPIs
            seriesDirString1="sess1-0004 sess1-0005 sess1-0007 sess1-0008"
        elif [ "$subjNum" == "172" ]; then
            # T1
            seriesDirString="sess1-0003"
            # EPIs
            seriesDirString1="sess1-0006 sess1-0007 sess1-0008 sess1-0009"
        elif [ "$subjNum" == "173" ]; then
            # T1
            seriesDirString="sess1-0002"
            # EPIs
            seriesDirString1="sess1-0006 sess1-0007 sess1-0008 sess1-0009"
        elif [ "$subjNum" == "177" ]; then
            # T1
            seriesDirString="sess1-0003"
            # EPIs
            seriesDirString1="sess1-0006 sess1-0007 sess1-0009 sess1-0010"
        elif [ "$subjNum" == "178" ]; then
            # T1
            seriesDirString="sess1-0002"
            # EPIs
            seriesDirString1="sess1-0005 sess1-0006 sess1-0007 sess1-0008"

        fi


        
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
	execute=0
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
		mri_convert ${subjRawDataDIR}/dcm/*${seriesDirString}*-00001.dcm --in_type siemens --out_type nii mprage.nii.gz
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
		@auto_tlrc -base ${atlasDIR}/MNI152_1mm_uni+tlrc -input anat_mprage_skullstripped+orig -no_ss
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
		execute=0
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
                            echo "1"
                        else
                            numTRs=$numTRsPRT
                            echo "2"
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
		execute=0
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
	execute=0
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
	execute=0
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}

		echo "-Run align_epi_anat.py to align EPIs to MPRAGE, motion correct, and Talairach transform EPIs (output in 333 space)-"
		echo "Make sure Python is version 2.6 or greater"
		python -V
		#Correcting for motion, aligning fMRI data to MPRAGE, and aligning fMRI data to Talairach template [applying all transformation at once reduces reslicing artifacts]
		#[You could alternatively analyze all of the data, then Talairach transform the statistics (though this would make extraction of time series based on Talairached ROIs difficult)]
		#Visit for more info: http://afni.nimh.nih.gov/pub/dist/doc/program_help/align_epi_anat.py.html
		align_epi_anat.py -overwrite -anat anat_mprage_skullstripped+orig -epi ${nextInputFilename}+orig -epi_base 10 -epi2anat -anat_has_skull no -AddEdge -epi_strip 3dSkullStrip -ex_mode quiet -volreg on -deoblique on -tshift off -tlrc_apar anat_mprage_skullstripped+tlrc -master_tlrc ${atlasDIR}/MNI_EPI_333+tlrc

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
	execute=0
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
	execute=0
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
	execute=0
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
	execute=0
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

            echo "Computing derivatives of WM and ventricle timeseries"
            1d_tool.py -overwrite -infile ${subjNum}_WM_timeseries_task.1D -derivative -write ${subjNum}_WM_timeseries_task_deriv.1D
            1d_tool.py -overwrite -infile ${subjNum}_ventricles_timeseries_task.1D -derivative -write ${subjNum}_ventricles_timeseries_task_deriv.1D

            # Split all task timings into 3 separate regressor 1D files
            # Opposite stimulus timings
            cat ../sdm/allTaskTimings.txt | awk '{ print $1 }' > ../sdm/OppTimes.1D
            cat ../sdm/allTaskTimings.txt | awk '{ print $2 }' > ../sdm/SameTimes.1D
            cat ../sdm/allTaskTimings.txt | awk '{ print $3 }' > ../sdm/LongAdaptTimes.1D

            echo "Run GLM"
            3dDeconvolve \
                    -input ${nextInputFilename}+tlrc \
                    -mask ${subjMaskDIR}/${subjNum}_gmMask_func_dil1vox+tlrc \
                    -concat "${concatString}" \
                    -polort A \
                    -num_stimts 13 \
                    -stim_file 1 ${subjfMRIDIR}/allruns_motion_params.1D'[0]' -stim_base 1 \
                    -stim_file 2 ${subjfMRIDIR}/allruns_motion_params.1D'[1]' -stim_base 2 \
                    -stim_file 3 ${subjfMRIDIR}/allruns_motion_params.1D'[2]' -stim_base 3 \
                    -stim_file 4 ${subjfMRIDIR}/allruns_motion_params.1D'[3]' -stim_base 4 \
                    -stim_file 5 ${subjfMRIDIR}/allruns_motion_params.1D'[4]' -stim_base 5 \
                    -stim_file 6 ${subjfMRIDIR}/allruns_motion_params.1D'[5]' -stim_base 6 \
                    -stim_file 7 ${subjfMRIDIR}/${subjNum}_ventricles_timeseries_task.1D -stim_base 7 \
                    -stim_file 8 ${subjfMRIDIR}/${subjNum}_ventricles_timeseries_task_deriv.1D -stim_base 8 \
                    -stim_file 9 ${subjfMRIDIR}/${subjNum}_WM_timeseries_task.1D -stim_base 9 \
                    -stim_file 10 ${subjfMRIDIR}/${subjNum}_WM_timeseries_task_deriv.1D -stim_base 10 \
                    -stim_file 11 ${subjfMRIDIR}/../sdm/OppTimes.1D -stim_label 11 Opp \
                    -stim_file 12 ${subjfMRIDIR}/../sdm/SameTimes.1D -stim_label 12 Same \
                    -stim_file 13 ${subjfMRIDIR}/../sdm/LongAdaptTimes.1D -stim_label 13 LongAdapt \
                    -xsave -x1D xmat_rall.x1D -xjpeg xmat_rall.jpg \
                    -fout -tout \
                    -jobs 4 -float \
                    -errts resid_${nextInputFilename} \
                    -cbucket glm_cbucket_${nextInputFilename} \
                    -bucket glm_outbucket_${nextInputFilename} -overwrite

            
            popd
	fi
        nextInputFilename=resid_${nextInputFilename}


	##Spatially smooth the data
	execute=0
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}

		echo "-Spatially smooth data-"
                3dBlurInMask -input ${nextInputFilename}+tlrc -FWHM $FWHMSmoothing -mask ${subjMaskDIR}/${subjNum}_gmMask_func_dil1vox+tlrc -float -prefix smInMask_${nextInputFilename}
		#rm -v ${nextInputFilename}.????

		popd
	fi
	nextInputFilename="smInMask_"${nextInputFilename}


        ## Split data into notACS and yestACS
        execute=1
        if [ $execute -eq 1 ]; then
            pushd ${subjfMRIDIR}

            echo "-Split timeseries according to yes/no tACS-"
            TRs1=`cat ../sdm/TaskTiming1.txt | wc -l`
            TRs2=`cat ../sdm/TaskTiming2.txt | wc -l`
            TRs3=`cat ../sdm/TaskTiming3.txt | wc -l`
            TRs4=`cat ../sdm/TaskTiming4.txt | wc -l`

            notACS=$(expr $TRs1 + $TRs2 - 1) # Since AFNI counts from 0
            yestACS=$(expr $notACS + 1)
            #yestACS=$(expr $TRs3 + $TRs4)

            3dcalc -a ${nextInputFilename}+tlrc[0..$notACS] -expr a -prefix notACS_${nextInputFilename}+tlrc -overwrite
            3dcalc -a ${nextInputFilename}+tlrc[$yestACS..$] -expr a -prefix yestACS_${nextInputFilename}+tlrc -overwrite

            popd
        fi

        execute=1
        if [ $execute -eq 1 ]; then
            pushd ${subjAnalysisDIR}

            echo "-Running seed-based correlation using right and left area MT on both stim and nostim conditions-"
            
            # First extract L and R area MT average timeseries on either condition
            # no stim first
            3dmaskave -quiet -mask ../../../areaMT_Left+tlrc ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc > nostim_leftAreaMT.1D
            3dmaskave -quiet -mask ../../../areaMT_Right+tlrc ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc > nostim_rightAreaMT.1D
            # stim
            3dmaskave -quiet -mask ../../../areaMT_Left+tlrc ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc > stim_leftAreaMT.1D
            3dmaskave -quiet -mask ../../../areaMT_Right+tlrc ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc > stim_rightAreaMT.1D

            #### 4 Output correlation maps
            # Polort = because of quadratic fitting of the baseline + drifting
            # First no stimulation left areaMT
            3dfim+ -input ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc \
                -polort 2 \
                -ideal_file nostim_leftAreaMT.1D \
                -out Correlation \
                -bucket nostim_leftareaMT_corrbucket -overwrite
            # no stimulation right areaMT
            3dfim+ -input ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc \
                -polort 2 \
                -ideal_file nostim_rightAreaMT.1D \
                -out Correlation \
                -bucket nostim_rightareaMT_corrbucket -overwrite

            # stimulation left areaMT
            3dfim+ -input ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc \
                -polort 2 \
                -ideal_file stim_leftAreaMT.1D \
                -out Correlation \
                -bucket stim_leftareaMT_corrbucket -overwrite
            # stimulation right areaMT
            3dfim+ -input ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc \
                -polort 2 \
                -ideal_file stim_rightAreaMT.1D \
                -out Correlation \
                -bucket stim_rightareaMT_corrbucket -overwrite
       
            echo "-Run FisherZ transform on data-"
            3dcalc -a nostim_leftareaMT_corrbucket+tlrc -expr 'atanh(a)' -prefix fz_nostim_leftareaMT_corrbucket -overwrite
            3dcalc -a nostim_rightareaMT_corrbucket+tlrc -expr 'atanh(a)' -prefix fz_nostim_rightareaMT_corrbucket -overwrite
            3dcalc -a stim_leftareaMT_corrbucket+tlrc -expr 'atanh(a)' -prefix fz_stim_leftareaMT_corrbucket -overwrite
            3dcalc -a stim_rightareaMT_corrbucket+tlrc -expr 'atanh(a)' -prefix fz_stim_rightareaMT_corrbucket -overwrite

            popd
        fi
        

done


execute=1
if [ $execute -eq 1 ]; then

    echo "-Running ANOVA-"
    
    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    pushd $datadir/results/

    3dANOVA3 -DAFNI_FLOATIZE=YES -type 4 -alevels 2 -blevels 2 -clevels 10 \
        -dset 1 1 1 ${datadir}/038/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
        -dset 1 2 1 ${datadir}/038/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
        -dset 2 1 1 ${datadir}/038/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
        -dset 2 2 1 ${datadir}/038/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
        -dset 1 1 2 ${datadir}/069/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
        -dset 1 2 2 ${datadir}/069/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
        -dset 2 1 2 ${datadir}/069/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
        -dset 2 2 2 ${datadir}/069/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
        -dset 1 1 3 ${datadir}/141/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
        -dset 1 2 3 ${datadir}/141/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
        -dset 2 1 3 ${datadir}/141/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
        -dset 2 2 3 ${datadir}/141/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
        -dset 1 1 4 ${datadir}/172/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
        -dset 1 2 4 ${datadir}/172/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
        -dset 2 1 4 ${datadir}/172/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
        -dset 2 2 4 ${datadir}/172/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
        -dset 1 1 5 ${datadir}/173/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
        -dset 1 2 5 ${datadir}/173/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
        -dset 2 1 5 ${datadir}/173/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
        -dset 2 2 5 ${datadir}/173/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
        -dset 1 1 6 ${datadir}/177/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
        -dset 1 2 6 ${datadir}/177/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
        -dset 2 1 6 ${datadir}/177/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
        -dset 2 2 6 ${datadir}/177/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
        -dset 1 1 7 ${datadir}/178/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
        -dset 1 2 7 ${datadir}/178/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
        -dset 2 1 7 ${datadir}/178/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
        -dset 2 2 7 ${datadir}/178/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
        -dset 1 1 8 ${datadir}/083/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
        -dset 1 2 8 ${datadir}/083/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
        -dset 2 1 8 ${datadir}/083/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
        -dset 2 2 8 ${datadir}/083/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
        -dset 1 1 9 ${datadir}/144/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
        -dset 1 2 9 ${datadir}/144/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
        -dset 2 1 9 ${datadir}/144/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
        -dset 2 2 9 ${datadir}/144/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
        -dset 1 1 10 ${datadir}/170/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
        -dset 1 2 10 ${datadir}/170/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
        -dset 2 1 10 ${datadir}/170/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
        -dset 2 2 10 ${datadir}/170/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
        -fa maineffect_stimulation -fb maineffect_leftright -fab interaction_stimLR \
        -amean 1 Stimulation -amean 2 NoStimulation -bmean 1 Left -bmean 2 Right \
        -adiff 1 2 diff_Stim_V_NoStim -bdiff 1 2 diff_Left_V_Right \
        -bucket ANOVA_stimulation_LR_subjN10 -overwrite

    popd
fi


execute=1
if [ $execute -eq 1 ]; then

    echo "-Running T-Tests-"
    
    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    pushd $datadir/results/

    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/069/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/141/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/172/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/173/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/177/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/178/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/083/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/144/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/178/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    -setB ${datadir}/038/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/069/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/141/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/172/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/173/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/177/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/178/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/083/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/144/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/170/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    -prefix stim_LeftVRight_ttest_n10 -paired -overwrite


    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/069/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/141/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/172/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/173/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/177/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/178/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/083/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/144/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/170/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    -setB ${datadir}/038/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/069/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/141/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/172/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/173/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/177/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/178/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/083/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/144/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/170/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    -prefix nostim_LeftVRight_ttest_n10 -paired -overwrite

    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/069/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/141/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/172/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/173/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/177/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/178/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/083/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/144/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    ${datadir}/170/${analdir}/fz_stim_leftareaMT_corrbucket+tlrc \
    -setB ${datadir}/038/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/069/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/141/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/172/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/173/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/177/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/178/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/083/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/144/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    ${datadir}/170/${analdir}/fz_nostim_leftareaMT_corrbucket+tlrc \
    -prefix Left_stimVnostim_ttest_n10 -paired -overwrite


    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/069/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/141/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/172/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/173/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/177/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/178/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/083/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/144/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    ${datadir}/170/${analdir}/fz_stim_rightareaMT_corrbucket+tlrc \
    -setB ${datadir}/038/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/069/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/141/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/172/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/173/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/177/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/178/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/083/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/144/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    ${datadir}/170/${analdir}/fz_nostim_rightareaMT_corrbucket+tlrc \
    -prefix Right_stimVnostim_ttest_n10 -paired -overwrite
fi

