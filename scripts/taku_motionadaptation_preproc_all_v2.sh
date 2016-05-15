#!/bin/bash
# Taku Ito
# MotionAdaptation TACS collaboration with Krekelberg lab


##List of AFNI help files: http://afni.nimh.nih.gov/afni/doc/program_help/index.html

##--Analysis parameters--
#Subject numbers: 5 6 7 8 9 10 11 12 14
listOfSubjects="038 069 141 172 173 177 178 083 144 170"
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
        runOnsets=""
	TRCount=0
        runNum=1
	for epirun in $seriesDirString1
	do
	
		echo "--Run ${runNum}--"
		
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
	        numTRs=$(expr $numTRs - 1)	
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
                runOnsets="${runOnsets} ${TRCount}"
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
	execute=0
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

        echo "${runOnsets}" > ${subjfMRIDIR}/../sdm/RunOnsets.txt
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
                    -num_stimts 14 \
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
                    -stim_file 14 ${subjfMRIDIR}/../sdm/all_tACS_binary_timings.txt -stim_label 14 tACS \
                    -xsave -x1D xmat_rall.x1D -xjpeg xmat_rall.jpg \
                    -fout -tout \
                    -jobs 4 -float \
                    -errts resid2_tACSreg_${nextInputFilename} \
                    -cbucket glm2_tACSreg_cbucket_${nextInputFilename} \
                    -bucket glm2_tACSreg_outbucket_${nextInputFilename} -overwrite

            
            popd
	fi
        nextInputFilename=resid2_tACSreg_${nextInputFilename}


	##Spatially smooth the data
	execute=0
	if [ $execute -eq 1 ]; then
		pushd ${subjfMRIDIR}

		echo "-Spatially smooth data-"
                3dBlurInMask -input ${nextInputFilename}+tlrc -FWHM $FWHMSmoothing -mask ${subjMaskDIR}/${subjNum}_gmMask_func_dil1vox+tlrc -float -prefix smInMask_${nextInputFilename} -overwrite
		#rm -v ${nextInputFilename}.????

		popd
	fi
	nextInputFilename="smInMask_"${nextInputFilename}

        pushd ${subjfMRIDIR}
     #   3dAFNItoNIFTI -prefix ${nextInputFilename} ${nextInputFilename}+tlrc
     #    3dAFNItoNIFTI -prefix yestACS_${nextInputFilename} yestACS_${nextInputFilename}+tlrc
     #    3dAFNItoNIFTI -prefix notACS${nextInputFilename} notACS_${nextInputFilename}+tlrc
         popd

        ## Split data into notACS and yestACS
        execute=0
        if [ $execute -eq 1 ]; then
            pushd ${subjfMRIDIR}

            echo "-Split timeseries according to yes/no tACS-"
            TRs1=`3dinfo -nv epi_short_r1+orig`
            TRs2=`3dinfo -nv epi_short_r2+orig`
            TRs3=`3dinfo -nv epi_short_r3+orig`
            TRs4=`3dinfo -nv epi_short_r4+orig`

            notACS=$(expr $TRs1 + $TRs2 - 1) # Since AFNI counts from 0
            yestACS=$(expr $notACS + 1)
            #yestACS=$(expr $TRs3 + $TRs4)

            3dcalc -a ${nextInputFilename}+tlrc[0..$notACS] -expr a -prefix notACS_${nextInputFilename}+tlrc -overwrite
            3dcalc -a ${nextInputFilename}+tlrc[$yestACS..$] -expr a -prefix yestACS_${nextInputFilename}+tlrc -overwrite

            popd
        fi


        execute=0
        if [ $execute -eq 1 ]; then
            pushd ${subjAnalysisDIR}

            echo "-Running seed-based correlation using right and left area MT on both stim and nostim conditions-"
            
            # First extract L and R area MT average timeseries on either condition
            # no stim first
            3dmaskave -quiet -mask ../../../areaMT_Left+tlrc ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc > nostim_leftAreaMT_v2.1D
            3dmaskave -quiet -mask ../../../areaMT_Right+tlrc ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc > nostim_rightAreaMT_v2.1D
            # stim
            3dmaskave -quiet -mask ../../../areaMT_Left+tlrc ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc > stim_leftAreaMT_v2.1D
            3dmaskave -quiet -mask ../../../areaMT_Right+tlrc ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc > stim_rightAreaMT_v2.1D

            #### 4 Output correlation maps
            # Polort = because of quadratic fitting of the baseline + drifting
            # First no stimulation left areaMT
#            3dfim+ -input ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc \
#                -polort 2 \
#                -ideal_file nostim_leftAreaMT_v2.1D \
#                -out Correlation \
#                -bucket nostim_leftareaMT_corrbucket_v2 -overwrite
#            # no stimulation right areaMT
#            3dfim+ -input ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc \
#                -polort 2 \
#                -ideal_file nostim_rightAreaMT_v2.1D \
#                -out Correlation \
#                -bucket nostim_rightareaMT_corrbucket_v2 -overwrite
#
#            # stimulation left areaMT
#            3dfim+ -input ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc \
#                -polort 2 \
#                -ideal_file stim_leftAreaMT_v2.1D \
#                -out Correlation \
#                -bucket stim_leftareaMT_corrbucket_v2 -overwrite
#            # stimulation right areaMT
#            3dfim+ -input ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc \
#                -polort 2 \
#                -ideal_file stim_rightAreaMT_v2.1D \
#                -out Correlation \
#                -bucket stim_rightareaMT_corrbucket_v2 -overwrite
       
            3dTcorr1D -pearson -overwrite -prefix nostim_leftareaMT_corrbucket_v2 ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc nostim_leftAreaMT_v2.1D 
            3dTcorr1D -pearson -overwrite -prefix nostim_rightareaMT_corrbucket_v2 ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc nostim_rightAreaMT_v2.1D 
            3dTcorr1D -pearson -overwrite -prefix stim_leftareaMT_corrbucket_v2 ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc stim_leftAreaMT_v2.1D 
            3dTcorr1D -pearson -overwrite -prefix stim_rightareaMT_corrbucket_v2 ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc stim_rightAreaMT_v2.1D 

            echo "-Run FisherZ transform on data-"
            3dcalc -a nostim_leftareaMT_corrbucket_v2+tlrc -expr 'atanh(a)' -prefix fz_nostim_leftareaMT_corrbucket_v2 -overwrite
            3dcalc -a nostim_rightareaMT_corrbucket_v2+tlrc -expr 'atanh(a)' -prefix fz_nostim_rightareaMT_corrbucket_v2 -overwrite
            3dcalc -a stim_leftareaMT_corrbucket_v2+tlrc -expr 'atanh(a)' -prefix fz_stim_leftareaMT_corrbucket_v2 -overwrite
            3dcalc -a stim_rightareaMT_corrbucket_v2+tlrc -expr 'atanh(a)' -prefix fz_stim_rightareaMT_corrbucket_v2 -overwrite

            popd
        fi
        

done


execute=0
if [ $execute -eq 1 ]; then

    echo "-Running ANOVA-"
    
    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    pushd $datadir/results/

    3dANOVA3 -DAFNI_FLOATIZE=YES -type 4 -alevels 2 -blevels 2 -clevels 10 \
        -dset 1 1 1 ${datadir}/038/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 1 ${datadir}/038/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 1 ${datadir}/038/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 1 ${datadir}/038/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 2 ${datadir}/069/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 2 ${datadir}/069/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 2 ${datadir}/069/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 2 ${datadir}/069/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 3 ${datadir}/141/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 3 ${datadir}/141/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 3 ${datadir}/141/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 3 ${datadir}/141/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 4 ${datadir}/172/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 4 ${datadir}/172/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 4 ${datadir}/172/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 4 ${datadir}/172/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 5 ${datadir}/173/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 5 ${datadir}/173/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 5 ${datadir}/173/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 5 ${datadir}/173/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 6 ${datadir}/177/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 6 ${datadir}/177/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 6 ${datadir}/177/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 6 ${datadir}/177/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 7 ${datadir}/178/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 7 ${datadir}/178/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 7 ${datadir}/178/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 7 ${datadir}/178/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 8 ${datadir}/083/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 8 ${datadir}/083/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 8 ${datadir}/083/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 8 ${datadir}/083/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 9 ${datadir}/144/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 9 ${datadir}/144/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 9 ${datadir}/144/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 9 ${datadir}/144/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 10 ${datadir}/170/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 10 ${datadir}/170/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 10 ${datadir}/170/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 10 ${datadir}/170/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -fa maineffect_stimulation -fb maineffect_leftright -fab interaction_stimLR \
        -amean 1 Stimulation -amean 2 NoStimulation -bmean 1 Left -bmean 2 Right \
        -adiff 1 2 diff_Stim_V_NoStim -bdiff 1 2 diff_Left_V_Right \
        -bucket ANOVA_stimulation_LR_subjN10_v2 -overwrite

    popd
fi


execute=0
if [ $execute -eq 1 ]; then

    echo "-Running T-Tests-"
    
    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    pushd $datadir/results/

    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    -setB ${datadir}/038/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    -prefix stim_LeftVRight_ttest_n10_v2 -paired -overwrite


    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    -setB ${datadir}/038/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    -prefix nostim_LeftVRight_ttest_n10_v2 -paired -overwrite

    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    -setB ${datadir}/038/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    -prefix Left_stimVnostim_ttest_n10_v2 -paired -overwrite


    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    -setB ${datadir}/038/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    -prefix Right_stimVnostim_ttest_n10_v2 -paired -overwrite
fi

# Run 3dClustSim
execute=0
if [ $execute -eq 1 ]; then
    
    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    pushd $datadir/results/
    
    3dClustSim -mask Left_stimVnostim_ttest_n10_v2+tlrc -fwhm 6 -prefix 3dclustsim_l_ttest_stimVnostim
    3dClustSim -mask ANOVA_stimulation_LR_subjN10_v2+tlrc -fwhm 6 -prefix 3dclustsim_ANOVA

fi



# Volume to Surface Mapping of Results
execute=0
if [ $execute -eq 1 ]; then

    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    pushd $datadir/results/
    surfaceatlas=/usr/local/workbench/atlases/Conte69_atlas-v2.LR.32k_fs_LR.wb/32k_ConteAtlas_v2/

    # left stimVnoStim
    wb_command -volume-to-surface-mapping Clust_l_stimVnoStim_v2_mask.nii ${surfaceatlas}/Conte69.L.midthickness.32k_fs_LR.surf.gii workbench_dir/Clust_l_stimVnoStim_v2_mask_L.shape.gii -trilinear 
    wb_command -volume-to-surface-mapping Clust_l_stimVnoStim_v2_mask.nii ${surfaceatlas}/Conte69.R.midthickness.32k_fs_LR.surf.gii workbench_dir/Clust_l_stimVnoStim_v2_mask_R.shape.gii -trilinear 
    # right stimVnoStim
    wb_command -volume-to-surface-mapping Clust_r_stimVnoStim_v2_mask.nii ${surfaceatlas}/Conte69.L.midthickness.32k_fs_LR.surf.gii workbench_dir/Clust_r_stimVnoStim_v2_mask_L.shape.gii -trilinear 
    wb_command -volume-to-surface-mapping Clust_r_stimVnoStim_v2_mask.nii ${surfaceatlas}/Conte69.R.midthickness.32k_fs_LR.surf.gii workbench_dir/Clust_r_stimVnoStim_v2_mask_R.shape.gii -trilinear 
   
    # stim left v right
    wb_command -volume-to-surface-mapping Clust_stim_leftVright_v2_mask.nii ${surfaceatlas}/Conte69.L.midthickness.32k_fs_LR.surf.gii workbench_dir/Clust_stim_leftVright_v2_mask_L.shape.gii -trilinear 
    wb_command -volume-to-surface-mapping Clust_stim_leftVright_v2_mask.nii ${surfaceatlas}/Conte69.R.midthickness.32k_fs_LR.surf.gii workbench_dir/Clust_stim_leftVright_v2_mask_R.shape.gii -trilinear 
    # noStim left v right
    wb_command -volume-to-surface-mapping Clust_nostim_leftVright_v2_mask.nii ${surfaceatlas}/Conte69.L.midthickness.32k_fs_LR.surf.gii workbench_dir/Clust_nostim_leftVright_v2_mask_L.shape.gii -trilinear 
    wb_command -volume-to-surface-mapping Clust_nostim_leftVright_v2_mask.nii ${surfaceatlas}/Conte69.R.midthickness.32k_fs_LR.surf.gii workbench_dir/Clust_nostim_leftVright_v2_mask_R.shape.gii -trilinear 
    
    
fi


execute=0
if [ $execute -eq 1 ]; then
    echo "-Running T-Tests against 0 for each connectivity map (to run union versus intersection analysis)-"
    
    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    pushd $datadir/results/

    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
    -prefix stim_leftMT_ttest_against0_n10 -overwrite


    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
    -prefix nostim_leftMT_ttest_against0_n10 -overwrite

    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
    -prefix stim_rightMT_ttest_against0_n10 -overwrite


    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/069/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/141/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/172/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/173/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/177/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/178/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/083/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/144/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    ${datadir}/170/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
    -prefix nostim_rightMT_ttest_against0_n10 -overwrite

fi

# Now take intersection of no stimulated connectivity maps between left and right MT seeds
execute=0
if [ $execute -eq 1 ]; then

    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    pushd $datadir/results/
    
    3dcalc -a nostim_leftMT_sigclusters_n10_mask+tlrc -b nostim_rightMT_sigclusters_n10_mask+tlrc -expr 'a*b' -prefix nostim_lrMT_sigclusters_intersection_n10+tlrc
    3dcalc -a nostim_leftMT_sigclusters_n10_mask+tlrc -b nostim_rightMT_sigclusters_n10_mask+tlrc -expr 'ispositive(a+b)' -prefix nostim_lrMT_sigclusters_union_n10+tlrc

fi

# Now get average timeseries for each of the clusters in the intersection mask
execute=0
if [ $execute -eq 1 ]; then
    
    datadir=${basedir}/data/
    pushd $datadir/results/
    
    for subjNum in $listOfSubjects
    do
	subjDir=${basedir}/data/${subjNum}/
	subjfMRIDIR=${subjDir}/fMRI/
	subjAnalysisDIR=${subjfMRIDIR}/${ANALYSISNAME}Analysis
        pushd ${subjfMRIDIR}

        echo "Running 3dmaskave on subject ${subjNum}"
        for i in {1..6}
        do
            3dmaskave -quiet -mask ${datadir}/results/nostim_lrMT_sigclusters_intersection_n10_v2_mask+tlrc -mrange ${i} ${i} \
                ${subjfMRIDIR}/notACS_smInMask_resid2_tACSreg_stc_epi_short_allruns_tlrc_al+tlrc > ${subjAnalysisDIR}/cluster${i}_notACS_timeseries.1D
            3dmaskave -quiet -mask ${datadir}/results/nostim_lrMT_sigclusters_intersection_n10_v2_mask+tlrc -mrange ${i} ${i} \
                ${subjfMRIDIR}/yestACS_smInMask_resid2_tACSreg_stc_epi_short_allruns_tlrc_al+tlrc > ${subjAnalysisDIR}/cluster${i}_yestACS_timeseries.1D

        done
    done
fi


# Run ANOVA on union
execute=0
if [ $execute -eq 1 ]; then

    echo "-Running ANOVA-"
    
    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    pushd $datadir/results/

    3dANOVA3 -DAFNI_FLOATIZE=YES -type 4 -alevels 2 -blevels 2 -clevels 10 \
        -mask ${datadir}/results/nostim_lrMT_sigclusters_union_n10+tlrc \
        -dset 1 1 1 ${datadir}/038/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 1 ${datadir}/038/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 1 ${datadir}/038/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 1 ${datadir}/038/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 2 ${datadir}/069/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 2 ${datadir}/069/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 2 ${datadir}/069/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 2 ${datadir}/069/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 3 ${datadir}/141/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 3 ${datadir}/141/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 3 ${datadir}/141/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 3 ${datadir}/141/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 4 ${datadir}/172/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 4 ${datadir}/172/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 4 ${datadir}/172/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 4 ${datadir}/172/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 5 ${datadir}/173/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 5 ${datadir}/173/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 5 ${datadir}/173/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 5 ${datadir}/173/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 6 ${datadir}/177/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 6 ${datadir}/177/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 6 ${datadir}/177/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 6 ${datadir}/177/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 7 ${datadir}/178/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 7 ${datadir}/178/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 7 ${datadir}/178/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 7 ${datadir}/178/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 8 ${datadir}/083/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 8 ${datadir}/083/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 8 ${datadir}/083/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 8 ${datadir}/083/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 9 ${datadir}/144/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 9 ${datadir}/144/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 9 ${datadir}/144/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 9 ${datadir}/144/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 10 ${datadir}/170/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 10 ${datadir}/170/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 10 ${datadir}/170/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 10 ${datadir}/170/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -fa maineffect_stimulation -fb maineffect_leftright -fab interaction_stimLR \
        -amean 1 Stimulation -amean 2 NoStimulation -bmean 1 Left -bmean 2 Right \
        -adiff 1 2 diff_Stim_V_NoStim -bdiff 1 2 diff_Left_V_Right \
        -bucket ANOVA_stimulation_LR_subjN10_lrMTunion_v2 -overwrite

    # Try intersection data set
    3dANOVA3 -DAFNI_FLOATIZE=YES -type 4 -alevels 2 -blevels 2 -clevels 10 \
        -mask ${datadir}/results/nostim_lrMT_sigclusters_intersection_n10+tlrc \
        -dset 1 1 1 ${datadir}/038/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 1 ${datadir}/038/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 1 ${datadir}/038/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 1 ${datadir}/038/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 2 ${datadir}/069/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 2 ${datadir}/069/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 2 ${datadir}/069/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 2 ${datadir}/069/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 3 ${datadir}/141/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 3 ${datadir}/141/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 3 ${datadir}/141/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 3 ${datadir}/141/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 4 ${datadir}/172/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 4 ${datadir}/172/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 4 ${datadir}/172/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 4 ${datadir}/172/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 5 ${datadir}/173/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 5 ${datadir}/173/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 5 ${datadir}/173/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 5 ${datadir}/173/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 6 ${datadir}/177/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 6 ${datadir}/177/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 6 ${datadir}/177/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 6 ${datadir}/177/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 7 ${datadir}/178/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 7 ${datadir}/178/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 7 ${datadir}/178/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 7 ${datadir}/178/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 8 ${datadir}/083/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 8 ${datadir}/083/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 8 ${datadir}/083/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 8 ${datadir}/083/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 9 ${datadir}/144/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 9 ${datadir}/144/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 9 ${datadir}/144/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 9 ${datadir}/144/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -dset 1 1 10 ${datadir}/170/${analdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc \
        -dset 1 2 10 ${datadir}/170/${analdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
        -dset 2 1 10 ${datadir}/170/${analdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
        -dset 2 2 10 ${datadir}/170/${analdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
        -fa maineffect_stimulation -fb maineffect_leftright -fab interaction_stimLR \
        -amean 1 Stimulation -amean 2 NoStimulation -bmean 1 Left -bmean 2 Right \
        -adiff 1 2 diff_Stim_V_NoStim -bdiff 1 2 diff_Left_V_Right \
        -bucket ANOVA_stimulation_LR_subjN10_lrMTintersection_v2 -overwrite
    
    
    3dClustSim -mask ANOVA_stimulation_LR_subjN10_lrMTintersection_v2+tlrc -fwhm 6 -prefix 3dclustsim_ANOVA_intersection
    popd
fi


#################################
## 11/6/15
## Taku Ito
## By request from Kohitij: 3-way ANOVA, 2 fixed-effects: stimulation X time (with lMT connectivity)
# Run ANOVA on stimulation v time
execute=0
if [ $execute -eq 1 ]; then

    for subjNum in $listOfSubjects
    do
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

        pushd ${subjAnalysisDIR}

        echo "-Running seed-based correlation using left area MT on both stim and nostim conditions X with runs 1 and 2 split up-"
        # Extract number of time points for each subjects runs
        pushd ${subjfMRIDIR}
        TRs1=`3dinfo -nv epi_short_r1+orig`
        TRs2=`3dinfo -nv epi_short_r2+orig`
        TRs3=`3dinfo -nv epi_short_r3+orig`
        TRs4=`3dinfo -nv epi_short_r4+orig`
        
        TRs1start=$(expr $TRs1 - 1) # Since AFNI counts from 0
        TRs3start=$(expr $TRs3 - 1)
        popd
        echo "$TRs1; $TRs1start"

        outdir=${basedir}/data/results/4WayANOVAMats_STIMxHEMIxTIMExSUBJ/STIMxTIMExSUBJANOVA_wholebrain/
        3dTcorr1D -pearson -overwrite -prefix ${outdir}/${subjNum}_nostim_lMT_corrbucket_run1_v1 ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc[0..$TRs1start] nostim_leftAreaMT_v2.1D{0..$TRs1start}
        3dTcorr1D -pearson -overwrite -prefix ${outdir}/${subjNum}_nostim_lMT_corrbucket_run2_v1 ${subjfMRIDIR}/notACS_${nextInputFilename}+tlrc[$TRs1..$] nostim_leftAreaMT_v2.1D{$TRs1..$}
        3dTcorr1D -pearson -overwrite -prefix ${outdir}/${subjNum}_stim_lMT_corrbucket_run1_v1 ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc[0..$TRs3start] stim_leftAreaMT_v2.1D{0..$TRs3start}
        3dTcorr1D -pearson -overwrite -prefix ${outdir}/${subjNum}_stim_lMT_corrbucket_run2_v1 ${subjfMRIDIR}/yestACS_${nextInputFilename}+tlrc[$TRs3..$] stim_leftAreaMT_v2.1D{$TRs3..$}

        echo "-Run FisherZ transform on data-"
        3dcalc -a ${outdir}/${subjNum}_nostim_lMT_corrbucket_run1_v1+tlrc -expr 'atanh(a)' -prefix ${outdir}/${subjNum}_fz_nostim_lMT_corrbucket_run1_v1 -overwrite
        3dcalc -a ${outdir}/${subjNum}_nostim_lMT_corrbucket_run2_v1+tlrc -expr 'atanh(a)' -prefix ${outdir}/${subjNum}_fz_nostim_lMT_corrbucket_run2_v1 -overwrite
        3dcalc -a ${outdir}/${subjNum}_stim_lMT_corrbucket_run1_v1+tlrc -expr 'atanh(a)' -prefix ${outdir}/${subjNum}_fz_stim_lMT_corrbucket_run1_v1 -overwrite
        3dcalc -a ${outdir}/${subjNum}_stim_lMT_corrbucket_run2_v1+tlrc -expr 'atanh(a)' -prefix ${outdir}/${subjNum}_fz_stim_lMT_corrbucket_run2_v1 -overwrite

        popd

    done
fi

execute=0
if [ $execute -eq 1 ]; then
    echo "-Running ANOVA-"
    
    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    outdir=${basedir}/data/results/4WayANOVAMats_STIMxHEMIxTIMExSUBJ/STIMxTIMExSUBJANOVA_wholebrain/
    pushd $outdir
    
    3dresample -orient LPI -rmode NN -overwrite -prefix Kohitij_adapted_voxels_resampled -inset Kohitij_adapted_voxels.nii

    3dANOVA3 -DAFNI_FLOATIZE=YES -type 4 -alevels 2 -blevels 2 -clevels 10 \
        -mask ${outdir}/Kohitij_adapted_voxels_qw_resampled+tlrc \
        -dset 1 1 1 ${outdir}/038_fz_stim_lMT_corrbucket_run1_v1+tlrc \
        -dset 1 2 1 ${outdir}/038_fz_stim_lMT_corrbucket_run2_v1+tlrc \
        -dset 2 1 1 ${outdir}/038_fz_nostim_lMT_corrbucket_run1_v1+tlrc \
        -dset 2 2 1 ${outdir}/038_fz_nostim_lMT_corrbucket_run2_v1+tlrc \
        -dset 1 1 2 ${outdir}/069_fz_stim_lMT_corrbucket_run1_v1+tlrc \
        -dset 1 2 2 ${outdir}/069_fz_stim_lMT_corrbucket_run2_v1+tlrc \
        -dset 2 1 2 ${outdir}/069_fz_nostim_lMT_corrbucket_run1_v1+tlrc \
        -dset 2 2 2 ${outdir}/069_fz_nostim_lMT_corrbucket_run2_v1+tlrc \
        -dset 1 1 3 ${outdir}/141_fz_stim_lMT_corrbucket_run1_v1+tlrc \
        -dset 1 2 3 ${outdir}/141_fz_stim_lMT_corrbucket_run2_v1+tlrc \
        -dset 2 1 3 ${outdir}/141_fz_nostim_lMT_corrbucket_run1_v1+tlrc \
        -dset 2 2 3 ${outdir}/141_fz_nostim_lMT_corrbucket_run2_v1+tlrc \
        -dset 1 1 4 ${outdir}/172_fz_stim_lMT_corrbucket_run1_v1+tlrc \
        -dset 1 2 4 ${outdir}/172_fz_stim_lMT_corrbucket_run2_v1+tlrc \
        -dset 2 1 4 ${outdir}/172_fz_nostim_lMT_corrbucket_run1_v1+tlrc \
        -dset 2 2 4 ${outdir}/172_fz_nostim_lMT_corrbucket_run2_v1+tlrc \
        -dset 1 1 5 ${outdir}/173_fz_stim_lMT_corrbucket_run1_v1+tlrc \
        -dset 1 2 5 ${outdir}/173_fz_stim_lMT_corrbucket_run2_v1+tlrc \
        -dset 2 1 5 ${outdir}/173_fz_nostim_lMT_corrbucket_run1_v1+tlrc \
        -dset 2 2 5 ${outdir}/173_fz_nostim_lMT_corrbucket_run2_v1+tlrc \
        -dset 1 1 6 ${outdir}/177_fz_stim_lMT_corrbucket_run1_v1+tlrc \
        -dset 1 2 6 ${outdir}/177_fz_stim_lMT_corrbucket_run2_v1+tlrc \
        -dset 2 1 6 ${outdir}/177_fz_nostim_lMT_corrbucket_run1_v1+tlrc \
        -dset 2 2 6 ${outdir}/177_fz_nostim_lMT_corrbucket_run2_v1+tlrc \
        -dset 1 1 7 ${outdir}/178_fz_stim_lMT_corrbucket_run1_v1+tlrc \
        -dset 1 2 7 ${outdir}/178_fz_stim_lMT_corrbucket_run2_v1+tlrc \
        -dset 2 1 7 ${outdir}/178_fz_nostim_lMT_corrbucket_run1_v1+tlrc \
        -dset 2 2 7 ${outdir}/178_fz_nostim_lMT_corrbucket_run2_v1+tlrc \
        -dset 1 1 8 ${outdir}/083_fz_stim_lMT_corrbucket_run1_v1+tlrc \
        -dset 1 2 8 ${outdir}/083_fz_stim_lMT_corrbucket_run2_v1+tlrc \
        -dset 2 1 8 ${outdir}/083_fz_nostim_lMT_corrbucket_run1_v1+tlrc \
        -dset 2 2 8 ${outdir}/083_fz_nostim_lMT_corrbucket_run2_v1+tlrc \
        -dset 1 1 9 ${outdir}/144_fz_stim_lMT_corrbucket_run1_v1+tlrc \
        -dset 1 2 9 ${outdir}/144_fz_stim_lMT_corrbucket_run2_v1+tlrc \
        -dset 2 1 9 ${outdir}/144_fz_nostim_lMT_corrbucket_run1_v1+tlrc \
        -dset 2 2 9 ${outdir}/144_fz_nostim_lMT_corrbucket_run2_v1+tlrc \
        -dset 1 1 10 ${outdir}/170_fz_stim_lMT_corrbucket_run1_v1+tlrc \
        -dset 1 2 10 ${outdir}/170_fz_stim_lMT_corrbucket_run2_v1+tlrc \
        -dset 2 1 10 ${outdir}/170_fz_nostim_lMT_corrbucket_run1_v1+tlrc \
        -dset 2 2 10 ${outdir}/170_fz_nostim_lMT_corrbucket_run2_v1+tlrc \
        -fa maineffect_stimulation -fb maineffect_run -fab interaction_stimXrun \
        -amean 1 Stimulation -amean 2 NoStimulation -bmean 1 Run1 -bmean 2 Run2 \
        -adiff 1 2 diff_Stim_V_NoStim -bdiff 1 2 diff_Run1_V_Run2 \
        -bucket ${outdir}/ANOVA_stimulationXrun_subjN10_lMTfc_masked_v2 -overwrite

fi

execute=0
if [ $execute -eq 1 ]; then
    
    outdir=${basedir}/data/results/4WayANOVAMats_STIMxHEMIxTIMExSUBJ/STIMxTIMExSUBJANOVA_wholebrain/
    # Analysis to plot out the different correlations of the significant cluster from the stimXtime interaction effect ANOVA
    for subjNum in $listOfSubjects
    do
        subjdir=${basedir}/data/${subjNum}/fMRI/
        pushd ${subjdir}
        TRs1=`3dinfo -nv ${subjdir}/epi_short_r1+orig`
        TRs2=`3dinfo -nv ${subjdir}/epi_short_r2+orig`
        TRs3=`3dinfo -nv ${subjdir}/epi_short_r3+orig`
        TRs4=`3dinfo -nv ${subjdir}/epi_short_r4+orig`
        
        TRs1start=$(expr $TRs1 - 1) # Since AFNI counts from 0
        TRs3start=$(expr $TRs3 - 1)
        popd
        echo "${subjNum}"
        echo "${TRs1}; ${TRs1start}; ${TRs3}; ${TRs3start}"
        pushd $outdir
        3dmaskave -quiet -mask Clust_stimXtime_interaction_clust_mask+tlrc ${subjdir}/notACS_smInMask_resid2_tACSreg_stc_epi_short_allruns_tlrc_al+tlrc[0..$TRs1start] > interactionEffect_sigClust_timeseries/${subjNum}_nostimXrun1_stimXtime_interactioneffect_clust.1D
        3dmaskave -quiet -mask Clust_stimXtime_interaction_clust_mask+tlrc ${subjdir}/notACS_smInMask_resid2_tACSreg_stc_epi_short_allruns_tlrc_al+tlrc[$TRs1..$] > interactionEffect_sigClust_timeseries/${subjNum}_nostimXrun2_stimXtime_interactioneffect_clust.1D
        3dmaskave -quiet -mask Clust_stimXtime_interaction_clust_mask+tlrc ${subjdir}/yestACS_smInMask_resid2_tACSreg_stc_epi_short_allruns_tlrc_al+tlrc[0..$TRs3start] > interactionEffect_sigClust_timeseries/${subjNum}_stimXrun1_stimXtime_interactioneffect_clust.1D
        3dmaskave -quiet -mask Clust_stimXtime_interaction_clust_mask+tlrc ${subjdir}/yestACS_smInMask_resid2_tACSreg_stc_epi_short_allruns_tlrc_al+tlrc[$TRs3..$] > interactionEffect_sigClust_timeseries/${subjNum}_stimXrun2_stimXtime_interactioneffect_clust.1D
        popd
 
    done
fi




# 1. Take the lMT stim v no-stim maps and subtract fz correlation values (for each subject)
# 2. Take the rMT stim v no-stim maps and subtract fz correlation values (for each subject).
# 3. Take the t-test of these differences lMT stimVnostim against rMT stimVnostim (across subjects).
# Hypothesis: if there is a significant difference with lMT > rMT, then it shows that the stimulation is 
# Preferentially affecting FC with lMT regions more than rMT regions
execute=0
if [ $execute -eq 1 ]; then

    for subjNum in $listOfSubjects
    do
	subjDir=${basedir}/data/${subjNum}/
	subjfMRIDIR=${subjDir}/fMRI/
	subjAnalysisDIR=${subjfMRIDIR}/${ANALYSISNAME}Analysis
        pushd ${subjfMRIDIR}
        3dcalc -a ${subjAnalysisDIR}/fz_stim_leftareaMT_corrbucket_v2+tlrc -b ${subjAnalysisDIR}/fz_nostim_leftareaMT_corrbucket_v2+tlrc \
            -expr 'a-b' -prefix ${subjAnalysisDIR}/fz_stim_minus_nostim_lMT_corrmap_v2

        3dcalc -a ${subjAnalysisDIR}/fz_stim_rightareaMT_corrbucket_v2+tlrc -b ${subjAnalysisDIR}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
            -expr 'a-b' -prefix ${subjAnalysisDIR}/fz_stim_minus_nostim_rMT_corrmap_v2

        3dcalc -a ${subjAnalysisDIR}/fz_stim_leftareaMT_corrbucket_v2+tlrc -b ${subjAnalysisDIR}/fz_stim_rightareaMT_corrbucket_v2+tlrc \
            -expr 'a-b' -prefix ${subjAnalysisDIR}/fz_lMT_minus_rMT_stim_corrmap_v2

        3dcalc -a ${subjAnalysisDIR}/fz_nostim_leftareaMT_corrbucket_v2+tlrc -b ${subjAnalysisDIR}/fz_nostim_rightareaMT_corrbucket_v2+tlrc \
            -expr 'a-b' -prefix ${subjAnalysisDIR}/fz_lMT_minus_rMT_nostim_corrmap_v2

    done
    
    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis

    # Now run ttest 
    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -setB ${datadir}/038/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -prefix ${datadir}/results/lVrMT_stim_minus_nostim_ttest_n10_v2 -paired -overwrite

    3dttest++ -DAFNI_FLOATIZE=YES -setA ${datadir}/038/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -setB ${datadir}/038/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -prefix ${datadir}/results/stimVnostim_lMT_minus_rMT_ttest_n10_v2 -paired -overwrite
fi

# Same analysis but on intersection mask
execute=0
if [ $execute -eq 1 ]; then
    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    # Now run ttest 
    3dttest++ -DAFNI_FLOATIZE=YES -mask ${datadir}/results/nostim_lrMT_sigclusters_intersection_n10_v2_mask+tlrc \
    -setA ${datadir}/038/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -setB ${datadir}/038/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -prefix ${datadir}/results/lVrMT_stim_minus_nostim_ttest_intersection_mask_n10_v2 -paired -overwrite

    3dttest++ -DAFNI_FLOATIZE=YES -mask ${datadir}/results/nostim_lrMT_sigclusters_union_n10+tlrc \
    -setA ${datadir}/038/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -setB ${datadir}/038/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -prefix ${datadir}/results/lVrMT_stim_minus_nostim_ttest_union_mask_n10_v2 -paired -overwrite
    
    
    3dttest++ -DAFNI_FLOATIZE=YES -mask ${datadir}/results/nostim_lrMT_sigclusters_intersection_n10+tlrc \
    -setA ${datadir}/038/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -setB ${datadir}/038/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -prefix ${datadir}/results/stimVnostim_lMT_minus_rMT_ttest_intersection_n10_v2 -paired -overwrite

    3dttest++ -DAFNI_FLOATIZE=YES -mask ${datadir}/results/nostim_lrMT_sigclusters_union_n10+tlrc \
    -setA ${datadir}/038/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -setB ${datadir}/038/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -prefix ${datadir}/results/stimVnostim_lMT_minus_rMT_ttest_union_n10_v2 -paired -overwrite

    3dttest++ -DAFNI_FLOATIZE=YES -mask ${datadir}/results/Kohitij_adapted_voxels_qw_resampled+tlrc \
    -setA ${datadir}/038/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -setB ${datadir}/038/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/069/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/141/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/172/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/173/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/177/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/178/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/083/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/144/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    ${datadir}/170/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -prefix ${datadir}/results/stimVnostim_lMT_minus_rMT_ttest_adaptedROIsKohitij_n10_v2 -paired -overwrite
fi


# Run Wilcoxon test on correlation difference maps
execute=0
if [ $execute -eq 1 ];then

    datadir=${basedir}/data/
    analdir=fMRI/tacs_motionadaptationAnalysis
    
    # Now run ttest 
    3dWilcoxon -dset 1 ${datadir}/038/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -dset 1 ${datadir}/069/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -dset 1 ${datadir}/141/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -dset 1 ${datadir}/172/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -dset 1 ${datadir}/173/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -dset 1 ${datadir}/177/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -dset 1 ${datadir}/178/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -dset 1 ${datadir}/083/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -dset 1 ${datadir}/144/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -dset 1 ${datadir}/170/${analdir}/fz_stim_minus_nostim_lMT_corrmap_v2+tlrc \
    -dset 2 ${datadir}/038/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -dset 2 ${datadir}/069/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -dset 2 ${datadir}/141/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -dset 2 ${datadir}/172/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -dset 2 ${datadir}/173/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -dset 2 ${datadir}/177/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -dset 2 ${datadir}/178/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -dset 2 ${datadir}/083/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -dset 2 ${datadir}/144/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -dset 2 ${datadir}/170/${analdir}/fz_stim_minus_nostim_rMT_corrmap_v2+tlrc \
    -out ${datadir}/results/lVrMT_stim_minus_nostim_wilcoxon_n10_v2 

    3dWilcoxon -dset 1 ${datadir}/038/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -dset 1 ${datadir}/069/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -dset 1 ${datadir}/141/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -dset 1 ${datadir}/172/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -dset 1 ${datadir}/173/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -dset 1 ${datadir}/177/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -dset 1 ${datadir}/178/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -dset 1 ${datadir}/083/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -dset 1 ${datadir}/144/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -dset 1 ${datadir}/170/${analdir}/fz_lMT_minus_rMT_stim_corrmap_v2+tlrc \
    -dset 2 ${datadir}/038/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -dset 2 ${datadir}/069/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -dset 2 ${datadir}/141/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -dset 2 ${datadir}/172/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -dset 2 ${datadir}/173/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -dset 2 ${datadir}/177/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -dset 2 ${datadir}/178/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -dset 2 ${datadir}/083/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -dset 2 ${datadir}/144/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -dset 2 ${datadir}/170/${analdir}/fz_lMT_minus_rMT_nostim_corrmap_v2+tlrc \
    -out ${datadir}/results/stimVnostim_lMT_minus_rMT_wilcoxon_n10_v2 
fi

# Taku Ito; 11/19/15
# Extract correlation values for all subjects 
# for lVrMT_stim_minus_nostim_ttest cluster (insula)
# also for left stim V no stim whole brain t-test
execute=0
if [ $execute -eq 1 ]; then
    
    outdir=${basedir}/data/results/
    
    # Remove files if they previously exist
    pushd $outdir/lVrMT_stim_minus_nostim
    rm -v lvrMT_stim_minus_nostim_corr_nostim_lMT.txt \
        lvrMT_stim_minus_nostim_corr_nostim_rMT.txt \
        lvrMT_stim_minus_nostim_corr_stim_lMT.txt \
        lvrMT_stim_minus_nostim_corr_stim_rMT.txt
    popd
    pushd $outdir/l_stimVnoStim_ttest_v2
    rm -v lMT_stimVnostim_corr_nostim_lMT_clust*.txt \
        lMT_stimVnostim_corr_nostim_rMT_clust*.txt \
        lMT_stimVnostim_corr_stim_lMT_clust*.txt \
        lMT_stimVnostim_corr_stim_rMT_clust*.txt
    popd

    for subjNum in $listOfSubjects
    do
        subjdir=${basedir}/data/${subjNum}/fMRI/tacs_motionadaptationAnalysis/
        
        echo "${subjNum}"
        pushd $outdir/lVrMT_stim_minus_nostim
        3dmaskave -quiet -mask Clust_lVrMT_stim_minus_nostim_ttest_union_mask+tlrc ${subjdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc >> lvrMT_stim_minus_nostim_corr_nostim_lMT.txt
        3dmaskave -quiet -mask Clust_lVrMT_stim_minus_nostim_ttest_union_mask+tlrc ${subjdir}/fz_nostim_rightareaMT_corrbucket_v2+tlrc >> lvrMT_stim_minus_nostim_corr_nostim_rMT.txt
        3dmaskave -quiet -mask Clust_lVrMT_stim_minus_nostim_ttest_union_mask+tlrc ${subjdir}/fz_stim_leftareaMT_corrbucket_v2+tlrc >> lvrMT_stim_minus_nostim_corr_stim_lMT.txt
        3dmaskave -quiet -mask Clust_lVrMT_stim_minus_nostim_ttest_union_mask+tlrc ${subjdir}/fz_stim_rightareaMT_corrbucket_v2+tlrc >> lvrMT_stim_minus_nostim_corr_stim_rMT.txt
        popd

        pushd $outdir/l_stimVnoStim_ttest_v2
        for i in {1..9}
        do
            3dmaskave -quiet -mask ../Clust_l_stimVnoStim_v2_mask+tlrc -mrange ${i} ${i} \
                ${subjdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc >> lMT_stimVnostim_corr_nostim_lMT_clust${i}.txt
            3dmaskave -quiet -mask ../Clust_l_stimVnoStim_v2_mask+tlrc -mrange ${i} ${i} \
                ${subjdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc >> lMT_stimVnostim_corr_nostim_rMT_clust${i}.txt
            3dmaskave -quiet -mask ../Clust_l_stimVnoStim_v2_mask+tlrc -mrange ${i} ${i} \
                ${subjdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc >> lMT_stimVnostim_corr_stim_lMT_clust${i}.txt
            3dmaskave -quiet -mask ../Clust_l_stimVnoStim_v2_mask+tlrc -mrange ${i} ${i} \
                ${subjdir}/fz_nostim_leftareaMT_corrbucket_v2+tlrc >> lMT_stimVnositm_corr_stim_rMT_clust${i}.txt
        done
        popd
    done
fi



execute=0
if [ $execute -eq 1 ]; then
    
    for subj in $listOfSubjects
    do
        basedir=/projects/Collaborations/KrekelbergCollaboration/MotionAdaptation_tACS_FC/
        subjfMRIDIR=${basedir}/data/${subj}/fMRI
        3dmaskave -quiet -mask ${basedir}/data/areaMT_Left_5mmrad+tlrc ${subjfMRIDIR}/${nextInputFilename}+tlrc > ${subjfMRIDIR}/leftAreaMT_5mmrad_v2.1D
        3dmaskave -quiet -mask ${basedir}/data/areaMT_Right_5mmrad+tlrc ${subjfMRIDIR}/${nextInputFilename}+tlrc > ${subjfMRIDIR}/rightAreaMT_5mmrad_v2.1D
    done
fi
