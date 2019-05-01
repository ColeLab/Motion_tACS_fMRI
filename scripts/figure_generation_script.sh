# Taku Ito
# 05/12/16

# Figure generation of surface visualizations of tACS findings
# Use workbench to generate figures

basedir=/projects/Collaborations/KrekelbergCollaboration/MotionAdaptation_tACS_FC/
datadir=${basedir}/data/

# Generate hMT+ to DAN regions surface visualization
execute=1
if [ $execute -eq 1 ]; then
    regions="251 252 256 258 259 260 261 263 264" # Region 257 and 262 is excluded since it overlaps with right hMT+
    
    pushd ${basedir}/data
    # Include hMT+s into single NIFTI
    3dcalc -overwrite -a areaMT_Right_5mmrad+tlrc -b areaMT_Left_5mmrad+tlrc -expr "(a+b)*-1" -prefix DAN_hMT_Regions.nii.gz

    for roi in $regions
    do
        3dcalc -overwrite -a PetersenLab264_MNI_WithLabels333.nii -b DAN_hMT_Regions.nii.gz -expr "equals(a,${roi})*100+b" -prefix DAN_hMT_Regions.nii.gz 
    done
    popd
fi

# Generate entire Power atlas for surface visualizations, highlighting DAN regions with certain number
execute=1
if [ $execute -eq 1 ]; then
    regions="251 252 256 258 259 260 261 263 264" # Region 257 and 262 is excluded since it overlaps with right hMT+
    
    pushd ${basedir}/data
    # Include hMT+s into single NIFTI
    3dcalc -overwrite -a areaMT_Right_5mmrad+tlrc -b areaMT_Left_5mmrad+tlrc -expr "(a+b)*-1" -prefix DAN_hMT_AllRegions.nii.gz
    
    for roi in {1..264};
    do
        if [ $roi -ne 262 ]; then
            if [ $roi -ne 257 ]; then
                3dcalc -overwrite -a PetersenLab264_MNI_WithLabels333.nii -b DAN_hMT_AllRegions.nii.gz -expr "equals(a,${roi})+b" -prefix DAN_hMT_AllRegions.nii.gz 
            fi
        fi
    done

    # Numerically label DAN regions
    for roi in $regions
    do
        3dcalc -overwrite -a DAN_hMT_AllRegions.nii.gz -expr "a+equals(a,${roi})*100" -prefix DAN_hMT_AllRegions.nii.gz
    done
fi

# Generate specific ROIs with p < 0.01
execute=1
if [ $execute -eq 1 ]; then
    regions="130 147 258"
    
    pushd ${basedir}/data
    # Include hMT+s into single NIFTI
    3dcalc -overwrite -a areaMT_Right_5mmrad+tlrc -b areaMT_Left_5mmrad+tlrc -expr "(a+b)*-1" -prefix hMT_IndividualRegions.nii.gz

    for roi in $regions
    do
        3dcalc -overwrite -a PetersenLab264_MNI_WithLabels333.nii -b hMT_IndividualRegions.nii.gz -expr "equals(a,${roi})*100+b" -prefix hMT_IndividualRegions.nii.gz 
    done
    popd
fi


# Volume to Surface Mapping of Results
execute=1
if [ $execute -eq 1 ]; then

    pushd $datadir
    surfaceatlas=/usr/local/workbench/atlases/Conte69_atlas-v2.LR.32k_fs_LR.wb/32k_ConteAtlas_v2/
    #surfaceatlas=/projects/IndivRITL/data/Parcels/
    #surfaceatlas=/usr/local/workbench/atlases/Conte69_atlas_164k_wb/
    #wb_command -volume-to-surface-mapping DAN_hMT_Regions.nii.gz ${surfaceatlas}/Conte69.R.midthickness.164k_fs_LR.surf.gii DAN_hMT_Regions_R.shape.gii -trilinear 
    #wb_command -volume-to-surface-mapping DAN_hMT_Regions.nii.gz ${surfaceatlas}/Conte69.L.midthickness.164k_fs_LR.surf.gii DAN_hMT_Regions_L.shape.gii -trilinear 
    # LR GIFTI generations of hMT+ to DAN interaction regions                                                                                                   
    wb_command -volume-to-surface-mapping DAN_hMT_Regions.nii.gz ${surfaceatlas}/Conte69.L.midthickness.32k_fs_LR.surf.gii DAN_hMT_Regions_L.shape.gii -trilinear  
    wb_command -volume-to-surface-mapping DAN_hMT_Regions.nii.gz ${surfaceatlas}/Conte69.R.midthickness.32k_fs_LR.surf.gii DAN_hMT_Regions_R.shape.gii -trilinear

    # Generate surface for all Power regions
    wb_command -volume-to-surface-mapping DAN_hMT_AllRegions.nii.gz ${surfaceatlas}/Conte69.L.midthickness.32k_fs_LR.surf.gii DAN_hMT_AllRegions_L.shape.gii -trilinear  
    wb_command -volume-to-surface-mapping DAN_hMT_AllRegions.nii.gz ${surfaceatlas}/Conte69.R.midthickness.32k_fs_LR.surf.gii DAN_hMT_AllRegions_R.shape.gii -trilinear
    
    # Generate surface for individual regions
    wb_command -volume-to-surface-mapping hMT_IndividualRegions.nii.gz ${surfaceatlas}/Conte69.L.midthickness.32k_fs_LR.surf.gii hMT_IndividualRegions_L.shape.gii -trilinear  
    wb_command -volume-to-surface-mapping hMT_IndividualRegions.nii.gz ${surfaceatlas}/Conte69.R.midthickness.32k_fs_LR.surf.gii hMT_IndividualRegions_R.shape.gii -trilinear
    popd
fi
