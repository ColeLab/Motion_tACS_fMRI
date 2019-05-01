function output = loadTaskVolumeData(subj)
% Taku Ito
% 3/11/16
%
% This script loads in the volume data using the Power et al. 2011 264 ROI volume segmentation scheme
% 
% Parameters: subj ( must be input with single quotations, i.e., as a string)

    %% First run for entire time series
    % Get data directory based 
    basedir = ['/projects/Collaborations/KrekelbergCollaboration/MotionAdaptation_tACS_FC/data/'];
    fmridir = [basedir '/' subj '/fMRI/'];
    taskdata = [fmridir '/smInMask_resid2_tACSreg_stc_epi_short_allruns_tlrc_al.nii'];
    numROIs = 264;

    % Each roi timeseries matrix is 264x581, so we will want to create an empty matrix of 264x4648 (8*582 = 4648)
    disp('Loading Petersen lab 264 ROI set')
    psetHDR3D = load_nifti('/projects/AnalysisTools/atlases/PetersenLab264_MNI_WithLabels333.nii');
    % Reshape to 1D
    psetvol = reshape(psetHDR3D.vol,[size(psetHDR3D.vol,1)*size(psetHDR3D.vol,2)*size(psetHDR3D.vol,3) 1]);

    %%%%%%% Task data
    disp(['Loading task data for subj ' num2str(subj) '...'])
    % Import subject's task data
    dset4D = load_nifti(taskdata);
    % Reshape to 2D (space x time)
    dsetvol_task = reshape(dset4D.vol,[size(dset4D.vol,1)*size(dset4D.vol,2)*size(dset4D.vol,3) size(dset4D.vol,4)]);
    clear dset4D
    
    % Get ROI task data
    ROIData_task = zeros(numROIs, size(dsetvol_task,2));
    for regionNum=1:numROIs
        % Get region's data
        ROIData_task(regionNum,:) = mean(dsetvol_task(psetvol==regionNum,:),1);
    end
    clear dsetvol_task
   
    leftMT5mm = [fmridir '/leftAreaMT_5mmrad_v2.1D'];
    rightMT5mm = [fmridir '/rightAreaMT_5mmrad_v2.1D'];
    
    lmt = importdata(leftMT5mm)';
    rmt = importdata(rightMT5mm)';
    
    ROIData_task = [ROIData_task; lmt; rmt];
    output.dtseries_task = ROIData_task;
    
    %% Run for tACS On and tACS Off timeseries separately
        %% First run for entire time series
    % Get data directory based 
    basedir = ['/projects/Collaborations/KrekelbergCollaboration/MotionAdaptation_tACS_FC/data/'];
    fmridir = [basedir '/' subj '/fMRI/'];
    taskdata1 = [fmridir '/notACS_smInMask_resid2_tACSreg_stc_epi_short_allruns_tlrc_al.nii'];
    taskdata2 = [fmridir '/yestACS_smInMask_resid2_tACSreg_stc_epi_short_allruns_tlrc_al.nii'];

    numROIs = 264;

    % Each roi timeseries matrix is 264x581, so we will want to create an empty matrix of 264x4648 (8*582 = 4648)
    disp('Loading Petersen lab 264 ROI set')
    psetHDR3D = load_nifti('/projects/AnalysisTools/atlases/PetersenLab264_MNI_WithLabels333.nii');
    % Reshape to 1D
    psetvol = reshape(psetHDR3D.vol,[size(psetHDR3D.vol,1)*size(psetHDR3D.vol,2)*size(psetHDR3D.vol,3) 1]);

    %% For tACS Off
    %%%%%%% Task data
    disp(['Loading task data for subj ' num2str(subj) '...'])
    % Import subject's task data
    dset4D = load_nifti(taskdata1);
    % Reshape to 2D (space x time)
    dsetvol_task = reshape(dset4D.vol,[size(dset4D.vol,1)*size(dset4D.vol,2)*size(dset4D.vol,3) size(dset4D.vol,4)]);
    clear dset4D
    
    % Get ROI task data
    ROIData_task1 = zeros(numROIs, size(dsetvol_task,2));
    for regionNum=1:numROIs
        % Get region's data
        ROIData_task1(regionNum,:) = mean(dsetvol_task(psetvol==regionNum,:),1);
    end
    clear dsetvol_task
   
    %% For tACS On
    % Import subject's task data
    dset4D = load_nifti(taskdata2);
    % Reshape to 2D (space x time)
    dsetvol_task = reshape(dset4D.vol,[size(dset4D.vol,1)*size(dset4D.vol,2)*size(dset4D.vol,3) size(dset4D.vol,4)]);
    clear dset4D
    
    % Get ROI task data
    ROIData_task2 = zeros(numROIs, size(dsetvol_task,2));
    for regionNum=1:numROIs
        % Get region's data
        ROIData_task2(regionNum,:) = mean(dsetvol_task(psetvol==regionNum,:),1);
    end
    clear dsetvol_task
    
    %% split up MT timeseries accordingly
    leftMT5mm = [fmridir '/leftAreaMT_5mmrad_v2.1D'];
    rightMT5mm = [fmridir '/rightAreaMT_5mmrad_v2.1D'];
    
    lmt = importdata(leftMT5mm)';
    rmt = importdata(rightMT5mm)';
    
    tacs_off_ind = size(ROIData_task1,2);
    lmt_off = lmt(1:tacs_off_ind);
    rmt_off = rmt(1:tacs_off_ind);
    lmt_on = lmt(tacs_off_ind+1:end);
    rmt_on = rmt(tacs_off_ind+1:end);
    
    ROIData_task1 = [ROIData_task1; lmt_off; rmt_off];
    ROIData_task2 = [ROIData_task2; lmt_on; rmt_on];
    output.dtseries_tacs_off = ROIData_task1;
    output.dtseries_tacs_on = ROIData_task2;
    
end
