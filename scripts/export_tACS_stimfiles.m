% Taku Ito
% 6/2/2015

% tACS Motion Adaptation project

% Takes in MATLAB stimulus files given by Kohitij and outputs them as txt files to be fed into AFNI's 3dDeconvolve

subjNums = {'038', '069', '083', '141', '144', '170', '172', '173', '177', '178'};
basedir = '/projects/Collaborations/KrekelbergCollaboration/MotionAdaptation_tACS_FC/data/';
inputdir = '/projects/Collaborations/KrekelbergCollaboration/MotionAdaptation_tACS_FC/data/stimDesign/';

for i=1:length(subjNums)
    disp(['Running subject number ' subjNums{i}])
    sdmdir = [basedir subjNums{i} '/sdm/'];
    inputMat1 = [inputdir 'design_' subjNums{i} '_MAE_fMRI_tES_stim_1'];
    inputMat2 = [inputdir 'design_' subjNums{i} '_MAE_fMRI_tES_stim_2'];
    
    numTRsAllRuns = importdata([sdmdir 'allTaskTimings.txt']);
    numTRsAllRuns = size(numTRsAllRuns,1);
    
    % Write first tACS stimulus run
    load(inputMat1)
    % Make sure correct number of TRs in each run
    tmp = importdata([sdmdir 'RunOnsets.txt']);
    tmp = tmp + 1;
    numTRsRun3 = tmp(4)-tmp(3) + 1;
    if numTRsRun3 > size(design,1)
        disp('resizing')
        design3 = zeros(numTRsRun3, 1);
        design3(1:size(design,1),1) = design;
        design = design3;
    else
        design3 = design(1:size(numTRsRun3,1),1);
        design = design3;
    end
    outfile1 = [basedir subjNums{i} '/sdm/tACS_timingfiles_stim1.txt'];
    dlmwrite(outfile1, design, 'delimiter', '\n');
    clear design

    % Write second tACS stimulus run
    load(inputMat2)
    % Make sure correct number of TRs in each run
    numTRsRun4 = numTRsAllRuns - tmp(4) - 2;
    if numTRsRun4 > size(design,1)
        design4 = zeros(numTRsRun4, 1);
        design4(1:size(design,1),1) = design;
        design = design4;
    else
        design4 = design(1:size(numTRsRun4,1),1);
        design = design4;
    end
    outfile2 = [basedir subjNums{i} '/sdm/tACS_timingfiles_stim2.txt'];
    dlmwrite(outfile2, design, 'delimiter', '\n');
    clear design

    % Now create a single array for all 4 runs (padding 0s for the first 2 runs)
    numTRsRun1 = tmp(2) -tmp(1) + 1; 

    numTRsRun2 = tmp(3) - tmp(2) + 1;

    first2runs = numTRsRun1 + numTRsRun2;


    % Count TR index for 3rd and 4th runs (i.e., stimulation runs)
    
    % 3rd Run tACS stimulus times
    trCount = 1;
    tACS_out = zeros(numTRsAllRuns,1);
    % 3rd run (1st stim run)
    trCount = trCount + first2runs;
    tACS_out(trCount:(trCount+numTRsRun3-1),1) = importdata(outfile1);
    
    % 4th Run tACS stimulus times
    trCount = trCount + numTRsRun3;
    tACS_out(trCount:end,1) = importdata(outfile2);

    % Write out concatenated tACS stimulus timing file
    dlmwrite([sdmdir 'all_tACS_binary_timings.txt'], tACS_out, 'delimiter', '\n');
    

end
