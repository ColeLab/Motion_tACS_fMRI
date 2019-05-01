% Taku Ito
% 4/11/16
% tACS Motion Adaptation Project

%% Subject parameters

subjNums = {'038', '069', '083', '141', '144', '170', '172', '173', ...
            '177', '178'};

%% Load in Subject data
subjcount = 1;
for subj=subjNums
    out = loadTaskVolumeData(subj{1});
    subj_mat{subjcount} = out.dtseries_task;
    subj_mat1{subjcount} = out.dtseries_tacs_off;
    subj_mat2{subjcount} = out.dtseries_tacs_on;
    subjcount = subjcount + 1;
end

%% Write subject data to CSV
outdir = '/projects/Collaborations/KrekelbergCollaboration/MotionAdaptation_tACS_FC/data/results/PowerAnalyses/';

subjcount = 1;
for subj=subjNums
    filename = [outdir subj{1} '_266Power_timeseries.csv'];
    csvwrite(filename, subj_mat{subjcount});
    % For tACS ON and OFF
    filename1 = [outdir subj{1} '_266Power_timeseries_tacsOff.csv'];
    csvwrite(filename1, subj_mat1{subjcount});
    filename2 = [outdir subj{1} '_266Power_timeseries_tacsOn.csv'];
    csvwrite(filename2, subj_mat2{subjcount});
    subjcount = subjcount + 1;
end