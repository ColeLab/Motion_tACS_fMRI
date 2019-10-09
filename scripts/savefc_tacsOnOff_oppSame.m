% Taku Ito
% 10/09/2019
% Save out FC data for Bart (resubmission, JNeurophys)

%% Load in data
N = maxNumCompThreads(10);
basedir = '/projects/Collaborations/KrekelbergCollaboration/MotionAdaptation_tACS_FC/';
datadir = [basedir 'data/results/PowerAnalyses/'];
subjNums = {'038', '069', '083', '141', '144', '170', '172', '173', ...
            '177', '178'};

% Load in Power Network Affiliation
power11 = importdata('Power11NewCoorOrder.txt');
power11 = power11(1:264,5);
power11 = [power11; -1; -1]; % Add hMTs to uncertain networks
%power11 = power11(:,5);
% Load data
subj_data = {};
subjcount = 1;
tacsontimes = {};
tacsofftimes = {};
opptimes = {};
sametimes = {};
for subj=subjNums
    subjfile = [datadir subj{1} '_266Power_timeseries.csv'];

    % Import entire timeseries for subject
    subj_data{subjcount} = zscore(csvread(subjfile),1,2);
    
    % Import the length of tacs off runs (2 runs total)
    tacsonsetfile = [datadir subj{1} '_266Power_timeseries_tacsOff.csv'];
    tmp = csvread(tacsonsetfile);
    tacsoffset = size(tmp,2); % Obtain length of timeseries only
    % Generate stimulus arrays
    tacsofftimes{subjcount} = zeros(size(subj_data{subjcount},2),1);
    tacsontimes{subjcount} = zeros(size(subj_data{subjcount},2),1);
    tacsofftimes{subjcount}(1:tacsoffset) = 1;
    tacsontimes{subjcount}((tacsoffset+1):end) = 1;
    clear tmp
    
    % Import convolved stim files for same condition (adapted condition)
    stimdir = [basedir 'data/' subj{1} '/sdm/'];
    same_tmp = importdata([stimdir 'SameTimes.1D']);
    opp_tmp = importdata([stimdir 'OppTimes.1D']);
    % Binarize convolved time series to obtain only portions in which HRF activation > 0.5
    sametimes{subjcount} = same_tmp > 0;
    opptimes{subjcount} = opp_tmp > 0;
    
    subjcount = subjcount + 1;
end

%% hMT+ GBC Analysis
% Organize matrices
gbc_vectors = zeros(8,10);
fc_tacs_on_opp = zeros(266,266,length(subjNums));
fc_tacs_off_opp = zeros(266,266,length(subjNums));
fc_tacs_on_same = zeros(266,266,length(subjNums));
fc_tacs_off_same = zeros(266,266,length(subjNums));

gbc_tacs_on_opp = zeros(266,length(subjNums));
gbc_tacs_off_opp = zeros(266,length(subjNums));
gbc_tacs_on_same = zeros(266,length(subjNums));
gbc_tacs_off_same = zeros(266,length(subjNums));

for subj=1:length(subjNums)
    
    % Exclude region 257 and 262 from all analyses and set timeseries to NaN
    % Overlaps with rMT ROI mask by 1 voxel
    subj_data{subj}(257,:) = nan;
    % Adjacent voxels with lMT ROI mask
    subj_data{subj}(262,:) = nan; 
    
    % tACS OFF | SAME condition
    tacsoff_same = sametimes{subj}.*tacsofftimes{subj};
    ind = find(tacsoff_same);
    fc = corrcoef(subj_data{subj}(:,ind)');
    fc_pos = fc > 0;
    fc = fc_pos.*fc;
    fc = atanh(fc);
    nandiag = diag(nan(266,1));
    fc = fc + nandiag;
    fc_tacs_off_same(:,:,subj) = fc;
    gbc_tacs_off_same(:,subj) = nanmean(fc,2);
    
    % tACS OFF | OPP condition
    tacsoff_opp = opptimes{subj}.*tacsofftimes{subj};
    ind = find(tacsoff_opp);
    fc = corrcoef(subj_data{subj}(:,ind)');
    fc_pos = fc > 0;
    fc = fc_pos.*fc;
    fc = atanh(fc);
    nandiag = diag(nan(266,1));
    fc = fc + nandiag;
    fc_tacs_off_opp(:,:,subj) = fc;
    gbc_tacs_off_opp(:,subj) = nanmean(fc,2);
    
    % tACS ON | SAME
    tacson_same = sametimes{subj}.*tacsontimes{subj};
    ind = find(tacson_same);
    fc = corrcoef(subj_data{subj}(:,ind)');
    fc_pos = fc > 0;
    fc = fc_pos.*fc;
    fc = atanh(fc);
    nandiag = diag(nan(266,1));
    fc = fc + nandiag;
    fc_tacs_on_same(:,:,subj) = fc;
    gbc_tacs_on_same(:,subj) = nanmean(fc,2);
    
    % tACS ON | OPP condition
    tacson_opp = opptimes{subj}.*tacsontimes{subj};
    ind = find(tacson_opp);
    fc = corrcoef(subj_data{subj}(:,ind)');
    fc_pos = fc > 0;
    fc = fc_pos.*fc;
    fc = atanh(fc);
    nandiag = diag(nan(266,1));
    fc = fc + nandiag;
    fc_tacs_on_opp(:,:,subj) = fc;
    gbc_tacs_on_opp(:,subj) = nanmean(fc,2);
end

% Save relvant variables
save('savefc_tacsOnOff_oppSame.mat', 'fc_tacs_on_opp', 'fc_tacs_on_same', 'fc_tacs_off_opp', 'fc_tacs_off_same', 'power11')
