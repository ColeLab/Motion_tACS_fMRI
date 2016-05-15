% Taku Ito
% 4/13/16
% GBC Analysis:
% 4-way interaction ANOVA analysis (3 fixed effects, 1 random effects)
% Fixed effects:
% 1. Adaptation (Opp v. Same)
% 2. Hemisphere (left v. right MT)
% 3. tACS (On v. Off)
% Random effects:
% 1. Subjects

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
    gbc_vectors(1,subj) = (nanmean(fc(265,:))); % left MT connectivity
    gbc_vectors(5,subj) = (nanmean(fc(266,:))); % right MT connectivity
    
    % tACS OFF | OPP condition
    tacsoff_opp = opptimes{subj}.*tacsofftimes{subj};
    ind = find(tacsoff_opp);
    fc = corrcoef(subj_data{subj}(:,ind)');
    fc_pos = fc > 0;
    fc = fc_pos.*fc;
    fc = atanh(fc);
    nandiag = diag(nan(266,1));
    fc = fc + nandiag;
    gbc_vectors(2,subj) = (nanmean(fc(265,:))); % left MT connectivity
    gbc_vectors(6, subj) = (nanmean(fc(266,:))); % right MT connectivity
    
    % tACS ON | SAME
    tacson_same = sametimes{subj}.*tacsontimes{subj};
    ind = find(tacson_same);
    fc = corrcoef(subj_data{subj}(:,ind)');
    fc_pos = fc > 0;
    fc = fc_pos.*fc;
    fc = atanh(fc);
    nandiag = diag(nan(266,1));
    fc = fc + nandiag;
    gbc_vectors(3,subj) = (nanmean(fc(265,:))); % left MT connectivity
    gbc_vectors(7,subj) = (nanmean(fc(266,:))); % right MT connectivity
    
    % tACS ON | OPP condition
    tacson_opp = opptimes{subj}.*tacsontimes{subj};
    ind = find(tacson_opp);
    fc = corrcoef(subj_data{subj}(:,ind)');
    fc_pos = fc > 0;
    fc = fc_pos.*fc;
    fc = atanh(fc);
    nandiag = diag(nan(266,1));
    fc = fc + nandiag;
    gbc_vectors(4,subj) = (nanmean(fc(265,:))); % left MT connectivity
    gbc_vectors(8,subj) = (nanmean(fc(266,:)));
end

% Run ANOVA
% Define group names
adapt = {'same';'opp';'same';'opp';'same';'opp';'same';'opp';};
adapt = repmat(adapt,10,1);
tacs = {'tacs_off';'tacs_off'; 'tacs_on'; 'tacs_on';'tacs_off';'tacs_off'; 'tacs_on'; 'tacs_on';};
tacs = repmat(tacs,10,1);
hemi = {'left';'left';'left';'left';'right';'right';'right';'right'};
hemi = repmat(hemi,10,1);
subjs = reshape(repmat(1:10',8,1),[80,1]);

indata = reshape(gbc_vectors,[80,1]);
[p,tbl] = anovan(indata,{adapt,tacs,hemi,subjs}, 'model', 'full','random',[4], 'varnames', {'adapt','tacs','hemi','subjs'});
% close(1)
gbc_interactionF = tbl{12,6};
gbc_interactionP = tbl{12,7};

% Print out results
disp('Interaction of GBC of MTs')
disp(['F = ' num2str(gbc_interactionF)])
disp(['p = ' num2str(gbc_interactionP)])

% Plot effect

figure
lmt_same = [mean(gbc_vectors(1,:)), mean(gbc_vectors(3,:))];
rmt_same = [mean(gbc_vectors(5,:)), mean(gbc_vectors(7,:))];
lmt_opp = [mean(gbc_vectors(2,:)), mean(gbc_vectors(4,:))];
rmt_opp = [mean(gbc_vectors(6,:)), mean(gbc_vectors(8,:))];
numsubjs = 10;
e_lmt_same = [std(gbc_vectors(1,:))/numsubjs, std(gbc_vectors(3,:))/numsubjs];
e_rmt_same = [std(gbc_vectors(5,:))/numsubjs, std(gbc_vectors(7,:))/numsubjs];
e_lmt_opp = [std(gbc_vectors(2,:))/numsubjs, std(gbc_vectors(4,:))/numsubjs];
e_rmt_opp = [std(gbc_vectors(6,:))/numsubjs, std(gbc_vectors(8,:))/numsubjs];

hold on
title('Weighted Degree Centrality for hMT+ to whole-brain')
errorbar(lmt_same,e_lmt_same,'b')
errorbar(lmt_opp,e_lmt_opp, 'c')
errorbar(rmt_same,e_rmt_same,'r')
errorbar(rmt_opp,e_rmt_opp, 'm')
xlim([0.9,2.1])
ylabel('Average Functional Connectivity')
set(gca, 'XTick', [1.0, 2.0])
set(gca, 'XTickLabel', {'tACS Off', 'tACS On'})
legend('left hMT+ adapt', 'left hMT+ non-adapt', 'right hMT+ adapt', 'right hMT+ non-adapt', 'Location', 'northwest')

%% Now run ANOVA using GBC within each of the different Power networks
% Organize matrices
gbc_vectors_network = {};
networks = unique(power11);
for netid=1:length(networks)
    netid_num = networks(netid);
    
    gbc_vectors_network{netid} = zeros(8,10);
    % Obtain indices for this network
    net_ind = power11(:)==netid_num;
%     nregions_net = sum(net_ind);
    
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
%         fc_pos = fc > 0;
        fc_pos = ones(266,1);
        ind265 = net_ind.*fc_pos(265,:)';
        ind266 = net_ind.*fc_pos(266,:)';
        fc = fc.*(fc>0);
        fc = atanh(fc);
        nandiag = diag(nan(266,1));
        fc = fc + nandiag;
        gbc_vectors_network{netid}(1,subj) = (nanmean(fc(265,ind265>0))); % left MT connectivity
        gbc_vectors_network{netid}(5,subj) = (nanmean(fc(266,ind266>0))); % right MT connectivity
        
        % tACS OFF | OPP condition
        tacsoff_opp = opptimes{subj}.*tacsofftimes{subj};
        ind = find(tacsoff_opp);
        fc = corrcoef(subj_data{subj}(:,ind)');
%         fc_pos = fc > 0;
        fc_pos = ones(266,1);
        ind265 = net_ind.*fc_pos(265,:)';
        ind266 = net_ind.*fc_pos(266,:)';
        fc = fc.*(fc>0);
        fc = atanh(fc);
        nandiag = diag(nan(266,1));
        fc = fc + nandiag;
        gbc_vectors_network{netid}(2,subj) = (nanmean(fc(265,ind265>0))); % left MT connectivity
        gbc_vectors_network{netid}(6, subj) = (nanmean(fc(266,ind266>0))); % right MT connectivity
        
        % tACS ON | SAME
        tacson_same = sametimes{subj}.*tacsontimes{subj};
        ind = find(tacson_same);
        fc = corrcoef(subj_data{subj}(:,ind)');
%         fc_pos = fc > 0;
        fc_pos = ones(266,1);
        ind265 = net_ind.*fc_pos(265,:)';
        ind266 = net_ind.*fc_pos(266,:)';   
        fc = fc.*(fc>0);
        fc = atanh(fc);
        nandiag = diag(nan(266,1));
        fc = fc + nandiag;
        gbc_vectors_network{netid}(3,subj) = (nanmean(fc(265,ind265>0))); % left MT connectivity
        gbc_vectors_network{netid}(7,subj) = (nanmean(fc(266,ind266>0))); % right MT connectivity
        
        % tACS ON | OPP condition
        tacson_opp = opptimes{subj}.*tacsontimes{subj};
        ind = find(tacson_opp);
        fc = corrcoef(subj_data{subj}(:,ind)');
%         fc_pos = fc > 0;
        fc_pos = ones(266,1);
        ind265 = net_ind.*fc_pos(265,:)';
        ind266 = net_ind.*fc_pos(266,:)';  
        fc = fc.*(fc>0);
        fc = atanh(fc);
        nandiag = diag(nan(266,1));
        fc = fc + nandiag;
        gbc_vectors_network{netid}(4,subj) = (nanmean(fc(265,ind265>0))); % left MT connectivity
        gbc_vectors_network{netid}(8,subj) = (nanmean(fc(266,ind266>0)));
    end
end

% Run ANOVA
% Define group names
adapt = {'same';'opp';'same';'opp';'same';'opp';'same';'opp';};
adapt = repmat(adapt,10,1);
tacs = {'tacs_off';'tacs_off'; 'tacs_on'; 'tacs_on';'tacs_off';'tacs_off'; 'tacs_on'; 'tacs_on';};
tacs = repmat(tacs,10,1);
hemi = {'left';'left';'left';'left';'right';'right';'right';'right'};
hemi = repmat(hemi,10,1);
subjs = reshape(repmat(1:10',8,1),[80,1]);


gbc_interactionF = zeros(length(networks),1);
gbc_interactionP = zeros(length(networks),1);
for netid=1:length(networks)
%     netid = networks(net)
    indata = reshape(gbc_vectors_network{netid},[80,1]);
    [p,tbl] = anovan(indata,{adapt,tacs,hemi,subjs}, 'model', 'full','random', [4], 'varnames', {'adapt','tacs','hemi','subjs'});
    close(1)
    gbc_interactionF(netid) = tbl{12,6};
    gbc_interactionP(netid) = tbl{12,7};
end

% Run FDR correction on gbc_interaction for all networks (except for -1)
[qvals_network_gbc] = mafdr(gbc_interactionP(1:end),'BHFDR', true);
% qvals_network_gbc = [nan; qvals_network_gbc]; % uncertain network is nan
bonferroni_network_gbc = gbc_interactionP * length(gbc_interactionP(1:end));
% Write p-values to csv for FDR correction in python
csvwrite('gbc_interactionP_networks.csv', gbc_interactionP);

% Display results
for net=1:length(networks)
    netid = networks(net);
    disp(['Results for GBC from MTs to network ' num2str(netid)])
    disp(['F = ' num2str(gbc_interactionF(net))])
    disp(['p = ' num2str(gbc_interactionP(net))])
    disp(['q = ' num2str(qvals_network_gbc(net))])
    disp(['Bonferonni = ' num2str(bonferroni_network_gbc(net))])
    disp([' '])
end

% Plot effect for DAN only

figure
lmt_same = [mean(gbc_vectors_network{13}(1,:)), mean(gbc_vectors_network{13}(3,:))];
rmt_same = [mean(gbc_vectors_network{13}(5,:)), mean(gbc_vectors_network{13}(7,:))];
lmt_opp = [mean(gbc_vectors_network{13}(2,:)), mean(gbc_vectors_network{13}(4,:))];
rmt_opp = [mean(gbc_vectors_network{13}(6,:)), mean(gbc_vectors_network{13}(8,:))];
numsubjs = 10;
e_lmt_same = [std(gbc_vectors_network{13}(1,:))/numsubjs, std(gbc_vectors_network{13}(3,:))/numsubjs];
e_rmt_same = [std(gbc_vectors_network{13}(5,:))/numsubjs, std(gbc_vectors_network{13}(7,:))/numsubjs];
e_lmt_opp = [std(gbc_vectors_network{13}(2,:))/numsubjs, std(gbc_vectors_network{13}(4,:))/numsubjs];
e_rmt_opp = [std(gbc_vectors_network{13}(6,:))/numsubjs, std(gbc_vectors_network{13}(8,:))/numsubjs];

hold on
title('Average FC from hMT+ to DAN')
errorbar(lmt_same,e_lmt_same,'b')
errorbar(lmt_opp,e_lmt_opp, 'c')
errorbar(rmt_same,e_rmt_same,'r')
errorbar(rmt_opp,e_rmt_opp, 'm')
xlim([0.9,2.1])
ylabel('Average Functional Connectivity')
set(gca, 'XTick', [1.0, 2.0])
set(gca, 'XTickLabel', {'tACS Off', 'tACS On'})
% ax.XTickLabel({'
legend('left hMT+ adapt', 'left hMT+ non-adapt', 'right hMT+ adapt', 'right hMT+ non-adapt', 'Location', 'northwest')



%% Now run ANOVA on whole brain for each ROI independently (no GBC measures)
% Organize matrices
fc_vectors = zeros(8,10,266);

for subj=1:length(subjNums)
    
    % tACS OFF | SAME condition
    tacsoff_same = sametimes{subj}.*tacsofftimes{subj};
    ind = find(tacsoff_same);
    fc = corrcoef(subj_data{subj}(:,ind)');
    fc_pos = fc > 0;
  fc = fc_pos.*fc;
    fc = atanh(fc);
    nandiag = diag(nan(266,1));
    fc = fc + nandiag;
    fc_vectors(1,subj,:) = fc(265,:); % left MT connectivity
    fc_vectors(5,subj,:) = fc(266,:); % right MT connectivity
    
    % tACS OFF | OPP condition
    tacsoff_opp = opptimes{subj}.*tacsofftimes{subj};
    ind = find(tacsoff_opp);
    fc = corrcoef(subj_data{subj}(:,ind)');
    fc_pos = fc > 0;
  fc = fc_pos.*fc;
    fc = atanh(fc);
    nandiag = diag(nan(266,1));
    fc = fc + nandiag;
    fc_vectors(2,subj,:) = fc(265,:); % left MT connectivity
    fc_vectors(6,subj,:) = fc(266,:); % right MT connectivity
    
    % tACS ON | SAME
    tacson_same = sametimes{subj}.*tacsontimes{subj};
    ind = find(tacson_same);
    fc = corrcoef(subj_data{subj}(:,ind)');
    fc_pos = fc > 0;
  fc = fc_pos.*fc;
    fc = atanh(fc);
    nandiag = diag(nan(266,1));
    fc = fc + nandiag;
    fc_vectors(3,subj,:) = fc(265,:); % left MT connectivity
    fc_vectors(7,subj,:) = fc(266,:); % right MT connectivity
    
    % tACS ON | OPP condition
    tacson_opp = opptimes{subj}.*tacsontimes{subj};
    ind = find(tacson_opp);
    fc = corrcoef(subj_data{subj}(:,ind)');
    fc_pos = fc > 0;
  fc = fc_pos.*fc;
    fc = atanh(fc);
    nandiag = diag(nan(266,1));
    fc = fc + nandiag;
    fc_vectors(4,subj,:) = fc(265,:); % left MT connectivity
    fc_vectors(8,subj,:) = fc(266,:);
end


% Define group names
adapt = {'same';'opp';'same';'opp';'same';'opp';'same';'opp';};
adapt = repmat(adapt,10,1);
tacs = {'tacs_off';'tacs_off'; 'tacs_on'; 'tacs_on';'tacs_off';'tacs_off'; 'tacs_on'; 'tacs_on';};
tacs = repmat(tacs,10,1);
hemi = {'left';'left';'left';'left';'right';'right';'right';'right'};
hemi = repmat(hemi,10,1);
subjs = reshape(repmat(1:10',8,1),[80,1]);
ps = cell(266,1);
interactionF = zeros(266,1);
interactionP = zeros(266,1);
for roi=1:266
    if roi==257
        continue
    end
    if roi==262
        continue
    end
    indata = reshape(fc_vectors(:,:,roi),[80,1]);
    [p,tbl] = anovan(indata,{adapt,tacs,hemi,subjs}, 'model', 'full','random',[4], 'varnames', {'adapt','tacs','hemi','subjs'});
    close(1)
    interactionF(roi) = tbl{12,6};
    interactionP(roi) = tbl{12,7};
end
% Write p-values to CSV for FDR correction
csvwrite('gbc_interactionP_266rois.csv', interactionP);

