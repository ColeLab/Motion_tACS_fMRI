% Taku Ito
% 10/2/15

% Run a 2 x 2 ANOVA on the correlations for stim v no_stim X left v right
% Import the ANOVA matrix as a csv for each cluster. 
% Originally was to do in python until I realized that python has no good 
% ANOVA implementation...

anovaMat = zeros(20,2,6);
for clust=1:6
    anovaMat(:,:,clust) = csvread(['cluster' num2str(clust) '_ANOVA_mat.csv']);

    
end

