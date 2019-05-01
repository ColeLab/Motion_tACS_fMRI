# Taku Ito
# 10/2/2015

import numpy as np

# Short script to take in time series for each of the 6 ROIs that intersected for significant FC for left and right MTs with in the 
# no stimulation condition. Want to compare the correlations with both MTs to these ROIs for the no stim and stim conditions across all subjects

listOfSubjects = ['038', '069', '141', '172', '173', '177', '178', '083', '144', '170']

basedir='/projects/Collaborations/KrekelbergCollaboration/MotionAdaptation_tACS_FC/data/'

corr_nostim_lMT = {}
corr_nostim_rMT = {}
corr_stim_lMT = {}
corr_stim_rMT = {}

for subj in listOfSubjects:
    corr_nostim_lMT[subj] = {}
    corr_nostim_rMT[subj] = {}
    corr_stim_lMT[subj] = {}
    corr_stim_rMT[subj] = {}

    subjdir = basedir + subj + '/fMRI/tacs_motionadaptationAnalysis/'

    nostim_lMT = np.loadtxt(subjdir + 'nostim_leftAreaMT_v2.1D', dtype='float')
    nostim_rMT = np.loadtxt(subjdir + 'nostim_rightAreaMT_v2.1D', dtype='float')
    stim_lMT = np.loadtxt(subjdir + 'stim_leftAreaMT_v2.1D', dtype='float') 
    stim_rMT = np.loadtxt(subjdir + 'stim_rightAreaMT_v2.1D', dtype='float') 

    for clust in range(1,7):

        # Run correlation for each subject's time series
        
        nostim_clust= np.loadtxt(subjdir + 'cluster' + str(clust) + '_notACS_timeseries.1D', dtype='float')
        stim_clust = np.loadtxt(subjdir + 'cluster' + str(clust) + '_yestACS_timeseries.1D', dtype='float')
        
        # Now get correlations of each of the 4 conditions
        corr_nostim_lMT[subj][clust] = np.corrcoef(nostim_clust, nostim_lMT)[0,1]
        corr_nostim_rMT[subj][clust] = np.corrcoef(nostim_clust, nostim_rMT)[0,1]
        corr_stim_lMT[subj][clust] = np.corrcoef(stim_clust, stim_lMT)[0,1]
        corr_stim_rMT[subj][clust] = np.corrcoef(stim_clust, stim_rMT)[0,1]


outMat = np.zeros(shape=(20,2,6)) # 2 x 2 ANOVA x 10 replications, by 6 clusters
# Inefficient but...
for clust in range(1,7):
    
    rowcount = 0
    # a levels are stimulation v no stim
    for alevel in range(1,3):
        for subj in listOfSubjects:
            if alevel==1: # Then nostim group
                outMat[rowcount,0,clust-1] = corr_nostim_lMT[subj][clust]
                outMat[rowcount,1,clust-1] = corr_nostim_rMT[subj][clust]
            elif alevel==2: # then stim group
                outMat[rowcount,0,clust-1] = corr_stim_lMT[subj][clust]
                outMat[rowcount,1,clust-1] = corr_stim_rMT[subj][clust]
        
            rowcount += 1 
        
    #write to csv
    np.savetxt('cluster' + str(clust) + '_ANOVA_mat.csv', outMat[:,:,clust-1], delimiter=',')
        




