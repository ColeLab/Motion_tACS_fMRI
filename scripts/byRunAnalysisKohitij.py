# Taku Ito
# This script splits up each of the runs within the tACS on v off condition
# Attempts to see if there is a main effect of time that may confound our results
# This will only output the correlation values of each run (4 runs) for each of the 6 clusters with either MTs
# Across all subjects
import scipy.stats as stats
import numpy as np

subjects = ['038', '069', '141', '172', '173', '177', '178', '083', '144', '170']

basedir = '/projects/Collaborations/KrekelbergCollaboration/MotionAdaptation_tACS_FC/data/'

tCSon = {}
tCSoff = {}
correlations = {}

# Output matrix
outmat = np.zeros(shape=(20,4,6))
subjcount = 0 
for subj in subjects:

    subjdir = basedir + subj + '/fMRI/tacs_motionadaptationAnalysis/'

    tCSon[subj] = {}
    tCSoff[subj] = {}

    # Instantiate correlation dict
    correlations[subj] = {}

    # Load tACS on timeseries by run and MT hemisphere
    seedon = {}
    seedon['left'] = {}
    seedon['right'] = {}
    seedon['left'][1] = np.loadtxt(subjdir + 'stim_leftAreaMT_v2.1D')[0:181]
    seedon['left'][2] = np.loadtxt(subjdir + 'stim_leftAreaMT_v2.1D')[181:]
    seedon['right'][1] = np.loadtxt(subjdir + 'stim_rightAreaMT_v2.1D')[0:181]
    seedon['right'][2] = np.loadtxt(subjdir + 'stim_rightAreaMT_v2.1D')[181:]

    # Load tACS off timeseries by run and MT hemisphere
    seedoff = {}
    seedoff['left'] = {}
    seedoff['right'] = {}
    seedoff['left'][1] = np.loadtxt(subjdir + 'nostim_leftAreaMT_v2.1D')[0:181]
    seedoff['left'][2] = np.loadtxt(subjdir + 'nostim_leftAreaMT_v2.1D')[181:]
    seedoff['right'][1] = np.loadtxt(subjdir + 'nostim_rightAreaMT_v2.1D')[0:181]
    seedoff['right'][2] = np.loadtxt(subjdir + 'nostim_rightAreaMT_v2.1D')[181:]


    for clust in range(1,7):
        # Instantiate dicts for each run separately (within condition)
        tCSon[subj][clust] = {}
        tCSoff[subj][clust] = {}
        correlations[subj][clust] = {}
        correlations[subj][clust]['left'] = {}
        correlations[subj][clust]['right'] = {}

        # load data
        tmpon = np.loadtxt(subjdir + 'cluster' + str(clust) + '_yestACS_timeseries.1D')
        tmpoff = np.loadtxt(subjdir + 'cluster' + str(clust) + '_notACS_timeseries.1D')
        
        tCSon[subj][clust][1] = tmpon[0:181]
        tCSon[subj][clust][2] = tmpon[181:]
        tCSoff[subj][clust][1] = tmpoff[0:181]
        tCSoff[subj][clust][2] = tmpoff[181:]

        # Compute actual correlations with lMT
        correlations[subj][clust]['left'][1], p = stats.pearsonr(tCSoff[subj][clust][1], seedoff['left'][1])
        correlations[subj][clust]['left'][2], p = stats.pearsonr(tCSoff[subj][clust][2], seedoff['left'][2])
        correlations[subj][clust]['left'][3], p = stats.pearsonr(tCSon[subj][clust][1], seedon['left'][1])
        correlations[subj][clust]['left'][4], p = stats.pearsonr(tCSon[subj][clust][2], seedon['left'][2])
        # Correlation with rMT
        correlations[subj][clust]['right'][1], p = stats.pearsonr(tCSoff[subj][clust][1], seedoff['right'][1])
        correlations[subj][clust]['right'][2], p = stats.pearsonr(tCSoff[subj][clust][2], seedoff['right'][2])
        correlations[subj][clust]['right'][3], p = stats.pearsonr(tCSon[subj][clust][1], seedon['right'][1])
        correlations[subj][clust]['right'][4], p = stats.pearsonr(tCSon[subj][clust][2], seedon['right'][2])

        outmat[subjcount,0,clust-1] = correlations[subj][clust]['left'][1]
        outmat[subjcount,1,clust-1] = correlations[subj][clust]['right'][1]
        outmat[subjcount,2,clust-1] = correlations[subj][clust]['left'][2]
        outmat[subjcount,3,clust-1] = correlations[subj][clust]['right'][2]
        
        outmat[subjcount+10,0,clust-1] = correlations[subj][clust]['left'][3]
        outmat[subjcount+10,1,clust-1] = correlations[subj][clust]['right'][3]
        outmat[subjcount+10,2,clust-1] = correlations[subj][clust]['left'][4]
        outmat[subjcount+10,3,clust-1] = correlations[subj][clust]['right'][4]


    subjcount += 1


# Now right out CSV files for ANOVA for each cluster
for clust in range(1,7):
    np.savetxt('4WayANOVA_STIMxLRxTIMExSUBJ_cluster' + str(clust) + '.csv', outmat[:,:,clust-1], delimiter=',')
