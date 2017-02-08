# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 14:49:42 2017

@author: Zoë
"""
# Use this script to try running each algorithm as one clustering and as multiple clusterings to confirm they are functioning as expected
#This script is not necessary for the functioning of the clustering

#from .utils import Atom, Residue, ActiveSite
from hw2skeleton import io
from hw2skeleton import cluster
import matplotlib.pyplot as plt
import numpy as np

active_sites = io.read_active_sites("C:\\Users\Zoë\Documents\GitHub\hw2-skeleton\data")
#site1 = active_sites[5]
#site2 = active_sites[7]
#print('site1: ', site1.categories)
#print('site2: ', site2.categories)
#sim = cluster.compute_similarity(site1,site2)

# Run for one clustering by kmeans
Pclusters, PmaxDistance = cluster.cluster_by_partitioning(active_sites)
##for i in clusters:
##    print(i.toStr())
io.write_clustering("clusterPk=10", Pclusters)


# Run for just one clustering by agglomerative clustering
Hclusters, distH = cluster.cluster_hierarchically(active_sites)
io.write_clustering("clusterHcutoff=0.3", Hclusters)


## Run for one clustering by agglomerative clustering
#Hclusters, HmaxDist, Hclusterings = cluster.cluster_hierarchically(active_sites)
#io.write_mult_clusterings("clusteringsH1", Hclusterings)

#%%
## Clusterings of multiple k values in kmeans
#kvals = [2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,120,130,136]
#kvalstest = [2,10,100]
#
#numtrials = 10
#allscores = np.zeros((numtrials, len(kvals)))
#for i in range(numtrials):        
#    Pclusterings = cluster.cluster_by_partitioning(active_sites)
#    #io.write_mult_clusterings("clusteringsP2", Pclusterings)
#    multclusterscoresP = []
#    for clustering in Pclusterings:
#        score = cluster.cluster_eval(clustering)
#        multclusterscoresP.append(score)
#    allscores[i] = multclusterscoresP
#avescoresP = np.mean(allscores,0)
#
#plt.plot(kvals, avescoresP, '.-')  
##plt.plot(kvals, multclusterscoresP, '.-')    
#plt.xlabel('k value (number of clusters')
#plt.ylabel('Cluster evaluation score')
#plt.title('Evaluation of clustering by k means')

    
#%%
# Clusterings of multiple cutoff values in agglomerative

#cutoffvalstest = [0.6, .2, .16]
#cutoffvals = [.6, .4, .3, .25, .2, .15, .1, .05, .01]
#Hclusterings = cluster.cluster_hierarchically(active_sites)
#io.write_mult_clusterings("clusteringsH3", Hclusterings)
#multclusterscoresH = []
#for clustering in Hclusterings:
#    score = cluster.cluster_eval(clustering)
#    multclusterscoresH.append(score)
# 
#plt.plot(cutoffvals, multclusterscoresH, '.-')  
##plt.plot(kvals, multclusterscoresP, '.-')    
#plt.xlabel('cutoff value (cophenetic distance)')
#plt.ylabel('Cluster evaluation score')
#plt.title('Evaluation of clustering by agglomerative hierarchical clustering')
    


