# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 14:49:42 2017

@author: Zoë
"""

#from .utils import Atom, Residue, ActiveSite
from hw2skeleton import io
from hw2skeleton import cluster

active_sites = io.read_active_sites("C:\\Users\Zoë\Documents\GitHub\hw2-skeleton\data")
#site1 = active_sites[5]
#site2 = active_sites[7]
#print('site1: ', site1.categories)
#print('site2: ', site2.categories)
#sim = cluster.compute_similarity(site1,site2)

#Pclusters, PmaxDistance = cluster.cluster_by_partitioning(active_sites)
#for i in clusters:
#    print(i.toStr())
#io.write_clustering("cluster6", Pclusters)

Hclusters, HmaxDistance = cluster.cluster_hierarchically(active_sites)
io.write_clustering("clusterH1", Hclusters)