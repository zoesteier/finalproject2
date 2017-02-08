from hw2skeleton import cluster
from hw2skeleton import io
import os

def test_similarity():
    filename_a = os.path.join("data", "14181.pdb")
    filename_b = os.path.join("data", "4629.pdb")

    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    # These two active sites should be identical
    assert cluster.compute_similarity(activesite_a, activesite_b) == 0.0

def test_partition_clustering():
    # tractable subset
    pdb_ids = [34047, 45127, 50362, 82238, 96099]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    clusters, dist = cluster.cluster_by_partitioning(active_sites)
    if len(clusters[0].points) == 3:
        clusterA = clusters[0]
        clusterB = clusters[1]
    else:
        clusterA = clusters[1]
        clusterB = clusters[0]
    
    # Expect one cluster will contain 34047, 45127, 50632 and a second cluster will contain 82238, 96099
    assert (active_sites[0] in clusterA.points)
    assert (active_sites[1] in clusterA.points)
    assert (active_sites[2] in clusterA.points)
    assert (active_sites[3] in clusterB.points)
    assert (active_sites[4] in clusterB.points)
    
    

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [34047, 45127, 50362, 82238, 96099]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # Expect one cluster will contain 34047, 45127, 50632 and a second cluster will contain 82238, 96099
    clusters, dist = cluster.cluster_hierarchically(active_sites)
    print(clusters[0].points)
    if len(clusters[0].points) == 3:
        clusterA = clusters[0]
        clusterB = clusters[1]
    else:
        clusterA = clusters[1]
        clusterB = clusters[0]
        
    # Expect one cluster will contain 34047, 45127, 50632 and a second cluster will contain 82238, 96099
    assert (active_sites[0] in clusterA.points)
    assert (active_sites[1] in clusterA.points)
    assert (active_sites[2] in clusterA.points)
    assert (active_sites[3] in clusterB.points)
    assert (active_sites[4] in clusterB.points)
    
    
    
