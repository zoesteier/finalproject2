from .utils import Atom, Residue, ActiveSite
import numpy as np
import random

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """

    """Similarity metric: Find the Euclidean distance between the categorized amino acids in each active site (e.g. difference in fraction of polar aa's, difference in fraction of hydrophobic aa's, etc.). Categories include acidic, basic, polar, nonpolar, and bulky hydrophobic residues. Residues are categorized in io.
    
    To calculate similarity, first find the difference between the fraction of aa's in each category. Then find the square root of the sum of squared differences. Multiply the square differences by a weight to adjust how important each category is in determining similarity.
    """
#   # From site.categories list of aa fractions, make a sorted list
#   # list = [(category1,value1), (category2,value2)]
#    alist = sorted(site_a.categories.items())
#    acatvals = np.array([])
#    for cat in alist:
#        acatvals = np.append(acatvals, cat[1]) # list of values associated with each category
#    blist = sorted(site_b.categories.items()) 
#    bcatvals = np.array([])
#    for cat in blist:
#        bcatvals = np.append(bcatvals, cat[1]) # list of values associated with each category
    
    # Find differences in fraction of aa's in each category                         
    differences = site_a.categories - site_b.categories
#    print('differences: ', differences)
    
    # Weight each category (weight = 0 means ignore this category, weight =  1 is max)
    weights = np.array([0.0]*5) # array of weights
    weights[0] = 1 # acidic
    weights[1] = 1 # basic
    weights[2] = 1 # hydrophobic
    weights[3] = 1 # nonpolar
    weights[4] = 1 # polar
    
    # Find the Euclidean distance as the similarity score
    similarity = (sum(weights*differences**2))**0.5    
    
#    print('similarity: ', similarity)
    
    return similarity

#%%    
### Cluster by partitioning with kmeans

# Make a class for clusters, where points are the active sites anc centroid is the mean of the the points within a cluster
class Cluster(object):
    def __init__(self, points):
        self.points = points
        self.centroid = self.computeCentroid()
    def getCentroid(self):
        return self.centroid
    def computeCentroid(self):
        # the centroid is of class ActiveSite and contains category values that are the mean of all of the points contained within this cluster
        # cluster across the number of dimensions included in the categories
        dimensions = len(self.points[0].categories) # currently 5 categories of aa's
        totVals = np.array([0.0]*dimensions) # sum of all points in the cluster
        for p in self.points:
            totVals += p.categories
        meanPoint = ActiveSite('mean')
        meanPoint.categories = totVals/float(len(self.points))                                    
        return meanPoint # an ActiveSite with mean category values
    def update(self, points): 
        # when given new points in a cluster, find a new centroid
        # return the similarity distance between the old and new centroid
        oldCentroid = self.centroid
        self.points = points
        if len(points) > 0:
            self.centroid = self.computeCentroid()
            return compute_similarity(oldCentroid, self.centroid)
        else:
            return 0.0
    def getPoints(self):
        return self.points
    def contains(self, name):
        for p in self.points:
            if p.getName() == name:
                return True
        return False
    def toStr(self):
        result = ''
        for p in self.points:
            result = result + repr(p) + ', '
        return result[:-2]
    def __str__(self):
        result = ''
        for p in self.points:
            result = result + p.repr() + ', '
        return result[:-2]
        
        
        
def kmeans(points, k, cutoff, minIters = 3, maxIters = 100, toPrint = False):
    """ Returns (Cluster list, max dist of any point to its cluster) """
    # k = number of clusters
    # cutoff: terminate if change in centroid distance between iterations is smaller than the cutoff
    # minIters = minimum iterations before termination
    # maxIters: terminate if reach maximum iterations even if centroids still changing
    
    #Uses random initial centroids
    print('k: ', k)
    initialCentroids = random.sample(points,k)
    clusters = []
    for p in initialCentroids: # put each initial centroid in its own cluster
        clusters.append(Cluster([p])) # clusters is a list of class Cluster
    numIters = 0 # keep count of number of iterations until termination
    biggestChange = cutoff
    while (biggestChange >= cutoff and numIters < maxIters) or numIters < minIters:
        print("Starting iteration " + str(numIters))
        newClusters = []
        for c in clusters:
            newClusters.append([]) # newClusters is a list of lists of points
        for p in points:
            smallestDistance = compute_similarity(p, clusters[0].getCentroid())
            # initialize smallestDistance to distance between current point and centroid of first cluster
            index = 0
            for i in range(len(clusters)): # find distance from p to each cluster's centroid
                distance = compute_similarity(p, clusters[i].getCentroid())
                if distance < smallestDistance:
                    smallestDistance = distance
                    index = i #index = index of the cluster with smallest distance between p and centroid
            newClusters[index].append(p)
        biggestChange = 0.0 # track the biggest difference between old and new cluster in order to reach a point where the biggest difference is minimized below cutoff
        for i in range(len(clusters)):
            change = clusters[i].update(newClusters[i]) # puts new points into clusters, updates centroid, change = distance from old to new centroid
            #print "Cluster " + str(i) + ": " + str(len(clusters[i].points))
            biggestChange = max(biggestChange, change)
        numIters += 1
        if toPrint:
            print('Iteration count =', numIters)
    maxDist = 0.0
    for c in clusters:
        for p in c.getPoints():
            pDist = compute_similarity(p, c.getCentroid())
            if pDist > maxDist:
                maxDist = pDist
    print('Total Number of iterations =', numIters, 'Max Diameter =', maxDist)
    print(biggestChange)
    return clusters, maxDist
        

def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
            
    Clustering method: K-means (minimizes distance from sites to mean of cluster)
    Start with K clusters. Assign K random sites as initial centroids. Assign each object to nearest cluster. Calculate updated centroid. Reassign sites until centroids don't move or have iterated some max number of times
    """
    k = 10
    cutoff = 0.00001
    clusters, maxDist = kmeans(active_sites, k, cutoff, minIters = 3, maxIters = 100, toPrint = False)
    
    return clusters, maxDist

#%%
### Cluster hierarchically with agglomerative clustering

def agglomerativecluster(points, cutoff = 0.01, minIters = 3, maxIters = 135, toPrint = False):
    """ Returns (Cluster list, max dist of any point to its cluster) """
    clusters = [] # store clusters in a list
    # Initialize by putting each point in its own cluster
    for p in points:
        clusters.append(Cluster([p]))
    
    # Termination criteria: stop agglomerating when the distance between clusters to be merged reaches a cutoff (meaning, clustering should stop when the two clusters to be merged are very distant from each other)
    mindist = 0
    numIters = 0 # keep count of number of iterations until termination
    maxscore = 5**0.5 # the max similarity score (use for similarity between a point and itself so it isn't recognized as the min distance)
    
    # Check similarity between clusters, find nearest two and join
    # Linkage criteria: centroid linkage (merge two clusters with the nearest distance between their centroids).
    while (mindist < cutoff and numIters < maxIters) or numIters < minIters:
        # Make an nxn array of zeros where n = number of clusters to store distances between each cluster's centroid
        clusterdistances = np.full((len(clusters),len(clusters)), maxscore)
        for i in range(len(clusters)):
            for j in range(len(clusters)):
                if i == j:
                    # similarity between point and itself will be 0
                    # distance was initialized to max, so skip finding distance
                    # this should not be recognized as the min distance
                    continue
                distance = compute_similarity(clusters[i].getCentroid(), clusters[j].getCentroid()) # distance between two clusters
                clusterdistances[i][j] = distance
                
        # find index of minimum cluster distance, output (r,c)
        #print(clusterdistances)
        mindistIndex = np.unravel_index(clusterdistances.argmin(), clusterdistances.shape)
        firstclusterI = mindistIndex[0] # index of first cluster
        secondclusterI = mindistIndex[1] # index of second cluster
        
        # get the minimum cluster distance from the array
        mindist = clusterdistances[firstclusterI][secondclusterI]
        
        # merge first and second clusters
        newCluster = []
        newCluster.extend(clusters[firstclusterI].points)
        newCluster.extend(clusters[secondclusterI].points)
        clusters[firstclusterI].update(newCluster) # add new points from j to cluster i, updates centroid, change = distance from old to new centroid
        if toPrint:
            print('Iteration count =', numIters)
            print('clusters[i].points new', clusters[firstclusterI].points)
            print('clusters[j].points ', clusters[secondclusterI].points)
        del clusters[secondclusterI] # delete old cluster j
        
        numIters += 1
            
    maxDist = 0.0
    for c in clusters:
        for p in c.getPoints():
            pDist = compute_similarity(p, c.getCentroid())
            if pDist > maxDist:
                maxDist = pDist
    print('Total Number of iterations =', numIters, 'Max Diameter =', maxDist)
    return clusters, maxDist
    
    

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    cutoff = 0.00001
    clusters, maxDist = agglomerativecluster(active_sites, cutoff, minIters = 3, maxIters = 135, toPrint = True)
    
    return clusters, maxDist
    

