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
    print('Total Number of iterations =', numIters, 'Max Distance =', maxDist)
    print('biggestChange: ', biggestChange)
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
    # k = 10 # a good k value is 10
    cutoff = 0.00001 # biggest change in centroids is typically 0 or 10^-16 (very small)
    
#    # for multiple k values
#    #kvals = [2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,120,130,136]
#    kvalstest = [10]
#    clusterings = []
#    for k in kvalstest:
#        clusters, maxDist = kmeans(active_sites, k, cutoff, minIters = 3, maxIters = 100, toPrint = False)
#        clusterings.append(clusters)
#    
#    return clusterings
    
    # for one k value (for testing)
    k = 2
    clusters, maxDist = kmeans(active_sites, k, cutoff, minIters = 3, maxIters = 100, toPrint = False)
    
    return clusters, maxDist

#%%
### Cluster hierarchically with agglomerative clustering

def agglomerativecluster(points, cutoff = 0.01, minIters = 3, maxIters = 136, toPrint = False):
    """ Returns (Cluster list, max dist of any point to its cluster) """
    """Note: it would be more efficient to collect the clusters as a clustering each time two clusters were merged, but I wasn't able to implement this without going all the way through to one final cluster.
        """
    #clusterings = [] # store each set of clusters in a list, to use for writing multiple clusterings at each level of the hierarchy
    clusters = [] # store clusters in a list
    # Initialize by putting each point in its own cluster
    for p in points:
        clusters.append(Cluster([p]))
    #newclusters = clusters[:]
    #clusterings.append(newclusters) # add initial clustering to the clusterings list
    
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
            print('Min dist: ', mindist)
            #print('clusters[i].points new', clusters[firstclusterI].points)
            #print('clusters[j].points ', clusters[secondclusterI].points)
        

        del clusters[secondclusterI] # delete old cluster j
        
        #problem with del: it deletes the original cluster so new clustering won't be appended, but I couldn't implement other methods of removing the old cluster
        #clusters = clusters[:secondclusterI] + clusters[secondclusterI+1 :]
        #clusterings.append(tempclusters)
        #clusters = tempclusters[:] # copy clusters
        
        numIters += 1
            
    maxDist = 0.0
    for c in clusters:
        for p in c.getPoints():
            pDist = compute_similarity(p, c.getCentroid())
            if pDist > maxDist:
                maxDist = pDist
    print('Total Number of iterations =', numIters, 'Max distance =', maxDist)
    
    return clusters, maxDist
    
    

def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    #cutoff = 0.00001 # cutoff results in 39 clusters (97 loops)
#    cutoff = 5 # use high cutoff (5) to get all clusterings for mult_cluterings
    #cutoff = 0.29 # results in 10 clusters
    cutoff = .03
    # for testing
    clusters, maxDist = agglomerativecluster(active_sites, cutoff, minIters = 2, maxIters = 3, toPrint = True)
    
    # for running normally
    #clusters, maxDist = agglomerativecluster(active_sites, cutoff, minIters = 3, maxIters = 135, toPrint = True)
#    
    return clusters, maxDist
    
##    # for multiple clusterings, use a range of cutoff scores
##    exponents = range(17)
#    cutoffvalstest = [0.6, .2, .16] # store a range of cutoff scores from highest to lowest (highest cutoff score means fewest number of clusters)
#    cutoffvals = [.6, .4, .3, .25, .2, .15, .1, .05, .01]
##    for e in exponents:
##        cutoffvals.append(10**-e)
##        
#    # find multiple clusterings for a range of cutoff values
#    clusterings = []
#    for cutoff in cutoffvals:
#        print('cutoff', cutoff)
#        clusters, maxDist = agglomerativecluster(active_sites, cutoff, minIters = 3, maxIters = 136, toPrint = False)
#        clusterings.append(clusters)
#    
#    return clusterings
    

#%%

# Cluster evaluation using a metric similar to the silhouette score
def cluster_eval(clusters):
    """ Evaluate clusters by the average of the distance between each point and its own cluster centroid compared with each point and its nearest neighbor centroid (similar to the silhouette score comparing dissimilarity within a cluster to dissimilarity between clusters).
    
    Score for one point is s = (b-a)/max(a,b) where b is intercluster dissimilarity and a is intracluster dissimilarity
    
    Return the average score over all points which gives the evaluation score of the clustering
    """
    
    scores = [] # store scores at each point
    for c in clusters:
        for p in c.getPoints():
            # calculate a score for every point (intracluster similarity)
            a = compute_similarity(p, c.getCentroid())
            
            # calculate b score for every point (intercluster similarity)
            b = 100 # initialize to high number, find min
            for nearestc in clusters: # try to find nearest cluster to p
                if nearestc == c: # skip if comparing point to its own cluster
                    continue
                interclustdist = compute_similarity(p, nearestc.getCentroid())
                if interclustdist < b and interclustdist > 0:
                    # interclustdist will be 0 if cluster is empty
                    b = interclustdist # set similarity to the nearest cluster
            
            # calculate score at each point
            score = (b - a)/(max(a,b))
            scores.append(score)
            
    # average the scores over all points
    scoresarray = np.array(scores)
    evalscore = np.mean(scoresarray)
    
    return evalscore