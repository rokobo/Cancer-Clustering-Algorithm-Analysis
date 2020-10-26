"""
Example code for creating and visualizing
cluster of county-based cancer risk data
"""

import math
import matplotlib.pyplot as plt
import csv
import random
import timeit

#############################################################
# Define colors for clusters.  Display a max of 16 clusters.
COLORS = ['Aqua', 'Yellow', 'Blue', 'Fuchsia', 'Black', 
          'Green', 'Lime', 'Maroon', 'Navy', 'Olive', 
          'Orange', 'Purple', 'Red', 'Brown', 'Teal']

#############################################################
# Clustering class

class Cluster:
    """
    Class for creating and merging clusters of counties
    """    
    def __init__(self, fips_codes, horiz_pos, vert_pos, population, risk):
        """
        Create a cluster based the models a set of counties' data
        """
        self._fips_codes = fips_codes
        self._horiz_center = horiz_pos
        self._vert_center = vert_pos
        self._total_population = population
        self._averaged_risk = risk        
        
    def __repr__(self):
        """
        String representation assuming the module is "alg_cluster".
        """
        rep = "Cluster("
        rep += str(self._fips_codes) + ", "
        rep += str(self._horiz_center) + ", "
        rep += str(self._vert_center) + ", "
        rep += str(self._total_population) + ", "
        rep += str(self._averaged_risk) + ")"
        return rep

    def fips_codes(self):
        """
        Get the cluster's set of FIPS codes
        """
        return self._fips_codes
    
    def horiz_center(self):
        """
        Get the averged horizontal center of cluster
        """
        return self._horiz_center
    
    def vert_center(self):
        """
        Get the averaged vertical center of the cluster
        """
        return self._vert_center
    
    def total_population(self):
        """
        Get the total population for the cluster
        """
        return self._total_population
    
    def averaged_risk(self):
        """
        Get the averaged risk for the cluster
        """
        return self._averaged_risk
         
    def copy(self):
        """
        Return a copy of a cluster
        """
        copy_cluster = Cluster(set(self._fips_codes), self._horiz_center, self._vert_center,
                               self._total_population, self._averaged_risk)
        return copy_cluster

    def distance(self, other_cluster):
        """
        Compute the Euclidean distance between two clusters
        """
        vert_dist = self._vert_center - other_cluster.vert_center()
        horiz_dist = self._horiz_center - other_cluster.horiz_center()
        return math.sqrt(vert_dist ** 2 + horiz_dist ** 2)
        
    def merge_clusters(self, other_cluster):
        """
        Merge one cluster into another
        The merge uses the relatively populations of each
        cluster in computing a new center and risk
        
        Note that this method mutates self
        """
        if len(other_cluster.fips_codes()) == 0:
            return self
        else:
            self._fips_codes.update(set(other_cluster.fips_codes()))
 
            # compute weights for averaging
            self_weight = float(self._total_population)                        
            other_weight = float(other_cluster.total_population())
            self._total_population = self._total_population + other_cluster.total_population()
            self_weight /= self._total_population
            other_weight /= self._total_population
                    
            # update center and risk using weights
            self._vert_center = self_weight * self._vert_center + other_weight * other_cluster.vert_center()
            self._horiz_center = self_weight * self._horiz_center + other_weight * other_cluster.horiz_center()
            self._averaged_risk = self_weight * self._averaged_risk + other_weight * other_cluster.averaged_risk()
            return self

    def cluster_error(self, data_table):
        """
        Input: data_table is the original table of cancer data used in creating the cluster.
        
        Output: The error as the sum of the square of the distance from each county
        in the cluster to the cluster center (weighted by its population)
        """
        # Build hash table to accelerate error computation
        fips_to_line = {}
        for line_idx in range(len(data_table)):
            line = data_table[line_idx]
            fips_to_line[line[0]] = line_idx
        
        # compute error as weighted squared distance from counties to cluster center
        total_error = 0
        counties = self.fips_codes()
        for county in counties:
            line = data_table[fips_to_line[county]]
            singleton_cluster = Cluster(set([line[0]]), line[1], line[2], line[3], line[4])
            singleton_distance = self.distance(singleton_cluster)
            total_error += (singleton_distance ** 2) * singleton_cluster.total_population()
        return total_error

#############################################################
# Plotting functions

def circle_area(pop):
    """
    Compute area of circle proportional to population
    """
    return math.pi * pop / (200.0 ** 2)

def plot_clusters(data_table, cluster_list, draw_centers = False):
    """
    Create a plot of clusters of counties
    """

    fips_to_line = {}
    for line_idx in range(len(data_table)):
        fips_to_line[data_table[line_idx][0]] = line_idx
     
    # Load map image
    map_img = plt.imread("Data\\USA_Counties.png")

    # Scale plot to get size similar to CodeSkulptor version
    ypixels, xpixels, bands = map_img.shape
    DPI = 60.0                  # adjust this constant to resize your plot
    xinch = xpixels / DPI
    yinch = ypixels / DPI
    plt.figure(figsize=(xinch,yinch))
    implot = plt.imshow(map_img)
   
    # draw the counties colored by cluster on the map
    if not draw_centers:
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            for fips_code in cluster.fips_codes():
                line = data_table[fips_to_line[fips_code]]
                plt.scatter(x = [line[1]], y = [line[2]], s =  circle_area(line[3]), lw = 1,
                            facecolors = cluster_color, edgecolors = cluster_color)

    # add cluster centers and lines from center to counties            
    else:
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            for fips_code in cluster.fips_codes():
                line = data_table[fips_to_line[fips_code]]
                plt.scatter(x = [line[1]], y = [line[2]], s =  circle_area(line[3]), lw = 1,
                            facecolors = cluster_color, edgecolors = cluster_color, zorder = 1)
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            cluster_center = (cluster.horiz_center(), cluster.vert_center())
            for fips_code in cluster.fips_codes():
                line = data_table[fips_to_line[fips_code]]
                plt.plot( [cluster_center[0], line[1]],[cluster_center[1], line[2]], cluster_color, lw=1, zorder = 2)
        for cluster_idx in range(len(cluster_list)):
            cluster = cluster_list[cluster_idx]
            cluster_color = COLORS[cluster_idx % len(COLORS)]
            cluster_center = (cluster.horiz_center(), cluster.vert_center())
            cluster_pop = cluster.total_population()
            plt.scatter(x = [cluster_center[0]], y = [cluster_center[1]], s =  circle_area(cluster_pop), lw = 2,
                        facecolors = "none", edgecolors = "black", zorder = 3)

    plt.show()

#############################################################
# Code to load data tables

def load_data(choice):
    """
    code for loading data
    """
    data = open("Data\\unifiedCancerData_" + str(choice) + ".csv")
    dataReader = csv.reader(data)
    dataData = list(dataReader)
    
    for line in range(len(dataData)):
        dataData[line][0] = str(dataData[line][0])
        dataData[line][1] = float(dataData[line][1])
        dataData[line][2] = float(dataData[line][2])
        dataData[line][3] = int(dataData[line][3])
        dataData[line][4] = float(dataData[line][4])  
    return dataData

#############################################################
# Clustering Algorithms

def sequential_clustering(singleton_list, num_clusters):
    """
    Take a data table and create a list of clusters
    by partitioning the table into clusters based on its ordering
    
    Note that method may return num_clusters or num_clusters + 1 final clusters
    """  
    cluster_list = []
    cluster_idx = 0
    total_clusters = len(singleton_list)
    cluster_size = float(total_clusters)  / num_clusters
    
    for cluster_idx in range(len(singleton_list)):
        new_cluster = singleton_list[cluster_idx]
        if math.floor(cluster_idx / cluster_size) != \
           math.floor((cluster_idx - 1) / cluster_size):
            cluster_list.append(new_cluster)
        else:
            cluster_list[-1] = cluster_list[-1].merge_clusters(new_cluster)           
    return cluster_list

def slow_closest_pair(cluster_list):
    """
    Compute the distance between the closest pair of clusters in a list (slow)

    Input: cluster_list is the list of clusters
    
    Output: tuple of the form (dist, idx1, idx2) where the centers of the 
    clusters cluster_list[idx1] and cluster_list[idx2] have minimum distance dist.       
    """         
    min_dist = float('inf')
    idx = (-1, -1)
    cluster_range = range(len(cluster_list))
    for point1 in cluster_range:
        for point2 in cluster_range:
            if point1 != point2:
                distance = cluster_list[point1].distance(cluster_list[point2])
                if distance < min_dist:
                    min_dist = distance
                    idx = (point1, point2)
    return (min_dist, min(idx), max(idx))     

def closest_pair_strip(cluster_list, horiz_center, half_width):
    """
    Helper function to compute the closest pair of clusters in a vertical strip
    
    Input: cluster_list is a list of clusters produced by fast_closest_pair
    horiz_center is the horizontal position of the strip's vertical center line
    half_width is the half the width of the strip (i.e; the maximum horizontal distance
    that a cluster can lie from the center line)

    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] lie in the strip and have minimum distance dist.       
    """      
    clusters = []
    for idx in range(len(cluster_list)):
        if abs(cluster_list[idx].horiz_center() - horiz_center) < half_width:
            clusters.append(idx)
            
    clusters.sort(key = lambda cluster: cluster_list[cluster].vert_center())
    cluster_len = len(clusters)
    pair_dist = (float("inf"), -1, -1)

    if cluster_len == 2:
        distance = cluster_list[clusters[0]].distance(cluster_list[clusters[1]])
        return (distance, min(clusters[0], clusters[1]), 
                max(clusters[0], clusters[1]))
    
    for idx1 in range(cluster_len - 1):
        for idx2 in range(idx1 + 1, min((idx1 + 4, cluster_len)), 1):
            distance = cluster_list[clusters[idx1]].distance(cluster_list[clusters[idx2]])
            if pair_dist[0] > distance:
                pair_dist = (distance, min(clusters[idx1], clusters[idx2]),
                             max(clusters[idx1], clusters[idx2]))    
    return pair_dist

def fast_closest_pair(cluster_list):
    """
    Compute the distance between the closest pair of clusters in a list (fast)

    Input: cluster_list is list of clusters SORTED such that horizontal positions of their
    centers are in ascending order
    
    Output: tuple of the form (dist, idx1, idx2) where the centers of the clusters
    cluster_list[idx1] and cluster_list[idx2] have minimum distance dist.       
    """
    num_points = len(cluster_list)
    
    #base case
    if num_points <= 3:
        pair_dist = slow_closest_pair(cluster_list)
        
    #recursive case
    else:
        half_len = int(num_points/2)
        set_one = cluster_list[:half_len]
        set_two = cluster_list[half_len:]
        pair_dist1 = fast_closest_pair(set_one)
        pair_dist2 = fast_closest_pair(set_two)
        pair_dist2 = (pair_dist2[0], 
                      pair_dist2[1] + half_len, pair_dist2[2] + half_len)
        pair_dist3 = min((pair_dist1, pair_dist2))
        middle = 0.5 * (cluster_list[half_len - 1].horiz_center() + 
                        cluster_list[half_len].horiz_center())
        pair_dist4 = closest_pair_strip(cluster_list, middle, pair_dist3[0])
        pair_dist = min(pair_dist3, pair_dist4)
    return pair_dist

def hierarchical_clustering(cluster_list, num_clusters):
    """
    Compute a hierarchical clustering of a set of clusters
    Note: the function may mutate cluster_list
    
    Input: List of clusters, integer number of clusters
    Output: List of clusters whose length is num_clusters
    """
    clusters = []
    for cluster in cluster_list:
        clusters.append(cluster)
        
    while len(clusters) > num_clusters:
        clusters.sort(key=lambda cluster: cluster.horiz_center())
        closest = fast_closest_pair(clusters)
        cluster2 = clusters[closest[2]]
        clusters[closest[1]].merge_clusters(cluster2)
        clusters.remove(cluster2)
    return clusters
    
def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    """
    Compute the k-means clustering of a set of clusters
    Note: the function may not mutate cluster_list
    
    Input: List of clusters, integers number of clusters and number of iterations
    Output: List of clusters whose length is num_clusters
    """
    # position initial clusters at the location of clusters with largest populations
    cluster = list(cluster_list)
    cluster.sort(key=lambda clstr: clstr.total_population(), reverse=True)
    big_clusters = cluster[:num_clusters]    
    centers = []
    for cluster in big_clusters:
        centers.append(Cluster(set(),cluster.horiz_center(), 
                                        cluster.vert_center(), 0, 0))
        
    for dummy_iterations in range(1, num_iterations + 1):
        clusters = []
        for dummy_cluster in range(num_clusters):
            clusters.append(Cluster(set(), 0, 0, 0, 0))
        
        for cluster in cluster_list:
            closest = min(range(len(centers)), 
                          key=lambda index: centers[index].distance(cluster))
            clusters[closest].merge_clusters(cluster)
                      
        for index in range(num_clusters):
            new_center = Cluster(set(), 
                                             clusters[index].horiz_center(),
                                             clusters[index].vert_center(), 
                                             0, 0)
            centers[index] = new_center                
    return clusters

#############################################################
#testing and analysis

def gen_random_clusters(num_clusters):
    """
    creates a list of random clusters positioned inside a 1x1 square(for tests)
    """
    clusters = []
    for dummy_cluster in range(num_clusters):
        new = Cluster(set(),random.uniform(-1, 1), 
                                  random.uniform(-1, 1), 0, 0)
        clusters.append(new)
    return clusters
  
def calculate_time(start, end):
    """
    function for calculating time elapsed for closest clusters algorithms
    performs end-start operations then plots the data
    """
    fast = []
    slow = []
    
    for num in range(start, end + 1):
        clusters = gen_random_clusters(num)
        
        start_time = timeit.default_timer()
        fast_closest_pair(clusters)
        fast.append(timeit.default_timer() - start_time)
        
        start_time = timeit.default_timer()
        slow_closest_pair(clusters)
        slow.append(timeit.default_timer() - start_time)
     
    #graph details
    plt.plot(slow, "-b", label="Slow closest pair")
    plt.plot(fast, "-r", label="Fast closest pair")  
    plt.legend(loc="upper left")
    plt.ylabel("Running time")
    plt.xlabel("Clusters")   
    plt.title("Running time analysis") 
    plt.show()

def calculate_time2(start, end):
    """
    function for calculating time elapsed for clustering algorithms
    performs end-start operations then plots the data
    """
    kmeans = []
    hierarchical = []
    
    for num in range(start, end + 1):
        clusters = gen_random_clusters(num)
        
        start_time = timeit.default_timer()
        hierarchical_clustering(clusters, 15)
        hierarchical.append(timeit.default_timer() - start_time)
        
        start_time = timeit.default_timer()
        kmeans_clustering(clusters, 15, 5)
        kmeans.append(timeit.default_timer() - start_time)
     
    #graph details
    plt.plot(kmeans, "-b", label="K-means clustering")
    plt.plot(hierarchical, "-r", label="Hierarchical clustering")  
    plt.legend(loc="upper left")
    plt.ylabel("Running time")
    plt.xlabel("Clusters")   
    plt.title("Running time analysis") 
    plt.show()

def compute_distortion(choice, start, end):
    """
    calculates the error margin of a given clustering algorithm
    """
    data_table = load_data(choice)
    original_clusters = []
    distortion_x = []
    h_distortion_y = []
    k_distortion_y = []
    
    for line in data_table:
        original_clusters.append(Cluster(set([line[0]]), line[1], 
                                                     line[2], line[3], line[4]))
    
    for num_clusters in range(start, end + 1):
        clusters = list(original_clusters)
        cluster_list2 = kmeans_clustering(clusters, num_clusters, 5)
        cluster_list = hierarchical_clustering(clusters, num_clusters)
        k_distortion = 0
        h_distortion = 0
        distortion_x.append(num_clusters)
        
        for center in cluster_list2:
            k_distortion += center.cluster_error(data_table) 
        k_distortion_y.append(k_distortion)
        
        for center in cluster_list:
            h_distortion += center.cluster_error(data_table)    
        h_distortion_y.append(h_distortion)

    #graph details
    plt.plot(distortion_x, k_distortion_y, "-b", label="K-means clustering")
    plt.plot(distortion_x, h_distortion_y, "-r", label="Hierarchical clustering")  
    plt.legend(loc="upper left")
    plt.ylabel("Distortion")
    plt.xlabel("Clusters")   
    plt.title("Clustering distortion with " + str(choice) + " counties") 
    plt.show()

# Code to load cancer data, compute a clustering and visualize the results
def run_example(display, choice, number_of_clusters, kmeans_iterations=None):
    """
    Load a data table, compute a list of clusters and 
    plot a list of clusters
    """  
    data_table = load_data(choice)
    singleton_list = []
    for line in data_table:
        singleton_list.append(Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))
    
    if display == "S":
        cluster_list = sequential_clustering(singleton_list, 15)	
        print("Displaying", len(cluster_list), "sequential clusters")
        plot_clusters(data_table, cluster_list, False)
        #add cluster centers
        plot_clusters(data_table, cluster_list, True)  
    
    elif display == "H":  
        cluster_list = hierarchical_clustering(singleton_list, number_of_clusters)
        print("Displaying", len(cluster_list), "hierarchical clusters")
        plot_clusters(data_table, cluster_list, False)
        #add cluster centers
        plot_clusters(data_table, cluster_list, True)  
    
    elif display == "K":   
        cluster_list = kmeans_clustering(singleton_list, number_of_clusters, kmeans_iterations)	
        print("Displaying", len(cluster_list), "k-means clusters")
        plot_clusters(data_table, cluster_list, False)
        #add cluster centers
        plot_clusters(data_table, cluster_list, True)  
    
    else:
        print("Invalid display option")

calculate_time(0, 200)
print("Q.1 - Closest pair running time analysis plot finished\n")
   
run_example("H", 3108, 15)
print("Q.2 - Hierarchical clustering, 3108 counties, 15 clusters plot finished\n")

run_example("K", 3108, 15, 5)
print("Q.2 - K-Means clustering, 3108 counties, 15 clusters, 5 iterations plot finished\n")

print("Q.4 - K-Means is faster since it has a linear time of O(n)\n\
Hierarchical clustering has an exponential time of O(n^2)\n")

run_example("H", 111, 9)
print("Q.5 - Hierarchical clustering, 111 counties, 9 clusters plot finished\n")

run_example("K", 111, 9, 5)
print("Q.6 - K-Means clustering, 111 counties, 9 clusters, 5 iterations plot finished\n")

print("Q.7 - K-Means clustering distortion: 2.712542*(10^11)\n\
Hierarchical clustering distortion: 1.751639*(10^11)\n")

print("Q.8 - The higher distortion for K-Means can be explained by the three out of nine \n\
initial clusters for K-Means in southern California. This indicates the clusters were more \n\
spaced out when comparing it with the hierarchical clusters. Therefore, since K-Means chooses\n\
the largest clusters as the initial centers, it could lead to small clusters being merged\n\
with very far large centers in case there are none close\n")

print("Q.9 - Hierarchical clustering requires less human supervision\n")

compute_distortion(111, 6, 20)
print("Q.10.1 - Clustering distortion analysis on 111 counties plot finished")
compute_distortion(896, 6, 20)
print("Q.10.2 - Clustering distortion analysis on 896 counties plot finished\n")

print("Q.11 - On all datasets the hierarchical clustering method produced lower distortions")