"""
Testing code for clustering methods 

hierarchical_clustering(cluster_list, num_clusters)
kmeans_clustering(cluster_list, num_clusters, num_iterations)
"""

import csv
import math

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
# Cluster algorithms for testing

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
# Tests

class TestSuite:
    """
    Create a suite of tests similar to unittest
    """
    def __init__(self):
        """
        Creates a test suite object
        """
        self.total_tests = 0
        self.failures = 0   
    
    def run_test(self, computed, expected, message = ""):
        """
        Compare computed and expected
        If not equal, print message, computed, expected
        """
        self.total_tests += 1
        if computed != expected:
            msg = message + " Computed: " + str(computed)
            msg += " Expected: " + str(expected)
            print(msg)
            self.failures += 1
    
    def report_results(self):
        """
        Report back summary of successes and failures
        from run_test()
        """
        msg = "Ran " + str(self.total_tests) + " tests. "
        msg += str(self.failures) + " failures."
        print (msg)

# Load data tables locally
DATA_3108 = "Data\\unifiedCancerData_3108.csv"
DATA_896 = "Data\\unifiedCancerData_896.csv"
DATA_290 = "Data\\unifiedCancerData_290.csv"
DATA_111 = "Data\\unifiedCancerData_111.csv"
DATA_24 = "Data\\unifiedCancerData_24.csv"

def load_data_table(data):
    """
    Import a table of county-based cancer risk data
    from a csv format file
    """
    data_file = open(data)
    dataReader = csv.reader(data_file)
    dataData = list(dataReader)
    
    for line in range(len(dataData)):
        dataData[line][0] = str(dataData[line][0])
        dataData[line][1] = float(dataData[line][1])
        dataData[line][2] = float(dataData[line][2])
        dataData[line][3] = int(dataData[line][3])
        dataData[line][4] = float(dataData[line][4])
        
    return dataData

# Helper function for converting list of clusters to set of county tuples

def set_of_county_tuples(cluster_list):
    """
    Input: A list of Cluster objects
    Output: Set of sorted tuple of counties corresponds to counties in each cluster
    """
    set_of_clusters = set([])
    for cluster in cluster_list:
        counties_in_cluster = cluster.fips_codes()
        
        # convert to immutable representation before adding to set
        county_tuple = tuple(sorted(list(counties_in_cluster)))
        set_of_clusters.add(county_tuple)
    return set_of_clusters

# Testing code

def test_hierarchical():
    """
    Test for hierarchical clustering
    Note that hierarchical_clustering mutates cluster_list
    """ 
    # load small data table
    print()
    print ("Testing hierarchical_clustering on 24 county set")
    data_24_table = load_data_table(DATA_24) 
    
    # test data of the form [size of output cluster, sets of county tuples]
    hierdata_24 = [[23, set([('11001', '51013'), ('01073',), ('06059',), ('06037',), ('06029',), ('06071',), ('06075',), ('08031',), ('24510',), ('34013',), ('34039',), ('34017',), ('36061',), ('36005',), ('36047',), ('36059',), ('36081',), ('41051',), ('41067',), ('51840',), ('51760',), ('55079',), ('54009',)])],
                   [22, set([('11001', '51013'), ('36047', '36081'), ('01073',), ('06059',), ('06037',), ('06029',), ('06071',), ('06075',), ('08031',), ('24510',), ('34013',), ('34039',), ('34017',), ('36061',), ('36005',), ('36059',), ('41051',), ('41067',), ('51840',), ('51760',), ('55079',), ('54009',)])],
                   [21, set([('11001', '51013'), ('36005', '36061'), ('36047', '36081'), ('01073',), ('06059',), ('06037',), ('06029',), ('06071',), ('06075',), ('08031',), ('24510',), ('34013',), ('34039',), ('34017',), ('36059',), ('41051',), ('41067',), ('51840',), ('51760',), ('55079',), ('54009',)])],
                   [20, set([('11001', '51013'), ('36005', '36061'), ('36047', '36081'), ('01073',), ('06059',), ('06037',), ('06029',), ('06071',), ('06075',), ('08031',), ('24510',), ('34039',), ('34013', '34017'), ('36059',), ('41051',), ('41067',), ('51840',), ('51760',), ('55079',), ('54009',)])],
                   [19, set([('34013', '34017', '34039'), ('11001', '51013'), ('36005', '36061'), ('36047', '36081'), ('01073',), ('06059',), ('06037',), ('06029',), ('06071',), ('06075',), ('08031',), ('24510',), ('36059',), ('41051',), ('41067',), ('51840',), ('51760',), ('55079',), ('54009',)])],
                   [18, set([('34013', '34017', '34039'), ('11001', '51013'), ('01073',), ('06059',), ('06037',), ('06029',), ('06071',), ('06075',), ('08031',), ('24510',), ('36059',), ('36005', '36047', '36061', '36081'), ('41051',), ('41067',), ('51840',), ('51760',), ('55079',), ('54009',)])],
                   [17, set([('11001', '51013'), ('01073',), ('06059',), ('06037',), ('06029',), ('06071',), ('06075',), ('08031',), ('24510',), ('36059',), ('34013', '34017', '34039', '36005', '36047', '36061', '36081'), ('41051',), ('41067',), ('51840',), ('51760',), ('55079',), ('54009',)])],
                   [16, set([('11001', '51013'), ('01073',), ('06059',), ('06037',), ('06029',), ('06071',), ('06075',), ('08031',), ('24510',), ('34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081'), ('41051',), ('41067',), ('51840',), ('51760',), ('55079',), ('54009',)])],
                   [15, set([('11001', '51013'), ('01073',), ('06059',), ('06037',), ('06029',), ('06071',), ('06075',), ('08031',), ('24510',), ('34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081'), ('41051', '41067'), ('51840',), ('51760',), ('55079',), ('54009',)])],
                   [14, set([('01073',), ('06059',), ('06037',), ('06029',), ('06071',), ('06075',), ('08031',), ('34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081'), ('41051', '41067'), ('51840',), ('51760',), ('55079',), ('54009',), ('11001', '24510', '51013')])],
                   [13, set([('06037', '06059'), ('01073',), ('06029',), ('06071',), ('06075',), ('08031',), ('34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081'), ('41051', '41067'), ('51840',), ('51760',), ('55079',), ('54009',), ('11001', '24510', '51013')])],
                   [12, set([('06037', '06059'), ('01073',), ('06029',), ('06071',), ('06075',), ('08031',), ('34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081'), ('41051', '41067'), ('51760',), ('55079',), ('54009',), ('11001', '24510', '51013', '51840')])],
                   [11, set([('06029', '06037', '06059'), ('01073',), ('06071',), ('06075',), ('08031',), ('34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081'), ('41051', '41067'), ('51760',), ('55079',), ('54009',), ('11001', '24510', '51013', '51840')])],
                   [10, set([('06029', '06037', '06059'), ('01073',), ('06071',), ('06075',), ('08031',), ('34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081'), ('41051', '41067'), ('55079',), ('54009',), ('11001', '24510', '51013', '51760', '51840')])],
                   [9, set([('01073',), ('06029', '06037', '06059', '06071'), ('06075',), ('08031',), ('34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081'), ('41051', '41067'), ('55079',), ('54009',), ('11001', '24510', '51013', '51760', '51840')])],
                   [8, set([('01073',), ('06029', '06037', '06059', '06071'), ('06075',), ('08031',), ('41051', '41067'), ('55079',), ('54009',), ('11001', '24510', '34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081', '51013', '51760', '51840')])],
                   [7, set([('01073',), ('06029', '06037', '06059', '06071'), ('06075',), ('08031',), ('41051', '41067'), ('55079',), ('11001', '24510', '34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081', '51013', '51760', '51840', '54009')])],
                   [6, set([('06029', '06037', '06059', '06071', '06075'), ('01073',), ('08031',), ('41051', '41067'), ('55079',), ('11001', '24510', '34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081', '51013', '51760', '51840', '54009')])],
                   [5, set([('06029', '06037', '06059', '06071', '06075'), ('08031',), ('41051', '41067'), ('01073', '55079'), ('11001', '24510', '34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081', '51013', '51760', '51840', '54009')])],
                   [4, set([('06029', '06037', '06059', '06071', '06075'), ('01073', '11001', '24510', '34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081', '51013', '51760', '51840', '54009', '55079'), ('08031',), ('41051', '41067')])],
                   [3, set([('06029', '06037', '06059', '06071', '06075', '41051', '41067'), ('01073', '11001', '24510', '34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081', '51013', '51760', '51840', '54009', '55079'), ('08031',)])],
                   [2, set([('01073', '11001', '24510', '34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081', '51013', '51760', '51840', '54009', '55079'), ('06029', '06037', '06059', '06071', '06075', '08031', '41051', '41067')])],
                   ]
   
    suite = TestSuite()
    
    for num_clusters, expected_county_tuple in hierdata_24:
        
        # build initial list of clusters for each test since mutation is allowed
        cluster_list = []
        for idx in range(len(data_24_table)):
            line = data_24_table[idx]
            cluster_list.append(Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))

        # compute answer
        Cluster_data = hierarchical_clustering(cluster_list, num_clusters)
        County_tuples = set_of_county_tuples(Cluster_data)
        
        # Prepare test
        error_message = "Testing hierarchical_clustering on 24 county table, num_clusters = " + str(num_clusters)
        error_message += "\nCounty tuples: " + str(County_tuples)
        error_message += "\nExpected county tuples: " + str(expected_county_tuple)
        suite.run_test(County_tuples == expected_county_tuple, True, error_message)

    suite.report_results()

def test_kmeans():
    """
    Test for k-means clustering
    kmeans_clustering should not mutate cluster_list, but make a new copy of each test anyways
    """
    
    # load small data table
    print()
    print ("Testing kmeans_clustering on 24 county set")
    data_24_table = load_data_table(DATA_24)
        
    kmeansdata_24 = [[15, 1, set([('34017', '36061'), ('06037',), ('06059',), ('36047',), ('36081',), ('06071', '08031'), ('36059',), ('36005',), ('55079',), ('34013', '34039'), ('06075',), ('01073',), ('06029',), ('41051', '41067'), ('11001', '24510', '51013', '51760', '51840', '54009')])], 
                     [15, 3, set([('34017', '36061'), ('06037', '06059'), ('06071',), ('36047',), ('36081',), ('08031',), ('36059',), ('36005',), ('55079',), ('34013', '34039'), ('06075',), ('01073',), ('06029',), ('41051', '41067'), ('11001', '24510', '51013', '51760', '51840', '54009')])],
                     [15, 5, set([('34017', '36061'), ('06037', '06059'), ('06071',), ('36047',), ('36081',), ('08031',), ('36059',), ('36005',), ('55079',), ('34013', '34039'), ('06075',), ('01073',), ('06029',), ('41051', '41067'), ('11001', '24510', '51013', '51760', '51840', '54009')])],
                     [10, 1, set([('34017', '36061'), ('06029', '06037', '06075'), ('11001', '24510', '34013', '34039', '51013', '51760', '51840', '54009'), ('06059',), ('36047',), ('36081',), ('06071', '08031', '41051', '41067'), ('36059',), ('36005',), ('01073', '55079')])],
                     [10, 3, set([('34013', '34017', '36061'), ('06029', '06037', '06075'), ('08031', '41051', '41067'), ('06059', '06071'), ('34039', '36047'), ('36081',), ('36059',), ('36005',), ('01073', '55079'), ('11001', '24510', '51013', '51760', '51840', '54009')])],
                     [10, 5, set([('34013', '34017', '36061'), ('06029', '06037', '06075'), ('08031', '41051', '41067'), ('06059', '06071'), ('34039', '36047'), ('36081',), ('36059',), ('36005',), ('01073', '55079'), ('11001', '24510', '51013', '51760', '51840', '54009')])],
                     [5, 1, set([('06029', '06037', '06075'), ('01073', '11001', '24510', '34013', '34017', '34039', '36047', '51013', '51760', '51840', '54009', '55079'), ('06059',), ('36005', '36059', '36061', '36081'), ('06071', '08031', '41051', '41067')])],
                     [5, 3, set([('06029', '06037', '06075'), ('11001', '24510', '34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081', '51013'), ('08031', '41051', '41067'), ('06059', '06071'), ('01073', '51760', '51840', '54009', '55079')])],
                     [5, 5, set([('06029', '06037', '06075'), ('08031', '41051', '41067'), ('06059', '06071'), ('01073', '55079'), ('11001', '24510', '34013', '34017', '34039', '36005', '36047', '36059', '36061', '36081', '51013', '51760', '51840', '54009')])]]    
        
    suite = TestSuite()    
    
    for num_clusters, num_iterations, expected_county_tuple in kmeansdata_24:
        
        # build initial list of clusters for each test since mutation is allowed
        cluster_list = []
        for idx in range(len(data_24_table)):
            line = data_24_table[idx]
            cluster_list.append(Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))

        # compute answer
        Cluster_data = kmeans_clustering(cluster_list, num_clusters, num_iterations)
        County_tuples = set_of_county_tuples(Cluster_data)
        
        # Prepare test
        error_message = "Testing kmeans_custering on 24 county table, num_clusters = " + str(num_clusters)
        error_message += " num_iterations = " + str(num_iterations)
        error_message += "\nCounty tuples: " + str(County_tuples)
        error_message += "\nExpected county tuples: " + str(expected_county_tuple)
        suite.run_test(County_tuples == expected_county_tuple, True, error_message)   

    suite.report_results()
    
#tests
test_hierarchical()
test_kmeans()