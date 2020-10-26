# Cancer cluster analysis in the US
In this project, there will be the implementation and evaluation of two methods for clustering data, two methods for computing closest pairs and two methods for clustering data. Then these two clustering methods will be compared in terms of efficiency, automation, and quality.

## Provided data
As part of the analysis, clustering methods will be applied to several sets of 2D data that include information about lifetime cancer risk from air toxics. Each entry in the data set corresponds to a county in the USA (identified by a unique 5 digit string called a FIPS county code) and includes information on the total population of the county and its per-capita lifetime cancer risk due to air toxics. To aid in visualizing this data, the county-level data includes the (x, y) position of each county when overlaid on this map of the USA:
![USA map](https://github.com/rokobo/Cancer-Clustering-Algorithm-Analysis/blob/master/Data/USA_Counties.png?raw=true)
 The program draws this map of the USA and overlays a collection of circles whose radius and color represent the total population and lifetime cancer risk of the corresponding county. Choosing a threshold for the lifetime cancer risk (multiplied by 10^-5) and eliminating those counties whose cancer risk is below that threshold yields smaller data sets with 896, 290, and 111 counties, using 3.5, 4.5, and 5.5 tresholds respectively. These four data sets will be the primary test data for the clustering methods.

## The Cluster class
The class initializer `Cluster(FIPS_codes, horiz_center, vert_center, total_population, average_risk)` takes as input a set of county codes, the horizontal and vertical position of the cluster's center, the total population and averaged cancer risks for the cluster. The class definition supports methods for extracting these attributes of the cluster. The Cluster class also implements two methods:

* `distance(other_cluster)` - Computes the Euclidean distance between the centers of the clusters self and other_cluster.
* `merge_clusters(other_cluster)` - Mutates the cluster self to contain the union of the counties in self and other_cluster. The method updates the center of the mutated cluster using a weighted average of the centers of the two clusters based on their respective total populations. An updated cancer risk for the merged cluster is computed in a similar manner.

## Closest pair functions
The functions for computing closest pairs should work on lists of Cluster objects and compute distances between clusters using the distance method. The three functions are:

* `slow_closest_pair(cluster_list)` - Takes a list of Cluster objects and returns a closest pair where the pair is represented by the tuple `(dist, idx1, idx2)` with idx1 < idx2 where dist is the distance between the closest pair `cluster_list[idx1]` and `cluster_list[idx2]`. This function should be implemented from the brute-force closest pair method:
![SCP](https://github.com/rokobo/Cancer-Clustering-Algorithm-Analysis/blob/master/Data/Slow_closest_pair.png?raw=true)
* `fast_closest_pair(cluster_list)` - Takes a list of Cluster objects and returns a closest pair where the pair is represented by the tuple `(dist, idx1, idx2)` with idx1 < idx2 where dist is the distance between the closest pair `cluster_list[idx1]` and `cluster_list[idx2]`. This function should be implemented from the divide-and-conquer closest pair method:
![FCP](https://github.com/rokobo/Cancer-Clustering-Algorithm-Analysis/blob/master/Data/Fast_closest_pair.png?raw=true)
* `closest_pair_strip(cluster_list, horiz_center, half_width)` - Takes a list of Cluster objects and two floats horiz_center and half_width. horiz_center specifies the horizontal position of the center line for a vertical strip and half_width specifies the maximal distance of any point in the strip from the center line. This function should return a tuple corresponding to the closest pair of clusters that lie in the specified strip and should be implemented from the helper function:
![FCP](https://github.com/rokobo/Cancer-Clustering-Algorithm-Analysis/blob/master/Data/Closest_pair_strip.png?raw=true)

Sorting a list of clusters by the vertical (or horizontal) positions of the cluster centers was at some point needed. A sort by vertical position was done in a single line of Python using the sort method for lists by providing a key argument of the form:
`cluster_list.sort(key=lambda cluster: cluster.vert_center())`
## Tests
Clustering_method_tests.py contains tests used to verify all functions were correct before doing any analysis on the data.

## Clustering functions

* `hierarchical_clustering(cluster_list, num_clusters)` - Takes a list of Cluster objects and applies hierarchical clustering as described in the following pseudo-code: 
![HC](https://github.com/rokobo/Cancer-Clustering-Algorithm-Analysis/blob/master/Data/Hierarchical_clustering.jpg?raw=true)
This clustering process should proceed until num_clusters clusters remain. The function then returns this list of clusters.

* `kmeans_clustering(cluster_list, num_clusters, num_iterations)` - Takes a list of Cluster objects and applies k-means clustering as described in the following pseudo-code:
![KMC](https://github.com/rokobo/Cancer-Clustering-Algorithm-Analysis/blob/master/Data/K-Means_clustering.jpg?raw=true)
 This function should compute an initial list of clusters (line 2 in the pseudo-code) with the property that each cluster consists of a single county chosen from the set of the num_cluster counties with the largest populations. The function should then compute num_iterations of k-means clustering and return this resulting list of clusters.
# Analysis
The analysis part of this project will analyze the performance of the clustering methods on various subsets of our county-level cancer risk data set. In particular, we will compare these methods in three areas:
* Efficiency - Which method computes clusterings more efficiently?
* Automation - Which method requires less human supervision to generate reasonable clusterings?
* Quality - Which method generates clusterings with less error?

## Efficiency

### Question 1 

Write a function `gen_random_clusters(num_clusters)` that creates a list of clusters where each cluster in this list corresponds to one randomly generated point in the square with corners (±1,±1). Use this function and a Python timing library to compute the running times of the functions `slow_closest_pair` and `fast_closest_pair` for lists of clusters of size 2 to 200. Once you have computed the running times for both functions, plot the result as two curves combined in a single line plot. The horizontal axis for your plot should be the the number of initial clusters while the vertical axis should be the running time of the function in seconds.

### Question 2
Create an image of the 15 clusters generated by applying hierarchical clustering to the 3108 county cancer risk data set. 

### Question 3 
Create an image of the 15 clusters generated by applying 5 iterations of k-means clustering to the 3108 county cancer risk data set.

### Question 4 
Which clustering method is faster when the number of output clusters is either a small fixed number or a small fraction of the number of input clusters? Provide a short explanation in terms of the asymptotic running times of both methods. You should assume that `hierarchical_clustering` uses `fast_closest_pair` and that `k-means clustering` always uses a small fixed number of iterations.

## Automation

### Question 5 
Create an image of the 9 clusters generated by applying hierarchical clustering to the 111 county cancer risk data set.

### Question 6 
Create an image of the 9 clusters generated by applying 5 iterations of k-means clustering to the 111 county cancer risk data set. 

### Question 7 
The clusterings that you computed in Questions 5 and 6 illustrate that not all clusterings are equal. In particular, some clusterings are better than others. One way to make this concept more precise is to formulate a mathematical measure of the error associated with a cluster. Given a cluster *C*, its error is the sum of the squares of the distances from each county in the cluster to the cluster's center, weighted by each county's population. If <img src="https://render.githubusercontent.com/render/math?math=p_i"> is the position of the county and <img src="https://render.githubusercontent.com/render/math?math=w_i"> is its population, then the cluster's error is:

<img src="https://render.githubusercontent.com/render/math?math=error(C) = \sum_{p_i \in C} w_i (d_{p_{i}c})^2">

Where *c* is the center of the cluster *C*. The Cluster class includes a method `cluster_error(data_table)` that takes a Cluster object and the original data table associated with the counties in the cluster and computes the error associated with a given cluster. Given a list of clusters L, the distortion of the clustering is the sum of the errors associated with its clusters.

<img src="https://render.githubusercontent.com/render/math?math=distortion(L) = \sum_{C \in L} error(C)">


`compute_distortion(cluster_list)` -  takes a list of clusters and uses cluster_error to compute its distortion. Compute the distortions of the two clusterings in questions 5 and 6. Find values for the distortions (with at least four significant digits) for these two clusterings.

As a check on the correctness, the distortions associated with the 16 output clusters produced by hierarchical clustering and k-means clustering (with 5 iterations) on the 290 county data set are approximately 2.575×10^11 and 2.323×10^11, respectively.

### Question 8
Examine the clusterings generated in Questions 5 and 6. In particular, focus your attention on the number and shape of the clusters located on the west coast of the USA. Describe the difference between the shapes of the clusters produced by these two methods on the west coast of the USA. What caused one method to produce a clustering with a much higher distortion? 

### Question 9 
Based on your answer to Question 8, which method (hierarchical clustering or k-means clustering) requires less human supervision to produce clusterings with relatively low distortion?

## Quality

### Question 10 
Compute the distortion of the list of clusters produced by hierarchical clustering and k-means clustering (using 5 iterations) on the 111, 290, and 896 county data sets, respectively, where the number of output clusters ranges from 6 to 20 (inclusive). Once you have computed these distortions for both clustering methods, create three separate plots (one for each data set) that compare the distortion of the clusterings produced by both methods. Each plot should include two curves drawn as line plots. The horizontal axis for each plot should indicate the number of output clusters while the vertical axis should indicate the distortion associated with each output clustering. For each plot, include a title that indicates the data set used in creating the plots and a legend that distinguishes the two curves.

### Question 11 
For each data set (111, 290, and 896 counties), does one clustering method consistently produce lower distortion clusterings when the number of output clusters is in the range 6 to 20? Is so, indicate on which data set(s) one method is superior to the other.