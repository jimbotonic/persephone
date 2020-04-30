var2clustersall is a 846 x 3 matrix.
Each row corresponds to a different ID (0 to 845) in the bt_symmetric.csv dataset.
Each column corresponds to a different set of clusters.
Column 1 = weeknights (Monday to Friday, 0000 to 0600)
Column 2 = weekdays (Monday to Thursday, 0700 to 1700; Friday 0700 to 1200)
Column 3 = weekend evenings (Friday to Saturday, 1600 to midnight)

Each value is either 0 or positive.
0 = this ID was not put into a cluster in this set.
positive value = number of the cluster within that set.
(Cluster numbers are NOT numbered consecutively.)
IDs with the same cluster number in the same column are in the same cluster.
Cluster numbers in different columns refer to different clusters, even if the numbers are equal.  
For example, cluster 307 in column 1 is NOT the same as cluster 307 in column 2.

Summary statistics:
Weeknights: 45 clusters with 133 IDs
Weekdays:   56 clusters with 334 IDs
Weekends:   22 clusters with  49 IDs
