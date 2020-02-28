# Implementation of the K-Means algorithm
This is an implementation of the k-means algorithm, for clustering n-dimensional points (where n >= 1) into k clusters (k >= 1). This algorithm is a "proof of concept" from first principles. 

The `RandomPoints.java` file is used to create x n-dimensional points originating from k clusters (parameters can be set as required). The output of this script is stored in the `points.txt` file - an example of which is provided in the repository. The `PointsReader.java` file implements the k-means algorithm, reading points from the `points.txt` file.
