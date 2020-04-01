import java.io.*;
import java.util.*;

// Define a class to hold n-dimensional sets of coordinates
final class Point {
  private final double[] coords;

  public Point(int size) {
    coords = new double[size];
  }

  public Point(double[] coordsIn) {
    coords = coordsIn;
  }

  public int dimensions() {
    return coords.length;
  }

  public double[] getCoords() {
    return coords;
  }

  public double getCoords(int i) {
    return coords[i];
  }

  public double distance(Point altPoint) {
    if (this.dimensions() != altPoint.dimensions()) {
      System.out.println("Dimension mismatch between points: returning 0");
      return 0;
    }
    else {
      double cumSquareDist = 0;
      for (int i = 0; i < this.dimensions(); i++) {
        double dimensionDist = altPoint.getCoords(i) - this.getCoords(i);
        cumSquareDist += Math.pow(dimensionDist, 2);
      }
      return Math.sqrt(cumSquareDist);
    }
  }
}

// Define a class to hold a Point together with an int
class Pair {
  private Point point;
  private int cluster;

  public Pair(Point pointIn, int clusterIn) {
    point = pointIn;
    cluster = clusterIn;
  }

  public Point getPoint() {
    return point;
  }

  public int getCluster() {
    return cluster;
  }

  public void setCluster(int newCluster) {
    cluster = newCluster;
  }

  public void setCluster(Point[] clusterCentres) {
    int numClusters = clusterCentres.length;
    double minDistance = point.distance(clusterCentres[0]);
    int minDistanceIndex = 0;

    for (int i = 1; i < numClusters; i++) {
      double distance = point.distance(clusterCentres[i]);

      if (distance < minDistance) {
        minDistance = distance;
        minDistanceIndex = i;
      }
    }
    cluster = minDistanceIndex + 1;
  }

  public void printClusterInfo() {
  System.out.println(Arrays.toString(point.getCoords()) + " -> Cluster "
  + cluster);
  }
}

/* Define a class that can randomly initialise cluster centres, and then
minimise the distance for points in that cluster to the centre to iterate
the values of those cluster centres */
class ClusterCentres {
  // k refers to the chosen number of clusters in the k-means algorithm
  private int k;
  private Point[] centres = new Point[k];
  private double lastIterDist = 100;

  // Randomly initialise cluster centres to points in the set
  public ClusterCentres(ArrayList<Pair> pairsIn, int kIn) {
    this.k = kIn;

    ArrayList<Integer> takenIndexes = new ArrayList<>();
    int maxN = pairsIn.size();
    Point[] initPointCentres = new Point[k];

    Random rand = new Random();

    while (takenIndexes.size() < k) {
      int candIndex = rand.nextInt(maxN);
      if (takenIndexes.contains(candIndex)) {
        continue;
      }
      else {
        takenIndexes.add(candIndex);
      }
    }

    for (int i = 0; i < takenIndexes.size(); i++) {
      int index = takenIndexes.get(i);
      Pair initPointPair = pairsIn.get(index);
      initPointCentres[i] = initPointPair.getPoint();
    }

    this.centres = initPointCentres;
  }

  public int getK() {
    return k;
  }

  public Point[] getCentres() {
    return centres;
  }

  public double getLastIterDist() {
    return lastIterDist;
  }

  public void printClusterCentreInfo() {
    System.out.println("Cluster centres:");
    for (Point centre : centres) {
      System.out.println(Arrays.toString(centre.getCoords()));
    }
    System.out.println();
  }

  /* Find the expected value of the points in each cluster, and update the
  cluster centre to that value */
  public void updateClusters(ArrayList<Pair> pairsIn) {
    int pointDim = pairsIn.get(0).getPoint().getCoords().length;
    Point[] newCentres = new Point[this.k];

    for (int i = 0; i < this.k; i++) {
      int noPointsCount = 0;
      double[] clusterMean = new double[pointDim];

      for (Pair pair : pairsIn) {
        if (pair.getCluster() == i + 1) {
          noPointsCount += 1;

          double[] pointCoords = pair.getPoint().getCoords();
          Arrays.parallelSetAll(clusterMean, index -> clusterMean[index]
          + pointCoords[index]);
        }
      }

      for (int ind = 0; ind < clusterMean.length; ind++) {
        clusterMean[ind] /= noPointsCount;
      }

      Point clusterMeanPoint = new Point(clusterMean);
      newCentres[i] = clusterMeanPoint;
    }

    // Find the distance the new centre has moved from the old
    double oldNewSquareDist = 0;
    for (int centreNo = 0; centreNo < this.k; centreNo++) {
      double oldNewDist = centres[centreNo].distance(newCentres[centreNo]);
      oldNewSquareDist += Math.pow(oldNewDist, 2);
    }

    this.lastIterDist = Math.sqrt(oldNewSquareDist);

    // Update the cluster centres to the new values
    this.centres = newCentres;
  }
}

public class PointsReader {
  // Define an array to hold the set of points used in the clustering
  public static ArrayList<Pair> pointClusters = new ArrayList<>();
  public static int k = 3;

  /* Method to convert a string containing comma-separated coordinates
  into a pair, with cluster initialised at 0 */
  static Pair stringToPoint(String line) {
    String[] tokens = line.split(",");
    double[] coordinates = new double[tokens.length];

    for (int i = 0; i < tokens.length; i++) {
      double coordDim = Double.parseDouble(tokens[i]);
      coordinates[i] = coordDim;
    }

    Point point = new Point(coordinates);
    Pair pointClusterInit = new Pair(point, 0);
    return pointClusterInit;
  }


  // Run the algorithm
  public static void main(String[] args) {
    // Set distance between iterations at which to stop the algorithm
    double STOP_ITER_DISTANCE = 0;

    try {
      // Read in the points from file
      FileReader pointsListFR = new FileReader("points.txt");
      BufferedReader pointsList = new BufferedReader(pointsListFR);

      for (String line; (line = pointsList.readLine()) != null; ) {
        Pair pointInit = stringToPoint(line);
        pointClusters.add(pointInit);
      }

      // Initialise the centre points of each cluster
      ClusterCentres centrePoints = new ClusterCentres(pointClusters, k);
      centrePoints.printClusterCentreInfo();

      // Iterate the algorithm
      int iterations = 1;

      while (centrePoints.getLastIterDist() > STOP_ITER_DISTANCE) {
        for (int i = 0; i < pointClusters.size(); i++) {
          Pair pointCluster = pointClusters.get(i);
          pointCluster.setCluster(centrePoints.getCentres());
          pointClusters.set(i, pointCluster);
        }

        centrePoints.updateClusters(pointClusters);
        centrePoints.printClusterCentreInfo();

        iterations += 1;
      }

      // Print out the final clusters corresponding to each of the points
      System.out.println("Number of iterations: " + iterations + "\n");
      for (Pair pointPair: pointClusters) {
        pointPair.printClusterInfo();
      }
    }
    catch (IOException ioe) {
      ioe.printStackTrace();
    }
  }
}
