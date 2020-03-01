import java.util.*;
import java.io.*;

class RandomMVNParams {
  // Set the ranges in which to calculate the multivariate normal parameters
  private int dimensions;
  private int muRange;
  private double muMin;
  private double sigma2Max;
  private boolean posDef = true;

  private double[] mu;
  private double[][] sigma;
  private double[][] lowerT;

  RandomMVNParams(int dimIn, int muRangeIn, double muMinIn, double
  sigma2MaxIn) {
    this.dimensions = dimIn;
    this.muRange = Math.abs(muRangeIn);
    this.muMin = muMinIn;
    this.sigma2Max = Math.abs(sigma2MaxIn);

    this.mu = new double[dimIn];
    this.sigma = new double[dimIn][dimIn];

    // Randomly generate the mean and variances of the multivariate normal
    Random rand = new Random();

    for (int dim = 0; dim < this.dimensions; dim++) {
      mu[dim] = muMin + rand.nextInt(muRange) + rand.nextDouble();
      sigma[dim][dim] = sigma2Max * rand.nextDouble();
    }

    // Generate covariances for the off-diagonals of the covariance matrix
    for (int rows = 1; rows < this.dimensions; rows++) {
      for (int cols = 0; cols < rows; cols++) {

        // Generate correlations between -0.7 and 0.7
        double rho = rand.nextDouble() - 0.7;
        double covariance = rho * Math.sqrt(sigma[rows][rows]) *
          Math.sqrt(sigma[cols][cols]);

        sigma[rows][cols] = covariance;
        sigma[cols][rows] = covariance;
      }
    }

    /* Use the method of Cholesky decomposition to generate a lower triangular
    martix from sigma, which can then be used to produce properly correlated
    points from the MVN distribution */
    this.lowerT = new double[dimIn][dimIn];

    for (int rows = 0; rows < this.dimensions; rows++) {
      for (int cols = 0; cols <= rows; cols++) {

        // Input the diagonal entries into lowerT
        if (rows == cols) {
          double offDiagSqrSum = 0;

          for (int iterRow = 0; iterRow < rows; iterRow++) {
            offDiagSqrSum += Math.pow(lowerT[rows][iterRow], 2);
          }

          lowerT[rows][rows] = Math.sqrt(sigma[rows][rows] - offDiagSqrSum);

          /* Flag if the randomly generated covariance matrix is not positive
          definite, in which case the Cholesky decomposition will not exist */
          if (Double.isNaN(lowerT[rows][rows])) {
            posDef = false;
          }
        }

        // Input the off-diagonal entries into lowerT
        else {
          double sum = 0;

          for (int iterRow = 0; iterRow < cols; iterRow++) {
            sum += lowerT[rows][iterRow] * lowerT[cols][iterRow];
          }

          lowerT[rows][cols] = (sigma[rows][cols] - sum) / lowerT[cols][cols];
        }
      }
    }
  }

  /* Method to generate multivariate normal points from the distribution as
  specified at construction of the object */
  public double[] generateMVN() {
    double[] MVNPoint = new double[dimensions];
    double[] standardGaussians = new double[dimensions];

    Random rand = new Random();

    for (int elem = 0; elem < dimensions; elem++) {
      standardGaussians[elem] = rand.nextGaussian();
    }

    for (int dim = 0; dim < dimensions; dim++) {
      double muComp = mu[dim];
      double varComp = 0;
      for (int elem = 0; elem < dimensions; elem++) {
        varComp += lowerT[dim][elem] * standardGaussians[elem];
      }

      MVNPoint[dim] = muComp + varComp;
    }
    return MVNPoint;
  }

  // Define methods to print the generated matrices to screen
  public double[] getMu() {
    return mu;
  }

  public void printMu() {
    System.out.println("mu:");
    for (int elem = 0; elem < mu.length; elem++) {
      System.out.printf("%8.2f %n", mu[elem]);
    }
    System.out.println();
  }

  public double[][] getSigma() {
    return sigma;
  }

  public void printSigma() {
    System.out.println("sigma:");
    for (int row = 0; row < sigma.length; row++) {
      for (int col = 0; col < sigma.length; col++) {
        System.out.printf("%8.2f", sigma[row][col]);
      }
      System.out.println();
    }
    System.out.println();
  }

  public double[][] getLowerT() {
    return lowerT;
  }

  public void printLowerT() {
    System.out.println("L:");
    for (int row = 0; row < lowerT.length; row++) {
      for (int col = 0; col < lowerT.length; col++) {
        System.out.printf("%8.2f", lowerT[row][col]);
      }
      System.out.println();
    }
    System.out.println();
  }

  public boolean isPosDef() {
    return posDef;
  }
}


/* Wrapper class for RandomMVNParams that ensures all covariance matrices are
positive definite */
class RandomMVN {
  public RandomMVNParams params;

  public RandomMVN(int dimIn, int muRangeIn, double muMinIn, double
  sigma2MaxIn) {
    RandomMVNParams initParams = new RandomMVNParams(dimIn, muRangeIn, muMinIn,
    sigma2MaxIn);

    while (initParams.isPosDef() == false) {
      initParams = new RandomMVNParams(dimIn, muRangeIn, muMinIn, sigma2MaxIn);
    }

    params = initParams;
  }

  public double[] generateMVN() {
    return params.generateMVN();
  }
}

public class MVNGenerator {
  // Define static constants
  static int NUM_CLUSTERS = 3;
  static int CLUSTER_DIMENSION = 4;
  static int CLUSTER_CENTRE_RANGE = 15;
  static int CLUSTER_CENTRE_MIN = -5;
  static int MAX_VARIANCE = 4;
  static int NUM_POINTS = 100;

  public static void main(String[] args) {
    // Generate random cluster centres
    RandomMVN[] clusterCentres = new RandomMVN[NUM_CLUSTERS];

    for (int i = 0; i < NUM_CLUSTERS; i++) {
      RandomMVN centre = new RandomMVN(CLUSTER_DIMENSION, CLUSTER_CENTRE_RANGE,
      CLUSTER_CENTRE_MIN, MAX_VARIANCE);
      clusterCentres[i] = centre;

      System.out.printf("Cluster %d %n%n", (i + 1));
      centre.params.printMu();
      centre.params.printSigma();
      centre.params.printLowerT();
      System.out.println();
    }

    // Generate n points randomly from one of the clusters
    double[][] randomPoints = new double[NUM_POINTS][CLUSTER_DIMENSION];

    Random rand = new Random();
    for (int i = 0; i < NUM_POINTS; i++) {
      int clustNum = rand.nextInt(NUM_CLUSTERS);
      double[] randPoint = clusterCentres[clustNum].generateMVN();
      randomPoints[i] = randPoint;
    }

    // Write data to a file as comma-separated values
    try {
      BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
        new FileOutputStream("points.txt"), "utf-8"));

      for (int j = 0; j < NUM_POINTS; j++) {
        double[] pointToWrite = randomPoints[j];

        writer.write(String.valueOf(pointToWrite[0]));

        for (int k = 1; k < CLUSTER_DIMENSION; k++) {
          writer.write(", ");
          writer.write(String.valueOf(pointToWrite[k]));
        }

        writer.newLine();
      }

      writer.flush();
      writer.close();
    }

    catch (IOException ioe) {
      ioe.printStackTrace();
    }
  }
}
