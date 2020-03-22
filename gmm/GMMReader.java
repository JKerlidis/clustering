import java.util.*;
import java.io.*;

// Define a class to calculate multivariate normal densities
class MVN {
  // Store the parameters of the multivariate normal
  int dim;
  double[] mu;
  double[][] sigma;
  /* In addition, store the lower triangular matrix that arises from a Cholesky
  decomposition of the covariance matrix, sigma */
  double[][] L;

  MVN(double[] muIn, double[][] sigmaIn) {
    dim = muIn.length;
    mu = muIn;
    sigma = sigmaIn;

    // Use Cholesky decomposition to generate the value of L
    L = new double[dim][dim];

    for (int row = 0; row < dim; row++) {
      for (int col = 0; col <= row; col++) {
        double sum = 0;

        for (int iterCol = 0; iterCol < col; iterCol++) {
          sum += L[row][iterCol] * L[col][iterCol];
        }

        if (row == col) {
          L[row][row] = Math.sqrt(sigma[row][row] - sum);
        }

        else {
          L[row][col] = (sigma[row][col] - sum) / L[col][col];
        }
      }
    }
  }

  /* Calculate the determinant of sigma as the sum of the square of the
  diagonal entries of L */
  double determinant() {
    double det = 0;

    for (int row = 0; row < this.dim; row++) {
      det += Math.pow(L[row][row], 2);
    }

    return det;
  }

  // Calculate the inverse of sigma by first calculating the inverse of L
  double[][] invSigma() {
    /* Separate L into a strictly lower triangular matrix and diagonal matrix,
    denoted as Ls and invD respectively */
    double[][] invD = new double[dim][dim];
    double[][] Ls = new double[dim][dim];

    for (int row = 0; row < this.dim; row++) {
      for (int col = 0; col < this.dim; col++) {

        if (row == col) {
          double element = this.L[row][row];
          invD[row][row] = 1/element;
        }

        else {
          Ls[row][col] = this.L[row][col];
        }
      }
    }

    /* We calculate the determinant based on the triangluar matrix identity:
    L^(-1) = (sum_{j=0}^{dim-1} (- D^(-1) %*% Ls)^j ) %*%  D^(-1) */
    double[][] sumPowersInvDLs = new double[dim][dim];
    double[][] powerInvDLs = new double[dim][dim];
    double[][] prevPowerInvDLs = new double[dim][dim];
    double[][] invDLs = new double[dim][dim];

    // Calculate invD %*% Ls, and initialise sumPowersInvDLs as I - invD %*% Ls
    for (int row = 0; row < this.dim; row++) {
      for(int col = 0; col < this.dim; col++) {

        if (row == col) {
          sumPowersInvDLs[row][row] = 1;
        }

        else if (row > col) {
          invDLs[row][col] = invD[row][row] * Ls[row][col];
          powerInvDLs[row][col] = invDLs[row][col];
          sumPowersInvDLs[row][col] = -invDLs[row][col];
        }
      }
    }

    // For each power of invD %*% Ls, add/subtract from sumPowersInvDLs
    for (int power = 2; power < this.dim; power++) {

      // Move the entries of powerInvDLs into prevPowerInvDLs
      for (int row = 0; row < this.dim; row++) {
        for (int col = 0; col < this.dim; col++) {
          prevPowerInvDLs[row][col] = powerInvDLs[row][col];
          powerInvDLs[row][col] = 0;
        }
      }

      // Calculate the next power of invDLs
      for (int row = 0; row < this.dim; row++) {
        for (int col = 0; col < this.dim; col++) {
          for (int iter = 0; iter < this.dim; iter++) {
            powerInvDLs[row][col] += prevPowerInvDLs[row][iter]
              * invDLs[iter][col];
          }
        }
      }

      for (int row = 0; row < this.dim; row++) {
        for (int col = 0; col <= row; col++) {
          if (power % 2 == 0) {
            sumPowersInvDLs[row][col] += powerInvDLs[row][col];
          }
          else {
            sumPowersInvDLs[row][col] -= powerInvDLs[row][col];
          }
        }
      }
    }

    // Multiply sumPowersInvDLs with invD to obtain invL
    double[][] invL = new double[dim][dim];

    for (int row = 0; row < this.dim; row++) {
      for (int col = 0; col < this.dim; col++) {
        for (int iter = 0; iter < this.dim; iter++) {
          invL[row][col] += sumPowersInvDLs[row][iter] * invD[iter][col];
        }
      }
    }

    // find invSigma using the identity invSigma = (invL)^T %*% invL
    double[][] invSig = new double[dim][dim];

    for (int row = 0; row < this.dim; row++) {
      for (int col = 0; col < this.dim; col++) {
        for (int iter = 0; iter < this.dim; iter++) {
          invSig[row][col] += invL[iter][row] * invL[iter][col];
        }
      }
    }

    return invSig;
  }

  // Calculate the density of a point under the specified MVN distribution
  double density(double[] point) {

    // Raise a (soft) error if the dimensions are wrong
    if (point.length != this.dim) {
      System.out.println("Dimension mismatch");
      return 0;
    }

    else {
      // Calculate invSigma using the invSigma method
      double[][] invSigma = this.invSigma();

      // Calculate the squared Mahalanobis distance
      double[] xLessMu = new double[dim];
      for (int iter = 0; iter < dim; iter++) {
        xLessMu[iter] = point[iter] - mu[iter];
      }

      double[] xLessMuInvSigma = new double[dim];
      for (int col = 0; col < dim; col++) {
        for (int iter = 0; iter < dim; iter++) {
          xLessMuInvSigma[col] += xLessMu[iter] * invSigma[iter][col];
        }
      }

      double mahalanobisSquare = 0;
      for (int iter = 0; iter < dim; iter++) {
        mahalanobisSquare += xLessMuInvSigma[iter] * xLessMu[iter];
      }

      // Use this to calculate the exponential term of the density function
      double densExp = Math.exp(-mahalanobisSquare / 2);

      // Calculate the normalising factor of the density
      double detSigma = this.determinant();
      double normFactor = Math.pow(2 * Math.PI, this.dim) * detSigma;
      normFactor = 1 / Math.sqrt(normFactor);

      // Calculate the density and return its value
      double dens = normFactor * densExp;
      return dens;
    }
  }
}


public class GMMReader {
  public static void main(String[] args) {
    double[] mu = {1, 1, 1};
    double[][] sigma = {{3, 0.5, 0.2}, {0.5, 2, -0.1}, {0.2, -0.1, 4}};

    MVN mvnDist = new MVN(mu, sigma);

    double[] x = {0, 1, 1};
    System.out.println(mvnDist.density(x));
  }
}
