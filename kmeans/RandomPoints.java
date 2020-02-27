import java.util.*;
import java.io.*;

public class RandomPoints {
  static ArrayList<Double[]> numbers = new ArrayList<>();

  static int NUM_CLUSTERS = 3;
  static int CLUSTER_DIMENSION = 3;
  static int CLUSTER_CENTRE_RANGE = 15;
  static int CLUSTER_CENTRE_MIN = -5;
  static int NUM_POINTS = 60;

  public static void main(String[] args) {
    Random rand = new Random();
    Double[][] clusterCentres = new Double[NUM_CLUSTERS][CLUSTER_DIMENSION];

    for (int clustNum = 0; clustNum < NUM_CLUSTERS; clustNum++) {
      for (int clustDim = 0; clustDim < CLUSTER_DIMENSION; clustDim++) {
        clusterCentres[clustNum][clustDim] = rand.nextInt(CLUSTER_CENTRE_RANGE)
          + CLUSTER_CENTRE_MIN + rand.nextDouble();
      }
    }

    for (int i = 0; i < NUM_POINTS; i++) {
      Double[] pointValues = new Double[CLUSTER_DIMENSION];
      int clusterNum = rand.nextInt(NUM_CLUSTERS);

      for (int j = 0; j < CLUSTER_DIMENSION; j++) {
        Double pointValue = clusterCentres[clusterNum][j]
          + rand.nextGaussian();
        pointValues[j] = pointValue;
      }

      numbers.add(pointValues);
    }

    try {
      BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
        new FileOutputStream("points.txt"), "utf-8"));

      for (int k = 0; k < numbers.size(); k++) {
        Double[] pointToWrite = numbers.get(k);

        writer.write(String.valueOf(pointToWrite[0]));

        for (int l = 1; l < CLUSTER_DIMENSION; l++) {
          writer.write(", ");
          writer.write(String.valueOf(pointToWrite[l]));
        }

        writer.newLine();
      }

      writer.flush();
      writer.close();

    } catch (IOException ioe) {
        ioe.printStackTrace();
    }
  }
}
