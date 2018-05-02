package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataWriter;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.*;
import edu.pitt.csb.KCI.KCI_Ind;

import java.io.*;

/**
 * Created by vinee_000 on 2/22/2016.
 */
public class testKCI {
    public static void main(String [] args) throws Exception {
        for (int fe = 0; fe < 5; fe++) {


            PrintStream out = new PrintStream("temp.txt");
            String filename = "DAG_" + fe + "_data_java";
            String filename2 = "DAG_" + fe + "_graph.txt";
            String extension = ".txt";
            Graph trueGraph = GraphUtils.loadGraphTxt(new File(filename2));
            DataSet data = MixedUtils.loadDataSet2(filename + extension);
            DataSet removeStuff = data.copy();
            //if i==0 do nothing
            int iii = 2;
            if (iii == 0) {

            } else if (iii == 1) {
                removeStuff.permuteRows();
                int[] removal = new int[4000];
                for (int las = 0; las < 4000; las++)
                    removal[las] = las;
                removeStuff.removeRows(removal);

            } else if (iii == 2) {
                removeStuff.permuteRows();
                int[] removal = new int[4500];
                for (int las = 0; las < 4500; las++)
                    removal[las] = las;
                removeStuff.removeRows(removal);
            } else {
                removeStuff.permuteRows();
                int[] removal = new int[4900];
                for (int las = 0; las < 4900; las++)
                    removal[las] = las;
            }
                    int labMult = 2;
                discretize(removeStuff);
                double[][] ges = new double[8][6];
                double[][] mges = new double[8][6];
                double[][] pcs = new double[8][6];
                double[][] mpcs = new double[8][6];
                PrintWriter pp = new PrintWriter(new File("curr_data.txt"));
                DataWriter.writeRectangularData(removeStuff, pp, '\t');
                DataSet d = MixedUtils.loadDataSet2("temp_data.txt");

                //System.out.println("True Graph " + p);
                // System.out.println("Data: " + data);

                double lab = .1 + .05*labMult;
                //Need to ask about this when you don't remove any variables
                double[] lambda = {lab, lab, lab};
                double tolerance = 1e-7; //convergeance tolerance
                int iterLimit = 500; //iteration limit
                MGM model = new MGM(removeStuff, lambda);
                long mgmTime = 0;
                mgmTime = System.currentTimeMillis();
                model.learnEdges(iterLimit);

                Graph mgmGraph = model.graphFromMGM();
                mgmTime = System.currentTimeMillis() - mgmTime;
                PrintStream testStuff = new PrintStream("DAG_test_info_" + fe + "_" + labMult + ".txt");
                IndependenceTest indTest = new IndTestMultinomialLogisticRegression(removeStuff, 0.05);
                IndependenceTest kci = new kci_matlab(removeStuff, 0.05, fe, trueGraph, testStuff);
                //Fgs2 g = new Fgs2(d);
                // Graph g1 = g.search();
                // PerformanceTests.graphComparison(g1, trueGraph, out);
                IndTestMultinomialLogisticRegression indind = (IndTestMultinomialLogisticRegression) indTest;
               // g.setInitialGraph(mgmGraph);
               // Graph g2 = g.search();
                //PerformanceTests.graphComparison(g2, trueGraph, out);

                PcStable p = new PcStable(indTest);
                p.setInitialGraph(mgmGraph);
                Graph g3 = p.search();
                //PerformanceTests.graphComparison(g3, trueGraph, out);
                KCI_Ind k2 = new KCI_Ind(removeStuff,removeStuff,.05);


                PcStable p2 = new PcStable(k2);
                long fullkciTime = System.currentTimeMillis();
                p2.search();
                fullkciTime = System.currentTimeMillis()-fullkciTime;
                PcStable p_prime = new PcStable(k2);
                long kciTime = System.currentTimeMillis();
                p_prime.setInitialGraph(mgmGraph);
                Graph g5 = p_prime.search();
                kciTime = System.currentTimeMillis()-kciTime;
              //  PerformanceTests.graphComparison(g5, trueGraph, out);
                System.out.println("Ind Test Called: " + k2.getTotal());
                testStuff.println(k2.getTotal());
                testStuff.println(kciTime);
                testStuff.println(fullkciTime);
                //p.setInitialGraph(mgmGraph);
                //Graph g4 = p.search();
                //PerformanceTests.graphComparison(g4, trueGraph, out);
                out.flush();
                out.close();
                testStuff.flush();
                testStuff.close();
                //run pc stable with kci independence by using matlab call

                BufferedReader b2 = new BufferedReader(new FileReader("temp.txt"));
                int count = 0;
                double[][] currConf = new double[8][6];
                while (b2.ready()) {
                    String line = b2.readLine();
                    if (line.contains("No Edge")) {
                        for (int k = 0; k < 8; k++) {
                            String[] row = b2.readLine().split(" ");
                            // System.out.println(Arrays.toString(row));
                            int currIndex = 0;
                            //6, 10,14,18,22,26,34
                            for (int m = 0; m < 6; m++) {

                                A:
                                while (true) {
                                    try {
                                        row[currIndex] = row[currIndex].trim();
                                        Double.parseDouble(row[currIndex]);
                                        break A;
                                    } catch (NumberFormatException e) {
                                        if (row[currIndex].equals("*"))
                                            break A;
                                        currIndex++;
                                    }
                                }
                                if (row[currIndex].equals("*"))
                                    row[currIndex] = "0";
                                currConf[k][m] = Double.parseDouble(row[currIndex]);
                                currIndex++;
                            }

                        }
                        if (count == 0) {
                            for (int ii = 0; ii < 8; ii++) {
                                for (int jj = 0; jj < 6; jj++) {
                                    ges[ii][jj] += currConf[ii][jj];
                                }
                            }


                        } else if (count == 1) {
                            for (int ii = 0; ii < 8; ii++) {
                                for (int jj = 0; jj < 6; jj++) {
                                    mges[ii][jj] += currConf[ii][jj];
                                }
                            }
                        } else if (count == 2) {
                            for (int ii = 0; ii < 8; ii++) {
                                for (int jj = 0; jj < 6; jj++) {
                                    mpcs[ii][jj] += currConf[ii][jj];
                                }
                            }
                        }
                        //save relevant statistics from array
                        count++;
                    }
                }//End whileb.ready
                b2.close();
                PrintStream output = new PrintStream("output" + fe + "_" + labMult+ "_DAG.txt");
               // output.println("Times Called: " + k2.timesCalled);
               // output.println("Changes: " + k2.changes);
               // output.println("Correct Changes: " + k2.correct);
                output.flush();
                output.println("MGM-GES");
                output.println();

                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 6; j++) {
                        output.print(ges[i][j] + "\t");

                    }
                    output.println();
                }
                output.println();
                output.println();
                output.println("MGM-Pcs-Regular");
                output.println();
                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 6; j++) {
                        output.print(mges[i][j] + "\t");

                    }
                    output.println();
                }
                output.println();
                output.println();
                output.println("MGM-PCS-KCI");
                output.println();
                output.println();
                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 6; j++) {
                        output.print(mpcs[i][j] + "\t");

                    }
                    output.println();
                }
                output.println();
                output.println();
                output.flush();
                output.close();

        }
    }
    public static double mean(double[] curr) {

        double sum = 0;
        int num = curr.length;
        for (int i = 0; i < num; i++)
            sum += curr[i];
        return sum / num;
    }

    public static double std(double[] curr, double mean) {
        double sum = 0;
        for (int i = 0; i < curr.length; i++) {
            sum += Math.pow((curr[i] - mean), 2);
        }
        sum = sum / curr.length;
        return Math.sqrt(sum);

    }

    public static void discretize(DataSet ds) throws Exception {
        PrintStream cu = new PrintStream("temp_data.txt");
        int nc = ds.getNumColumns();
        cu.println("/variables");
        for(int i = 0; i < nc;i++)
        {
            cu.println("X" + (i + 1) + ":0.0000 1.0000 2.0000");
        }
        cu.println("/data");
        for(int i = 0; i < nc;i++)
        {
            cu.print("X" + (i+1));
            if(i!=nc-1)
                cu.print("\t");
        }
        cu.println();
        double [][] fullSet =  ds.getDoubleData().toArray();

        int [][] data = new int[ds.getNumRows()][ds.getNumColumns()];
        for (int i = 0; i < ds.getNumColumns(); i++) {
            double b = fullSet[0][i];
            String b2 = Double.toString(b);
            boolean isDiscrete = false;
            if(b2.length() < 4)
                isDiscrete = true;

            if(isDiscrete)
            {
                for(int ii = 0; ii < ds.getNumRows();ii++)
                {

                    data[ii][i] = (int)fullSet[ii][i];

                }
            }
            else {
                double[] curr = new double[ds.getNumRows()];
                for (int j = 0; j < ds.getNumRows(); j++) {
                    curr[j] = fullSet[j] [i];
                }
                double mean = mean(curr);
                double std = std(curr, mean);
                //.67448*std + mean, mean - .67448*std = range
                int[] newCurr = new int[curr.length];
                for (int j = 0; j < curr.length; j++) {
                    if (curr[j] < mean - 1.5 * std)
                        data[j][i] = 0;
                    else if (curr[j] > mean + 1.5 * std) {
                        data[j][i] = 2;
                    } else
                        data[j][i] = 1;
                }
            }
        }

        for(int i = 0; i < ds.getNumRows();i++)
        {
            for(int j = 0; j < nc; j++)
            {
                if(data[i][j] == 0)
                    cu.print("0.0000");
                else if(data[i][j] == 1)
                    cu.print("1.0000");
                else
                    cu.print("2.0000");

                if(j!=nc-1)
                    cu.print("\t");
            }
            cu.println();
        }
        cu.flush();
        cu.close();
    }
}
