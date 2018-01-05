package edu.pitt.csb.mgm;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Arrays;

/**
 * Created by vinee_000 on 10/21/2015.
 */
//gathers results for original experiment with varying hubs and timing
    //and produces precision and recall based upon the method of least specific = negative
public class gatherResults {
    public static void main(String[] args) throws Exception {

        String fullPath = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/Sedgewick_CausalMGM_data/Sedgewick_CausalMGM_data/simulated_data/data";
        double[][][][][] MGM = new double[6][5][50][3][14];        //Lambda, Hubs, Graph Num, Edge Type, Score type
        double[][][][] FCI = new double[5][50][3][14];
        double[][][][] MAX = new double[5][50][3][14];

        // 0 = lambda 0.15, Hubs straightforward, Graph num straightforward, Edge Type 0 = CC, 1 = DC, 2 = DD, Score Type 0 = TP, 1 = FP, 2 = FN, 3 = TN, 0 = ADJ, 1 = Orientation


        long[][][] MGMTime = new long[5][6][50];
        long[][][] MGMFciTime = new long[5][6][50];
        long[][] FCITime = new long[5][50];
        long[][] MaxTime = new long[5][50];
        for (int i = 0; i < 5; i++) //hubs
        {
            for (int j = 0; j < 6; j++) //lambda
            {
                double lambda = 0.1 + j * .05;
                for (int k = 0; k < 50; k++) //for each data/graph
                {

                    String filename = "SF_" + k + "_data_output_" + i + "_" + new DecimalFormat("##.##").format(lambda) + ".txt";
                    System.out.println(filename);
                    BufferedReader b;
                    try {
                        b = new BufferedReader(new FileReader(fullPath + "/" + filename));
                    } catch (Exception e) {
                        continue;
                    }
                    int count = 0;

                    while (b.ready()) {
                        String[] currLine = b.readLine().split(" ");
                        System.out.println(Arrays.toString(currLine));

                        if (count == 0) {
                            count++;
                            for (int eT = 0; eT < 3; eT++) {

                                String[] line = b.readLine().split("\t");
                                //System.out.println(Arrays.toString(line));
                                for (int curr = 1; curr < 9; curr++) {
                                        MGM[j][i][k][eT][curr - 1] = Double.parseDouble(line[curr]);
                                }
                                //(F6*I6-(G6*H6))/(SQRT((F6+G6)*(F6+H6)*(I6+G6)*(I6+H6)))
                                MGM[j][i][k][eT][8] = MGM[j][i][k][eT][4] / (MGM[j][i][k][eT][5] + MGM[j][i][k][eT][4]);
                                MGM[j][i][k][eT][9] = MGM[j][i][k][eT][4] / (MGM[j][i][k][eT][4] + MGM[j][i][k][eT][6]);
                                MGM[j][i][k][eT][10] = (MGM[j][i][k][eT][4]*MGM[j][i][k][eT][7] - MGM[j][i][k][eT][5]*MGM[j][i][k][eT][6]) / Math.sqrt((MGM[j][i][k][eT][4]+MGM[j][i][k][eT][5])*(MGM[j][i][k][eT][4]+MGM[j][i][k][eT][6])*(MGM[j][i][k][eT][6]+MGM[j][i][k][eT][7])*(MGM[j][i][k][eT][5]+MGM[j][i][k][eT][7]));
                                MGM[j][i][k][eT][11] = MGM[j][i][k][eT][0] / (MGM[j][i][k][eT][1] + MGM[j][i][k][eT][0]);
                                MGM[j][i][k][eT][12] = MGM[j][i][k][eT][0] / (MGM[j][i][k][eT][0] + MGM[j][i][k][eT][2]);
                                MGM[j][i][k][eT][13] = (MGM[j][i][k][eT][0]*MGM[j][i][k][eT][3] - MGM[j][i][k][eT][1]*MGM[j][i][k][eT][2]) / Math.sqrt((MGM[j][i][k][eT][0]+MGM[j][i][k][eT][1])*(MGM[j][i][k][eT][0]+MGM[j][i][k][eT][2])*(MGM[j][i][k][eT][2]+MGM[j][i][k][eT][3])*(MGM[j][i][k][eT][1]+MGM[j][i][k][eT][3]));
                                //Copy these for FCI and MAX
                                if(Double.isNaN(MGM[j][i][k][eT][10]))
                                    MGM[j][i][k][eT][10] = 0;
                                if(Double.isNaN(MGM[j][i][k][eT][13]))
                                    MGM[j][i][k][eT][13] = 0;
                            }
                            String[] timeLine = b.readLine().split(" ");
                            MGMTime[i][j][k] = Long.parseLong(timeLine[1]);
                            continue;
                        }
                        if (currLine[0].equals("MGM-FCI:")) {
                            MGMFciTime[i][j][k] = Long.parseLong(currLine[1]);
                        }
                        if (count == 1 && currLine[0].contains("Edge_Type")) {

                            count++;
                            for (int eT = 0; eT < 3; eT++) {
                                String[] line = b.readLine().split("\t");
                                for (int curr = 1; curr < 9; curr++) {
                                        FCI[i][k][eT][curr - 1]= Double.parseDouble(line[curr]);

                                }
                                FCI[i][k][eT][8] = FCI[i][k][eT][4] / (FCI[i][k][eT][5] + FCI[i][k][eT][4]);
                                FCI[i][k][eT][9] = FCI[i][k][eT][4] / (FCI[i][k][eT][4] + FCI[i][k][eT][6]);
                                FCI[i][k][eT][10] = (FCI[i][k][eT][4]*FCI[i][k][eT][7] - FCI[i][k][eT][5]*FCI[i][k][eT][6]) / Math.sqrt((FCI[i][k][eT][4]+FCI[i][k][eT][5])*(FCI[i][k][eT][4]+FCI[i][k][eT][6])*(FCI[i][k][eT][6]+FCI[i][k][eT][7])*(FCI[i][k][eT][5]+FCI[i][k][eT][7]));
                                FCI[i][k][eT][11] = FCI[i][k][eT][0] / (FCI[i][k][eT][1] + FCI[i][k][eT][0]);
                                FCI[i][k][eT][12] = FCI[i][k][eT][0] / (FCI[i][k][eT][0] + FCI[i][k][eT][2]);
                                FCI[i][k][eT][13] = (FCI[i][k][eT][0]*FCI[i][k][eT][3] - FCI[i][k][eT][1]*FCI[i][k][eT][2]) / Math.sqrt((FCI[i][k][eT][0]+FCI[i][k][eT][1])*(FCI[i][k][eT][0]+FCI[i][k][eT][2])*(FCI[i][k][eT][2]+FCI[i][k][eT][3])*(FCI[i][k][eT][1]+FCI[i][k][eT][3]));
                                //Copy these for FCI and MAX
                                if(Double.isNaN(FCI[i][k][eT][10]))
                                    FCI[i][k][eT][10] = 0;
                                if(Double.isNaN(FCI[i][k][eT][13]))
                                    FCI[i][k][eT][13] = 0;
                            }
                            FCITime[i][k] = Long.parseLong(b.readLine().split(" ")[1]);

                        }
                        else if (count == 2 && currLine[0].contains("Edge_Type")) {
                            count++;
                            for (int eT = 0; eT < 3; eT++) {
                                String[] line = b.readLine().split("\t");
                                for (int curr = 1; curr < 9; curr++) {

                                        MAX[i][k][eT][curr - 1] = Double.parseDouble(line[curr]);

                                }
                                MAX[i][k][eT][8] = MAX[i][k][eT][4] / (MAX[i][k][eT][5] + MAX[i][k][eT][4]);
                                MAX[i][k][eT][9] = MAX[i][k][eT][4] / (MAX[i][k][eT][4] + MAX[i][k][eT][6]);
                                MAX[i][k][eT][10] = (MAX[i][k][eT][4]*MAX[i][k][eT][7] - MAX[i][k][eT][5]*MAX[i][k][eT][6]) / Math.sqrt((MAX[i][k][eT][4]+MAX[i][k][eT][5])*(MAX[i][k][eT][4]+MAX[i][k][eT][6])*(MAX[i][k][eT][6]+MAX[i][k][eT][7])*(MAX[i][k][eT][5]+MAX[i][k][eT][7]));
                                MAX[i][k][eT][11] = MAX[i][k][eT][0] / (MAX[i][k][eT][1] + MAX[i][k][eT][0]);
                                MAX[i][k][eT][12] = MAX[i][k][eT][0] / (MAX[i][k][eT][0] + MAX[i][k][eT][2]);
                                MAX[i][k][eT][13] = (MAX[i][k][eT][0]*MAX[i][k][eT][3] - MAX[i][k][eT][1]*MAX[i][k][eT][2]) / Math.sqrt((MAX[i][k][eT][0]+MAX[i][k][eT][1])*(MAX[i][k][eT][0]+MAX[i][k][eT][2])*(MAX[i][k][eT][2]+MAX[i][k][eT][3])*(MAX[i][k][eT][1]+MAX[i][k][eT][3]));
                                //Copy these for FCI and MAX
                                if(Double.isNaN(MAX[i][k][eT][10]))
                                    MAX[i][k][eT][10] = 0;
                                if(Double.isNaN(MAX[i][k][eT][13]))
                                  MAX[i][k][eT][13] = 0;
                            }

                            MaxTime[i][k] = Long.parseLong(b.readLine().split(" ")[1]);

                        }
                    }
                }
            }
        }

        //need to write out the Matrices to a file
        for (int i = 0; i < 5; i++) {
            String outputFile = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/outputs/hubs_" + i + ".txt";
            String outputFile3 = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/outputs/hubs_" + i + "_Time.txt";
            PrintWriter p1 = new PrintWriter(new BufferedWriter(new FileWriter(outputFile)));
            PrintWriter p3 = new PrintWriter(new BufferedWriter(new FileWriter(outputFile3)));
            for (int type = 0; type < 3; type++) {
                if (type == 0) {
                    p1.println("MGM-FCI");
                    p3.println("MGM-FCI");
                }
                if (type == 1) {

                    p1.println("FCI");
                    p3.println("FCI");
                } else {
                    p1.println("FCI-MAX");
                    p3.println("FCI-MAX");
                }
                if (type != 0) {

                    for (int edge = 0; edge < 3; edge++) {

                        if (edge == 0) {
                            p1.println("Continuous-Continuous");

                        } else if (edge == 1) {
                            p1.println("Continuous-Discrete");

                        } else {
                            p1.println("Discrete-Discrete");

                        }
                        p1.println("Trial\tTPU\tFPU\tFNU\tTNU\tTPD\tFPD\tFND\tTND\tPrecision_Directed\tRecall_Directed\tMCC_Directed\tPrecision_Undirected\tRecall_Undirected\tMCC_Undirected");
                        p3.println("Trial\tTime");
                        double []fciSum = new double[14];
                        double []maxSum = new double[14];


                        for (int k = 0; k < 50; k++) {
                            p1.print(k+1 + "\t");
                            p3.print(k+1 + "\t");
                            for (int scoreType = 0; scoreType < 14; scoreType++) {
                                if (type == 1) {
                                    fciSum[scoreType] = fciSum[scoreType] + FCI[i][k][edge][scoreType];
                                    p1.print(FCI[i][k][edge][scoreType] + "\t");
                                } else {
                                    maxSum[scoreType] = maxSum[scoreType] + MAX[i][k][edge][scoreType];
                                    p1.print(MAX[i][k][edge][scoreType] + "\t");
                                }

                            }
                            if (type == 1) {
                                p3.println(FCITime[i][k]);
                            } else if (type == 2) {
                                p3.println(MaxTime[i][k]);
                            }
                            p1.println();


                        }
                        p1.print("Average\t");

                        for(int s = 0; s < 14; s++)
                        {
                            fciSum[s] = fciSum[s]/50;
                            maxSum[s] = maxSum[s]/50;
                            if(type==1)
                                p1.print(fciSum[s]+"\t");
                            else
                                p1.print(maxSum[s]+"\t");
                        }
                        p1.println();
                        p1.print("Standard Error\t");
                        double [] stdSum = new double[14];
                        if(type==1)
                        {

                            for(int s = 0; s < 14; s++) {
                                for (int k = 0; k < 50; k++) {

                                    stdSum[s] = stdSum[s] + Math.pow((FCI[i][k][edge][s] - fciSum[s]),2);

                                }
                                stdSum[s] = Math.sqrt(stdSum[s]/49);
                                p1.print(stdSum[s]+"\t");
                            }
                            p1.println();
                            p1.println();
                        }
                        else
                        {
                            for(int s = 0; s < 14; s++) {
                                for (int k = 0; k < 50; k++) {

                                    stdSum[s] = stdSum[s] + Math.pow((MAX[i][k][edge][s] - maxSum[s]),2);

                                }
                                stdSum[s] = Math.sqrt(stdSum[s]/49);
                                p1.print(stdSum[s]/Math.sqrt(50)+"\t");
                            }
                            p1.println();
                            p1.println();
                        }

                    }
                } else {
                    for (int lam = 0; lam < 6; lam++) {
                        double lambda = 0.1 + .05 * lam;
                        p1.println("Lambda = " + lambda);
                        p3.println("Lambda = " + lambda);
                        for (int edge = 0; edge < 3; edge++) {

                            if (edge == 0) {
                                p1.println("Continuous-Continuous");
                            } else if (edge == 1) {
                                p1.println("Continuous-Discrete");

                            } else {
                                p1.println("Discrete-Discrete");

                            }
                            p1.println("Trial\tTPU\tFPU\tFNU\tTNU\tTPD\tFPD\tFND\tTND\tPrecision_Directed\tRecall_Directed\tMCC_Directed\tPrecision_Undirected\tRecall_Undirected\tMCC_Undirected");
                            p3.println("Trial\tPlain MGM Time\tPlain FCI Time");
                            double[] sum = new double[14];
                            for (int k = 0; k < 50; k++) {
                                p1.print(k + "\t");
                                p3.print(k + "\t");
                                for (int scoreType = 0; scoreType < 14; scoreType++) {
                                    //l, hubs, k, edget, scoret
                                    p1.print(MGM[lam][i][k][edge][scoreType] + "\t");
                                    sum[scoreType] = sum[scoreType] + MGM[lam][i][k][edge][scoreType];

                                }

                                p1.println();
                                p3.print(MGMTime[i][lam][k] + "\t");
                                p3.println(MGMFciTime[i][lam][k]);

                            }
                            p1.print("Average\t");
                            for(int s = 0; s < 14; s++)
                            {
                                sum[s] = sum[s]/50;
                                p1.print(sum[s]+"\t");
                            }
                            p1.println();
                            p1.print("Standard Error\t");
                            double stdSum = 0;
                            for(int s = 0; s < 14; s++)
                            {
                                for(int k = 0; k < 50; k++)
                                {
                                    stdSum = stdSum + Math.pow((MGM[lam][i][k][edge][s]-sum[s]),2);
                                }
                                stdSum = Math.sqrt(stdSum/49);
                                stdSum = stdSum/Math.sqrt(50);
                                p1.print(stdSum + "\t");
                                stdSum = 0;
                            }
                            p1.println();
                            p1.println();
                            p1.println();
                        }
                    }
                }
            }
            p1.flush();
            p3.flush();
            p1.close();
            p3.close();
        }
    }
}
