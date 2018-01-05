package edu.pitt.csb.mgm;

import java.io.*;
import java.text.DecimalFormat;
import java.util.Arrays;

/**
 * Created by vinee_000 on 11/17/2015.
 */
//Gathers results for confusion matrix experiment with non varying sample size
public class gatherResults2 {
    public static void main(String[] args) throws Exception {

        String fullPath2 = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/Sedgewick_CausalMGM_data/Sedgewick_CausalMGM_data/simulated_data/data";
        String fullPath = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/sampleSizeOutput/";
        for (int ss = 0; ss < 4; ss++) {


            double[][][][] mgm1 = new double[5][8][6][50];
            double[][][][] mgm2 = new double[5][8][6][50];
            double[][][][] fci = new double[5][8][6][50];
            double[][][][] max = new double[5][8][6][50];
            double[][][][] mgmax = new double[5][8][6][50];
            for (int i = 0; i < 50; i++) {
                for (int j = 0; j < 4; j++) {
                    System.out.println(i + " " + j);
                    String filename2 = "SF_" + i + "_data_output_" + j + "_" + 0.2 + "_pt.txt";
                    String filename = "SF_" + i + "_graph_output_" + j + "_" + ss + "_ss.txt";
                    BufferedReader b;
                    if(ss==2)
                        b = new BufferedReader(new FileReader(fullPath2 + "/" + filename2));
                    else
                        b = new BufferedReader(new FileReader(fullPath + "/" + filename));
                    int count = 0;
                    while (b.ready()) {
                        String line = b.readLine();
                        if (line.contains("No Edge")) {
                            for (int k = 0; k < 8; k++) {
                                String[] row = b.readLine().split(" ");
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
                                    if (count == 0)
                                        mgm1[j][k][m][i] = Double.parseDouble(row[currIndex]);
                                    else if (count == 1)
                                        fci[j][k][m][i] = Double.parseDouble(row[currIndex]);
                                    else if (count == 2)
                                        max[j][k][m][i] = Double.parseDouble(row[currIndex]);
                                    else
                                        mgmax[j][k][m][i] = Double.parseDouble(row[currIndex]);
                                    currIndex++;
                                }

                            }
                            count++;
                        }
                    }
                    b.close();


                }

            }
            DecimalFormat df = new DecimalFormat("##.###");
            String topLine = "\t---\to-o\to->\t-->\t<->\tNo Edge";
            String[] stuff = {"---", "o-o", "o->", "<-o", "-->", "<--", "<->", "No Edge"};
            for (int j = 0; j < 5; j++) {
                System.out.println("Computing hubs " + j);
                int tempss = 0;
                if(ss == 0)
                    tempss = 100;
                if (ss == 1)
                    tempss = 1000;
                if(ss==2)
                    tempss = 500;
                if(ss ==3 )
                    tempss = 100;
                PrintWriter p = new PrintWriter(new BufferedWriter(new FileWriter("C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/outputsSS" + tempss + "/hubs_" + j + ".txt")));
                p.println("MGM Lambda = 0.15");
                p.println();
                p.println();
                p.println(topLine);
                for (int i = 0; i < 8; i++) {
                    p.print(stuff[i] + "\t");
                    for (int k = 0; k < 6; k++) {
                        double mea = mean(mgm2,j,i,k);
                        double stddev = std(mgm2,j,i,k,mea);
                        p.print(df.format(mea) + "(" + df.format(stddev) + ")" + "\t");
                    }
                    p.println();
                }
                precRecall(mgm2,p,j);
                p.println();
                p.println("MGM Lambda = 0.2");
                p.println();
                p.println();
                p.println(topLine);
                for (int i = 0; i < 8; i++) {
                    p.print(stuff[i] + "\t");
                    for (int k = 0; k < 6; k++) {
                        double mea = mean(mgm1,j,i,k);
                        double stddev = std(mgm1,j,i,k,mea);
                        p.print(df.format(mea) + "(" + df.format(stddev) + ")" + "\t");
                    }
                    p.println();
                }
                precRecall(mgm1,p,j);
                p.println("FCI");
                p.println();
                p.println();
                p.println(topLine);
                for (int i = 0; i < 8; i++) {
                    p.print(stuff[i] + "\t");
                    for (int k = 0; k < 6; k++) {
                        double mea = mean(fci,j,i,k);
                        double stddev = std(fci,j,i,k,mea);
                        p.print(df.format(mea) + "(" + df.format(stddev) + ")" + "\t");
                    }
                    p.println();
                }
                precRecall(fci,p,j);
                p.println("FCI-MAX");
                p.println();
                p.println();
                p.println(topLine);
                for (int i = 0; i < 8; i++) {
                    p.print(stuff[i] + "\t");
                    for (int k = 0; k < 6; k++) {
                        double mea = mean(max,j,i,k);
                        double stddev = std(max,j,i,k,mea);
                        p.print(df.format(mea) + "(" + df.format(stddev) + ")" + "\t");
                    }
                    p.println();
                }
                precRecall(max,p,j);
                p.println("MGM-FCI-MAX");
                p.println();
                p.println();
                p.println(topLine);
                for (int i = 0; i < 8; i++) {
                    p.print(stuff[i] + "\t");
                    for (int k = 0; k < 6; k++) {
                        double mea = mean(mgmax,j,i,k);
                        double stddev = std(mgmax,j,i,k,mea);
                        p.print(df.format(mea) + "(" + df.format(stddev) + ")" + "\t");
                    }
                    p.println();
                }
                precRecall(mgmax,p,j);
                p.flush();
                p.close();
            }
        }
    }
    public static double mean(double[][][][] arr,int i,int j, int k ){
        double sum = 0;
        for(int m = 0; m < 50; m++)
        {
            sum = sum + arr[i][j][k][m];
        }
        return sum/50;
    }
    public static double std(double[][][][] arr,int i,int j, int k, double mean) {
            double num = 0;
        for(int m = 0; m < 50; m++)
        {
            num = num + Math.pow((arr[i][j][k][m] - mean),2);
        }
        num = num /49;
        return Math.sqrt(num);
    }
    public static double getMean(double[]arr)
    {
        double mean = 0;
        int count = arr.length;
        for(int i = 0; i < arr.length;i++) {
            if(arr[i] < 0)
                count--;
            else
                mean += arr[i];
        }
        return mean/count;
    }
    public static double getStd(double[]arr,double mea)
    {
        double num = 0;
        int count = arr.length;
        for(int m = 0; m < arr.length; m++)
        {
            if(arr[m] < 0)
                count--;
            else
             num+= Math.pow((arr[m]-mea),2);
        }
        num = num /(count-1);
        return Math.sqrt(num);
    }
    //TO INCLUDE ONLY EDGES THAT BOTH HAD IN COMMON, SIMPLY EXCLUDE THE 8TH ROW, and the 5th column from all loops ( j < 9, k < 6)
    public static void precRecall(double[][][][] arr, PrintWriter p,int i)
    {
        p.println();
        p.println("Recall Information");
        p.println("Edge Type\tRecall");
        String[] stuff = {"---", "o-o", "o->", "<-o", "-->", "<--", "<->", "No Edge"};

        double [][] edgeRecall = new double[4][50];
        double[] recall = new double[50];
        double [] acc = new double[50];
        //i = 0 -> 5 hubs, j = 0 -> 8 = rows, k = 0 -> 6 cols, m = 0 -> 50 trials
        for(int m = 0; m < 50; m++) {
            double[] probs = new double[8];
            double[] corr = new double[8];
            double totCorr = 0;
            double compTotal = 0;
            double recTot = 0;
            for(int j = 0; j < 8; j++)
            {
                recTot = recTot + arr[i][j][5][m];
            }
            for (int j = 0; j < 8; j++) {

                int curr = -1;
                if (j == 0 || j == 3 || j == 5 || j == 7) {
                    //for(int m = 0; m < 50; m++)
                    //{
                    for (int k = 0; k < 5; k++) {
                        probs[j] = probs[j] + arr[i][j][k][m];
                        if (j != 7)
                            compTotal += arr[i][j][k][m];
                    }
                    // }
                    continue;
                }
                if (j == 1) {
                    curr = 1;
                } else if (j == 2)
                    curr = 2;
                else if (j == 4)
                    curr = 3;
                else if (j == 6)
                    curr = 4;
                //p.print(stuff[j] + "\t");
                double correct = 0;
                double total = 0;
                //for(int m = 0; m < 50; m++)
                //{
                correct = correct + arr[i][j][curr][m];
                totCorr += correct;
                for (int k = 0; k < 5; k++)
                    total = total + arr[i][j][k][m];
                compTotal += total;
                //}
                probs[j] = total;
                corr[j] = correct;
                //p.println((correct / total));
                if(total==0)
                    edgeRecall[curr-1][m] = -1;
                else
                    edgeRecall[curr - 1][m] = correct / total;


            }

           /* double fullTotal = 0;
            for (int j = 0; j < probs.length - 1; j++) {
                fullTotal += probs[j];
            }
            double finalSum = 0;
            for (int j = 0; j < probs.length; j++) {
                finalSum = finalSum + corr[j] / fullTotal;

            }
            if(fullTotal == 0)
                recall[m] = -1;
            else
                recall[m] = finalSum;*/
            acc[m] = totCorr/compTotal;
            recall[m] = compTotal/(compTotal+recTot);
        }
        String [] otherStuff = {"o-o","o->","-->","<-->"};
        for(int ii = 0; ii < 4; ii++)
        {
            double [] temp = edgeRecall[ii];
            double mean = getMean(temp);
            double std = getStd(temp,mean);
            p.println(otherStuff[ii] + "\t" + mean + " (" + std + ")");

        }
        double fullMean = getMean(acc);
        double fullStd = getStd(acc, fullMean);
        double recMean = getMean(recall);
        double recStd = getStd(recall,recMean);
        p.println();
        p.println("Overall Accuracy: " + fullMean + " (" + fullStd + ")");
        p.println("Overall Recall: " + recMean + " (" + recStd + ") ");
        p.println();
        p.println("Precision Information");
        p.println("Edge Type\tPrecision");
        String[] stuff2 = {"---", "o-o", "o->", "-->", "<->", "No Edge"};
       double [] [] edgePrecision = new double[5][50];
        double[] precision = new double[50];

        for(int m = 0; m < 50; m++) {
            double[] probs = new double[6];
            double[] corr = new double[6];
            double precTotal = 0;
            double fullTotal = 0;
            int[] inds = {-1, 1, 2, 4, 6, 7};
            for (int k = 1; k < 6; k++) {
               // p.print(stuff2[k] + "\t");
                double correct = 0;
                double total = 0;
                // for(int m = 0; m < 50; m++)
                //{
                correct = correct + arr[i][inds[k]][k][m];
                for (int j = 0; j < 8; j++) {
                    total = total + arr[i][j][k][m];
                    if(k!=5)
                    {
                        fullTotal = fullTotal + arr[i][j][k][m];
                        if(j!=7)
                            precTotal = precTotal + arr[i][j][k][m];
                    }
                }
                //}
                probs[k] = total;
                corr[k] = correct;
                //p.println((correct / total));
                if(total==0)
                    edgePrecision[k-1][m] = -1;
                else
                    edgePrecision[k-1][m] = correct/total;

            }
           precision[m] = precTotal / fullTotal;
        }
        for(int ii = 0; ii <4;ii++ )
        {
            double[] temp = edgePrecision[ii];
            double mean = getMean(temp);
            double std = getStd(temp,mean);
            p.println(stuff2[ii+1]+"\t" + mean + " (" + std + ")");

        }
        double finalMean = getMean(precision);
        double finalStd = getStd(precision,finalMean);
        p.println();
        p.println("Overall Precision: " + finalMean + " (" + finalStd + ")");
        p.flush();

    }
}
