package edu.pitt.csb.mgm;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintStream;

/**
 * Created by vinee_000 on 7/19/2017.
 */
public class createLatexTablePriors {
    public static void main(String [] args)throws Exception
    {
        String [] algs = {"mgm_one_steps","oracle_one","STEPS","oracle","mgm_priors"};
        int numRuns = 20;
        double [] prior = {0.1,0.3,0.6};
        int [] numVars = {50,100};
        int [] samples = {500,1000,2000,3000};
        double [][][][][] result = new double[prior.length][numVars.length][samples.length][algs.length][4];
        for(int i = 0; i < prior.length;i++)
        {
            for(int j = 0; j < numVars.length;j++)
            {
                for(int k = 0; k < samples.length;k++)
                {
                    for(int m = 0; m < algs.length;m++) {
                        BufferedReader b = new BufferedReader(new FileReader(algs[m] + "_" + prior[i] + "_5_" + numVars[j] + "_" + samples[k] + "_10.txt"));
                        double [][] temp = new double[4][numRuns];
                        b.readLine();//eat the header
                        for(int r = 0; r < numRuns;r++)
                        {
                            String [] line = b.readLine().split("\t");
                            for(int type = 0; type < 4; type++)
                            {
                                temp[type][r] = Double.parseDouble(line[type+1]);
                            }
                            for(int type = 0; type < 4;type++) {
                                result[i][j][k][m][type] = nanmean(temp[type]);
                            }
                        }
                        b.close();
                    }
                }
            }
        }
        String [] tString = {"CC","CD","DD","All"};
        for(int type = 0; type < 4;type++) {
            PrintStream out = new PrintStream("results_"  + tString[type] + ".txt");
            //TODO Construct Latex Table from Data
        }

    }
    private static double nanmean(double [] x)
    {
        double sum = 0;
        int count = 0;
        for(int i = 0; i < x.length;i++)
        {
            if(!Double.isNaN(x[i]))
            {
                sum+= x[i];
                count++;
            }
        }
        if(count==0)
            return Double.NaN;
        else
            return sum/count;
    }
}
