package edu.pitt.csb.Priors;

import edu.cmu.tetrad.data.DataSet;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

/**
 * Created by vinee_000 on 3/21/2019.
 */
public class PriorUtils {

    public static double[][] loadPriorSif(File prior,DataSet data)
    {
        try{
            double [][] result = new double[data.getNumColumns()][data.getNumColumns()];
            BufferedReader b = new BufferedReader(new FileReader(prior));
            while(b.ready())
            {
                String [] line = b.readLine().split("\t");
                if(data.getVariable(line[0])!=null && data.getVariable(line[1])!=null)
                {
                    int x = data.getColumn(data.getVariable(line[0]));
                    int y = data.getColumn(data.getVariable(line[1]));
                    result[x][y] = Double.parseDouble(line[2]);
                    result[y][x] = Double.parseDouble(line[2]);
                }


            }
            b.close();
            return result;

        }catch(Exception e)
        {
            System.err.println("Exception in reading file: " + prior + " as an SIF file");
            return null;
        }
    }


    public static double [][] loadPrior(File data, int numVariables) throws Exception
    {
        double [][] temp = new double[numVariables][numVariables];
        BufferedReader b = new BufferedReader(new FileReader(data));
        for(int i = 0; i < numVariables;i++)
        {
            String [] line = b.readLine().split("\t");
            if(line.length < numVariables)
                line = b.readLine().split("\t"); //To fix issues with header vs no header
            for(int j = 0; j < numVariables;j++)
            {
                temp[i][j] = Double.parseDouble(line[j]);
            }
        }
        return temp;
    }
}
