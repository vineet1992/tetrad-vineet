package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DelimiterType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;

import java.io.PrintStream;
import java.util.ArrayList;

/**
 * Created by vinee_000 on 1/5/2018.
 */
public class runSteps {
    public static void main(String [] args)throws Exception
    {
        int index = 0;
        String directory = ".";
        String file = "Final_Data_Exercise.txt";
        String stabilityFile = "";
        String graphFile = "";
        int ns = 20;
        double g = 0.01;
        int numLambdas = 40;
        double lambdaLow = 0.05;
        double lambdaHigh = 0.95;
        boolean LOO = true;
        int maxCategoriesForDiscrete = 4;
        ArrayList<String> varsToRemove = new ArrayList<String>();
        while (index < args.length) {
            if (args[index].equals("-d"))
            {
                directory = args[index + 1];
                index += 2;
            }
            else if(args[index].equals("-f"))
            {
                file = args[index+1];
                index+=2;
            }
            else if(args[index].equals("-rv"))
            {
                while(index +1 < args.length && !args[index+1].contains("-"))
                {
                    varsToRemove.add(args[index+1]);
                    index++;
                }
                index++;
            }
            else if(args[index].equals("-stabOut"))
            {
                stabilityFile = args[index + 1];
                index+=2;
            }
            else if(args[index].equals("-graphOut"))
            {
                graphFile = args[index+1];
                index+=2;
            }
            else if(args[index].equals("-ns"))
            {
                ns = Integer.parseInt(args[index+1]);
                index+=2;

            }
            else if (args[index].equals("-g"))
            {
                g = Double.parseDouble(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-nl"))
            {
                numLambdas = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-lamLow"))
            {
                lambdaLow = Double.parseDouble(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-lamHigh"))
            {
                lambdaHigh = Double.parseDouble(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-LOO"))
            {
                LOO = Boolean.parseBoolean(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-mc"))
            {
                maxCategoriesForDiscrete = Integer.parseInt(args[index+1]);
                index+=2;
            }
        }
        DelimiterType d2 = DelimiterType.TAB;
        DataSet d = MixedUtils.loadDataSet2(directory + "/" + file,d2,maxCategoriesForDiscrete);
        for(int i = 0; i < varsToRemove.size();i++)
        {
            d.removeColumn(d.getVariable(varsToRemove.get(i)));
        }
        double [] lambda = new double[numLambdas];
        for(int i = 0; i < numLambdas;i++)
        {
            lambda[i] = lambdaLow+(lambdaHigh/numLambdas)*i;
        }
        if(stabilityFile.equals(""))
            stabilityFile = "Stability_" + file;
        if(graphFile.equals(""))
            graphFile = "Graph_" + file;

        for(String s: varsToRemove)
        {
            d.removeColumn(d.getVariable(s));
        }
        System.out.println(d);
        STEPS s = new STEPS(d,lambda,g,d.getNumRows(),LOO);
        s.setComputeStabs(true);
        Graph g2 = s.runStepsPar();
        double [][] stab = s.stabilities;
        PrintStream out = new PrintStream(directory + "/" + graphFile);
        out.println(g2);
        out.flush();
        out.close();
        IndependenceTest ii = new IndTestMultinomialAJ(d,0.05);
        FciMaxP p = new FciMaxP(ii);
        p.setInitialGraph(g2);
        System.out.println(p.search());
        out = new PrintStream(directory + "/" + stabilityFile);
        printStability(out,d,stab);
    }

    public static void printStability(PrintStream out, DataSet d, double [][] stab)
    {
        for(int i = 0; i < d.getNumColumns();i++)
        {
            out.print(d.getVariable(i).getName());
            if(i < d.getNumColumns()-1)
                out.print("\t");
            else
                out.println();
        }
        for(int i = 0; i < d.getNumColumns();i++)
        {
            out.print(d.getVariable(i).getName()+"\t");
            for(int j = 0; j < d.getNumColumns();j++)
            {
                out.print(stab[i][j]);
                if(j < d.getNumColumns()-1)
                    out.print("\t");
                else
                    out.println();
            }
        }

        out.flush();
        out.close();
    }
}
