package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;

import static edu.pitt.csb.Pref_Div.Functions.computeDistance;

/**
 * Created by vinee_000 on 7/20/2017.
 */

/*
Goal is to determine if dissimilarity scores between genes are related to edges showing up
We will examine the causal graphs produced from all methods, and loop through edges to see if dissim relates to edge appearence

 */
public class relateIntensityToEdges {
    public static void main(String [] args) throws Exception
    {
        int numSub = 10;
        String directory = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Current Projects/Pref_Div/";
        ArrayList<Gene> g = Functions.loadGeneData(directory,null,null,true);
        System.out.println("...Loaded Gene Data");
        DataSet d = MixedUtils.loadDataSet2(directory + "/Top_K_Subtype/Tetrad_Sets/dataset_1.txt");
        d.removeColumn(d.getVariable("Race"));
        double [] lambda = new double[40];
        for(int i = 0; i < lambda.length;i++)
            lambda[i] = i*(0.7/lambda.length) + 0.1;
        STEPS s = new STEPS(d,lambda,0.01,numSub);

        //For the optimal lambda, what percentage of subsampled graphs did the edge appear in?
        int[][]adj = s.runStepsPD();
        Graph g2 = s.pdGraph;
        System.out.println("...Completed STEPS");
        PrintStream out = new PrintStream("similarity_result_1.txt");
        out.println("Stability\tDD\tCD\tFamily\tCorr");
      /*  for(Edge e: g2.getEdges())
        {
            Gene one = find(g,e.getNode1().getName());
            Gene two = find(g,e.getNode2().getName());
            if(one==null || two==null)
                continue;
           double [] result = Functions.computeDistance(one,two);
            int x = d.getColumn(d.getVariable(e.getNode1().getName()));
            int y = d.getColumn(d.getVariable(e.getNode2().getName()));
            out.println(adj[x][y]/(double)numSub + "\t" + result[0] + "\t" + result[1] + "\t" + result[2] + "\t" + result[3]);
            out.flush();
        }*/
        int samples = 10000;
        int count = 0;
        Random rand = new Random();
       A:while(count < samples)
       {
           int x = rand.nextInt(d.getNumColumns());
           int y = rand.nextInt(d.getNumColumns());
          // if(g2.getEdge(g2.getNode(d.getVariable(x).getName()),g2.getNode(d.getVariable(y).getName()))==null)
           //{
               Gene one = find(g,d.getVariable(x).getName());
               Gene two = find(g,d.getVariable(y).getName());
               if(one==null||two==null)
                   continue A;
               double [] result = Functions.computeDistance(one,two);
               out.println(adj[x][y]/(double)numSub + "\t" + result[0] + "\t" + result[1] + "\t" + result[2] + "\t" + result[3]);
               out.flush();
               count++;
           //}


       }
        out.close();
    }
    private static Gene find(ArrayList<Gene> g, String x)
    {
        for(int i = 0; i < g.size();i++)
        {
            if(g.get(i).symbol.equals(x))
                return g.get(i);
        }
        return null;
    }
}
