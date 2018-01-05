package edu.pitt.csb.Olja_Cancer_Analysis;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DelimiterType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 7/27/2017.
 */
/*Runs steps on second batch of nanostring data, just plain right now,
Will also try averaging pooled biological repliacates, and will try incorporating priors on the suggested gene list

 */
public class STEPS_NS_2 {
    public static void main(String [] args)throws Exception
    {

        DelimiterType d2 = DelimiterType.TAB;
        DataSet d = MixedUtils.loadDataSet2("NPN_Normalized_data.txt",d2);
        System.out.println(d);
        d.removeColumn(d.getVariable("Sample"));
        d.removeColumn(d.getVariable("Response"));
        d.removeColumn(d.getVariable("Screen_IGG"));
        d.removeColumn(d.getVariable("Week12_IGG"));
      //  DataSet d3 = MixedUtils.completeCases(d);
        int ns = 20;
        double g = 0.01;
        int numLambdas = 40;
        double [] lambda = new double[numLambdas];
        for(int i = 0; i < numLambdas;i++)
        {
            lambda[i] = .1+(.7/numLambdas)*i;
        }
        STEPS s = new STEPS(d,lambda,g,d.getNumRows(),true);
        Graph g2 = s.runSteps();
        double [][] stab = s.stabilities;
        PrintStream out = new PrintStream("Cluster_Graph_With_Age_Sex_NPN.txt");
        out.println(g2);
        out.flush();
        out.close();
        out = new PrintStream("STEPS_Stabilities_Cluster.txt");
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
        IndependenceTest i = new IndTestMultinomialAJ(d,.2);
        FciMaxP f = new FciMaxP(i);
        f.setInitialGraph(g2);
        out = new PrintStream("Max_Cluster_Graph_With_Age_Sex_NPN.txt");
        out.println(f.search());
        out.flush();
        out.close();
        d = MixedUtils.loadDataSet2("NPN_Normalized_data.txt",d2);
        d.removeColumn(d.getVariable("Cluster_Group"));
        d.removeColumn(d.getVariable("Sample"));
        d.removeColumn(d.getVariable("Screen_IGG"));
        d.removeColumn(d.getVariable("Week12_IGG"));
        //d3 = MixedUtils.completeCases(d);
        s = new STEPS(d,lambda,g,d.getNumRows(),true);
        g2 = s.runSteps();
        out = new PrintStream("Response_Graph_With_Age_Sex_NPN.txt");
        out.println(g2);
        out.flush();
        out.close();
        out = new PrintStream("STEPS_Stabilities_Response.txt");
        stab = s.stabilities;
        for(int j = 0; j < d.getNumColumns();j++)
        {
            out.print(d.getVariable(j).getName());
            if(j < d.getNumColumns()-1)
                out.print("\t");
            else
                out.println();
        }
        for(int k = 0; k < d.getNumColumns();k++)
        {
            out.print(d.getVariable(k).getName()+"\t");
            for(int j = 0; j < d.getNumColumns();j++)
            {
                out.print(stab[k][j]);
                if(j < d.getNumColumns()-1)
                    out.print("\t");
                else
                    out.println();
            }
        }
        out.flush();
        out.close();
        i = new IndTestMultinomialAJ(d,.2);
        out = new PrintStream("Max_Response_Graph_With_Age_Sex_NPN.txt");
        f = new FciMaxP(i);
        f.setInitialGraph(g2);
        out.println(f.search());
        out.flush();
        out.close();

    }
}
