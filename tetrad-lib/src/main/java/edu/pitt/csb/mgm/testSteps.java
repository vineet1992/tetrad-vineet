package edu.pitt.csb.mgm;


import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.Knowledge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import jgpml.GaussianProcess;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 10/6/2016.
 */
public class testSteps {
    public static void main(String[]args) throws Exception
    {
        DataSet d = MixedUtils.loadDataSet2("genes_with_clinical.txt");
        PrintStream out = new PrintStream("mfm_statin_graph_separate_sparsity.txt");
        System.out.println(d);
        double [] m = {.55,.3,.3};
        MGM m2 = new MGM(d,m);
        System.out.println("Intialized");
        m2.learnEdges(1000);
        /*double [] lambda = {.6,.64,.55,.5,.4,.35};

        STEPS s = new STEPS(d,lambda,.01,3);
        Graph g = s.runSteps();*/
        Knowledge k = new Knowledge();
        k.addToTier(1,"Statin_Sensitivity");
        for(String x:d.getVariableNames())
        {
            if(!x.equals("Statin_Sensitivity"))
                k.addToTier(2,x);
        }
        Graph g= m2.graphFromMGM();
        IndependenceTest i = new IndTestMultinomialAJ(d,.05);
        FciMaxP f = new FciMaxP(i);
        f.setKnowledge(k);
        f.setInitialGraph(g);
       Graph output = f.search();
        out.println(output);
        out.flush();
        System.out.println(g.getAdjacentNodes(g.getNode("Statin_Sensitivity")));
        System.out.println(output.getAdjacentNodes(output.getNode("Statin_Sensitivity")));
        out.close();
    }
}
