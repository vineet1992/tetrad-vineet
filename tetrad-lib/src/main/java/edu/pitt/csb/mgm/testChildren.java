package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.Knowledge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by vinee_000 on 11/28/2016.
 */
public class testChildren {
    public static void main(String [] args) throws Exception {
        DataSet d = MixedUtils.loadDataSet2("just_data_with_clusters.txt");
        System.out.println(d);

        /*double [] lambda = {.5,.4,.3,.2};
        STEPS s = new STEPS(d,lambda,.01,15);
       Graph g = s.runSteps();*/
        double [] lambda = {.55,.3,.5};
        MGM m = new MGM(d,lambda);
        m.learnEdges(1000);
        Graph g = m.graphFromMGM();
        IndependenceTest i = new IndTestMultinomialAJ(d, .05);
        FciMaxP f = new FciMaxP(i);
        f.setInitialGraph(g);
        PrintStream out = new PrintStream("max_clusters.txt");
        out.println(f.search());
        out.flush();
        out.close();

    }
}
