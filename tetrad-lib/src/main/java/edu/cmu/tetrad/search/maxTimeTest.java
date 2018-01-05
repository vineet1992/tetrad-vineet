package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.performance.PerformanceTests;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.File;

/**
 * Created by vinee_000 on 8/13/2016.
 */
public class maxTimeTest {
    public static void main(String [] args) throws Exception
    {
        DataSet d = MixedUtils.loadDataSet2("SF_42_data.txt");
        d.removeColumn(4);
        Graph g3 = GraphUtils.loadGraphTxt(new File("SF_42_graph.txt"));
        g3.getNode("X5").setNodeType(NodeType.LATENT);
        DagToPag p = new DagToPag(g3);
        Graph pag = p.convert();
        IndependenceTest i = new IndTestMultinomialAJ(d,.05);
        double [] l = {.15,.15,.15};
        MGM m = new MGM(d,l);
       m.learnEdges(1000);
        FciMax f = new FciMax(i);
       f.setInitialGraph(m.graphFromMGM());
        Graph g = f.search();
        //System.out.println(f.search());
        i = new IndTestMultinomialAJ(d,.05);
        FciMaxP f2 = new FciMaxP(i);
        f2.setInitialGraph(m.graphFromMGM());
        //System.out.println(f2.search());
        Graph g2 = f2.search();
       System.out.println(PerformanceTests.endpointMisclassification(d.getVariables(), g, pag));
        System.out.println(PerformanceTests.endpointMisclassification(d.getVariables(), g2, pag));

System.out.println("Original Max: " + g);
        System.out.println("Modfiied Max: " + g2);
        System.out.println("True Graph: " + g3);
    }
}
