package edu.pitt.csb.CSI;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.IndTestGSquare;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;

/**
 * Created by vinee_000 on 8/26/2017.
 */
public class exploreCSI {
    public static void main(String [] args) throws Exception
    {
        int samplesContext = 100;
        int samplesOutContext = 100;
        Graph g= new EdgeListGraph();
        g.addNode(new ContinuousVariable("X"));
        g.addNode(new ContinuousVariable("Y"));
        g.addNode(new ContinuousVariable("Z"));
        g.addNode(new ContinuousVariable("A1"));
        g.addNode(new ContinuousVariable("A2"));
        g.addNode(new DiscreteVariable("G"));

        g.addDirectedEdge(g.getNode("Y"),g.getNode("Z"));
        g.addDirectedEdge(g.getNode("A1"),g.getNode("A2"));
        g.addDirectedEdge(g.getNode("X"),g.getNode("Y"));
        MixedLeeHastieSimulation m = new MixedLeeHastieSimulation();
        m.setTrueGraph(g);
        Parameters p = new Parameters();
                p.setValue("sampleSize",samplesContext);
        m.simulate(p);
        DataSet g_0 = m.getDataSet(0);
        System.out.println(g_0);
        g.removeEdge(g.getEdge(g.getNode("X"),g.getNode("Y")));
        m = new MixedLeeHastieSimulation();
        m.setTrueGraph(g);
        p.setValue("sampleSize",samplesOutContext);
        m.simulate(p);
        DataSet g_1 = m.getDataSet(0);
        for(int i = 0; i < g_0.getNumRows();i++)
        {
            g_0.setInt(i,g_0.getColumn(g_0.getVariable("G")),0);
        }
        for(int i = 0; i < g_1.getNumRows();i++)
        {
            g_1.setInt(i,g_1.getColumn(g_1.getVariable("G")),1);
        }


        DataSet both = DataUtils.concatenate(g_0,g_1);
        System.out.println(both);
        IndependenceTest i = new IndTestMultinomialAJ(both,.08);
        //     System.out.println(i.isIndependent(i.getVariable("G"),i.getVariable("X"),i.getVariable("Y")));
        double [] lambda = {.2,.2,.2};
        //       MGM m = new MGM(d,lambda);
        // m.learnEdges(1000);
        //   System.out.println(m.graphFromMGM());
        PcStable p2 = new PcStable(i);
        System.out.println(p2.search());






    /*    DataSet d = MixedUtils.loadDataSet2("csi_categorical.txt", DelimiterType.TAB);
       // if(!d.isMixed())
        //    System.exit(0);
        int numSamples = 800;
        int [] subset = new int[numSamples];
        for(int j = 0; j < numSamples;j++)
            subset[j] = j;
        d.permuteRows();
      d = d.subsetRows(subset);*/

    }
}
