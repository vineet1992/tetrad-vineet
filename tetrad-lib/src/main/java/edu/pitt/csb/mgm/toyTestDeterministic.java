package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.PcStable;

import java.io.File;

/**
 * Created by vinee_000 on 2/2/2017.
 */
public class toyTestDeterministic {

    public static void main(String [] args)
    {
        Graph g = GraphUtils.loadGraphTxt(new File("toy_example.txt"));
        Parameters p = new Parameters();
        p.setValue("sampleSize",200);
        MixedLeeHastieSimulation c = new MixedLeeHastieSimulation();
        c.setTrueGraph(g);
        c.simulate(p);
        PcStable p2 = new PcStable(new IndTestMultinomialAJ(c.getDataSet(0), .05));
        Graph est = p2.search();
        System.out.println(est);
    }
}

