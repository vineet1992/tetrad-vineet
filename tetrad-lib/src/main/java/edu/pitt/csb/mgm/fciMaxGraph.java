package edu.pitt.csb.mgm;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.performance.PerformanceTests;

import java.io.File;

/**
 * Created by vinee_000 on 10/23/2016.
 */
public class fciMaxGraph {
    public static void main(String [] args)
    {
        for(int i = 0; i < 10; i ++)
        {
            Graph g = GraphUtils.loadGraphTxt(new File("out_" + i + "_graph.txt"));
            Graph h = GraphUtils.loadGraphTxt(new File("out_" + i + "_graph_joe.txt"));
            System.out.println(PerformanceTests.endpointMisclassification(g.getNodes(),g,h));
        }
    }
}
