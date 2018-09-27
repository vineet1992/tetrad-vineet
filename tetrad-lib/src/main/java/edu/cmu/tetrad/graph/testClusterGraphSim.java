package edu.cmu.tetrad.graph;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 3/28/2018.
 */
public class testClusterGraphSim {
    public static void main(String [] args) throws Exception
    {
        int numNodes = 1000;
        int numComponents = 20;
        int numConnectedComponents = 10;
        int numEdges = 2000;
        boolean evenDistribution = false;
        boolean connected = true;
        int maxDegree = 10;
        int maxInDegree = 10;
        int maxOutDegree = 10;

        PrintStream out = new PrintStream("Graph_Test.txt");
        out.println(GraphUtils.randomGraphForwardEdgesClusters(numNodes,numComponents,numConnectedComponents,numEdges,evenDistribution,connected,maxDegree,maxInDegree,maxOutDegree));
        out.flush();
        out.close();
    }
}
