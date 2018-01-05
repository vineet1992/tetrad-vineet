package edu.pitt.csb.mgm;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;

import java.io.File;
import java.io.PrintStream;

/**
 * Created by vinee_000 on 5/23/2017.
 */
public class graphStatistics {
    public static void main(String [] args) throws Exception
    {
        String folder = "Graphs_50_Nodes";
        int numGraphs = 10;
        PrintStream out = new PrintStream(folder + "/Graph_Statistics.txt");
        out.println("Graph_ID\tNum_Edges\tMax_Node_Degree\tAvg_Node_Degree\tStandard Deviation Node Degree");

        for(int i = 0; i < numGraphs;i++)
        {
            Graph temp = GraphUtils.loadGraphTxt(new File(folder + "/DAG_" + i + ".txt"));

            out.println(i+"\t" + temp.getNumEdges()+"\t" + maxNode(temp) + "\t" + avgNode(temp) + "\t" + stdNode(temp));
            out.flush();
        }
        out.close();

    }
    public static double maxNode(Graph g)
    {
        double max = 0;
        for(Node n: g.getNodes())
        {
            if(g.getAdjacentNodes(n).size()>max)
                max = g.getAdjacentNodes(n).size();
        }
        return max;
    }
    public static double avgNode(Graph g)
    {
        int size = g.getNodes().size();
        double avg = 0;
        for(Node n:g.getNodes())
        {
            avg+= g.getAdjacentNodes(n).size();
        }
        return avg/size;
    }
    public static double stdNode(Graph g)
    {
        double avg = avgNode(g);
        double sum = 0;
        for(Node n:g.getNodes())
        {
            sum+=Math.pow((g.getAdjacentNodes(n).size()-avg),2);
        }
        sum= sum/(g.getNodes().size()-1);
        return Math.sqrt(sum);
    }
}
