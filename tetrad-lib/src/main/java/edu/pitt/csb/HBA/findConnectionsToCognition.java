package edu.pitt.csb.HBA;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import org.apache.commons.collections4.CollectionUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created by vinee_000 on 5/24/2018.
 */
public class findConnectionsToCognition {

    public static void main(String [] args)throws Exception
    {
        String graphName = "Imaging_Cognitive_Graph_No_DIFFSS_0.05";
        String graphFile = graphName + ".txt";
        String outFile = graphName + "_" + "Sub_Cognitive.txt";
        String outTableFile = graphName + "_Table.txt";
        PrintStream tableOut =new PrintStream(outTableFile);

        Graph g = GraphUtils.loadGraphTxt(new File(graphFile));
        Set<Edge> s = g.getEdges();

        //String [] cogVars = {"X_A1_5_Chg","X_ALDFR_Chg","X_ASDFR_Chg","FASTotal_Chg","CWI3CSSFinal_Chg","DTMT4_Chg","DIFFSS_Chg"};
        String [] cogVars = {"X_A1_5_Chg","X_ALDFR_Chg","X_ASDFR_Chg","FASTotal_Chg","CWI3CSSFinal_Chg","DTMT4_Chg"};
        //String [] cogVars = {"X_A1_5_Chg","X_ALDFR_Chg","X_ASDFR_Chg","FASTotal_Chg","CWI3CSSFinal_Chg","DIFFSS_Chg"};

        List<Node> temp = new ArrayList<Node>();
        for(String x:cogVars)
        {
            temp.add(g.getNode(x));
        }
        List<Node> adjNodes = g.getAdjacentNodes(g.getNode("group"));
        HashSet<Node> adjSet = new HashSet(adjNodes);
        List<Node> allNodes = new ArrayList<Node>(); //Nodes adjacent to group
        allNodes.addAll(adjNodes);
        allNodes.addAll(temp);
        allNodes.add(g.getNode("group"));
        //Check each cognitive variable neighborhood, and extract subnetworks of vars connected to both cog var and group
        for(int i = 0; i < cogVars.length;i++)
        {
            temp = g.getAdjacentNodes(g.getNode(cogVars[i]));
            HashSet<Node> tempSet = new HashSet(temp);
            ArrayList<Node> finalSet = (ArrayList<Node>)CollectionUtils.intersection(adjSet,tempSet);
            for(int j = 0; j < finalSet.size();j++)
                tableOut.println(cogVars[i] + "\t" + finalSet.get(j));
            allNodes.addAll(finalSet);
        }
        tableOut.flush();
        tableOut.close();
        try {
            PrintStream out = new PrintStream(outFile);
            out.println(g.subgraph(allNodes));
            out.flush();
            out.close();
        }
        catch(Exception e)
        {
            e.printStackTrace();
            System.exit(-1);
        }

    }
}
