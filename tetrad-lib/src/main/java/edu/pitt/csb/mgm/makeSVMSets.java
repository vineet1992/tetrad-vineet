package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by vinee_000 on 1/22/2017.
 */
public class makeSVMSets {
    public static void main(String [] args) throws Exception
    {
        String dir = "Highest_Variance";
        String interest = "Subtype";

        File f = new File(dir + "/Classification_Sets");
        f.mkdir();

        if(dir.equals("Random"))
        {
            for(int j = 0; j < 5; j++) {
                for (int i = 0; i < 3; i++) {

                    PrintStream out = new PrintStream(dir + "/Classification_Sets/output_" + i + "_" + j + ".txt");
                    PrintStream out2 = new PrintStream(dir + "/Classification_Sets/outputPC_" + i + "_" + j + ".txt");
                    Graph g = GraphUtils.loadGraphTxt(new File(dir + "/Causal_Graphs/mgm_" + i + "_" + j + ".txt"));
                    List<Node> side = g.getAdjacentNodes(g.getNode(interest));
                    side.add(g.getNode("Subtype"));
                    DataSet d = MixedUtils.loadDataSet2(dir + "/Tetrad_Sets/dataset_" + i + "_" + j +  ".txt");
                    System.out.println(side);
                    out.println(d.subsetColumns(side));
                    g = GraphUtils.loadGraphTxt(new File(dir + "/Causal_Graphs/graph_" + i + "_" + j + ".txt"));
                    side = g.getAdjacentNodes(g.getNode(interest));
                    for(Node n: g.getChildren(g.getNode(interest)))
                    {
                        List<Node> temp =  g.getParents(n);
                        temp.remove(g.getNode(interest));
                        for(Node x: temp)
                        {
                            if(!side.contains(x)) {
                                side.add(x);
                            }
                        }
                    }
                    side.add(g.getNode(interest));
                    out2.println(d.subsetColumns(side));
                    out2.flush();
                    out2.close();
                    out.flush();
                    out.close();

                }
            }
        }
        else
        {
            for(int i = 0; i < 3; i ++) {

                PrintStream out = new PrintStream(dir + "/Classification_Sets/output_" + i + ".txt");
                PrintStream out2 = new PrintStream(dir + "/Classification_Sets/outputPC_" + i + ".txt");

                Graph g = GraphUtils.loadGraphTxt(new File(dir + "/Causal_Graphs/mgm_" + i + ".txt"));
                List<Node> side = g.getAdjacentNodes(g.getNode(interest));
                side.add(g.getNode(interest));
                DataSet d = MixedUtils.loadDataSet2(dir + "/Tetrad_Sets/dataset_" + i + ".txt");
                out.println(d.subsetColumns(side));

                g = GraphUtils.loadGraphTxt(new File(dir + "/Causal_Graphs/graph_" + i + ".txt"));
                side = g.getAdjacentNodes(g.getNode(interest));
                System.out.println(side);
                for(Node n: g.getChildren(g.getNode(interest)))
                {
                    List<Node> temp =  g.getParents(n);
                    temp.remove(g.getNode(interest));
                    for(Node x: temp)
                    {
                        if(!side.contains(x)) {
                            System.out.println(x);
                            side.add(x);
                        }
                    }

                }
                side.add(g.getNode(interest));
                out2.println(d.subsetColumns(side));
                out2.flush();
                out2.close();
                out.flush();
                out.close();

            }
        }







    }
}
