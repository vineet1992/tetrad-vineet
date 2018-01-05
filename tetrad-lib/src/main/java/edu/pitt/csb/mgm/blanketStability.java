package edu.pitt.csb.mgm;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;

import java.io.File;
import java.io.PrintStream;
import java.util.List;

/**
 * Created by vinee_000 on 1/27/2017.
 */
public class blanketStability {
public static void main(String [] args)throws Exception {
    PrintStream out = new PrintStream("tumor_blanket_stability.txt");
    String [] type = {"Highest_Variance","Balanced_Pref_Div_A=06_intense=25","Top_K"};
    String var = "Subtype";
    int [] top = {50, 100, 500};
    out.println("Type\tPercent_100_Blanket\tPercent_500_Blanket\tPercent_100_All\tPercent_500_All");
    for (int k = 0; k < 3; k++) {
System.out.println(type[k]);

        Graph g1 = GraphUtils.loadGraphTxt(new File(type[k] + "/Causal_Graphs/graph_0.txt"));

        Graph g2 = GraphUtils.loadGraphTxt(new File(type[k] + "/Causal_Graphs/graph_1.txt"));
        Graph g3 = GraphUtils.loadGraphTxt(new File(type[k] + "/Causal_Graphs/graph_2.txt"));

        double [] x = graphDistance.blanketStable(g1,g2,g3,var);

        out.println(type[k] + "\t" + x[0] + "\t" + x[1] + "\t" + x[2] + "\t" + x[3]);
        out.flush();

    }
    out.close();
}
}
