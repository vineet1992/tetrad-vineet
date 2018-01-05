package edu.pitt.csb.mgm;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;

import java.io.File;
import java.io.PrintStream;

/**
 * Created by vinee_000 on 1/23/2017.
 */
public class testDistances {
    public static void main(String [] args) throws Exception {
        PrintStream out = new PrintStream("distance_output.txt");
        out.println("Method_1\tMethod_2\tK\tDistance");
        int [] top = {50,100,500};
        for (int k = 0; k < 3; k++) {


            Graph g1 = GraphUtils.loadGraphTxt(new File("Highest_Variance/Causal_Graphs/graph_" + k + ".txt"));
            Graph g2 = GraphUtils.loadGraphTxt(new File("Balanced_Pref_Div_A=06_intense=25/Causal_Graphs/graph_" + k + ".txt"));
            Graph g3 = GraphUtils.loadGraphTxt(new File("Top_K/Causal_Graphs/graph_" + k + ".txt"));



            out.println("Highest_Variance\tBalanced_Pref_Div\t" + top[k] + "\t" + graphDistance.distance(g1,g2));
            out.println("Top_K\tBalanced_Pref_Div\t" + top[k] + "\t" + graphDistance.distance(g2,g3));
            out.println("Highest_Variance\tTop_K\t" + top[k] + "\t" + graphDistance.distance(g1,g3));


            String target = "Subtype";
            File f = new File(target + "_Intensity_Distributions");
            f.mkdir();
            graphDistance.intensityDegree(target + "_Intensity_Distributions/intensity_distribution_highest_variance_" + k + ".txt",g1,target);
            graphDistance.intensityDegree(target + "_Intensity_Distributions/intensity_distribution_balanced_" + k + ".txt",g2,target);
            graphDistance.intensityDegree(target + "_Intensity_Distributions/intensity_distribution_top_k_" + k + ".txt",g3,target);
        }
    }
}
