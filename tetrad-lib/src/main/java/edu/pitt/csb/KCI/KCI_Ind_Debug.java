package edu.pitt.csb.KCI;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.PcStable;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;

/**
 * Created by vinee_000 on 3/15/2016.
 */
public class KCI_Ind_Debug {
    public static void main(String[] args) throws Exception {
        PrintStream out = new PrintStream("result_problem_edges.txt");
        for (int i = 0; i < 1; i++) {
            PrintStream out2 = new PrintStream("problem_graph_" + i + ".txt");

            Graph g = GraphUtils.loadGraphTxt(new File("C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/Sedgewick_CausalMGM_data/Sedgewick_CausalMGM_data/simulated_data/dataAndNetworks/DAG_" + i + "_graph.txt"));
            DataSet ds1 = MixedUtils.loadDataSet2("C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/Sedgewick_CausalMGM_data/Sedgewick_CausalMGM_data/simulated_data/updated_data/DAG_" + i + "_data_java.txt");
            DataSet ds2 = MixedUtils.loadDataSet2("C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/Sedgewick_CausalMGM_data/Sedgewick_CausalMGM_data/simulated_data/updated_data/DAG_" + i + "_data.txt");

            int[] selectedRows = new int[4500];
            for (int j = 0; j < 4500; j++)
                selectedRows[j] = j;
            ds1.removeRows(selectedRows);
            ds2.removeRows(selectedRows);
            // System.out.println(ds.getNumRows());
            double alpha2 = 0.05;
            double alpha = 0.05;
            // KCI_Ind k = new KCI_Ind(ds1,ds2,alpha,g,2);
            //KCI_Ind k2 = new KCI_Ind(ds1,ds2,alpha,g,3);
            KCI_Ind k3 = new KCI_Ind(ds1, ds2, alpha, alpha2, g, Integer.MAX_VALUE);

            IndTestMultinomialAJ k4 = new IndTestMultinomialAJ(ds1, .05);
            // double [] lab = {.15,.15,.15};
            // MGM m  = new MGM(ds1,lab);
            // m.learnEdges(1000);
            //PcStable p = new PcStable(k);
            //PcStable p2 = new PcStable(k2);
            PcStable p3 = new PcStable(k3);
            PcStable p4 = new PcStable(k4);

            long currtime = System.currentTimeMillis();
            Graph g2 = p3.search();
            Graph g3 = p4.search();
            DecimalFormat d = new DecimalFormat("##.####");
            out2.println("Variable Pair?\tCorrect Conditioning Sets\tIncorrect Removal");
            for (String s : k3.justNames.keySet()) {

                out2.println(s + "\t" + k3.corrRemovals.get(s) + "\t" + (k3.corrP.get(s)));


            }
          /*  for(String s: k3.corrRemovals.keySet())
            {
                String x = k3.corrRemovals.get(s);
                String ps = k3.corrP.get(s);
                String[] ps2 = ps.split(",");
                String [] x2 = x.split(",");
                if(ps2.length != x2.length){
                    System.out.println("Error");
                    System.exit(0);
                }
                for(int j = 0; j < x2.length;j++)
                {
                    if(x2[j].length() > 0)
                     out2.println("1\t" + x2[j] + "\t" + ps2[j]);
                }

            }*/
            out2.flush();
            out2.close();
            out.println("The graph difference between KCI and AJ is:");
            //PerformanceTests.graphComparison(g2, g3, out);

            currtime = System.currentTimeMillis() - currtime;
           /* long currtime2 = System.currentTimeMillis();
            Graph g3 = p2.search();
            currtime2 = System.currentTimeMillis()-currtime2;
            long currtime3 = System.currentTimeMillis();
            Graph g4 = p3.search();
            currtime3 = System.currentTimeMillis() - currtime3;
            Graph g5 = p4.search();
             out.println("The set size of the incorrect removals were");
            out.println(k3.setSize);
            out.println("Removals Correct?");
            out.println(k3.corrRemovals);
*/

            out.println("Time:\t" + currtime);
            out.println("Conditional Test Results");
            out.println("Changes-CC:\t" + k3.changesCC);
            out.println("Correct-CC:\t" + k3.correctCC);
            out.println("Times Called-CC:\t" + k3.timesCalledCC);

            out.println("Changes-CD:\t" + k3.changesCD);
            out.println("Correct-CD:\t" + k3.correctCD);
            out.println("Times Called-CD:\t" + k3.timesCalledCD);

            out.println("Changes-DD:\t" + k3.changesDD);
            out.println("Correct-DD:\t" + k3.correctDD);
            out.println("Times Called-DD:\t" + k3.timesCalledDD);

            out.println("Unconditional Results");
            out.println("Changes-CC:\t" + k3.changesCCZ);
            out.println("Correct-CC:\t" + k3.correctCCZ);
            out.println("Times Called-CC:\t" + k3.timesCalledCCZ);

            out.println("Changes-CD:\t" + k3.changesCDZ);
            out.println("Correct-CD:\t" + k3.correctCDZ);
            out.println("Times Called-CD:\t" + k3.timesCalledCDZ);

            out.println("Changes-DD:\t" + k3.changesDDZ);
            out.println("Correct-DD:\t" + k3.correctDDZ);
            out.println("Times Called-DD:\t" + k3.timesCalledDDZ);
            out.println("KCI vs. True Graph");
           // PerformanceTests.graphComparison(g2, g, out);
            out.println("AJ vs. True Graph");
            //PerformanceTests.graphComparison(g3, g, out);
            // out.println("Cutoff limit - 4 vs. True Graph");
            // PerformanceTests.graphComparison(g3, g, out);
            // out.println("No cutoff limit vs. True Graph");
            // PerformanceTests.graphComparison(g4, g, out);
            // out.println("AJ test vs. True Graph");
            // PerformanceTests.graphComparison(g5, g, out);
            // out.println("Cutoff-limit 4 vs no cutoff");
            // PerformanceTests.graphComparison(g3,g4,out);
            // out.println("Cutoff-limit 2 vs. no cutoff");
            //  PerformanceTests.graphComparison(g2,g4,out);

            out.flush();
            out.close();
        }
    }
}