package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;

import java.io.PrintStream;
import java.util.ArrayList;

/**
 * Created by vinee_000 on 1/15/2017.
 */
public class Pref_Div_Causal {
    public static void main(String [] args) throws Exception
    {
        boolean random = false;

        if(random) {
            for (int i = 2; i < 3; i++) {

                for (int j = 0; j < 5; j++) {
                    DataSet d = MixedUtils.loadDataSet2("Tetrad_Sets/dataset_" + i + "_" + j + ".txt");
                    ArrayList<Node> bad = new ArrayList<Node>();
                    for (int ii = 0; ii < d.getNumColumns(); ii++) {
                        boolean good = false;
                        for (int jj = 0; jj < d.getNumRows(); jj++) {
                            if (d.getDouble(jj, ii) != 0)
                                good = true;

                        }
                        if (!good)
                            bad.add(d.getVariable(ii));
                    }

                    for (Node x : bad) {
                        System.out.println(x);
                        d.removeColumn(x);
                    }
                    PrintStream out2 = new PrintStream("Causal_Graphs/graph_" + i + "_" + j + ".txt");
                    PrintStream out = new PrintStream("Causal_Graphs/mgm_" + i + "_" + j + ".txt");
                    double[] lambda = {0.55, 0.35, 0.35};
                    MGM m = new MGM(d, lambda);
                    System.out.println("Init_MGM");
                    m.learnEdges(1000);
                    System.out.println("Learned_MGM");
                    Graph g = m.graphFromMGM();
                    IndependenceTest ii = new IndTestMultinomialAJ(d, .05);
                    PcStable f = new PcStable(ii);
                    f.setDepth(1);
                    f.setInitialGraph(g);
                    Graph h = f.search();
                    out.println(g);
                    out2.println(h);
                    out.flush();
                    out2.flush();
                    out.close();
                    out2.close();
                }
            }
        }
        else {

            for (int i = 0; i < 3; i++) {
                DataSet d = MixedUtils.loadDataSet2("Tetrad_Sets/dataset_" + i + ".txt");
                ArrayList<Node> bad = new ArrayList<Node>();
                for (int ii = 0; ii < d.getNumColumns(); ii++) {
                    boolean good = false;
                    for (int jj = 0; jj < d.getNumRows(); jj++) {
                        if (d.getDouble(jj, ii) != 0)
                            good = true;

                    }
                    if (!good)
                        bad.add(d.getVariable(ii));
                }

                for (Node x : bad) {
                    System.out.println(x);
                    d.removeColumn(x);
                }
                PrintStream out2 = new PrintStream("Causal_Graphs/graph_" + i + ".txt");
                PrintStream out = new PrintStream("Causal_Graphs/mgm_" + i + ".txt");
                double[] lambda = {0.55, 0.35, 0.35};
                MGM m = new MGM(d, lambda);
                System.out.println("Init_MGM");
                m.learnEdges(1000);
                System.out.println("Learned_MGM");
                Graph g = m.graphFromMGM();
                IndependenceTest ii = new IndTestMultinomialAJ(d, .05);
                PcStable f = new PcStable(ii);
                f.setDepth(1);
                f.setInitialGraph(g);
                Graph h = f.search();
                out.println(g);
                out2.println(h);
                out.flush();
                out2.flush();
                out.close();
                out2.close();
            }
        }
    }
}
