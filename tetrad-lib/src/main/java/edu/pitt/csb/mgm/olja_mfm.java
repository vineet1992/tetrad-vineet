package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by vinee_000 on 1/9/2017.
 */
public class olja_mfm {
    public static void main(String [] args) throws Exception
    {
        DataSet d = MixedUtils.loadDataSet2("tetrad_data.txt");
        d.removeColumn(d.getVariable("NEG_A"));
        d.removeColumn(d.getVariable("NEG_B"));
        d.removeColumn(d.getVariable("NEG_C"));
        d.removeColumn(d.getVariable("NEG_D"));
        d.removeColumn(d.getVariable("NEG_E"));
        d.removeColumn(d.getVariable("NEG_F"));
        d.removeColumn(d.getVariable("NEG_G"));
        d.removeColumn(d.getVariable("NEG_H"));
        d.removeColumn(d.getVariable("POS_A"));
        d.removeColumn(d.getVariable("POS_B"));
        d.removeColumn(d.getVariable("POS_C"));
        d.removeColumn(d.getVariable("POS_D"));
        d.removeColumn(d.getVariable("POS_E"));
        d.removeColumn(d.getVariable("POS_F"));
        d.removeColumn(d.getVariable("OR2T11"));
        d.removeColumn(d.getVariable("IFNA"));
        System.out.println(d);
        double [] lambda = {0.6,0.2,0.2};
        MGM m = new MGM(d,lambda);
        System.out.println("Init MGM");
        m.learnEdges(1000);
        System.out.println("Learned MGM");
        Graph g = m.graphFromMGM();
        IndependenceTest i = new IndTestMultinomialAJ(d,.05);
        FciMaxP f = new FciMaxP(i);
        f.setInitialGraph(g);
        PrintStream out = new PrintStream("max.txt");

        out.println(f.search());
        out.flush();
        out.close();
        out = new PrintStream("mgm.txt");
        out.println(g);
        out.flush();
        out.close();
    }
}
