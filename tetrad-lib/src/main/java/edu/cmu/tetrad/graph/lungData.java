package edu.cmu.tetrad.graph;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 11/17/2016.
 */
public class lungData {
    public static void main(String [] args)throws Exception{
        DataSet data = MixedUtils.loadDataSet2("data.txt");
        data.removeColumn(data.getVariable("Nodule_location"));
        data.removeColumn(data.getVariable("Death_of_LungCA"));
        data.removeColumn(data.getVariable("Death"));
        System.out.println(data);
        double [] lambda = {.4,.3,.25,.2,.1};
        STEPS s = new STEPS(data,lambda,.01,5);
       Graph g =s.runSteps();
        PrintStream steps = new PrintStream("steps.txt");
        steps.println(g);
        IndependenceTest i = new IndTestMultinomialAJ(data,.1);
        PcStable p = new PcStable(i);
        p.setInitialGraph(g);
        FciMaxP f = new FciMaxP(i);
        f.setInitialGraph(g);
        PrintStream out = new PrintStream("pcstable.txt");
        out.println(p.search());
        out.flush();
        out.close();
        PrintStream out2 = new PrintStream("max.txt");
        out2.println(f.search());
        out2.flush();
        out2.close();
    }
}
