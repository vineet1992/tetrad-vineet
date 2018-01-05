package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 3/3/2017.
 */
public class ExampleScript {
    public static void main(String [] args) throws Exception
    {
        DataSet d = MixedUtils.loadData(".","FILENAME");
        double [] lambda = {0.5,0.5,0.5}; //specify the lambda parameters Continuous-Continuous, Continuous-Discrete, Discrete-Discrete
        MGM m = new MGM(d,lambda);
        m.learnEdges(1000);
        IndependenceTest i = new IndTestMultinomialAJ(d,.05); //Change Alpha parameter here
        PcStable p = new PcStable(i);
        p.setInitialGraph(m.graphFromMGM());
        PrintStream out = new PrintStream("output.txt");//Change output filename here
        out.println(p.search());
        out.flush();
        out.close();


    }
}
