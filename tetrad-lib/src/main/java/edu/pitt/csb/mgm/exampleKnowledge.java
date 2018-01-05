package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 4/9/2017.
 */
public class exampleKnowledge {
    public static void main(String [] args) throws Exception
    {

        PrintStream out = new PrintStream("OUTPUT FILE NAME");
        DataSet d = MixedUtils.loadDataSet2("INPUT FILE NAME");
        double [] lambda = {.2,.2,.2}; // CC, CD, DD Lambda values
        MGM m = new MGM(d,lambda);
        m.learnEdges(1000);
        IKnowledge k  = new Knowledge();

        k.addToTier(0,"VARIABLES WITH NO CAUSES");
        k.addToTier(1,"VARIABLE NAMES THAT CANNOT CAUSE TIER 0");
        IndependenceTest i = new IndTestMultinomialAJ(d,.05);
        PcStable p = new PcStable(i);
        p.setKnowledge(k);
        p.setInitialGraph(m.graphFromMGM());
       out.println( p.search());

    }

}
