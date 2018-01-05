package edu.pitt.csb.CSI;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DelimiterType;
import edu.cmu.tetrad.search.IndTestGSquare;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;

/**
 * Created by vinee_000 on 8/26/2017.
 */
public class exploreCSI {
    public static void main(String [] args) throws Exception
    {
        DataSet d = MixedUtils.loadDataSet2("csi_categorical.txt", DelimiterType.TAB);
       // if(!d.isMixed())
        //    System.exit(0);
        int numSamples = 800;
        int [] subset = new int[numSamples];
        for(int j = 0; j < numSamples;j++)
            subset[j] = j;
        d.permuteRows();
      d = d.subsetRows(subset);
        IndependenceTest i = new IndTestMultinomialAJ(d,.08);
   //     System.out.println(i.isIndependent(i.getVariable("G"),i.getVariable("X"),i.getVariable("Y")));
        double [] lambda = {.2,.2,.2};
 //       MGM m = new MGM(d,lambda);
       // m.learnEdges(1000);
     //   System.out.println(m.graphFromMGM());
        PcStable p = new PcStable(i);
        System.out.println(p.search());
    }
}
