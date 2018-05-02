package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;
import edu.pitt.csb.KCI.KCI_Ind;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 5/26/2016.
 */
public class kci_LGRC {
    public static void main(String [] args) throws Exception {
        DataSet d = MixedUtils.loadDataSet2("kci_data_final.txt");
        DataSet d2 = d.copy();
        IndependenceTest i = new IndTestMultinomialAJ(d,.05);
        PcStable p = new PcStable(i);
        IndependenceTest i2 = new KCI_Ind(d,d2,.05);
        PcStable p2 = new PcStable(i2);
        PrintStream out = new PrintStream("mult_out.txt");
        PrintStream out2 = new PrintStream("kci_out.txt");
        out.println(p.search());
        out2.println(p2.search());
        out.flush();
        out.close();
        out2.flush();
        out2.close();
    }
}
