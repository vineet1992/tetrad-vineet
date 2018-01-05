package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 4/12/2017.
 */
public class tcga_statin {
    public static void main(String [] args)throws Exception
    {
        DataSet d = MixedUtils.loadDataSet2("statin_dataset.txt");
        d.removeColumn(d.getVariable("ID"));
        d.removeColumn(d.getVariable("Statin_Sensitivity"));
        d.removeColumn(d.getVariable("Metastasis"));
        double [] lambda = {.55,.35,.35};
        MGM m = new MGM(d,lambda);
        m.learnEdges(1000);
        IndependenceTest i = new IndTestMultinomialAJ(d,.05);
        FciMaxP f = new FciMaxP(i);
        PrintStream out2 = new PrintStream("statin_brca_mgm_fixed_genes.txt");
        out2.println(m.graphFromMGM());
        out2.flush();
        out2.close();
        f.setInitialGraph(m.graphFromMGM());
        PrintStream  out = new PrintStream("statin_brca_max_fixed_genes.txt");
        out.println(f.search());
        out.flush();
        out.close();
    }
}
