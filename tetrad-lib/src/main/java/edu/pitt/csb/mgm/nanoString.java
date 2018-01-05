package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.Knowledge;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.*;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 2/27/2017.
 */
public class nanoString {
    public static void main(String [] args) throws Exception
    {
        DataSet d = MixedUtils.loadDataSet2("tetrad_set.txt");
        d.removeColumn(d.getVariable("IFNA"));
        d.removeColumn(d.getVariable("OR2T11"));
        System.out.println(d);
DataSet cont = MixedUtils.getContinousData(d);
        cont = DataUtils.getNonparanormalTransformed(cont);
        System.out.println(cont);
        for(int i = 0; i < cont.getNumColumns();i++)
        {
            for(int j = 0; j < cont.getNumRows();j++)
            {
                d.setDouble(j,i,cont.getDouble(j,i));
            }
        }

        Knowledge k = new Knowledge();
        for(Node n: d.getVariables())
        {
            k.addToTier(0,n.getName());
        }
        double [] lambda = {.2,.2,.2};
        MGM m = new MGM(d,lambda);
        m.learnEdges(1000);

        IndependenceTest i = new IndTestMultinomialAJ(d,.25);
        FciMaxP p = new FciMaxP(i);
        p.setKnowledge(k);
        p.setInitialGraph(m.graphFromMGM());
        PrintStream out = new PrintStream("causal_search.txt");
        out.println(p.search());
        out.flush();
        out.close();
        out = new PrintStream("mgm_search.txt");
        out.println(m.graphFromMGM());
        out.flush();
        out.close();

    }
}
