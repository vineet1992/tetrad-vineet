package edu.pitt.csb.Olja_Cancer_Analysis;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DelimiterType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 8/19/2017.
 */
public class STEPS_Microarray_Original {
    public static void main(String [] args)throws Exception
    {
        DelimiterType d2 = DelimiterType.TAB;
        DataSet d = MixedUtils.loadDataSet2("eset_selected_genes_only_TNFMUC.txt",d2);
        System.out.println(d);
        d.removeColumn(d.getVariable("Patient_ID"));
        d.removeColumn(d.getVariable("Screen IgG"));
        d.removeColumn(d.getVariable("Week 12 IgG"));
        DataSet d3 = MixedUtils.completeCases(d);
        System.out.println(d3);
        int ns = 4;
        double g = 0.01;
        int numLambdas = 40;
        double [] lambda = new double[numLambdas];
        for(int i = 0; i < numLambdas;i++)
        {
            lambda[i] = .1+(.7/numLambdas)*i;
        }
        STEPS s = new STEPS(d3,lambda,g,ns,true);
        Graph g2 = s.runSteps();
        double [][] stab = s.stabilities;
        PrintStream out = new PrintStream("STEPS_graph_TNFMUC.txt");
        out.println(g2);
        out.flush();
        out.close();
        out = new PrintStream("STEPS_Stabilities_TNFMUC.txt");
        for(int i = 0; i < d.getNumColumns();i++)
        {
            out.print(d.getVariable(i).getName());
            if(i < d.getNumColumns()-1)
                out.print("\t");
            else
                out.println();
        }
        for(int i = 0; i < d.getNumColumns();i++)
        {
            out.print(d.getVariable(i).getName()+"\t");
            for(int j = 0; j < d.getNumColumns();j++)
            {
                out.print(stab[i][j]);
                if(j < d.getNumColumns()-1)
                    out.print("\t");
                else
                    out.println();
            }
        }
        out.flush();
        out.close();
        IndependenceTest i = new IndTestMultinomialAJ(d3,.2);
        FciMaxP f = new FciMaxP(i);
        f.setDepth(1);
        f.setInitialGraph(g2);
        out = new PrintStream("max_graph_TNFMUC.txt");
        out.println(f.search());
        out.flush();
        out.close();

    }
}
