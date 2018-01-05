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
 * Created by vinee_000 on 10/13/2017.
 */
public class STEPS_RNA_SEQ {
    public static void main(String [] args)throws Exception
    {
        DelimiterType d2 = DelimiterType.TAB;
        DataSet d = MixedUtils.loadDataSet2(args[0],d2);
        d.removeColumn(d.getVariable("Response"));
        //d.removeColumn(d.getVariable("IgG_Week0"));
       // d.removeColumn(d.getVariable("IgG_Week12"));
        d.removeColumn(d.getVariable("Response_DP"));
       // d.removeColumn(d.getVariable("IgG_Ratio"));
        System.out.println(d);

       // d.removeColumn(d.getVariable("Gender"));
       // d.removeColumn(d.getVariable("Age"));
        //  DataSet d3 = MixedUtils.completeCases(d);
        int ns = 20;
        double g = 0.01;
        int numLambdas = 40;
        double [] lambda = new double[numLambdas];
        for(int i = 0; i < numLambdas;i++)
        {
            lambda[i] = .1+(.7/numLambdas)*i;
        }
        STEPS s = new STEPS(d,lambda,g,d.getNumRows(),true);
        Graph g2 = s.runStepsPar();
        double [][] stab = s.stabilities;
        PrintStream out = new PrintStream("Dp_Full_RPKM_IgG.txt");
        out.println(g2);
        out.flush();
        out.close();
        out = new PrintStream("Stabilities_Dp_Full_RPKM_IgG.txt");
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
    }
}
