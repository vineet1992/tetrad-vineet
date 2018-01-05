package edu.pitt.csb.Olja_Cancer_Analysis;

import cern.colt.matrix.DoubleMatrix2D;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DelimiterType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.SearchWrappers;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 9/1/2017.
 */
public class STEPS_Immuno_Seq {
    public static void main(String [] args)throws Exception
    {
        DelimiterType d2 = DelimiterType.TAB;
        DataSet d = MixedUtils.loadDataSet2("data_cleaned.txt",d2);
        System.out.println(d);

      //  d.removeColumn(d.getVariable("Sample"));
      //  d.removeColumn(d.getVariable("Response"));
        d.removeColumn(d.getVariable("Screen_IGG"));
        d.removeColumn(d.getVariable("Week12_IGG"));
       d.removeColumn(d.getVariable("IgG_Ratio"));
        //  DataSet d3 = MixedUtils.completeCases(d);
        int ns = 10;
        double g = 0.01;
        int numLambdas = 20;
        boolean loo = true;
        double low = 0.15;
        double high = 0.7;
        double [] lambda = new double[numLambdas];
        for(int i = 0; i < numLambdas;i++)
        {
            lambda[i] = low+(high/numLambdas)*i;
        }
        STEPS s = new STEPS(d,lambda,g,ns,loo);

        Graph g2 = s.runStepsPar();
       // double [] la = {0.28,0.25,0.2};
       //DoubleMatrix2D mat = StabilityUtils.StabilitySearchPar(d,new SearchWrappers.MGMWrapper(la));
      // double [][] stab = mat.toArray();
      // MGM m = new MGM(d,s.lastLambda);
      // m.learnEdges(1000);
      // Graph g2 = m.graphFromMGM();
        double [][] stab = s.stabilities;
        PrintStream out = new PrintStream("STEPS_IS.txt");
        out.println(g2);
        out.flush();
        out.close();
        out = new PrintStream("STEPS_Stabilities_IS.txt");
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
        IndependenceTest i = new IndTestMultinomialAJ(d,.2);
        FciMaxP f = new FciMaxP(i);
        f.setInitialGraph(g2);
        out = new PrintStream("Max_IS.txt");
        out.println(f.search());
        out.flush();
        out.close();


    }
}
