package edu.pitt.csb.mgm;

import cern.colt.matrix.DoubleMatrix2D;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.SearchWrappers;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * Created by vinee_000 on 4/6/2017.
 */
public class lung_analysis {
    public static void main(String [] args)throws Exception
    {
        int N = 10;
        double g = 0.05;
        DataSet d = MixedUtils.loadDataSet2("MGM_Full_Data.txt");
        System.out.println(d);
       // d.removeColumn(d.getVariable("Score"));
        d.removeColumn(d.getVariable("GOLD4"));
        d.removeColumn(d.getVariable("Smk_Status"));
        d.removeColumn(d.getVariable("Pack_Years"));
      //  d.removeColumn(d.getVariable("Smoke_Duration"));
      //  d.removeColumn(d.getVariable("Years_quit"));
        d.removeColumn(d.getVariable("Nodule.size.group"));

        File f3 = new File("Subsamples_Results");
        if(!f3.exists())
            f3.mkdir();

       // d.removeColumn(d.getVariable("Years_quit"));
        int [][] subs = StabilityUtils.generateSubsamples(N,d.getNumRows());

        DataSet [] subsamples = new DataSet[N];
        for(int i = 0; i < N;i++) {
            subsamples[i] = d.copy();
            Arrays.sort(subs[i]);
            subsamples[i].removeRows(subs[i]);
            PrintStream temp = new PrintStream("Subsamples_Results/Training_Set_" + i + ".txt");
            for(int j = 0; j < d.getNumRows();j++)
            {
                if(Arrays.binarySearch(subs[i],j)<0)
                    temp.println(j+1);
            }
            temp.flush();
            temp.close();
            temp = new PrintStream("Subsamples_Results/Test_Set_" + i + ".txt");
            for(int j = 0; j < subs[i].length;j++)
                temp.println(subs[i][j] + 1);
            temp.flush();
            temp.close();
        }
        int numLambdas = 40;
        double [] lambda = new double[numLambdas];
        for(int i = 0; i < numLambdas;i++)
        {
            lambda[i] = .1+(.7/numLambdas)*i;
        }
        IndependenceTest i  = new IndTestMultinomialAJ(d,.2);
        STEPS s = new STEPS(d,lambda,g,subsamples);
        Graph g2 = s.runStepsPar();
        double [][] stab = s.stabilities;

        double []params = {s.lastLambda[0],s.lastLambda[1],s.lastLambda[2],0.2};
        DataGraphSearch gs = new SearchWrappers.MFMWrapper(params);
        DoubleMatrix2D db = StabilityUtils.StabilitySearchPar(d,gs,subsamples);
        FciMaxP f = new FciMaxP(i);
        f.setInitialGraph(g2);
        PrintStream out = new PrintStream("FCI_MAX_Full_Data_Trimmed_Variables.txt");
        out.println(f.search());
        out = new PrintStream("MGM_Stabilities_Full_Data_Trimmed_Variables.txt");
        for(int j = 0; j < d.getNumColumns();j++)
        {
            out.print(d.getVariable(j).getName());
            if(j < d.getNumColumns()-1)
                out.print("\t");
            else
                out.println();
        }
        for(int k = 0; k < d.getNumColumns();k++)
        {
            out.print(d.getVariable(k).getName()+"\t");
            for(int j = 0; j < d.getNumColumns();j++)
            {
                out.print(stab[k][j]);
                if(j < d.getNumColumns()-1)
                    out.print("\t");
                else
                    out.println();
            }
        }
        out.flush();
        out.close();
        out = new PrintStream("MAX_Stabilities_Full_Data_Trimmed_Variables.txt");
        for(int j = 0; j < d.getNumColumns();j++)
        {
            out.print(d.getVariable(j).getName());
            if(j < d.getNumColumns()-1)
                out.print("\t");
            else
                out.println();
        }
        for(int k = 0; k < d.getNumColumns();k++)
        {
            out.print(d.getVariable(k).getName()+"\t");
            for(int j = 0; j < d.getNumColumns();j++)
            {
                out.print(db.get(k,j));
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
