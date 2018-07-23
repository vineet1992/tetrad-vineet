package edu.pitt.csb.Lung_Nodule_Analysis;

import cern.colt.matrix.DoubleMatrix2D;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.SearchWrappers;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by vinee_000 on 4/6/2017.
 */
public class lung_analysis {
    public static void main(String [] args)throws Exception
    {
        boolean fiveByTwo = true;
        int N = 10;
        double alpha = 0.05;
        double g = 0.05;
        String subDir = "Subsamples_Five_By_Two";

        DataSet d = MixedUtils.loadDataSet2("MGM_Full_Data.txt", DelimiterType.TAB);

       // d.removeColumn(d.getVariable("Score"));

      //  d.removeColumn(d.getVariable("Pack_Years"));

      //  d.removeColumn(d.getVariable("Years_quit"));

        /*Removed in the current analysis
        d.removeColumn(d.getVariable("Nodule.size.group"));
        d.removeColumn(d.getVariable("Cigs_Per_Day"));
        d.removeColumn(d.getVariable("Smoke_Duration"));
        d.removeColumn(d.getVariable("GOLD4"));
        d.removeColumn(d.getVariable("Smk_Status"));
        */
       // System.out.println(d);
        System.out.println(MixedUtils.getContinousData(d));
        System.out.println(MixedUtils.getDiscreteData(d));

        File f3 = new File(subDir);
        if(!f3.exists())
            f3.mkdir();

       // d.removeColumn(d.getVariable("Years_quit"));

            int[][] subs = StabilityUtils.generateSubsamples(N, d.getNumRows());
            if(fiveByTwo)
            {
                subs = StabilityUtils.subSampleNoReplacement(d.getNumRows(),d.getNumRows()/2,5);
            }

            DataSet[] subsamples = new DataSet[N];
            for (int i = 0; i < N; i++) {
                subsamples[i] = d.copy();
                Arrays.sort(subs[i]);
                subsamples[i].removeRows(subs[i]);
                PrintStream temp = new PrintStream(subDir + "/Training_Set_" + i + ".txt");
                for (int j = 0; j < d.getNumRows(); j++) {
                    if (Arrays.binarySearch(subs[i], j) < 0)
                        temp.println(j + 1);
                }
                temp.flush();
                temp.close();
                temp = new PrintStream(subDir + "/Test_Set_" + i + ".txt");
                for (int j = 0; j < subs[i].length; j++)
                    temp.println(subs[i][j] + 1);
                temp.flush();
                temp.close();
            }
            System.out.println(Arrays.toString(subsamples));
            int numLambdas = 40;
            double[] lambda = new double[numLambdas];
            for (int i = 0; i < numLambdas; i++) {
                lambda[i] = .1 + (.7 / numLambdas) * i;
            }
            IndependenceTest ind = new IndTestMultinomialAJ(d, alpha);
            STEPS s = new STEPS(d, lambda, g, subs);
            Graph g2 = s.runStepsPar();
            double[][] stab = s.stabilities;

            double[] params = {s.lastLambda[0], s.lastLambda[1], s.lastLambda[2], 0.2};
            DataGraphSearch gs = new SearchWrappers.MFMWrapper(params);
            DoubleMatrix2D db = StabilityUtils.StabilitySearchPar(d, gs, subs);
            FciMaxP f = new FciMaxP(ind);
            f.setInitialGraph(g2);
            IKnowledge k2 = new Knowledge2();//Nothing will cause sex, education, age, and Lung Cancer won't cause anything else
            for (Node n : d.getVariables()) {
                k2.addVariable(n.getName());
                if (n.getName().equals("Lung.cancer"))
                    k2.addToTier(2, n.getName());
                else if (n.getName().equals("Sex") || n.getName().equals("Education") || n.getName().equals("Age"))
                    k2.addToTier(0, n.getName());
                else
                    k2.addToTier(1, n.getName());

            }
            f.setKnowledge(k2);
            PrintStream out = new PrintStream("FCI_MAX_Graph.txt");
            out.println(f.search());
            out = new PrintStream("MGM_Stabilities.txt");
            for (int j = 0; j < d.getNumColumns(); j++) {
                out.print(d.getVariable(j).getName());
                if (j < d.getNumColumns() - 1)
                    out.print("\t");
                else
                    out.println();
            }
            for (int k = 0; k < d.getNumColumns(); k++) {
                out.print(d.getVariable(k).getName() + "\t");
                for (int j = 0; j < d.getNumColumns(); j++) {
                    out.print(stab[k][j]);
                    if (j < d.getNumColumns() - 1)
                        out.print("\t");
                    else
                        out.println();
                }
            }
            out.flush();
            out.close();
            out = new PrintStream("MAX_Stabilities.txt");
            for (int j = 0; j < d.getNumColumns(); j++) {
                out.print(d.getVariable(j).getName());
                if (j < d.getNumColumns() - 1)
                    out.print("\t");
                else
                    out.println();
            }
            for (int k = 0; k < d.getNumColumns(); k++) {
                out.print(d.getVariable(k).getName() + "\t");
                for (int j = 0; j < d.getNumColumns(); j++) {
                    out.print(db.get(k, j));
                    if (j < d.getNumColumns() - 1)
                        out.print("\t");
                    else
                        out.println();
                }
            }
            out.flush();
            out.close();

            for (int i = 0; i < N; i++) {

                DataSet train = subsamples[i];
                double[] lam = {s.lastLambda[0], s.lastLambda[1], s.lastLambda[2]};
                MGM m = new MGM(train, lam);
                m.learnEdges(1000);
                Graph g3 = m.graphFromMGM();

                IndependenceTest ii = new IndTestMultinomialAJ(train, alpha);
                FciMaxP f2 = new FciMaxP(ii);
                f2.setInitialGraph(g3);
                Graph fin = f2.search();

                PrintStream p = new PrintStream(subDir + "/Neighbors_" + i + ".txt");

                for (Node n : fin.getAdjacentNodes(fin.getNode("Lung.cancer")))
                    p.println(n.getName());

                p.flush();
                p.close();

            }
    }

}
