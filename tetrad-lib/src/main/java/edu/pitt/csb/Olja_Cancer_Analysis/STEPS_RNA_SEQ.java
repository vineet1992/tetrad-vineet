package edu.pitt.csb.Olja_Cancer_Analysis;

import cern.colt.matrix.DoubleMatrix2D;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DelimiterType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.SearchWrappers;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * Created by vinee_000 on 10/13/2017.
 */
public class STEPS_RNA_SEQ {
    public static void main(String [] args)throws Exception
    {
        DelimiterType d2 = DelimiterType.TAB;
        DataSet d = MixedUtils.loadDataSet2(args[0],d2);
        String graphOut = args[1];
        String stabOut = args[2];
        String subDir = args[3];
        String target = args[4];
        for(int i = 5; i < args.length;i++)
        {
            d.removeColumn(d.getVariable(args[i]));
        }

        File f = new File(subDir);
        if(!f.exists())
            f.mkdir();

       /* d.removeColumn(d.getVariable("Response"));
        d.removeColumn(d.getVariable("IgG_Week0"));
        d.removeColumn(d.getVariable("IgG_Week12"));
        d.removeColumn(d.getVariable("Response_DP"));*/
       // d.removeColumn(d.getVariable("IgG_Ratio"));
       // d.removeColumn(d.getVariable("Gender"));
       // d.removeColumn(d.getVariable("Age"));
        //  DataSet d3 = MixedUtils.completeCases(d);


        //ns set to number of samples since we're doing leave-one-out cross validation
        int ns = d.getNumRows();
        double g = 0.01;
        int numLambdas = 40;
        double [] lambda = new double[numLambdas];
        for(int i = 0; i < numLambdas;i++)
        {
            lambda[i] = .1+(.7/numLambdas)*i;
        }
        STEPS s = new STEPS(d,lambda,g,d.getNumRows(),true);
        System.out.print("Running StEPS...");
        Graph g2 = s.runStepsPar();
        System.out.println("Done");
        double [][] stab = s.stabilities;
        double [] l = s.lastLambda;

        System.out.println("Lambdas Chosen: " + Arrays.toString(s.lastLambda));


        //This first segment tells us what the model will look like using the full dataset
        double [] params = {s.lastLambda[0],s.lastLambda[1],s.lastLambda[2]};
        DataGraphSearch gs = new SearchWrappers.MGMWrapper(params);
        int [][] samps = StabilityUtils.generateSubsamples(d.getNumRows());
        DataSet [] subs = new DataSet[ns];
        for(int i = 0; i < ns;i++) {
            subs[i] = d.copy();
            Arrays.sort(samps[i]);
            subs[i] = subs[i].subsetRows(samps[i]);
            PrintStream temp = new PrintStream(subDir + "/Test_Set_" + i + ".txt");
            for(int j = 0; j < d.getNumRows();j++)
            {
                if(Arrays.binarySearch(samps[i],j)<0)
                    temp.println(j+1);
            }
            temp.flush();
            temp.close();
            temp = new PrintStream(subDir + "/Training_Set_" + i + ".txt");
            for(int j = 0; j < samps[i].length;j++)
                temp.println(samps[i][j] + 1);
            temp.flush();
            temp.close();
        }
        System.out.print("Finding Full Data Stabilities...");
        DoubleMatrix2D edges = StabilityUtils.StabilitySearchPar(d,gs,subs);
        System.out.println("Done");
        PrintStream out = new PrintStream(args[1]);
        out.println(g2);
        out.flush();
        out.close();
        out = new PrintStream(args[2]);
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

        //This portion tells us how good the modeling procedure is on this dataset (cross-validation of the full procedure)

        for(int i = 0; i < subs.length;i++)
        {

            System.out.println("Cross Validation for Subsample " + i + " out of " + subs.length + "...");
            //TODO fix this up and possibly repeat STEPS procedure on each subsample? Not sure what to do here exactly
            //Only select the features that are stable across subsamples, and that are connected when running MGM on full train
            DataSet train = subs[i];
            s = new STEPS(train,lambda,g,ns,true);
            Graph gOut = s.runStepsPar();
            double [][] stabs = s.stabilities;
            double [] lam = {s.lastLambda[0],s.lastLambda[1],s.lastLambda[2]};
            PrintStream p = new PrintStream(subDir + "/Neighbors_" + i + ".txt");
            PrintStream p2 = new PrintStream(subDir + "/Neighbors_PCS_" + i + ".txt");
            PrintStream p3 = new PrintStream(subDir + "/Neighbors_PCS_Causes_" + i + ".txt");


            System.out.print("Running PCS...");
            IndependenceTest indy = new IndTestMultinomialAJ(train,0.05);
            PcStable pc = new PcStable(indy);
            pc.setInitialGraph(gOut);
            Graph gOut2 = pc.search();
            System.out.println("Done");

            //MGM Adjacencies
            for(Node n: gOut.getAdjacentNodes(gOut.getNode(target)))
            {
                int x = d.getColumn(d.getVariable(n.getName()));
                int y = d.getColumn(d.getVariable(target));

                if(stabs[x][y] > 0.5 || stabs[y][x] > 0.5)
                    p.println(n.getName() + "\t" + stabs[x][y]);
            }

//PCS Adjacencies
            for(Node n: gOut2.getAdjacentNodes(gOut2.getNode(target)))
            {
                int x = d.getColumn(d.getVariable(n.getName()));
                int y = d.getColumn(d.getVariable(target));

                if(stabs[x][y] > 0.5 || stabs[y][x] > 0.5)
                    p2.println(n.getName() + "\t" + stabs[x][y]);
            }


            //PCS Direct Causes Only
            for(Node n: gOut2.getParents(gOut2.getNode(target)))
            {
                int x = d.getColumn(d.getVariable(n.getName()));
                int y = d.getColumn(d.getVariable(target));

                if(stabs[x][y] > 0.5 || stabs[y][x] > 0.5)
                    p3.println(n.getName() + "\t" + stabs[x][y]);
            }


            p.flush();
            p.close();

            System.out.println("Done");
        }
    }
}
