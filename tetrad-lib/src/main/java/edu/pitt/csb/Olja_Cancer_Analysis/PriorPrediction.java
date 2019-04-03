package edu.pitt.csb.Olja_Cancer_Analysis;

import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.LogisticRegression;
import edu.cmu.tetrad.regression.Regression;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.search.CpcStable;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.pitt.csb.Priors.PriorUtils;
import edu.pitt.csb.Priors.mgmPriors;
import edu.pitt.csb.Priors.realDataPriorTest;
import edu.pitt.csb.Priors.runPriors;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by vinee_000 on 7/26/2018.
 */
public class PriorPrediction {
    public static void main(String [] args)throws Exception
    {


        String target = "IgG_Ratio";
        String priorDir = "Priors_All";
        String dataFile = "Highest_Variance_Prediction_Log.txt";
        int numLambdas = 10;
        boolean loocv = true;
        int ns = 5;
        double lamLow = 0.1;
        double lamHigh = 0.9;
        int start = 0;
        int end = 5;
        boolean steps = false;
        double g = 0.05;

        int index = 0;
        while(index < args.length) {
            if (args[index].equals("-ns")) {
                ns = Integer.parseInt(args[index + 1]);
                index += 2;
            } else if (args[index].equals("-nl")) {
                numLambdas = Integer.parseInt(args[index + 1]);
                index += 2;
            } else if (args[index].equals("-llow")) {
                lamLow = Double.parseDouble(args[index + 1]);
                index += 2;
            } else if (args[index].equals("-lhigh")) {
                lamHigh = Double.parseDouble(args[index + 1]);
                index += 2;
            } else if (args[index].equals("-target")) {
                target = args[index + 1];
                index += 2;
            } else if (args[index].equals("-priorDir")) {
                priorDir = args[index+1];
                index+=2;
            } else if (args[index].equals("-dataFile")) {
                dataFile = args[index + 1];
                index += 2;
            }
            else if (args[index].equals("-start")) {
                start = Integer.parseInt(args[index + 1]);
                index += 2;
            }
            else if (args[index].equals("-end")) {
                end = Integer.parseInt(args[index + 1]);
                index += 2;
            }
            else if(args[index].equals("-steps"))
            {
                steps = true;
                index++;
            }
            else if(args[index].equals("-g"))
            {
                g = Double.parseDouble(args[index+1]);
                index+=2;
            }
            else
            {
                System.err.println("Unrecognized Command Line Option: " + args[index]);
                System.exit(-1);
            }
        }


        double [] initLambdas = new double[numLambdas];
        for(int i = 0; i < initLambdas.length;i++)
        {
            initLambdas[i] = lamLow + (lamHigh)*i/(double)numLambdas;
        }

        String runName = target + "_" + priorDir;
        DataSet data = MixedUtils.loadDataSet2(dataFile);

        if(loocv)
            ns = data.getNumRows();



        File f = new File("Predictions_" + runName);
        if(!f.exists())
            f.mkdir();
        File pd = new File(priorDir);

        PrintStream pi1 = new PrintStream(f.getName() + "/Predictions.txt");
        PrintStream pi2 = new PrintStream(f.getName() + "/Features.txt");
        pi1.println("Run\tPrediction_piMGM\tPrediction_CPC_MB\tPrediction_CPC\tActual");
        pi2.println("Run\tFeatures_piMGM\tFeatures_CPC_MB\tFeatures_CPC");
        if(end==-1)
            end = data.getNumRows();

        SparseDoubleMatrix2D[] priors = new SparseDoubleMatrix2D[1];
        String[] priorNames = new String[1];

        if(!steps) {
            if (!pd.isDirectory()) {
                System.out.println("Prior Directory doesn't exist...exiting");
                System.exit(-1);
            }
             priors = new SparseDoubleMatrix2D[pd.listFiles().length];
            priorNames = new String[pd.listFiles().length];

            int ct = 0;
            for (File ff : pd.listFiles()) {
                priorNames[ct] = ff.getName();
                BufferedReader temp = new BufferedReader(new FileReader(ff));

                //Assume this prior knowledge file is in sif format
                if (temp.readLine().split("\t").length < data.getNumColumns()) {
                    temp.close();
                    priors[ct] = new SparseDoubleMatrix2D(PriorUtils.loadPriorSif(ff, data));
                }
                //Otherwise assume the full matrix is given
                else {
                    temp.close();
                    priors[ct] = new SparseDoubleMatrix2D(PriorUtils.loadPrior(ff, data.getNumColumns()));
                }
                ct++;
            }

        }
        //Two Output Files for both piMGM alone and piMGM with causal discovery algorithm downstream
        //First file is list of features selected in each run
        //Second file is True values and prediction made for each run

        for(int i = start; i < end;i++)
        {
            System.out.println("Running CV set " + i + " out of " + data.getNumRows());
            int [] rows = new int[data.getNumRows()-1];
            int count = 0;
            for(int j = 0; j < data.getNumRows();j++)
            {
                if(i!=j)
                {
                    rows[count] = j;
                    count++;
                }
            }
            DataSet temp = data.subsetRows(rows);
            int [][] samps = runPriors.genSubs(temp,ns,loocv);
            ns = samps.length;
            Graph out;
            mgmPriors p;

            if(steps)
            {
                System.out.print("Running STEPS...");
                STEPS s = new STEPS(temp,initLambdas,g,samps);
                s.runStepsArrayPar();
                out = s.lastGraph;
                System.out.println("Done");
            }
            else
            {
                p = new mgmPriors(ns,initLambdas,temp,priors,samps,true);
                System.out.println("Running piMGM");
                out = p.runPriors();
                System.out.println("Priors: " + Arrays.toString(priorNames));
                System.out.println("Weights: " + Arrays.toString(p.normalizedExpertWeights));
                System.out.println("P-Values; " + Arrays.toString(p.uncorrectedPValues));
            }
            IndependenceTest ind = new IndTestMultinomialAJ(temp,0.05);
            System.out.println("Running CPC Stable");
            CpcStable cpc = new CpcStable(ind);
            cpc.setInitialGraph(out);
            Graph cpcOut = cpc.search();

            List<Node> neighbors = out.getAdjacentNodes(out.getNode(target));
            System.out.println("Adjacent To Target: " + neighbors);
            double pred = getRegressionResult(neighbors,temp,data,i,target);
            double real = -1;
            if(data.getVariable(target)instanceof DiscreteVariable)
            {
                real = data.getInt(i,data.getColumn(data.getVariable(target)));
            }else
            {
                real = data.getDouble(i,data.getColumn(data.getVariable(target)));
            }
            System.out.println("piMGM Prediction: " + pred);
            pi1.print(i + "\t" + pred +"\t");
            pi2.print(i + "\t" + neighbors + "\t");

            neighbors = markovBlanket(cpcOut,target);
            if(neighbors.size()==0) {
                neighbors = out.getAdjacentNodes(out.getNode(target));
            }
            pred = getRegressionResult(neighbors,temp,data,i,target);
            System.out.println("CPC-MB Neighbors: " + neighbors);
            System.out.println("CPC-MB Prediction: " + pred);
            pi1.print(pred + "\t");
            pi2.print(neighbors + "\t");
            System.out.println("Actual value: " + real);

            neighbors = cpcOut.getAdjacentNodes(cpcOut.getNode(target));
            if(neighbors.size()==0)
                neighbors = out.getAdjacentNodes(out.getNode(target));
            pred = getRegressionResult(neighbors,temp,data,i,target);
            System.out.println("CPC Neighbors: " + neighbors);
            System.out.println("CPC Prediction: " + pred);
            pi1.println(pred + "\t" + real);
            pi2.println(neighbors);

            pi1.flush();
            pi2.flush();



        }
        pi1.close();
        pi2.close();




    }

    public static double getRegressionResult(List<Node> neighbors,DataSet small, DataSet full, int row, String target)
    {
        //DO NOT ALLOW DISCRETE VARIABLES (MAYO-FINN)
        List<Node> temp = new ArrayList<Node>();
        for (Node n : neighbors) {
            if(!(small.getVariable(n.getName())instanceof DiscreteVariable))
            temp.add(small.getVariable(n.getName()));
        }


        if(small.getVariable(target)instanceof DiscreteVariable)
        {
            LogisticRegression lr  = new LogisticRegression(small);
            LogisticRegression.Result r = lr.regress((DiscreteVariable)small.getVariable(target),temp);
            double sum = 0;
            for (int x = 0; x < temp.size(); x++) {
                sum+= r.getCoefs()[x] * full.getDouble(row, full.getColumn(temp.get(x)));
            }
            sum = sum + r.getIntercept();
            return 1/(1 + Math.exp(-1*sum));

        }
        else {


            RegressionDataset rd = new RegressionDataset(small);
            RegressionResult res = rd.regress(small.getVariable(target), temp);
            double[] testVec = new double[temp.size()];
            for (int x = 0; x < temp.size(); x++) {
                testVec[x] = full.getDouble(row, full.getColumn(temp.get(x)));
            }
            return res.getPredictedValue(testVec);
        }
    }
    public static List<Node> markovBlanket(Graph out, String target)
    {
        List<Node> neighbors = out.getAdjacentNodes(out.getNode(target));
        for(Node n: out.getChildren(out.getNode(target)))
        {
            List<Node> parents = out.getParents(n);
            for(Node p: parents)
            {
                if(!neighbors.contains(p))
                    neighbors.add(p);
            }
        }
        return neighbors;
    }
}
