package edu.pitt.csb.Priors;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.MGM_Priors;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * Created by vinee_000 on 12/6/2017.
 */
public class realDataPriorTest {

    public static void main(String [] args) throws Exception
    {
        boolean runNoPrior = false;
        boolean runIrrelevantPrior = false;
        boolean runRelevantPrior = false;
        boolean runBoth = true;
        boolean computeMetrics = false;
        int numLambda = 40;
        int numPriors = 9;
        double g = 0.01;
        double [] lambda = new double[numLambda];
        double low = .05;
        double high = .9;
        double inc = (high-low)/(numLambda-1);
        for(int i = 0; i < numLambda;i++)
        {
            lambda[i] = i*inc + low;
        }

        DataSet data = MixedUtils.loadDataSet2("genes_with_clinical.txt");
        //System.out.println(data);
        data = MixedUtils.completeCases(data);
        //System.out.println(data);
        File f = new File("Subsamples");
        if(!f.exists())
            f.mkdir();
        int numSub = 20;
        int b = (int) Math.floor(10 * Math.sqrt(data.getNumRows()));
        if (b > data.getNumRows())
            b = data.getNumRows() / 2;
        int[][] samps = StabilityUtils.subSampleNoReplacement(data.getNumRows(), b, numSub);
        DataSet [] subsamples = new DataSet[numSub];
        for(int i = 0; i< numSub;i++)
        {
            //File temp = new File("Subsamples/Subsample_" + i + ".txt");
            File temp = new File("Subsample_" + i + ".txt");
            if(!temp.exists())
            {
                DataSet tData = data.subsetRows(samps[i]);
                PrintStream out2 = new PrintStream(temp.getAbsolutePath());
                out2.println(tData);
                subsamples[i] = tData;
            }
            else
            {
                subsamples[i] = MixedUtils.loadDataSet2(temp.getAbsolutePath());
            }
        }

        //Test no Prior situation
        if(runNoPrior) {
            STEPS s = new STEPS(data, lambda, g, subsamples);
            Graph graph = s.runStepsPar();
            double[][] stab = s.stabilities;
            PrintStream out = new PrintStream("Stabilities_No_Prior.txt");
            for (int i = 0; i < data.getNumColumns(); i++) {
                if (i == data.getNumColumns() - 1)
                    out.println(data.getVariable(i).getName());
                else
                    out.print(data.getVariable(i).getName() + "\t");
            }
            for (int i = 0; i < stab.length; i++) {
                out.print(data.getVariable(i).getName() + "\t");
                for (int j = 0; j < stab[i].length; j++) {
                    if (j == stab[i].length - 1)
                        out.println(stab[i][j]);
                    else
                        out.print(stab[i][j] + "\t");
                }
            }
            out.flush();
            out.close();
            out = new PrintStream("Graph_No_Prior.txt");
            out.println(graph);
            out.flush();
            out.close();
        }

        //Test Irrelevant Prior Situation

        if(runIrrelevantPrior)
        {

            TetradMatrix [] priors = new TetradMatrix[numPriors];
            for(int i = 0; i < numPriors;i++) {
                //priors[i] = new TetradMatrix(loadPrior(new File("prior_sources/Irr_Prior_" + i + ".txt"),data.getNumColumns()));
                priors[i] = new TetradMatrix(loadPrior(new File("Irr_Prior_" + i + ".txt"),data.getNumColumns()));
                System.out.println(priors[i].rows() + "," + priors[i].columns());
            }
            System.out.println("Constructing lambdas...");
            mgmPriors m = new mgmPriors(subsamples.length,lambda,data,priors,subsamples);
            PrintStream lb = new PrintStream("Irrelevant_Lambdas.txt");
            for(int i = 0; i < m.getLambdas().length;i++)
            {
                lb.println(Arrays.toString(m.getLambdas()[i]));
            }
            lb.flush();
            lb.close();
            System.out.println("Running Priors...");
            Graph g2 = m.runPriors();
            double [][] stab = m.edgeScores;
            double [] weights = m.expertWeights;
            PrintStream w = new PrintStream("Weights_Irrelevant_Prior.txt");
            for(int i = 0; i < weights.length;i++)
            {
                w.println(i + "\t" + weights[i]);
            }
            w.flush();
            w.close();
            PrintStream out = new PrintStream("Stabilities_Irrelevant_Prior.txt");
            for (int i = 0; i < data.getNumColumns(); i++) {
                if (i == data.getNumColumns() - 1)
                    out.println(data.getVariable(i).getName());
                else
                    out.print(data.getVariable(i).getName() + "\t");
            }
            for (int i = 0; i < stab.length; i++) {
                out.print(data.getVariable(i).getName() + "\t");
                for (int j = 0; j < stab[i].length; j++) {
                    if (j == stab[i].length - 1)
                        out.println(stab[i][j]);
                    else
                        out.print(stab[i][j] + "\t");
                }
            }
            out.flush();
            out.close();
            out = new PrintStream("Graph_Irrelevant_Prior.txt");
            out.println(g2);
            out.flush();
            out.close();


        }




        //Test Relevant Pathway Situation
        if(runRelevantPrior)
        {

            TetradMatrix [] priors = new TetradMatrix[numPriors+1];
            for(int i = 0; i < numPriors;i++) {
                priors[i] = new TetradMatrix(loadPrior(new File("Prior_" + i + ".txt"),data.getNumColumns()));

            }
            priors[priors.length-1] = new TetradMatrix(loadPrior(new File("Prior_PAM50.txt"),data.getNumColumns()));

            mgmPriors m = new mgmPriors(subsamples.length,lambda,data,priors,subsamples);
            PrintStream lb = new PrintStream("Relevant_Lambdas.txt");
            for(int i = 0; i < m.getLambdas().length;i++)
            {
                lb.println(Arrays.toString(m.getLambdas()[i]));
            }
            lb.flush();
            lb.close();
            System.out.println("Running Priors...");
            Graph g2 = m.runPriors();
            double [][] stab = m.edgeScores;
            double [] weights = m.expertWeights;
            PrintStream w = new PrintStream("Weights_Relevant_Prior.txt");
            for(int i = 0; i < weights.length;i++)
            {
                w.println(i + "\t" + weights[i]);
            }
            w.flush();
            w.close();
            PrintStream out = new PrintStream("Stabilities_Relevant_Prior.txt");
            for (int i = 0; i < data.getNumColumns(); i++) {
                if (i == data.getNumColumns() - 1)
                    out.println(data.getVariable(i).getName());
                else
                    out.print(data.getVariable(i).getName() + "\t");
            }
            for (int i = 0; i < stab.length; i++) {
                out.print(data.getVariable(i).getName() + "\t");
                for (int j = 0; j < stab[i].length; j++) {
                    if (j == stab[i].length - 1)
                        out.println(stab[i][j]);
                    else
                        out.print(stab[i][j] + "\t");
                }
            }
            out.flush();
            out.close();
            out = new PrintStream("Graph_Relevant_Prior.txt");
            out.println(g2);
            out.flush();
            out.close();


        }

        if(runBoth)
        {

            TetradMatrix [] priors = new TetradMatrix[numPriors*2+1];
            for(int i = 0; i < numPriors;i++) {
                priors[i] = new TetradMatrix(loadPrior(new File("Prior_" + i + ".txt"),data.getNumColumns()));
                priors[numPriors+i] = new TetradMatrix(loadPrior(new File("Irr_Prior_" + i + ".txt"),data.getNumColumns()));
            }
            priors[priors.length-1] = new TetradMatrix(loadPrior(new File("Prior_PAM50.txt"),data.getNumColumns()));
            mgmPriors m = new mgmPriors(subsamples.length,lambda,data,priors,subsamples);
            PrintStream lb = new PrintStream("Both_Lambdas.txt");
            for(int i = 0; i < m.getLambdas().length;i++)
            {
                lb.println(Arrays.toString(m.getLambdas()[i]));
            }
            lb.flush();
            lb.close();
            System.out.println("Running Priors...");
            Graph g2 = m.runPriors();
            double [][] stab = m.edgeScores;
            double [] weights = m.expertWeights;
            PrintStream w = new PrintStream("Weights_Both_Prior.txt");
            for(int i = 0; i < weights.length;i++)
            {
                w.println(i + "\t" + weights[i]);
            }
            w.flush();
            w.close();
            PrintStream out = new PrintStream("Stabilities_Both_Prior.txt");
            for (int i = 0; i < data.getNumColumns(); i++) {
                if (i == data.getNumColumns() - 1)
                    out.println(data.getVariable(i).getName());
                else
                    out.print(data.getVariable(i).getName() + "\t");
            }
            for (int i = 0; i < stab.length; i++) {
                out.print(data.getVariable(i).getName() + "\t");
                for (int j = 0; j < stab[i].length; j++) {
                    if (j == stab[i].length - 1)
                        out.println(stab[i][j]);
                    else
                        out.print(stab[i][j] + "\t");
                }
            }
            out.flush();
            out.close();
            out = new PrintStream("Graph_Both_Prior.txt");
            out.println(g2);
            out.flush();
            out.close();


        }


        //Output metrics, precision and recall of PAM50 Genes (Check the neighborhood of the Subtype variable)
        //, Recovery of relevant pathway connections
    }
    public static double [][] loadPrior(File data, int numVariables) throws Exception
    {
        double [][] temp = new double[numVariables][numVariables];
        BufferedReader b = new BufferedReader(new FileReader(data));
        b.readLine();
        for(int i = 0; i < numVariables;i++)
        {
            String [] line = b.readLine().split("\t");
            for(int j = 0; j < numVariables;j++)
            {
                temp[i][j] = Double.parseDouble(line[j]);
            }
        }
        return temp;
    }
}
