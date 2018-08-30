package edu.pitt.csb.Priors;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MGM_Priors;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;
import edu.pitt.csb.stability.CrossValidationSets;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Created by vinee_000 on 12/6/2017.
 */
public class realDataNew {

    public static void main(String [] args) throws Exception
    {

        //Three Groups, No Prior = run with the data alone (just to compute the PAM50 enrichment and prediction accuracy)
        //Irrelevant Prior = run with the Irrelevant PAM50 Genesets alone (also to compute PAM50 enrichment and prediction accuracy)
        //Relevant Prior = run with Irrelevant and actual PAM50 to compare weighting scores, and to compute Prediction accuracy of this network
        //Pathway Prior = run with just all of the pathways to compute weighting

        String priorDirectory = "";
        String dataFile = "";
        String runDir = "";
        int type = -1;
        // boolean computeMetrics = false;
        boolean useStabilities = false;

        int numLambda = 15;
        int numPriors = 5;
        double g = 0.05;

        double low = .05;
        double high = .95;
        int k = 10;
        int numSub = 10;



        System.out.print("Parsing Arguments...");

        int index = 0;
        while(index < args.length)
        {
            if(args[index].equals("-nl"))
            {
                numLambda = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-priors"))
            {
                priorDirectory = args[index+1];
                index+=2;
            }
            else if(args[index].equals("-runName"))
            {
                runDir = args[index+1];
                index+=2;
            }
            else if(args[index].equals("-data"))
            {
                dataFile = args[index+1];
                index+=2;
            }
            else if(args[index].equals("-g"))
            {
                g = Double.parseDouble(args[index+1]);
                index+=2;
            }
            else if (args[index].equals("-us"))
            {
                useStabilities = true;
                index++;
            }
            else if(args[index].equals("-type"))
            {
                type = Integer.parseInt(args[index + 1]);
                index+=2;
            }
            else if(args[index].equals("-kfold"))
            {
                k = Integer.parseInt(args[index+1]);
                index+=2;
            }

        }

        System.out.println("Done");

        double [] lambda = new double[numLambda];
        double inc = (high-low)/(numLambda-1);
        for(int i = 0; i < numLambda;i++)
        {
            lambda[i] = i*inc + low;
        }
        DataSet data = null;


        System.out.print("Reading Subsamples...");
        data = MixedUtils.loadDataSet2(dataFile);
        //System.out.println(data);
        data = MixedUtils.completeCases(data);

        data = MixedUtils.getNonparanormalTransform(data);
        System.out.println(data);
        /*HashMap<Integer,Integer> dist = new HashMap<Integer,Integer>();
        int col = data.getColumn(data.getVariable("Subtype"));
        for(int j = 0; j < data.getNumRows();j++)
        {
            if(dist.get(data.getInt(j,col))==null)
                dist.put(data.getInt(j,col),1);
            else
                dist.put(data.getInt(j,col),dist.get(data.getInt(j,col))+1);

        }

        System.out.println(dist);*/
        //System.out.println(data);
        File f = new File(runDir);
        if(!f.exists())
            f.mkdir();
        f = new File(runDir + "/Subsamples");
        if(!f.exists())
            f.mkdir();

        int b = (int) Math.floor(10 * Math.sqrt(data.getNumRows()));
        if (b > data.getNumRows())
            b = 3*data.getNumRows() / 4;
        int [][] samps = null;
        boolean exit = false;
        while(!exit) {
            samps = StabilityUtils.subSampleNoReplacement(data.getNumRows(), b, numSub);
            exit = true;
            for(int i = 0; i < samps.length;i++)
            {
                if(runPriors.checkForVariance(data.subsetRows(samps[i]),data)!=-1)
                    exit = false;
            }
        }
        File temp = new File(runDir + "/Subsamples/Subsamples.txt");
        if(!temp.exists())
        {
            PrintStream out2 = new PrintStream(temp.getAbsolutePath());
            for(int i = 0; i < samps.length;i++) {

                for (int j = 0; j < samps[i].length; j++) {
                    if(j==samps[i].length-1)
                        out2.println(samps[i][j]);
                    else
                        out2.print(samps[i][j] + "\t");
                }
            }
            out2.flush();
            out2.close();
        }
        else
        {
            BufferedReader b2 = new BufferedReader(new FileReader(temp));
            for(int j = 0; j < numSub;j++)
            {
                String [] line = b2.readLine().split("\t");
                samps[j] = new int[line.length];
                for(int l = 0; l< line.length;l++)
                    samps[j][l] = Integer.parseInt(line[l]);
            }
        }

        boolean runNoPrior = false;
        if(priorDirectory.equals(""))
            runNoPrior=true;
        System.out.println("Done");
        //Test no Prior situation
        if(runNoPrior) {
            STEPS s = new STEPS(data, lambda, g, samps);
            Graph graph = s.runStepsPar();
            double[][] stab = s.stabilities;
            System.out.println("Cross Validating");
            CrossValidationSets cv = new CrossValidationSets(data,s.lastLambda,runDir+"/Subsamples",runDir,"Subtype",k,runDir);
            cv.crossValidate();
            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            PrintStream out;
            out = new PrintStream(runDir + "/Stabilities.txt");
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
            out = new PrintStream(runDir + "/Graph.txt");
            out.println(graph);
            out.flush();
            out.close();
        }
        else
        {
            File f2 = new File(priorDirectory);
            File[] priorFiles = f2.listFiles();
            numPriors = priorFiles.length;
            String [] priorNames = new String[numPriors];
            SparseDoubleMatrix2D[] priors = new SparseDoubleMatrix2D[numPriors];

            for(int j =0; j < priorFiles.length;j++)
            {
                priorNames[j] = priorFiles[j].getName();
            }
            for(int i = 0; i < numPriors;i++) {
                //priors[i] = new TetradMatrix(loadPrior(new File("prior_sources/Irr_Prior_" + i + ".txt"),data.getNumColumns()));
                priors[i] = new SparseDoubleMatrix2D(loadPrior(priorFiles[i],data.getNumColumns()));
                System.out.println(priors[i].rows() + "," + priors[i].columns());
            }
            System.out.println("Constructing lambdas...");
            mgmPriors m = new mgmPriors(samps.length,lambda,data,priors,samps);
            PrintStream lb = new PrintStream(runDir + "/Lambdas.txt");
            for(int i = 0; i < m.getLambdas().length;i++)
            {
                lb.println(Arrays.toString(m.getLambdas()[i]));
            }
            lb.flush();
            lb.close();
            System.out.println("Running Priors...");
            Graph g2 = m.runPriors();

            //TODO SAME AS ABOVE SHOULD BE PARALLELIZED BUT TIME CRUNCH
            System.out.println("Cross Validating");
            CrossValidationSets cv = new CrossValidationSets(data,m.lastNPLambda,m.lastWPLambda,runDir + "/Subsamples",runDir,"Subtype",m.lastHavePrior,k,runDir);
            cv.setPriors(priors);
            cv.crossValidate();
            //////////////////////////////////////////////////////////////


            /////////////////////
            double [][] stab = m.edgeScores;
            double [] weights = m.normalizedExpertWeights;
            double [] pValues = m.pValues;

            PrintStream w;
            w = new PrintStream(runDir + "/Weights.txt");
            for(int i = 0; i < weights.length;i++)
            {
                w.println("" + i  + "\t" + weights[i] + "\t" + pValues[i]);
            }
            w.flush();
            w.close();
            PrintStream out;
            out = new PrintStream(runDir + "/Stabilities.txt");
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
            out = new PrintStream(runDir + "/Graph.txt");
            out.println(g2);
            out.flush();
            out.close();


        }



        //Output metrics, precision and recall of PAM50 Genes (Check the neighborhood of the Subtype variable)
        //, Recovery of relevant pathway connections
    }

    public static double [][] loadPAM50(File data, int numVariables) throws Exception
    {
        double [][] temp = new double[numVariables][numVariables];
        BufferedReader b = new BufferedReader(new FileReader(data));
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
    public static double [][] loadStability(String data, int numVariables) throws Exception
    {

        double [][] temp = new double[numVariables][numVariables];
        BufferedReader b = new BufferedReader(new FileReader(new File(data)));
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
