package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.HashMap;

/**
 * The Purpose of this class is to produce the first figure for the experiments section of the Pref-Div Paper
 * How well can we evaluate the prior knowledge we are getting from both intensity and dissimilarity files?
 */
public class evaluatePriorRecovery {

    public static void main(String [] args) {


        //TODO Implement ability to switch experiments
        int experiment = 0; //0 -> prior evaluation, 1 -> Accuracy of chosen parameters vs Optimal, 2-> Comparing different summarization methods for clustered genes
        //3 -> Evaluating how good cross validation is to select alpha parameters, 4-> Sensitivity analysis on number of folds, number of subsamples, number of parameters, ss, numvars, etc.


        int numRuns = 20;
        int numGenes = 200;
        int sampleSize = 1000;
        double amountPrior = 0.6;//percentage of edges to have prior knowledge for
        boolean boot = false; //Should we use bootstrap samples for PiPrefDiv
        boolean loocv = false; //Should we use leave-one-out CV for PiPrefDiv
        boolean useCausalGraph = false; //Are we selecting genes to then use a causal graph to find connections or direct selection?
        int numSamples = 20; //Number of bootstrap/sub-sampled samples to use
        int numParams = 10;//Number of parameters to sweep over



        int numPriors = 10; //Number of prior knowledge sources
        int numReliable = 5; //Number of reliable sources
        int numComponents = 10; //How many components do we have for cluster simulation?
        int minTargetParents = 5; //How many true parents of the target are there?
        boolean amountRandom = false; //Should the priors have a random amount of prior knowledge?
        boolean targetContinuous = true; //Is the target variable continuous?
        boolean evenDistribution = true; //Is the distribution of nodes in each cluster even?
        int numCategories = 4; //number of categories for discrete variables


        File rFile = new File("Results");
        if (!rFile.isDirectory())
            rFile.mkdir();

        try {
            PrintStream out = new PrintStream("Results/Prior_Evaluation_" + numGenes + "_" + sampleSize + "_" + numRuns + "_" + numPriors + "_" + amountPrior + "_" + amountRandom + "_" + boot + "_" + loocv + "_"  + evenDistribution + "_" + targetContinuous + "_" + useCausalGraph +  "_" + numReliable + "_" + numComponents + "_" + minTargetParents + "_" + numParams + ".txt");

            out.println("Run\tAmount_Prior\tReliable?\tIntensity?\tPredicted_Weight\tActual_Reliability");
            for (int j = 0; j < numRuns; j++) {
                HashMap<String, Integer> clusters = new HashMap<String, Integer>();

                System.out.print("Loading graph and data for run " + j + "...");
                Graph g = null;
                DataSet d = null;
                File graphFile = new File("Graphs/Graph_" + numGenes + "_" + minTargetParents + "_" + numComponents + "_" + evenDistribution + "_" + j + ".txt");
                if (graphFile.exists())
                    g = GraphUtils.loadGraphTxt(graphFile);
                else {
                    System.err.println("Couldn't find graph file with these parameters, please double check");
                }
                File dataFile = null;
                dataFile = new File("Data/Dataset_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + j + ".txt");
                try {
                    if (dataFile.exists())
                        d = MixedUtils.loadDataSet2(dataFile.getAbsolutePath());
                    else {
                        System.err.println("Could not load data with specified parameters, please double check");
                        System.exit(-1);
                    }
                } catch (Exception e) {
                    System.err.println("Unable to load dataset");
                    System.exit(-1);
                }
                System.out.println("Done");


                System.out.print("Loading clusters for run " + j + "...");
                File clusterFile = null;
                clusterFile = new File("Clusters/Cluster_" + numGenes + "_" + minTargetParents + "_" + numComponents + "_" + evenDistribution + "_" + j + ".txt");
                if (clusterFile.exists()) {
                    clusters = new HashMap<String, Integer>();
                    try {
                        BufferedReader bTemp = new BufferedReader(new FileReader(clusterFile));
                        while (bTemp.ready()) {
                            String[] line = bTemp.readLine().split("\t");
                            clusters.put(line[0], Integer.parseInt(line[1]));
                        }
                    } catch (Exception e) {
                        System.err.println("Couldn't load clusters from file");
                        e.printStackTrace();
                        System.exit(-1);
                    }

                } else {
                    System.err.println("Could not load cluster file with specified parameters");
                    System.exit(-1);
                }
                System.out.println("Done");


                File subsFile = null;
                subsFile = new File("Subsamples/Subsample_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + boot + "_" + numSamples + "_" + loocv + "_" + j + ".txt");
                int[][] subs = null;
                if (subsFile.exists()) {
                    try {
                        BufferedReader b = new BufferedReader(new FileReader(subsFile));
                        int lineCount = 0;
                        while (b.ready()) {
                            b.readLine();
                            lineCount++;
                        }
                        subs = new int[lineCount][];
                        b = new BufferedReader(new FileReader(subsFile));
                        lineCount = 0;
                        while (b.ready()) {
                            String[] line = b.readLine().split("\t");
                            subs[lineCount] = new int[line.length];
                            for (int x = 0; x < line.length; x++) {
                                subs[lineCount][x] = Integer.parseInt(line[x]);
                            }
                            lineCount++;
                        }
                        b.close();
                    } catch (Exception e) {
                        System.err.println("Unable to load subsamples");
                        e.printStackTrace();
                        System.exit(-1);
                    }
                } else {
                    System.err.println("Could not load subsamples file with specifieid parameters, please double check");
                    System.exit(-1);
                }


                /***LOAD PRIOR KNOWLEDGE FILES***/
                String priorIntensity = "Priors/Prior_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + numPriors + "_"  + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" +  j + "_intensity.txt";
                String[] dFile = new String[numPriors];
                for (int k = 0; k < numPriors; k++) {
                    dFile[k] = "Priors/Prior_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + k + "_" + j +  "_dissimilarity.txt";
                }

                /***LOAD RELIABILITY SCORES***/
                File rd = new File("Reliabilities_D_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + j + ".txt");
                File ri = new File("Reliabilities_I_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + j + ".txt");
                double [] reliDis = new double[numPriors];
                double [] reliInt = new double[numPriors];
                BufferedReader b = new BufferedReader(new FileReader(rd));
                BufferedReader b2 = new BufferedReader(new FileReader(ri));

                for(int i = 0; i < numPriors;i++)
                {
                    reliDis[i] = Double.parseDouble(b.readLine());
                    reliInt[i] = Double.parseDouble(b2.readLine());
                }
                b.close();
                b2.close();


                /***RUN PI PREF DIV***/
                PiPrefDiv p = new PiPrefDiv(d, "Target", minTargetParents,numParams);
                p.setSubsamples(subs);
                p.setVerbose();
                double[][] weights = p.evaluatePriors(boot, numSamples, priorIntensity, dFile, useCausalGraph);

                /***PRINT RESULTS TO FILE***/
                for(int i = 0; i < numPriors;i++)
                {
                    out.println(j + "\t" + getAmountPrior(priorIntensity,i)/((double)g.getNumNodes()-1) + "\t" + (i<numReliable) + "\tT\t" + weights[0][i] + "\t" + reliInt[i]);
                }
                for(int i = 0; i < numPriors;i++)
                {
                    out.println(j + "\t" + (double)(getAmountPrior(dFile[i]))/(g.getNumNodes()*(g.getNumNodes()-1)/2) + "\t" + (i<numReliable) + "\tF\t" + weights[1][i] + "\t" + reliDis[i]);
                }

                out.flush();

            }
            out.close();
        }
        catch(Exception e)
        {
            System.err.println("Error with printstream");
            e.printStackTrace();
            System.exit(-1);
        }
    }


    //Dissimilarity (in multiple files), this loads from a single one and you call it many times
    public static int getAmountPrior(String file)
    {
        try{
            int result = 0;
            BufferedReader b = new BufferedReader(new FileReader(new File(file)));
            int count = 1;
            while(b.ready())
            {
                String [] line = b.readLine().split("\t");
                for(int i = count;i<line.length;i++)
                {
                    if(Double.parseDouble(line[i])!=-1)
                        result++;
                }
                count++;
            }
            return result;
        }
        catch(Exception e)
        {
            System.err.println("Couldn't compute amount of prior information");
            e.printStackTrace();
            System.exit(-1);
        }
        return -1;

    }

    //Intensity (all in one file)
    public static int getAmountPrior(String file,int k)
    {
        try{
            int result = 0;
            BufferedReader b = new BufferedReader(new FileReader(new File(file)));
            b.readLine();//Eat the header
            while(b.ready())
            {
                if(Double.parseDouble(b.readLine().split("\t")[k+1])!=-1)
                    result++;
            }
            return result;
        }
        catch(Exception e)
        {
            System.err.println("Couldn't compute amount of prior information");
            e.printStackTrace();
            System.exit(-1);
        }
        return -1;
    }

}
