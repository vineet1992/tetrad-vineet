package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.data.DataSet;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * The purpose of this class is to evaluate the runtime ability of the PiPrefDiv method
 * and to ensure that the
 */
public class PPD_Runtime {


        static int numRuns = 1;
        static int numGenes = 10000;

        static boolean boot = false; //Should we use bootstrap samples for PiPrefDiv
        static boolean loocv = false; //Should we use leave-one-out CV for PiPrefDiv
        static boolean useCausalGraph = false; //Are we selecting genes to then use a causal graph to find connections or direct selection?
        static int numSamples = 20; //Number of bootstrap/sub-sampled samples to use
        static int numParams = 15;//Number of parameters to sweep over
        static boolean noiseRandom = true;//Is the reliability range of reliable and unreliable priors set randomly?
        static int numFolds = 5; //Number of Folds for CV to pick alpha


        static int numPriors = 10; //Number of prior knowledge sources
        static int numReliable = 5; //Number of reliable sources
        static int numComponents = 500; //How many components do we have for cluster simulation?
        static int minTargetParents = 100; //How many true parents of the target are there?
        static boolean amountRandom = false; //Should the priors have a random amount of prior knowledge?
        static boolean targetContinuous = true; //Is the target variable continuous?
        static boolean evenDistribution = true; //Is the distribution of nodes in each cluster even?
        static int numCategories = 4; //number of categories for discrete variables
        static boolean stabilitySelection = false; //Should stability selection be used within Pref-Div?
        static int sampleSize = 200;
        static double amountPrior = 0.3;
        static int[][] subs;
        static boolean parallel = false; //Should we run the experiment with parallel processing?
        static boolean partialCorr = false;
        static boolean pdStability = false;
        public static void main(String [] args) throws Exception
        {
            System.out.println("Sample Size: " + sampleSize + ", Amount Prior: " + amountPrior + ", Number Reliable: " + numReliable);
            File dataFile = null;
            DataSet d = null;
            dataFile = new File("Data/Dataset_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + 0 + ".txt");
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

            ArrayList<Gene> temp = PiPrefDiv2.createGenes(d,"Target");

            long time = System.nanoTime();
            float [] one = Functions.computeAllCorrelations(temp,d,false,false,false,0.05);
            double result = (System.nanoTime()-time)/Math.pow(10,9);

            time = System.nanoTime();
            float [] two = Functions.computeAllCorrelationsTest(temp,d,false,false,false,0.05);
            double result2 = (System.nanoTime()-time)/Math.pow(10,9);

            System.out.println("Time OG: " + result + ", Time New: " + result2);

            int diffs = 0;
            for(int i = 0; i < one.length;i++)
            {
                if(Math.abs(one[i]-two[i])>0.0001)
                    diffs++;
            }
            System.out.println("Differences: " + diffs);
            System.exit(-1);



            /***LOAD PRIOR KNOWLEDGE FILES***/
            String priorIntensity = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + 0 + "_intensity.txt";
            String[] dFile = new String[numPriors];
            for (int k = 0; k < numPriors; k++) {
                dFile[k] = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + k + "_" + 0 + "_dissimilarity.txt";
            }

            PiPrefDiv2 p = new PiPrefDiv2(d, "Target", minTargetParents, numParams);
            p.setSubsamples(subs);
            //p.setVerbose();
            p.setUseStabilitySelection(stabilitySelection);
            p.setParallel(false);
            p.setPartialCorrs(partialCorr);
            p.setPdStability(pdStability);
            ArrayList<Gene> selected = p.selectGenes(boot, numSamples, priorIntensity, dFile, useCausalGraph);




        }
}
