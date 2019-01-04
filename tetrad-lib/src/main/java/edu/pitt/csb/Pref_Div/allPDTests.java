package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

/**
 * The Purpose of this class is to produce the first figure for the experiments section of the Pref-Div Paper
 * How well can we evaluate the prior knowledge we are getting from both intensity and dissimilarity files?
 */
public class allPDTests {



    static int numRuns = 10;
    static int numGenes = 2000;

    static boolean boot = false; //Should we use bootstrap samples for PiPrefDiv
    static boolean loocv = false; //Should we use leave-one-out CV for PiPrefDiv
    static boolean useCausalGraph = false; //Are we selecting genes to then use a causal graph to find connections or direct selection?
    static int numSamples = 20; //Number of bootstrap/sub-sampled samples to use
    static int numParams = 15;//Number of parameters to sweep over
    static boolean noiseRandom = true;//Is the reliability range of reliable and unreliable priors set randomly?
    static int numFolds = 5; //Number of Folds for CV to pick alpha


    static int numPriors = 10; //Number of prior knowledge sources
    static int numReliable = 5; //Number of reliable sources
    static int numComponents = 200; //How many components do we have for cluster simulation?
    static int minTargetParents = 50; //How many true parents of the target are there?
    static boolean amountRandom = false; //Should the priors have a random amount of prior knowledge?
    static boolean targetContinuous = true; //Is the target variable continuous?
    static boolean evenDistribution = true; //Is the distribution of nodes in each cluster even?
    static int numCategories = 4; //number of categories for discrete variables
    static boolean stabilitySelection = false; //Should stability selection be used within Pref-Div?
    static int sampleSize = 200;
    static double amountPrior = 0.3;
    static int[][] subs;
    static boolean parallel = true; //Should we run the experiment with parallel processing?
    static boolean partialCorr = false;
    static boolean pdStability = false;
    public static void main(String [] args) {


        double [] ap = new double[]{0.1,0.3,0.6};
        int [] nr = new int[]{1,5,10};
        int [] ss = new int[]{200};

        for(int ii = 0; ii < ap.length;ii++) {
            amountPrior = ap[ii];
            for (int jj = 0; jj < nr.length; jj++) {
                numReliable = nr[jj];
                for (int kk = 0; kk < ss.length; kk++) {
                    sampleSize = ss[kk];



                    int experiment = 0; //0 -> prior evaluation, 1 -> Accuracy of chosen parameters vs Optimal,
                    //2-> Comparing separate prior parameters vs not, 3-> Comparing different summarization methods for clustered genes
                    //4 -> Runtime profiling


                    String[] types = new String[]{"FS", "Prediction", "Cluster"}; //Help to create file header for experiment 1
                    String[] starts = new String[]{"Best_Radius", "Best_Threshold", "Best_Accuracy", "Pred_Accuracy"}; //More help
                    RunPrefDiv.ClusterType[] allTypes = RunPrefDiv.ClusterType.values();


                    /***Create Results directory***/
                    File rFile = new File("Results");
                    if (!rFile.isDirectory())
                        rFile.mkdir();

                    /***Create special files for this particular experiment***/
                    try {
                        String start = "Prior_Evaluation_";
                        if (experiment == 1) {
                            start = "Parameter_Accuracy_";
                        } else if (experiment == 2) {
                            start = "Separate_Param_Comparison_";
                        } else if (experiment == 3) {
                            start = "Summarization_Comparison_";
                        }

                        PrintStream out = new PrintStream("Results/" + start + numGenes + "_" + sampleSize + "_" + numRuns + "_" + numPriors + "_" + amountPrior + "_" + amountRandom + "_" + boot + "_" + loocv + "_" + evenDistribution + "_" + targetContinuous + "_" + useCausalGraph + "_" + numReliable + "_" + numComponents + "_" + minTargetParents + "_" + numParams + "_" + stabilitySelection + ".txt");



                        if (experiment == 0) {
                            out.println("Run\tAmount_Prior\tReliable?\tIntensity?\tPredicted_Weight\tActual_Reliability");
                        } else if (experiment == 1) {
                            out.print("Run\tPredicted_Radius_NP\tPredicted_Radius_WP\tPredicted_Threshold_NP\tPredicted_Threshold_WP\tBest_Radius_NP\tBest_Threshold_NP\t");

                            for (int i = 0; i < types.length; i++) {
                                for (int j = 0; j < starts.length; j++) {
                                    out.print(starts[j] + "_" + types[i] + "\t");
                                }
                            }
                            out.println();
                        } else if (experiment == 2) {
                            out.println("Run\tRadius_Both\tThreshold_Both\tBoth_Feat_Accuracy\tBoth_Pred_Accuracy\tBoth_Clust_Accuracy\tRadius_WP\tRadius_NP\tThreshold_WP\tThreshold_NP\tSeparate_Feat_Accuracy\tSeparate_Pred_Accuracy\tSeparate_Clust_Accuracy");
                        } else if (experiment == 3) {
                            out.print("Run\t");
                            for (int i = 0; i < allTypes.length; i++)
                                out.print(allTypes[i] + "\t");
                            out.println();
                        }


                        /***Main Loop going through numRuns times***/
                        for (int j = 0; j < numRuns; j++) {
                            HashMap<String, Integer> clusters = new HashMap<String, Integer>();


                            /***Load Graph and Data from the files***/
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


                            /***Load Clusters from file***/
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


                            /***Load subsamples from file***/
                            File subsFile = null;
                            if (experiment == 0)
                                subsFile = new File("Subsamples/Subsample_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + boot + "_" + numSamples + "_" + loocv + "_" + j + ".txt");
                            else
                                subsFile = new File("Subsamples/Subsample_Train_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + boot + "_" + numSamples + "_" + loocv + "_" + j + ".txt");

                            subs = null;
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
                            } else if (experiment == 0) {
                                System.err.println("Could not load subsamples file with specified parameters, please double check");
                            } else {
                                System.out.println("Couldn't load " + subsFile);
                            }


                            /***LOAD PRIOR KNOWLEDGE FILES***/
                            String priorIntensity = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + "_intensity.txt";
                            String[] dFile = new String[numPriors];
                            for (int k = 0; k < numPriors; k++) {
                                dFile[k] = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + k + "_" + j + "_dissimilarity.txt";
                            }


                            /***RUN PI PREF DIV***/
                            DataSet toRun = d;
                            DataSet test = d;
                            if (experiment >= 1) {
                                /***Load Train-test split***/
                                int trainLength = (int) (d.getNumRows() * 0.9);
                                int[] trainInds = new int[trainLength];
                                for (int i = 0; i < trainLength; i++) {
                                    trainInds[i] = i;
                                }
                                toRun = d.subsetRows(trainInds);
                                test.removeRows(trainInds);
                            }

                            PiPrefDiv2 p = new PiPrefDiv2(toRun, "Target", minTargetParents, numParams);
                            p.setSubsamples(subs);
                            //p.setVerbose();
                            p.setUseStabilitySelection(stabilitySelection);
                            p.setParallel(false);
                            p.setPartialCorrs(partialCorr);
                            p.setPdStability(pdStability);


                            String dir = "Results/Detailed_Score_Evaluation/";
                            File f = new File(dir);
                            if (experiment == 1) {
                                if (!f.exists())
                                    f.mkdir();
                                p.setOutputScores(new PrintStream(dir + "Parameter_Scores_" + j + "_" + numGenes + "_" + sampleSize + "_" + numRuns + "_" + numPriors + "_" + amountPrior + "_" + amountRandom + "_" + boot + "_" + loocv + "_" + evenDistribution + "_" + targetContinuous + "_" + useCausalGraph + "_" + numReliable + "_" + numComponents + "_" + minTargetParents + "_" + numParams + "_" + stabilitySelection + ".txt"));
                            }


                            if (experiment < 0) //Current debugging test
                            {
                                ArrayList<Gene> selected = p.selectGenes(boot, numSamples, priorIntensity, dFile, useCausalGraph);
                                Map<Gene, List<Gene>> lastCluster = p.getLastCluster();

                            }
                            if (experiment == 0) //Evaluate prior knowledge
                            {
                                /***LOAD RELIABILITY SCORES***/
                                File rd = new File("Reliabilities_D_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");
                                File ri = new File("Reliabilities_I_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");
                                double[] reliDis = new double[numPriors];
                                double[] reliInt = new double[numPriors];
                                BufferedReader b = new BufferedReader(new FileReader(rd));
                                BufferedReader b2 = new BufferedReader(new FileReader(ri));

                                for (int i = 0; i < numPriors; i++) {
                                    reliDis[i] = Double.parseDouble(b.readLine());
                                    reliInt[i] = Double.parseDouble(b2.readLine());
                                }
                                b.close();
                                b2.close();
                                long time = System.nanoTime();
                                p.selectGenes(boot, numSamples, priorIntensity, dFile, useCausalGraph);
                                double[][] weights = new double[2][numPriors];
                                weights[0] = p.getLastIntensityWeights();
                                weights[1] = p.getLastSimilarityWeights();

                                time = System.nanoTime() - time;
                                System.out.println("Took " + time / Math.pow(10, 9) + " seconds");
                                /***PRINT RESULTS TO FILE***/
                                for (int i = 0; i < numPriors; i++) {
                                    out.println(j + "\t" + getAmountPrior(priorIntensity, i) / ((double) g.getNumNodes() - 1) + "\t" + (i < numReliable) + "\tT\t" + weights[0][i] + "\t" + reliInt[i]);
                                }
                                for (int i = 0; i < numPriors; i++) {
                                    out.println(j + "\t" + (double) (getAmountPrior(dFile[i])) / (g.getNumNodes() * (g.getNumNodes() - 1) / 2) + "\t" + (i < numReliable) + "\tF\t" + weights[1][i] + "\t" + reliDis[i]);
                                }
                            }
                            /***Evaluate Parameter Accuracy, by sweeping over potential parameter choices***/
                            else if (experiment == 1) //Evaluate Parameter accuracy
                            {

                                ArrayList<Gene> selected = p.selectGenes(boot, numSamples, priorIntensity, dFile, useCausalGraph);
                                System.out.println("Selected Genes: " + selected);
                                System.out.println(g.getAdjacentNodes(g.getNode("Target")));
                                Map<Gene, List<Gene>> lastCluster = p.getLastCluster();
                                System.out.println("Cluster: " + lastCluster);
                                double clustAccuracy = getClusterAccuracy(lastCluster, clusters, numComponents);
                                double featAccuracy = -1;
                                if (minTargetParents == g.getAdjacentNodes(g.getNode("Target")).size()) {
                                    //Correct amount of selections so use accuracy
                                    featAccuracy = getFeatAccuracy(selected, g, "Target", "ACC");

                                } else {
                                    //Incorrect amount of selected variables so use F1
                                    featAccuracy = getFeatAccuracy(selected, g, "Target", "F1");
                                }
                                double predAccuracy = getAccuracy(toRun, test, selected);
                                System.out.println("Run #" + j + ": Prediction=" + predAccuracy + ", Features=" + featAccuracy + ", Cluster=" + clustAccuracy);

                                double[] radii = p.getTestedRadii();
                                double[] thresholds = p.getTestedThresholds();


                                ArrayList<Gene> meanGenes = PiPrefDiv.createGenes(d, "Target");
                                meanGenes = Functions.computeAllIntensitiesWithP(meanGenes, 1, d, "Target", false);

                                //0 is correlation, 1 is p-value
                                float[][] meanDis = Functions.computeAllCorrelationsWithP(meanGenes, d, false, false, false, 1);
                                //Use the weights from PiPrefDiv, we are interested in determining how good param selection is given that weights are good
                                //Because experiment 1 demonstrated that weights are good


                                double[][][] results;
                                double[] npParams;
                                //for every radii and threshold combination, give the three accuracies

                                if (!parallel) {
                                    npParams = getBestNPParams(radii, thresholds, numComponents, meanDis, meanGenes, p.getDissimilarityPriors(), clusters, g, toRun, test); //Feature Select, Cluster Acc, Prediction Acc
                                    results = getBestParams(npParams[0], npParams[1], radii, thresholds, numComponents, meanDis, meanGenes, p.getIntensityPriors(), p.getDissimilarityPriors(), clusters, g, toRun, test);
                                } else {
                                    npParams = getBestNPParamsPar(radii, thresholds, numComponents, meanDis, meanGenes, p.getDissimilarityPriors(), clusters, g, toRun, test);
                                    results = getBestParamsPar(npParams[0], npParams[1], radii, thresholds, numComponents, meanDis, meanGenes, p.getIntensityPriors(), p.getDissimilarityPriors(), clusters, g, toRun, test);
                                }


                                out.print(j + "\t" + p.getLastRadius()[0] + "\t" + p.getLastRadius()[1] + "\t" + p.getLastThreshold()[0] + "\t" + p.getLastThreshold()[1] + "\t");

                                //print out results to file
                                printResults(results, out, radii, thresholds, new double[]{featAccuracy, predAccuracy, clustAccuracy}, npParams[0], npParams[1]);

                                printDetailedScores(results, new PrintStream("Results/Detailed_Results_" + j + "_" + numGenes + "_" + sampleSize + "_" + numRuns + "_" + numPriors + "_" + amountPrior + "_" + amountRandom + "_" + boot + "_" + loocv + "_" + evenDistribution + "_" + targetContinuous + "_" + useCausalGraph + "_" + numReliable + "_" + numComponents + "_" + minTargetParents + "_" + numParams + "_" + stabilitySelection + ".txt"));


                                //Roll through all tested radii,threshold combos and RunPrefDiv directly to get accuracies
                            }

                            /***Evaluate using priors to choose parameters vs. choosing separate parameters for with/without prior information***/
                            else if (experiment == 2) {
                                out.print(j + "\t");
                                ArrayList<Gene> selected = p.selectGenes(boot, numSamples, priorIntensity, dFile, useCausalGraph);
                                System.out.println("Selected Genes Separate: " + selected);
                                Map<Gene, List<Gene>> lastCluster = p.getLastCluster();
                                System.out.println("Cluster: " + lastCluster);
                                double clustAccuracy = getClusterAccuracy(lastCluster, clusters, numComponents);
                                double featAccuracy = -1;
                                if (minTargetParents == g.getAdjacentNodes(g.getNode("Target")).size()) {
                                    //Correct amount of selections so use accuracy
                                    featAccuracy = getFeatAccuracy(selected, g, "Target", "ACC");

                                } else {
                                    //Incorrect amount of selected variables so use F1
                                    featAccuracy = getFeatAccuracy(selected, g, "Target", "F1");
                                }
                                double predAccuracy = getAccuracy(toRun, test, selected);

                                out.print(p.getLastRadius() + "\t" + p.getLastThreshold() + "\t");
                                out.print(featAccuracy + "\t" + predAccuracy + "\t" + clustAccuracy + "\t");
                                System.out.println("Unified Params...Run #" + j + ": Prediction=" + predAccuracy + ", Features=" + featAccuracy + ", Cluster=" + clustAccuracy);

                                System.out.println("Running PiPrefDiv with Separate params...");
                                selected = p.selectGenes(boot, numSamples, priorIntensity, dFile, useCausalGraph);
                                System.out.println("Selected Genes Separate: " + selected);
                                lastCluster = p.getLastCluster();
                                System.out.println("Cluster: " + lastCluster);
                                clustAccuracy = getClusterAccuracy(lastCluster, clusters, numComponents);
                                featAccuracy = -1;
                                if (minTargetParents == g.getAdjacentNodes(g.getNode("Target")).size()) {
                                    //Correct amount of selections so use accuracy
                                    featAccuracy = getFeatAccuracy(selected, g, "Target", "ACC");

                                } else {
                                    //Incorrect amount of selected variables so use F1
                                    featAccuracy = getFeatAccuracy(selected, g, "Target", "F1");
                                }
                                predAccuracy = getAccuracy(toRun, test, selected);


                                System.out.println("Separate Params...Run #" + j + ": Prediction=" + predAccuracy + ", Features=" + featAccuracy + ", Cluster=" + clustAccuracy);
                                out.print(p.getLastRadiusWP() + "\t" + p.getLastRadius() + "\t" + p.getLastThresholdWP() + "\t" + p.getLastThreshold() + "\t");
                                out.println(featAccuracy + "\t" + predAccuracy + "\t" + clustAccuracy);

                                out.flush();

                            }
                            /***Evaluate different summarization methods***/
                            else if (experiment == 3) {
                                out.print(j + "\t");
                                ArrayList<Gene> selected = p.selectGenes(boot, numSamples, priorIntensity, dFile, useCausalGraph);
                                Map<Gene, List<Gene>> lastCluster = p.getLastCluster();

                                for (int i = 0; i < allTypes.length; i++) {

                                    DataSet summarized = RunPrefDiv.summarizeData(toRun, selected, lastCluster, allTypes[i]);
                                    DataSet summarizedTest = RunPrefDiv.summarizeData(test, selected, lastCluster, allTypes[i]);
                                    double predAccuracy = getAccuracy(summarized, summarizedTest, null);

                                    out.print(predAccuracy + "\t");


                                }
                                out.println();
                            }

                            out.flush();

                        }
                        out.close();

                    } catch (Exception e) {
                        System.err.println("Error with printstream");
                        e.printStackTrace();
                        System.exit(-1);
                    }
                }
            }
        }


    }

    public static void printDetailedScores(double[][][]results, PrintStream out)
    {
        for(int k = 0; k < 3;k++) {
            for (int i = 0; i < results.length; i++) {
                for (int j = 0; j < results[i].length; j++) {
                    if(j==results[i].length-1)
                        out.println(results[i][j][k]);
                    else
                        out.print(results[i][j][k] + "\t");
                }
                out.println();
            }
        }
    }

    public static void printResults(double[][][]results, PrintStream out, double [] radii, double[] thresholds, double [] predAcc, double radiusNP, double threshNP)
    {
        out.print(radiusNP + "\t" + threshNP + "\t");

        /**Feature Selection Results Printing**/
        double[] bestRadii = new double[3];
        double []maxValue = new double[3];
        maxValue[1] = Double.MAX_VALUE;
        double [] bestThreshold = new double[3];
        for(int i = 0; i < results.length;i++)
        {
            for(int j = 0; j < results[i].length;j++)
            {
                for(int k = 0; k < 3;k++)
                {
                    if((k!=1 && results[i][j][k]>maxValue[k]) || (k==1 && results[i][j][k] < maxValue[k]))
                    {
                        maxValue[k] = results[i][j][k];
                        bestRadii[k] = radii[i];
                        bestThreshold[k] = thresholds[j];
                    }
                }
            }
        }
        for(int k = 0; k < 3;k++)
        {
            if(k==2)
                out.println(radii[k] + "\t" + thresholds[k] + "\t" + maxValue[k] + "\t" + predAcc[k]);
            else
                out.print(radii[k] + "\t" + thresholds[k] + "\t" + maxValue[k] + "\t" + predAcc[k] + "\t");
        }
    }



    public static double[][][] getBestParamsPar(final double radiiNP, final double  thresholdNP, final double[] radii, final double[] threshold, final int numComponents,final float[][] corrs, final ArrayList<Gene> meanGenes, final boolean[] iPriors,final boolean[] dPriors, final Map<String,Integer>clusters, final Graph g, final DataSet toRun,final DataSet test)
    {
        final double[][][] result = new double[radii.length][threshold.length][3];
        final ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

        class StabilityAction extends RecursiveAction {
            private int chunk;
            private int from;
            private int to;

            public StabilityAction(int chunk, int from, int to){
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }


            private synchronized void addAccuracy(List<Gene>selected,Map<Gene,List<Gene>> lastCluster,int i,int j)
            {
                result[i][j][2] = getClusterAccuracy(lastCluster,clusters,numComponents);

                String type = "F1";
                if(selected.size()==g.getAdjacentNodes(g.getNode("Target")).size())
                    type = "ACC";
                result[i][j][0] = getFeatAccuracy(selected,g,"Target",type);
                result[i][j][1] = getAccuracy(toRun,test,selected);
            }

            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        try {

                            int i = s / radii.length;
                            int j = s % threshold.length;
                            ArrayList<Gene> currGenes = shrinkByThreshold(meanGenes,iPriors,(float)Math.exp(thresholdNP),(float)Math.exp(threshold[j]));
                            float [] currCorr = shrinkByThreshold(corrs,dPriors,(float)Math.exp(thresholdNP),(float)Math.exp(threshold[j]));
                            Collections.sort(currGenes,Gene.IntensityComparator);
                            PrefDiv pd = new PrefDiv(currGenes,minTargetParents,0,radiiNP,(float)radii[i],dPriors,currCorr);
                            pd.setCluster(true);
                            List<Gene> selected = pd.diverset();
                            Map<Gene,List<Gene>> lastCluster = pd.clusters;
                            addAccuracy(selected,lastCluster,i,j);


                        }catch(Exception e)
                        {
                            e.printStackTrace();
                            System.exit(-2);
                        }
                    }

                    return;
                } else {
                    List<StabilityAction> tasks = new ArrayList<>();

                    final int mid = (to + from) / 2;

                    tasks.add(new StabilityAction(chunk, from, mid));
                    tasks.add(new StabilityAction(chunk, mid, to));

                    invokeAll(tasks);

                    return;
                }
            }
        }

        final int chunk = 5;
        StabilityAction sa = new StabilityAction(chunk,0, radii.length*threshold.length);
        pool.invoke(sa);
        return result;
    }



    private static double[] getBestNPParamsPar(final double [] radii,final double [] threshold, final int numComponents,final float[][]corrs, final ArrayList<Gene> meanGenes, final boolean [] dPriors, final  Map<String,Integer>clusters, final Graph g, final DataSet toRun,final DataSet test)
    {
        final double[][][] result = new double[radii.length][threshold.length][3];
        final ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();



        class StabilityAction extends RecursiveAction {
            private int chunk;
            private int from;
            private int to;

            public StabilityAction(int chunk, int from, int to){
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }


            private synchronized void addAccuracy(List<Gene>selected,Map<Gene,List<Gene>> lastCluster,int i,int j)
            {
                result[i][j][2] = getClusterAccuracy(lastCluster,clusters,numComponents);

                String type = "F1";
                if(selected.size()==g.getAdjacentNodes(g.getNode("Target")).size())
                    type = "ACC";
                result[i][j][0] = getFeatAccuracy(selected,g,"Target",type);
                result[i][j][1] = getAccuracy(toRun,test,selected);
            }

            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        try {

                            int i = s / radii.length;
                            int j = s % threshold.length;

                            ArrayList<Gene> currGenes = shrinkByThreshold(meanGenes,(float)Math.exp(threshold[j]));
                            float [] currCorr = shrinkByThreshold(corrs, (float)Math.exp(threshold[j]));
                            Collections.sort(currGenes,Gene.IntensityComparator);


                            PrefDiv pd = new PrefDiv(currGenes,minTargetParents,0,(float)radii[i],(float)radii[i],dPriors,currCorr);
                            pd.setCluster(true);
                            List<Gene> selected = pd.diverset();
                            Map<Gene,List<Gene>> lastCluster = pd.clusters;
                            addAccuracy(selected,lastCluster,i,j);


                        }catch(Exception e)
                        {
                            e.printStackTrace();
                            System.exit(-2);
                        }
                    }

                    return;
                } else {
                    List<StabilityAction> tasks = new ArrayList<>();

                    final int mid = (to + from) / 2;

                    tasks.add(new StabilityAction(chunk, from, mid));
                    tasks.add(new StabilityAction(chunk, mid, to));

                    invokeAll(tasks);

                    return;
                }
            }
        }

        final int chunk = 5;
        StabilityAction sa = new StabilityAction(chunk,0, radii.length*threshold.length);
        pool.invoke(sa);



        return bestAcc(result,radii,threshold);
    }


    private static double [] bestAcc(double[][][]result, double [] radii, double[] threshold)
    {
        double bestAcc = 10;
        double bestRadii = -1;
        double bestThreshold = -1;
        for(int i = 0; i < radii.length;i++)
        {
            for(int j = 0; j < threshold.length;j++)
            {
                if(result[i][j][1] < bestAcc)
                {
                    bestAcc = result[i][j][1];
                    bestRadii = radii[i];
                    bestThreshold = threshold[j];
                }
            }
        }
        return new double[]{bestRadii,bestThreshold};
    }

    private static double[] getBestNPParams(double [] radii, double [] threshold, int numComponents,float[][]corrs, ArrayList<Gene> genes, boolean[]dPrior, Map<String,Integer> clusters,Graph g, DataSet toRun,DataSet test)
    {
        double[][][] result = new double[radii.length][threshold.length][3];

        //Find optimal parameters for no priors, and then loop again for with priors
            for(int j =0 ; j < threshold.length;j++)
            {
                float [] currCorr = shrinkByThreshold(corrs,(float)Math.exp(threshold[j]));
                for (int i = 0; i < radii.length; i++)
                {
                    ArrayList<Gene> currGenes = shrinkByThreshold(genes,(float)Math.exp(threshold[j]));

                        Collections.sort(currGenes,Gene.IntensityComparator);
                        PrefDiv pd = new PrefDiv(currGenes,minTargetParents,0,(float)radii[i],(float)radii[i],dPrior,currCorr);
                        pd.setCluster(true);
                        List<Gene> selected = pd.diverset();
                        Map<Gene, List<Gene>> lastCluster = pd.clusters;
                        result[i][j][2] = getClusterAccuracy(lastCluster, clusters, numComponents);

                        String type = "F1";
                        if (selected.size() == g.getAdjacentNodes(g.getNode("Target")).size())
                            type = "ACC";
                        result[i][j][0] = getFeatAccuracy(selected, g, "Target", type);
                        result[i][j][1] = getAccuracy(toRun, test, selected);

              }
            }
        return bestAcc(result,radii,threshold);
    }


    private static double[][][] getBestParams(double radiiNP, double thresholdNP, double [] radii, double [] threshold, int numComponents,float[][]corrs, ArrayList<Gene> genes, boolean [] iPrior, boolean[]dPrior, Map<String,Integer> clusters,Graph g, DataSet toRun,DataSet test)
    {
        double[][][] result = new double[radii.length][threshold.length][3];

        //Find optimal parameters for no priors, and then loop again for with priors
        for(int j =0 ; j < threshold.length;j++)
        {
            float [] currCorr = shrinkByThreshold(corrs,dPrior,(float)Math.exp(thresholdNP),(float)Math.exp(threshold[j]));
            for (int i = 0; i < radii.length; i++)
            {
                ArrayList<Gene> currGenes = shrinkByThreshold(genes,iPrior,(float)Math.exp(thresholdNP),(float)Math.exp(threshold[j]));
                Collections.sort(currGenes,Gene.IntensityComparator);
                PrefDiv pd = new PrefDiv(currGenes,minTargetParents,0,(float)radiiNP,(float)radii[i],dPrior,currCorr);
                pd.setCluster(true);
                List<Gene> selected = pd.diverset();
                Map<Gene, List<Gene>> lastCluster = pd.clusters;
                result[i][j][2] = getClusterAccuracy(lastCluster, clusters, numComponents);

                String type = "F1";
                if (selected.size() == g.getAdjacentNodes(g.getNode("Target")).size())
                    type = "ACC";
                result[i][j][0] = getFeatAccuracy(selected, g, "Target", type);
                result[i][j][1] = getAccuracy(toRun, test, selected);

               // System.out.println(selected);
                //System.out.println("Radius: " + radii[i] + ", Threshold: " + threshold[j] + ", Feature: " + result[i][j][0] + ", Prediction: " + result[i][j][1] + ", Cluster: " + result[i][j][2]);
            }
        }
        return result;
    }


    //TODO Incorporate discrete target variables
    public static double getAccuracy(DataSet train, DataSet test, List<Gene> selected)
    {
        RegressionDataset rd = new RegressionDataset(train);
        List<Node> regressors = new ArrayList<Node>();
        if(selected!=null) {
            for (int i = 0; i < selected.size(); i++)
                regressors.add(train.getVariable(selected.get(i).symbol));
        }
        else
        {
            for (int i = 0; i < train.getNumColumns(); i++) {
                if(train.getVariable(i).getName().equals("Target"))
                    continue;
                regressors.add(train.getVariable(i));
            }
        }

        RegressionResult res = rd.regress(train.getVariable("Target"),regressors);
        int [] cols = new int[regressors.size()];

        for(int i = 0; i < regressors.size();i++)
        {
            cols[i] = test.getColumn(test.getVariable(regressors.get(i).getName()));
        }
        double[]pred = new double[test.getNumRows()];
        double [] actual = new double[test.getNumRows()];
        for(int i = 0; i < test.getNumRows();i++)
        {
            double [] x = new double[cols.length];
            for(int j = 0; j < cols.length;j++) {
                x[j] = test.getDouble(i, cols[j]);

            }
            pred[i] = res.getPredictedValue(x);
            actual[i] = test.getDouble(i,test.getColumn(test.getVariable("Target")));


        }

        return RMSE(pred,actual);

        //Prediction accuracy via linear regression
    }

    public static double RMSE(double [] pred, double[] actual)
    {
        double err = 0;
        for(int i = 0; i < pred.length;i++)
        {
            err += Math.pow(pred[i]-actual[i],2);
        }
        return Math.sqrt(err/pred.length);
    }

    //Computes cluster similarity using the Hungarian Algorithm approach
    public static double getClusterAccuracy(Map<Gene,List<Gene>> est, Map<String,Integer> actual, int numClusters)
    {
        List<List<Gene>> temp = new ArrayList<List<Gene>>();
        for(int i = 0; i < numClusters;i++)
            temp.add(new ArrayList<Gene>());
        int x = 0;
        for(String s: actual.keySet())
        {
            Gene g = new Gene(x);
            g.symbol = s;
            temp.get(actual.get(s)).add(g);
            x++;
        }
        HashMap<Gene,List<Gene>> real = new HashMap<Gene,List<Gene>>();
        for(int i = 0; i < temp.size();i++)
        {
            List<Gene> curr = temp.get(i);
            Gene key = curr.get(0);
            curr.remove(0);
            real.put(key,curr);
        }
        return RunPrefDiv.clusterSim(est,real);
    }


    private static synchronized float[] shrinkByThreshold(float[][]corrs, boolean [] dPrior, float pThreshNP, float pThresh)
    {
        float [] corr = new float[corrs.length];
        for(int i = 0; i < corrs.length;i++)
        {
            if((corrs[i][1] < pThreshNP && !dPrior[i]) ||(corrs[i][1] < pThresh  && dPrior[i]))
                corr[i] = corrs[i][0];

        }
        return corr;
    }
    private static synchronized float[] shrinkByThreshold(float[][] corrs, float thresh)
       {
           float [] corr = new float[corrs.length];
           for(int i = 0; i < corrs.length;i++)
           {
               if(corrs[i][1] < thresh)
                   corr[i] = corrs[i][0];
           }
           return corr;
       }


    private static synchronized ArrayList<Gene> shrinkByThreshold(ArrayList<Gene> genes, boolean [] iPrior, float threshNP, float thresh)
    {
        ArrayList<Gene> curr = new ArrayList<Gene>();
        for(int i = 0; i < genes.size();i++)
        {
            Gene t = new Gene(i);
            t.symbol = genes.get(i).symbol;
            if((iPrior[i] && genes.get(i).intensityP<thresh) || (!iPrior[i] && genes.get(i).intensityP < threshNP))
            {
                t.intensityValue =genes.get(i).intensityValue;
                t.intensityP = genes.get(i).intensityP;
            }
            else
            {
                t.intensityValue = 0;
            }
            curr.add(t);
        }
        return curr;
    }
    private static synchronized ArrayList<Gene> shrinkByThreshold(ArrayList<Gene> genes, float thresh)
        {
            ArrayList<Gene> curr = new ArrayList<Gene>();
            for(int i = 0; i < genes.size();i++)
            {
                Gene t = new Gene(i);
                t.symbol = genes.get(i).symbol;
                if(genes.get(i).intensityP<thresh)
                {
                    t.intensityValue =genes.get(i).intensityValue;
                    t.intensityP = genes.get(i).intensityP;
                }
                else
                {
                    t.intensityValue = 0;
                }
                curr.add(t);
            }
            return curr;
        }
    public static double getFeatAccuracy(List<Gene> selected, Graph g, String target, String type)
    {
        double acc = 0;
        double tp = 0;
        double fp = 0;
        double fn = 0;
        List<Node> correct = g.getAdjacentNodes(g.getNode(target));
        for(int i = 0; i < selected.size();i++)
        {
            boolean found = false;
            for(Node n: correct)
            {
                if(n.getName().equals(selected.get(i).symbol))
                {
                    acc++;
                    tp++;
                    found = true;
                }

            }
            if(!found)
            {
                fp++;
            }
        }
        if(type.equals("ACC"))
            return acc/correct.size();

        for(int i = 0; i < correct.size();i++)
        {
            Node curr = correct.get(i);
            boolean found = false;

            for(int j = 0; j < selected.size();j++)
            {
                if(selected.get(j).symbol.equals(curr.getName()))
                    found = true;
            }
            if(!found)
                fn++;
        }

        double prec = tp/(tp+fp);
        double rec = tp/(tp+fn);

        return (2*(prec*rec)/(prec+rec));

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