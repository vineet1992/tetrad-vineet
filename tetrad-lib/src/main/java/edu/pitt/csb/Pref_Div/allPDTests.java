package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.pitt.csb.Priors.runPriors;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;
import edu.pitt.dbmi.data.Dataset;

import java.io.*;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

/**
 * The Purpose of this class is to produce the first figure for the experiments section of the Pref-Div Paper
 * How well can we evaluate the prior knowledge we are getting from both intensity and dissimilarity files?
 */
public class allPDTests {



    static int numRuns = 25;
    static int numGenes = 300;

    static boolean boot = false; //Should we use bootstrap samples for PiPrefDiv
    static boolean loocv = false; //Should we use leave-one-out CV for PiPrefDiv
    static boolean useCausalGraph = false; //Are we selecting genes to then use a causal graph to find connections or direct selection?
    static int numSamples = 20; //Number of bootstrap/sub-sampled samples to use
    static int numParams = 30;//Number of parameters to sweep over
    static boolean noiseRandom = false;//Is the reliability range of reliable and unreliable priors set randomly?
    static int numFolds = 5; //Number of Folds for CV to pick alpha


    static int numPriors = 5; //Number of prior knowledge sources
    static int numReliable = 5; //Number of reliable sources
    static int numComponents = 20; //How many components do we have for cluster simulation?
    static int minTargetParents = 10; //How many true parents of the target are there?
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
    static boolean simulateTogether = true; //Is all prior knowledge contained in a single file?

    public static void main(String [] args) {


        double [] ap = new double[]{0.01,0.1,0.5,1.0};
        //double [] ap = new double[]{1.0};
        int [] nr = new int[]{0,1,3,5};
        //int [] nr = new int[]{5};
        int [] ss = new int[]{200};

        for(int ii = 0; ii < ap.length;ii++) {
            amountPrior = ap[ii];
            for (int jj = 0; jj < nr.length; jj++) {
                numReliable = nr[jj];
                for (int kk = 0; kk < ss.length; kk++) {
                    sampleSize = ss[kk];

                    ArrayList<Integer> experiment = new ArrayList<Integer>();
                    experiment.add(1);
                    experiment.add(2);
                    /*experiment.add(3);
                    experiment.add(4);
                    experiment.add(5);*/
                    //0 ->prior evaluation, 1 -> Accuracy of chosen parameters vs Optimal
                    //2->With priors vs no Priors 3-> Comparing different summarization methods for clustered genes
                    //4->Realistic pathway tests where cluster determination is the key
                    //5 ->Summarization with graphical modeling

                    String[] types = new String[]{"FS", "Prediction", "Cluster"}; //Help to create file header for experiment 1
                    String[] starts = new String[]{"Best_Radius", "Best_Accuracy", "Pred_Accuracy"}; //More help

                    String [] startsNP = new String[]{"Pred_Accuracy","Pred_Accuracy_NP"};

                    RunPrefDiv.ClusterType[] allTypes = RunPrefDiv.ClusterType.values();


                    /***Create Results directory***/
                    File rFile = new File("Results");
                    if (!rFile.isDirectory())
                        rFile.mkdir();

                    /***Create special files for this particular experiment***/
                    try {
                        List<PrintStream> out = new ArrayList<PrintStream>();
                        int curr = 0;
                        for(Integer x: experiment) {
                            String start = "Prior Evaluation Experiment/Results/Prior_Evaluation_";
                            if (x== 1) {
                                start = "Parameter Selection Experiment/Results/Parameter_Accuracy_";
                            } else if(x==2)
                            {
                                start = "No Prior Comparison/Results/Comparison_";
                            }else if (x==3) {
                                start = "Summarization Experiment/Results/Summarization_Comparison_";
                            } else if(x==4){
                                start = "Real Clusters/Results/Cluster_Identification_";
                            } else if(x==5){
                                start = "Summarization Causal/Results/Comparison_";
                            } 

                            out.add(new PrintStream(start + numGenes + "_" + sampleSize + "_" + numRuns + "_" + numPriors + "_" + amountPrior + "_" + amountRandom + "_" + boot + "_" + loocv + "_" + evenDistribution + "_" + targetContinuous + "_" + useCausalGraph + "_" + numReliable + "_" + numComponents + "_" + minTargetParents + "_" + numParams + "_" + stabilitySelection + ".txt"));


                            if (x == 0) {
                                out.get(curr).println("Run\tAmount_Prior\tReliable?\tIntensity?\tPredicted_Weight\tActual_Reliability");
                            } else if (x == 1) {
                                out.get(curr).println("Run\tPredicted_Radius_NP\tPredicted_Radius_WP\tPredicted_Radius\tBest_Radius_NP\tBest_Radius_WP\tBest_Radius\tPredicted_Accuracy_WP\tPredicted_Accuracy_NP\tBest_Accuracy_WP\tBest_Accuracy_NP");
                            } else if(x==4)  {

                                out.get(curr).print("Run\tPredicted_Radius_NP\tPredicted_Radius_WP\tPredicted_Threshold_NP\tPredicted_Threshold_WP\tPredicted_Radius\tPredicted_Threshold\t");
                                out.get(curr).println("Clusters_WP\tClusters_NP");
                            }


                            else if (x == 3) {
                                out.get(curr).print("Run\t");
                                for (int i = 0; i < allTypes.length; i++)
                                    out.get(curr).print(allTypes[i] + "\t" + (allTypes[i]+"_NP") + "\t");
                                out.get(curr).println();
                            }
                            else if (x == 5) {
                                out.get(curr).print("Run\t");
                                for (int i = 0; i < allTypes.length; i++)
                                    out.get(curr).print(allTypes[i] + "\t" + (allTypes[i]+"_NP") + "\t");
                                out.get(curr).println();
                            }

                            else if(x==2)
                            {
                                out.get(curr).print("Run\t");
                                for (int i = 0; i < types.length; i++) {
                                    for (int j = 0; j < startsNP.length; j++) {
                                        out.get(curr).print(startsNP[j] + "_" + types[i] + "\t");
                                    }
                                }
                                out.get(curr).println("Radius_NP\tThreshold_NP\tRadius_WP\tThreshold_WP\tRadius\tThreshold");


                            }
                            curr++;
                        }


                        /***Main Loop going through numRuns times***/
                        for (int j = 0; j < numRuns; j++) {
                            HashMap<String, List<Integer>> clusters = new HashMap<String, List<Integer>>();


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
                                clusters = new HashMap<String, List<Integer>>();
                                try {
                                    BufferedReader bTemp = new BufferedReader(new FileReader(clusterFile));
                                    while (bTemp.ready()) {
                                        String[] line = bTemp.readLine().split("\t");
                                        List<Integer> result = new ArrayList<Integer>();
                                        for(int p = 1; p < line.length;p++)
                                            result.add(Integer.parseInt(line[p]));
                                        clusters.put(line[0], result);
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
                            File subsFile = new File("Subsamples/Subsample_Train_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + boot + "_" + numSamples + "_" + loocv + "_" + j + ".txt");

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
                            } else {
                                System.out.println("Couldn't load " + subsFile);
                            }


                            /***LOAD PRIOR KNOWLEDGE FILES***/
                            String priorIntensity = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + "_intensity.txt";
                            String[] dFile = new String[numPriors];
                            for (int k = 0; k < numPriors; k++) {
                                if(simulateTogether)
                                    dFile[k] = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + k + "_" + j + ".txt";
                                else
                                    dFile[k] = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + k + "_" + j + "_dissimilarity.txt";
                            }


                            /***RUN PI PREF DIV***/

                            DataSet toRun;
                            DataSet test = d;
                                /***Load Train-test split***/
                                int trainLength = (int) (d.getNumRows() * 0.9);
                                int[] trainInds = new int[trainLength];
                                for (int i = 0; i < trainLength; i++) {
                                    trainInds[i] = i;
                                }
                                toRun = d.subsetRows(trainInds);
                                test.removeRows(trainInds);

                            PiPrefDiv4 p = new PiPrefDiv4(toRun, "Target", minTargetParents, numParams);
                            p.setSubsamples(subs);
                            p.setUseStabilitySelection(stabilitySelection);
                            p.setParallel(false);
                            p.setPartialCorrs(partialCorr);
                            p.setPdStability(pdStability);
                            p.setVerbose();

                            PiPrefDiv4 noPrior = new PiPrefDiv4(toRun,"Target",minTargetParents,numParams);
                            noPrior.setSubsamples(subs);
                            noPrior.setUseStabilitySelection(stabilitySelection);
                            noPrior.setParallel(false);
                            noPrior.setPartialCorrs(partialCorr);
                            noPrior.setPdStability(pdStability);
                            noPrior.setVerbose();

                            //TODO Switch back and include no priors as an option
                            ArrayList<Gene> selected = p.selectGenes(boot, numSamples, dFile);
                            ArrayList<Gene> npGenes = noPrior.selectGenes(boot,numSamples);

                            //ArrayList<Gene> selected = p.selectGenes(boot,numSamples,useCausalGraph);


                            curr = 0;

                            for(Integer x: experiment) {

                                /***Comparison of selection, cluster, prediction accuracy vs. no prior information at all***/
                                if(x==2)
                                {
                                    System.out.println("Selected Genes: " + selected);
                                    System.out.println("Selected Genes, No Prior " + npGenes);
                                    System.out.println("Truth: " + g.getAdjacentNodes(g.getNode("Target")));
                                    Map<Gene, List<Gene>> lastCluster = p.getLastCluster();
                                    Map<Gene,List<Gene>> npCluster = noPrior.getLastCluster();

                                    System.out.println("Cluster: " + lastCluster);
                                    System.out.println("Cluster: " + npCluster);


                                    /***Do we correctly identify the clusters of variables***/
                                    double clustAccuracy = getClusterAccuracy(lastCluster, clusters, numComponents);
                                    double npClustAccuracy = getClusterAccuracy(npCluster,clusters,numComponents);
                                    double featAccuracy = -1;
                                    double npFeatAccuracy = -1;

                                    /***Do we pick the correct causes of the target variable***/
                                    if (minTargetParents == g.getAdjacentNodes(g.getNode("Target")).size()) {
                                        //Correct amount of selections so use accuracy
                                        featAccuracy = getFeatAccuracy(selected, g, "Target", "ACC");
                                        npFeatAccuracy = getFeatAccuracy(npGenes,g,"Target","ACC");

                                    } else {
                                        //Incorrect amount of selected variables so use F1
                                        featAccuracy = getFeatAccuracy(selected, g, "Target", "F1");
                                        npFeatAccuracy = getFeatAccuracy(npGenes,g,"Target","F1");
                                    }

                                    /***Are all the correct clusters represented in the selected features***/
                                    double repAccuracy = getRepresentAccuracy(selected,g,clusters,"Target");
                                    double npRepAccuracy = getRepresentAccuracy(npGenes,g,clusters,"Target");

                                    /**** Can we accurately predict the target variable on the test set***/
                                    double predAccuracy = getAccuracy(toRun, test, selected);
                                    double npPredAccuracy = getAccuracy(toRun,test,npGenes);


                                    System.out.println("Run #" +j + " With " + amountPrior + " priors, " + numReliable + " reliable, and " + sampleSize + " samples \n : Prediction=" + predAccuracy + ", Features=" + featAccuracy + ", Cluster=" + clustAccuracy + ", Representation=" + repAccuracy);
                                    System.out.println("No Priors: Prediction=" + npPredAccuracy + ", Features=" + npFeatAccuracy + ", Cluster=" + npClustAccuracy + ", Representation=" + npRepAccuracy);

                                    out.get(curr).print(j + "\t" +repAccuracy + "\t" + npRepAccuracy + "\t" +  predAccuracy + "\t" + npPredAccuracy + "\t" + clustAccuracy + "\t" + npClustAccuracy + "\t");
                                    out.get(curr).println(p.getLastRadius()[0] + "\t" + p.getLastThreshold()[0] + "\t" + p.getLastRadiusWP() + "\t" + p.getLastThresholdWP() + "\t" + noPrior.getLastRadius()[0] + "\t" + noPrior.getLastThreshold()[0]);


                                }
                                if (x == 0) //Evaluate prior knowledge
                                {
                                    if(!simulateTogether)
                                    {
                                        File ri = new File("Reliabilities_I_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");
                                        double[] reliInt = new double[numPriors];
                                        BufferedReader b2 = new BufferedReader(new FileReader(ri));
                                        for (int i = 0; i < numPriors; i++) {
                                            reliInt[i] = Double.parseDouble(b2.readLine());
                                        }
                                        b2.close();
                                        double []weights = p.getLastIntensityWeights();

                                        for (int i = 0; i < numPriors; i++) {
                                            out.get(curr).println(j + "\t" + getAmountPrior(priorIntensity, i) / ((double) g.getNumNodes() - 1) + "\t" + (i < numReliable) + "\tT\t" + weights[i] + "\t" + reliInt[i]);
                                        }
                                    }

                                    /***LOAD RELIABILITY SCORES***/
                                    File rd = new File("Reliabilities_D_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");
                                    if(simulateTogether)
                                        rd = new File("Reliabilities_D_All_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");

                                    double[] reliDis = new double[numPriors];
                                    BufferedReader b = new BufferedReader(new FileReader(rd));

                                    for (int i = 0; i < numPriors; i++) {
                                        reliDis[i] = Double.parseDouble(b.readLine());
                                    }
                                    b.close();
                                    double [] weights = p.getLastSimilarityWeights();


                                    /***PRINT RESULTS TO FILE***/

                                    for (int i = 0; i < numPriors; i++) {
                                        out.get(curr).println(j + "\t" + (double) (getAmountPrior(dFile[i])) / (g.getNumNodes() * (g.getNumNodes() - 1) / 2) + "\t" + (i < numReliable) + "\tF\t" + weights[i] + "\t" + reliDis[i]);
                                    }
                                }
                                /***Evaluate Parameter Accuracy, by sweeping over potential parameter choices***/
                                else if (x == 1) //Evaluate Parameter accuracy
                                {
                                    System.out.println("Selected Genes: " + selected);
                                    System.out.println(g.getAdjacentNodes(g.getNode("Target")));
                                    Map<Gene, List<Gene>> lastCluster = p.getLastCluster();
                                    System.out.println("Cluster: " + lastCluster);


                                    /***Do we correctly identify the clusters of variables***/
                                    double clustAccuracy = getClusterAccuracy(lastCluster, clusters, numComponents);
                                    double featAccuracy = -1;

                                    /***Do we pick the correct causes of the target variable***/
                                    if (minTargetParents == g.getAdjacentNodes(g.getNode("Target")).size()) {
                                        //Correct amount of selections so use accuracy
                                        featAccuracy = getFeatAccuracy(selected, g, "Target", "ACC");

                                    } else {
                                        //Incorrect amount of selected variables so use F1
                                        featAccuracy = getFeatAccuracy(selected, g, "Target", "F1");
                                    }

                                    /***Are all the correct clusters represented in the selected features***/
                                    double repAccuracy = getRepresentAccuracy(selected,g,clusters,"Target");

                                    /**** Can we accurately predict the target variable on the test set***/
                                    double predAccuracy = getAccuracy(toRun, test, selected);
                                    System.out.println("Run #" + j + ": Prediction=" + predAccuracy + ", Features=" + featAccuracy + ", Cluster=" + clustAccuracy + ", Representation=" + repAccuracy);

                                    double[] radii = p.getTestedRadii();

                                    ArrayList<Gene> meanGenes = PiPrefDiv4.createGenes(toRun, "Target",false);
                                    meanGenes = Functions.computeAllIntensities(meanGenes,1,toRun,"Target",false,false,false,1);

                                    //0 is correlation, 1 is p-value
                                    float[] meanDis = Functions.computeAllCorrelations(meanGenes, toRun, false, false, false, 1);
                                    //Use the weights from PiPrefDiv, we are interested in determining how good param selection is given that weights are good
                                    //Because experiment 1 demonstrated that weights are good


                                    double[] results;
                                    double[]npParam;
                                    //for every radii and threshold combination, give the three accuracies

                                    if (!parallel) {
                                        npParam = getBestNPParams(radii, numComponents, meanDis, meanGenes, p.getDissimilarityPriors(), clusters, g, toRun, test); //Feature Select, Cluster Acc, Prediction Acc
                                        results = getBestParams(radii, radii, numComponents, meanDis, meanGenes, p.getIntensityPriors(), p.getDissimilarityPriors(), clusters, g, toRun, test);
                                    } else {
                                        npParam = getBestNPParamsPar(radii, numComponents, meanDis, meanGenes, p.getDissimilarityPriors(), clusters, g, toRun, test);
                                        results = getBestParamsPar(radii, radii, numComponents, meanDis, meanGenes, p.getIntensityPriors(), p.getDissimilarityPriors(), clusters, g, toRun, test);
                                    }

                                    out.get(curr).println(j + "\t" + p.getLastRadius()[0] + "\t" + p.getLastRadius()[1] + "\t" + noPrior.getLastRadius()[0] + "\t" +
                                    results[0] + "\t" + results[1] + "\t" + npParam[0] + "\t"  + clustAccuracy + "\t" +  getClusterAccuracy(noPrior.getLastCluster(), clusters, numComponents)  + "\t" + results[2] + "\t" + npParam[1]);


                                    //Roll through all tested radii,threshold combos and RunPrefDiv directly to get accuracies
                                }

                                /***Real pathway cluster identification**/ 
                                else if (x==4)
                                {
                                    System.out.println("Selected Genes: " + selected);
                                    System.out.println("Selected Genes, No Prior " + npGenes);
                                    System.out.println("Truth: " + g.getAdjacentNodes(g.getNode("Target")));
                                    Map<Gene, List<Gene>> lastCluster = p.getLastCluster();
                                    Map<Gene,List<Gene>> npCluster = noPrior.getLastCluster();

                                    System.out.println("Cluster: " + lastCluster);
                                    System.out.println("Cluster: " + npCluster);


                                    /***Do we correctly identify the clusters of variables***/
                                    double clustAccuracy = getClusterAccuracy(lastCluster, clusters, numComponents);
                                    double npClustAccuracy = getClusterAccuracy(npCluster,clusters,numComponents);


                                    out.get(curr).print(j + "\t" + p.getLastRadius()[0] + "\t" + p.getLastThreshold()[0] + "\t" + p.getLastRadiusWP() + "\t" + p.getLastThresholdWP() + "\t" + noPrior.getLastRadius()[0] + "\t" + noPrior.getLastThreshold()[0]);
                                    out.get(curr).println("\t" + clustAccuracy + "\t" + npClustAccuracy);


                                }
                                /***Summarization Method with Causal modeling***/
                                else if(x==5)
                                {
                                    out.get(curr).print(j + "\t");
                                    Map<Gene, List<Gene>> lastCluster = p.getLastCluster();

                                    Map<Gene,List<Gene>> npCluster = noPrior.getLastCluster();
                                    for (int i = 0; i < allTypes.length; i++) {
                                        System.out.print("Running " + allTypes[i] + " with priors");

                                        DataSet summarized = RunPrefDiv.summarizeData(toRun, selected, lastCluster, allTypes[i]);
                                        DataSet causes = getCausalFeatures(summarized,"Target");

                                        DataSet summarizedTest = RunPrefDiv.summarizeData(test, selected, lastCluster, allTypes[i]);
                                        DataSet causesTest = subsetByCols(causes,summarizedTest);
                                        double predAccuracy = -1;

                                        try {
                                            predAccuracy = getAccuracy(causes, causesTest, null);
                                        }catch(Exception e)
                                        {
                                            predAccuracy = Double.NaN;
                                            System.out.println(causes);
                                        }

                                        out.get(curr).print(predAccuracy + "\t");
                                        System.out.println("Done");

                                        System.out.print("Running " + allTypes[i] + " for no Priors...");
                                        summarized = RunPrefDiv.summarizeData(toRun, npGenes, npCluster, allTypes[i]);
                                        causes = getCausalFeatures(summarized,"Target");

                                        summarizedTest = RunPrefDiv.summarizeData(test, npGenes, npCluster, allTypes[i]);
                                        causesTest = subsetByCols(causes,summarizedTest);

                                        try {
                                            predAccuracy = getAccuracy(causes, causesTest, null);
                                        }catch(Exception e)
                                        {
                                            predAccuracy = Double.NaN;
                                            System.out.println(causes);
                                        }
                                        System.out.println("Done");

                                        out.get(curr).print(predAccuracy + "\t");


                                    }
                                    out.get(curr).println();
                                }
                                /***Summarization Method Tests***/
                                else if (x == 3) {
                                    out.get(curr).print(j + "\t");
                                    Map<Gene, List<Gene>> lastCluster = p.getLastCluster();
                                    Map<Gene,List<Gene>> npCluster = noPrior.getLastCluster();

                                    for (int i = 0; i < allTypes.length; i++) {

                                        DataSet summarized = RunPrefDiv.summarizeData(toRun, selected, lastCluster, allTypes[i]);
                                        DataSet summarizedTest = RunPrefDiv.summarizeData(test, selected, lastCluster, allTypes[i]);
                                        double predAccuracy = -1;
                                        try {
                                            predAccuracy = getAccuracy(summarized, summarizedTest, null);
                                        }
                                        catch(Exception e)
                                        {
                                            System.out.println(summarized);
                                            predAccuracy = Double.NaN;
                                        }

                                        out.get(curr).print(predAccuracy + "\t");

                                        summarized = RunPrefDiv.summarizeData(toRun, npGenes, npCluster, allTypes[i]);
                                        summarizedTest = RunPrefDiv.summarizeData(test, npGenes, npCluster, allTypes[i]);
                                        try {
                                            predAccuracy = getAccuracy(summarized, summarizedTest, null);
                                        }
                                        catch(Exception e)
                                        {
                                            predAccuracy = Double.NaN;
                                            System.out.println(summarized);

                                        }

                                        out.get(curr).print(predAccuracy + "\t");


                                    }
                                    out.get(curr).println();
                                }

                                out.get(curr).flush();
                                curr++;
                            }
                        }
                        for(PrintStream temp:out)
                        {
                            temp.flush();
                            temp.close();
                        }

                    } catch (Exception e) {
                        System.err.println("Error with printstream");
                        e.printStackTrace();
                        System.exit(-1);
                    }
                }
            }
        }


    }


    public static DataSet subsetByCols(DataSet train, DataSet test)
    {
        List<Node> vars = new ArrayList<Node>();
        for(String s: train.getVariableNames())
        {
            if(test.getVariable(s)==null)
            {
                System.err.println(train + "\n and test:\n" + test);
                System.err.println(s);
                System.exit(-1);
            }
            vars.add(test.getVariable(s));
        }
        return test.subsetColumns(vars);
    }

    /***
     *
     * @param results 2-D array of numRadii x 3, Col 0 = FS, Col 1 = Prediction RMSE, Col 2 = Cluster Acc
     * @param out PrintStream to write results to
     * @param radii List of tested radii
     * @param predAcc The accuracy that our method got (with prior information)
     * @param radiusNP Best radius for the entire dataset as one
     */
    public static void printResults(double[][]results, PrintStream out, double [] radii, double [] predAcc, double radiusNP)
    {
        out.print(radiusNP +  "\t");

        /**Feature Selection Results Printing**/
        double[] bestRadii = new double[3];
        double []maxValue = new double[3];
        maxValue[1] = Double.MAX_VALUE;
        for(int i = 0; i < results.length;i++)
        {
                for(int k = 0; k < 3;k++)
                {
                    if((k!=1 && results[i][k]>maxValue[k]) || (k==1 && results[i][k] < maxValue[k]))
                    {
                        maxValue[k] = results[i][k];
                        bestRadii[k] = radii[i];
                    }
                }
        }
        for(int k = 0; k < 3;k++)
        {
            out.print(radii[k] + "\t" + maxValue[k] + "\t" + predAcc[k] + "\t");
        }
        out.println(predAcc[3]);/**To include cluster representation score***/
    }


    /*** Run graphical modeling algorithm and get the features related to the target**/
    public static DataSet getCausalFeatures(DataSet input,String target)
    {
        try{
            runPriors.addDummy(input);
            int numLambda = 40;
            double low = 0.05;
            double high = 0.9;
            double [] lambdas = new double[numLambda];
            for(int i = 0; i < numLambda;i++)
            {
                lambdas[i] = low + i*(high-low)/numLambda;
            }

            STEPS s= new STEPS(input,lambdas,0.05,numSamples);
            Graph graph = s.runStepsPar();
            /***Columns to keep**/
            List<Node> cols = new ArrayList<Node>();
            cols.add(input.getVariable(target));

            for(Node n: graph.getAdjacentNodes(graph.getNode(target)))
            {
                if(!n.getName().equals("Dummy"))
                    cols.add(n);
            }

            /***If no variables are selected, use all except the dummy**/
            if(cols.size()==1)
            {
                for(Node r: input.getVariables())
                {
                    if(!r.getName().equals("Dummy"))
                        cols.add(r);
                }
            }

            return input.subsetColumns(cols);

        }catch(Exception e)
        {
            e.printStackTrace();
            return null;
        }
    }

    /****
     *
     * @param radiiNP List of tested radii for No Priors (NP)
     * @param radii List of tested radii for With Priors (WP)
     * @param numComponents Number of clusters in the data generating graph
     * @param corrs Float [] correlation between each pair of variables
     * @param meanGenes List of genes
     * @param iPriors Which relationships between genes and targets do we have prior information for?
     * @param dPriors Which relationships do we have prior information for?
     * @param clusters True clusters in the data generating graph
     * @param g Data generating graph
     * @param toRun Training Data
     * @param test Testing Data
     * @return 1-D Array with two entries: Best Radii NP, Best Radii -> Both are with respect to cluster accuracy
     */
    public static double[] getBestParamsPar(final double [] radiiNP, final double[] radii, final int numComponents,final float[] corrs, final ArrayList<Gene> meanGenes, final boolean[] iPriors, final boolean[] dPriors, final Map<String,List<Integer>>clusters, final Graph g, final DataSet toRun,final DataSet test)
    {
        final double[][] result = new double[radii.length][radiiNP.length];
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
                result[i][j] = getClusterAccuracy(lastCluster,clusters,numComponents);

            }

            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        //RadiiWP is the columns, RadiiNP is the rows
                        int j = s % radii.length; //Column number (radiiWP)
                        int i = s / radii.length; //Row Number (radiiNP)
                        try {
                            ArrayList<Gene> currGenes = shrinkByThreshold(meanGenes,iPriors,(float)radiiNP[i],(float)radii[j]);
                            Collections.sort(currGenes,Gene.IntensityComparator);
                            PrefDiv pd = new PrefDiv(currGenes,minTargetParents,0,radiiNP[i],(float)radii[j],dPriors,corrs);
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
        StabilityAction sa = new StabilityAction(chunk,0, radii.length * radiiNP.length);
        pool.invoke(sa);

        double [] res = new double[3];
        for(int i = 0; i < result.length;i++)
        {
            for(int j = 0; j < result[i].length;j++)
            {
                if(result[i][j] > res[2])
                {
                    res[2] = result[i][j];
                    res[0] = radiiNP[i];
                    res[1] = radii[j];
                }
            }
        }
        return res;
    }


    /****
     *
     * @param radii Tested radii values
     * @param numComponents Number of components in the simulated graph
     * @param corrs Correlation between every pair of variables
     * @param meanGenes List of genes
     * @param dPriors Whether or not each relationship had prior information
     * @param clusters The true clusters in the data generating graph
     * @param g The data generating graph
     * @param toRun Training dataset
     * @param test Testing dataset
     * @return The radii corresponding to the best cluster accuracy
     */
    private static double [] getBestNPParamsPar(final double [] radii,final int numComponents,final float[]corrs, final ArrayList<Gene> meanGenes, final boolean [] dPriors, final  Map<String,List<Integer>>clusters, final Graph g, final DataSet toRun,final DataSet test)
    {
        final double[][] result = new double[radii.length][3];
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


            private synchronized void addAccuracy(List<Gene>selected,Map<Gene,List<Gene>> lastCluster,int i)
            {
                result[i][2] = getClusterAccuracy(lastCluster,clusters,numComponents);

                String type = "F1";
                if(selected.size()==g.getAdjacentNodes(g.getNode("Target")).size())
                    type = "ACC";
                result[i][0] = getFeatAccuracy(selected,g,"Target",type);
                result[i][1] = getAccuracy(toRun,test,selected);
            }

            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        try {

                            ArrayList<Gene> currGenes = shrinkByThreshold(meanGenes,(float)radii[s]);
                            Collections.sort(currGenes,Gene.IntensityComparator);


                            PrefDiv pd = new PrefDiv(currGenes,minTargetParents,0,radii[s],radii[s],dPriors,corrs);
                            pd.setCluster(true);
                            List<Gene> selected = pd.diverset();
                            Map<Gene,List<Gene>> lastCluster = pd.clusters;
                            addAccuracy(selected,lastCluster,s);


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
        StabilityAction sa = new StabilityAction(chunk,0, radii.length);
        pool.invoke(sa);



        return bestAcc(result,radii);
    }


    /***
     *
     * @param result 2-D Array of size numRadii x 3, Column 0 is Feature selection accuracy, Column 1 is prediction accuracy, Column 2 is cluster accuracy
     * @param radii List of tested radii
     * @return The radii that gives the best cluster accuracy and the accuracy itself
     */
    private static double [] bestAcc(double[][]result, double [] radii)
    {
        double bestAcc = -1;
        double bestRadii = -1;
        for(int i = 0; i < radii.length;i++)
        {
            if(result[i][2] > bestAcc) {
                bestAcc = result[i][2];
                bestRadii = radii[i];
            }

        }
        return new double[] {bestRadii,bestAcc};
    }

    /***
     *
     *  Serial version of the best NP Params function above
     */
    private static double[] getBestNPParams(double [] radii, int numComponents,float[]corrs, ArrayList<Gene> genes, boolean[]dPrior, Map<String,List<Integer>> clusters,Graph g, DataSet toRun,DataSet test)
    {
        double[][] result = new double[radii.length][3];
                for (int i = 0; i < radii.length; i++)
                {
                    ArrayList<Gene> currGenes = shrinkByThreshold(genes,(float)radii[i]);
                    Collections.shuffle(currGenes);
                        Collections.sort(currGenes,Gene.IntensityComparator);
                        PrefDiv pd = new PrefDiv(currGenes,minTargetParents,0,radii[i],radii[i],dPrior,corrs);
                        pd.setCluster(true);
                        List<Gene> selected = pd.diverset();
                        Map<Gene, List<Gene>> lastCluster = pd.clusters;
                        result[i][2] = getClusterAccuracy(lastCluster, clusters, numComponents);

                        String type = "F1";
                        if (selected.size() == g.getAdjacentNodes(g.getNode("Target")).size())
                            type = "ACC";
                        result[i][0] = getFeatAccuracy(selected, g, "Target", type);
                        result[i][1] = getAccuracy(toRun, test, selected);

              }
        return bestAcc(result,radii);
    }

    /*** Serial implementation of the above parallel method***/
    private static double[] getBestParams(double [] radiiNP, double [] radii, int numComponents,float[]corrs, ArrayList<Gene> genes, boolean [] iPrior, boolean[]dPrior, Map<String,List<Integer>> clusters,Graph g, DataSet toRun,DataSet test)
    {
        double[] result = new double[3];

        //Find optimal parameters for no priors, and then loop again for with priors
        for (int i = 0; i < radii.length; i++) {
            for (int j = 0; j < radiiNP.length; j++) {
                ArrayList<Gene> currGenes = shrinkByThreshold(genes, iPrior, (float) radiiNP[j], (float) radii[i]);
                Collections.shuffle(currGenes);
                Collections.sort(currGenes, Gene.IntensityComparator);
                PrefDiv pd = new PrefDiv(currGenes, minTargetParents, 0, (float) radiiNP[j], (float) radii[i], dPrior, corrs);
                pd.setCluster(true);

                List<Gene> selected = pd.diverset();
                Map<Gene, List<Gene>> lastCluster = pd.clusters;
                if(getClusterAccuracy(lastCluster,clusters,numComponents)> result[2])
                {
                     result[2] = getClusterAccuracy(lastCluster, clusters, numComponents);
                    result[0] = radii[i];
                    result[1] = radiiNP[j];
                }
             }
        }
        return result;
    }


    //TODO Incorporate discrete target variables
    public static double getAccuracy(DataSet train, DataSet test, List<Gene> selected)
    {
        test = DataUtils.standardizeData(test);
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

    /***
     *
     * @param est Estimated clusters (lists of genes)
     * @param actual Actual cluster membership map
     * @param numClusters Total number of clusters
     * @return cluster accuracy metric via Hungarian algorithm
     */
    public static double getClusterAccuracy(Map<Gene,List<Gene>> est, Map<String,List<Integer>> actual, int numClusters)
    {

        /***Create list of genes representing each actual pathway***/
        List<List<Gene>> temp = new ArrayList<List<Gene>>();
        for(int i = 0; i < numClusters;i++)
            temp.add(new ArrayList<Gene>());
        int x = 0;
        for(String s: actual.keySet())
        {
            Gene g = new Gene(x);
            g.symbol = s;
            List<Integer> currPaths = actual.get(s);
            for(int i = 0; i < currPaths.size();i++)
            {
                temp.get(currPaths.get(i)).add(g);
            }
            x++;
        }

        /***Create list of genes representing estimated pathways***/
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


    /***
     *
     * @param selected List of genes that were  by the feature selection method
     * @param g  True data generating graph
     * @param clusts True cluster assignments for each gene
     * @param targName The target variable's name as a String
     * @return What percent of the true causal clusters are represented by the result
     * NOTE THIS IS NOT REPRESENTATIVE OF GENES IN MULTIPLE CLUSTERS (ARBITRARILY CHOOSES THE FIRST ONE)
     */
    public static double getRepresentAccuracy(List<Gene> selected, Graph g, Map <String,List<Integer>> clusts, String targName)
    {

        int trueLength = 0; /***How many clusts are represented in the truth?***/
        ArrayList<Integer> trueClusts = new ArrayList<Integer>();


        /***Get true cluster assignments***/
        Node target = g.getNode(targName);
        for(Node n: g.getAdjacentNodes(target))
        {
            trueClusts.add(clusts.get(n.getName()).get(0));
        }

        trueLength = trueClusts.size();
        ArrayList<Integer> estClusts = new ArrayList<Integer>();

        /***Get estimated cluster assignments***/
        for(Gene n: selected)
        {
            estClusts.add(clusts.get(n.symbol).get(0));
        }

        /***Compute how many clusters in the truth are represented in the est***/

        int count = 0;
        for(int i = 0; i < estClusts.size();i++)
        {
            if(trueClusts.contains(estClusts.get(i))) {
                trueClusts.remove(estClusts.get(i));
                count++;
            }
        }
        return count/(double)trueLength;

    }


    /***
     *
     * @param genes List of genes with intensity values filled
     * @param iPrior Boolean [] specifying if we had prior information for gene i
     * @param threshNP Threshold for genes that we did not have prior information for
     * @param threshWP Threshold for genes that we had prior information for
     * @return A new list of genes with intensities under the threshold shrunk to zero
     */
    private static synchronized ArrayList<Gene> shrinkByThreshold(ArrayList<Gene> genes, boolean [] iPrior, float threshNP, float threshWP)
    {
        ArrayList<Gene> curr = new ArrayList<Gene>();
        for(int i = 0; i < genes.size();i++)
        {
            Gene t = new Gene(i);
            t.symbol = genes.get(i).symbol;
            if((iPrior[i] && genes.get(i).intensityValue>=threshWP) || (!iPrior[i] && genes.get(i).intensityValue >= threshNP))
            {
                t.intensityValue =genes.get(i).intensityValue;
            }
            else
            {
                t.intensityValue = 0;
            }
            curr.add(t);
        }
        return curr;
    }

    /***
     *
     * @param genes List of genes with correlations to the target variable stored in the intensityValue
     * @param thresh An absolute threshold to shrink correlations to 0 (Must be greater than the threshold to be non-zero)
     * @return A List of genes with shrunk correlations
     */
    private static synchronized ArrayList<Gene> shrinkByThreshold(ArrayList<Gene> genes, float thresh)
        {
            ArrayList<Gene> curr = new ArrayList<Gene>();
            for(int i = 0; i < genes.size();i++)
            {
                Gene t = new Gene(i);
                t.symbol = genes.get(i).symbol;
                if(genes.get(i).intensityValue>=(1-thresh))
                {
                    t.intensityValue =genes.get(i).intensityValue;
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
