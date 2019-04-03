package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.Pref_Div.Comparisons.ClusterSim;
import edu.pitt.csb.Pref_Div.Comparisons.PrefDivComparator;
import edu.pitt.csb.Priors.runPriors;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;

import java.io.*;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

/**
 * The Purpose of this class is to produce the first figure for the experiments section of the Pref-Div Paper
 * How well can we evaluate the prior knowledge we are getting from both intensity and dissimilarity files?
 */
public class allPDTests {


    //TODO Deal with Causal experiments



    static int numRuns = 10;
    static int numGenes = 3000;

    static boolean boot = false; //Should we use bootstrap samples for PiPrefDiv
    static boolean loocv = false; //Should we use leave-one-out CV for PiPrefDiv
    static boolean useCausalGraph = false; //Are we selecting genes to then use a causal graph to find connections or direct selection?
    static int numSamples = 20; //Number of bootstrap/sub-sampled samples to use
    static int numParams = 30;//Number of parameters to sweep over
    static boolean noiseRandom = true;//Is the reliability range of reliable and unreliable priors set randomly?
    static int numFolds = 5; //Number of Folds for CV to pick alpha


    static int numPriors = 5; //Number of prior knowledge sources
    static int numReliable = 5; //Number of reliable sources
    static int numComponents = 300; //How many components do we have for cluster simulation?
    static int minTargetParents = 75; //How many true parents of the target are there?
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

    /***For experiment 7 only***/
    static double wpHigh = 0.95;
    static double wpLow = 0.5;
    static int numWp = 15;

    /***For experiment 6 only***/
    static int numFeatsHigh = 75;
    static int numFeatsLow = 5;
    static int numFeatsTotal = 10;

    public static void main(String [] args) {


        double [] ap = new double[]{0.5,0.75};
        //double [] ap = new double[]{1.0,0.5};
        int [] nr = new int[]{1,3,5,0};
        //int [] nr = new int[]{5};
        int [] ss = new int[]{100,200};

        double [] wpCutoff = new double[numWp];
        int [] numFeats = new int[numFeatsTotal];

        for(int i = 0; i < wpCutoff.length;i++)
        {
            wpCutoff[i] = wpLow + (wpHigh-wpLow)*i/wpCutoff.length;
        }

        for(int i = 0; i < numFeats.length;i++)
        {
            numFeats[i] = numFeatsLow + (numFeatsHigh - numFeatsLow)*i/numFeatsTotal;
        }


        for(int ii = 0; ii < ap.length;ii++) {
            amountPrior = ap[ii];
            for (int jj = 0; jj < nr.length; jj++) {
                numReliable = nr[jj];
                for (int kk = 0; kk < ss.length; kk++) {
                    sampleSize = ss[kk];

                    ArrayList<Integer> experiment = new ArrayList<Integer>();
                    experiment.add(0);
                    //experiment.add(1);
                    experiment.add(2);
                    //experiment.add(3);
                    //experiment.add(4);
                    experiment.add(5);
                    //experiment.add(6);
                    //experiment.add(7);
                    //0 ->prior evaluation, 1 -> Accuracy of chosen parameters vs Optimal
                    //2->With priors vs no Priors 3-> Comparing different summarization methods for clustered genes
                    //4->Realistic pathway tests where cluster determination is the key
                    //5 ->Summarization with graphical modeling

                    //6 -> # of selected clusters sensitivity analysis
                    //7 -> wpCutoff sensitivity analysis

                    String[] types = new String[]{"FS", "Prediction", "Cluster","Friendly_Cluster"}; //Help to create file header for experiment 1

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
                            else if(x==6)
                            {
                                start = "Number Features/Results/Sensitivity_";
                            }
                            else if(x==7)
                            {
                                start = "Cutoff Comparison/Results/Sensitivity_";
                            }

                            out.add(new PrintStream(start + numGenes + "_" + sampleSize + "_" + numRuns + "_" + numPriors + "_" + amountPrior + "_" + amountRandom + "_" + boot + "_" + loocv + "_" + evenDistribution + "_" + targetContinuous + "_" + useCausalGraph + "_" + numReliable + "_" + numComponents + "_" + minTargetParents + "_" + numParams + "_" + stabilitySelection + ".txt"));


                            if (x == 0) {
                                out.get(curr).println("Run\tAmount_Prior\tReliable?\tIntensity?\tPredicted_Weight\tActual_Reliability");
                            } else if (x == 1) {
                                out.get(curr).println("Run\tPredicted_Radius_NP\tPredicted_Radius_WP\tPredicted_Radius\tBest_Radius_NP\tBest_Radius_WP\tBest_Radius\tPredicted_Accuracy_WP\tPredicted_Accuracy_NP\tBest_Accuracy_WP\tBest_Accuracy_NP");
                            } else if(x==4)  {

                                out.get(curr).print("Run\tPredicted_Radius_NP\tPredicted_Radius_WP\tPredicted_Radius\t");
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
                                out.get(curr).print("Perfect_Acc");
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
                                out.get(curr).println("Radius_NP\tRadius_WP\tRadius");


                            }
                            else if(x==6)
                            {
                                //Run, Are Priors Included?, # of clusters selected, % of all clusters picked up, % of true clusters picked up, friendly cluster accuracy
                                out.get(curr).println("Run\tWith_Prior\tParam\tAll_Clusters\tTrue_Clusters\tFriendly_Cluster");
                            }
                            else if(x==7)
                            {
                                //Run, Are priors included?, Parameter, Representation Acc, Cluster Acc, Friendly Cluster Acc, Prediction Acc
                                out.get(curr).println("Run\tWith_Prior\tParam\tRep_Accuracy\tCluster_Accuracy\tFriendly_Cluster_Accuracy\tPrediction_Accuracy");
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
                            if (graphFile.exists()) {
                                g = GraphUtils.loadGraphTxt(graphFile);


                                /***Remove any nodes in the graph that are latents***/
                                List<Node> toRemove = new ArrayList<Node>();
                                for (Node n : g.getNodes()) {
                                    if (n.getName().startsWith("L"))
                                        toRemove.add(n);
                                }
                                g.removeNodes(toRemove);



                            }
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

                                /***Must be standardized separately using same parameters***/
                                standardizeData(toRun,test);


                            PiPrefDiv4 p = new PiPrefDiv4(toRun, "Target", minTargetParents, numParams);
                            p.setSubsamples(subs);
                            p.setParallel(false);
                            p.setPartialCorrs(partialCorr);
                            //p.setVerbose();

                            PiPrefDiv4 noPrior = new PiPrefDiv4(toRun,"Target",minTargetParents,numParams);
                            noPrior.setSubsamples(subs);
                            noPrior.setParallel(false);
                            noPrior.setPartialCorrs(partialCorr);
                           // noPrior.setVerbose();

                            //TODO Switch back and include no priors as an option
                            ArrayList<Gene> selected = p.selectGenes(boot, numSamples, dFile);

                            //selected = p.selectGenes(0,1,dFile);

                            ArrayList<Gene> npGenes = noPrior.selectGenes(boot,numSamples);

                            //ArrayList<Gene> selected = p.selectGenes(boot,numSamples,useCausalGraph);
                            PrefDivComparator compare = new PrefDivComparator(p,"Target",g,clusters);
                            PrefDivComparator npCompare = new PrefDivComparator(noPrior,"Target",g,clusters);
                            double clustAccuracy = compare.getClusterAccuracy(numComponents);
                            double npClustAccuracy = npCompare.getClusterAccuracy(numComponents);
                            double featAccuracy = -1;
                            double npFeatAccuracy = -1;


                            /***Are all the correct clusters represented in the selected features***/
                            double repAccuracy = compare.getRepresentAccuracy();
                            double npRepAccuracy = npCompare.getRepresentAccuracy();

                            /**** Can we accurately predict the target variable on the test set***/
                            double predAccuracy = compare.getPredictionAccuracy(toRun, test);
                            double npPredAccuracy = npCompare.getPredictionAccuracy(toRun,test);


                            /***Do we pick the correct causes of the target variable***/
                            if (minTargetParents == g.getAdjacentNodes(g.getNode("Target")).size()) {
                                //Correct amount of selections so use accuracy
                                featAccuracy = compare.getFeatAccuracy( "ACC");
                                npFeatAccuracy = npCompare.getFeatAccuracy("ACC");

                            } else {
                                //Incorrect amount of selected variables so use F1
                                featAccuracy = compare.getFeatAccuracy("F1");
                                npFeatAccuracy = npCompare.getFeatAccuracy("F1");
                            }

                            double friendlyClust = compare.getFriendlyCluster(numComponents);
                            double npFriendly = npCompare.getFriendlyCluster(numComponents);
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
                                    System.out.println("Cluster, No Prior: " + npCluster);



                                    /***Do we correctly identify the clusters of variables***/





                                    System.out.println("Run #" +j + " With " + amountPrior + " priors, " + numReliable + " reliable, and " + sampleSize + " samples \n : Prediction=" + predAccuracy + ", Features=" + featAccuracy + ", Cluster=" + clustAccuracy + ", Representation=" + repAccuracy);
                                    System.out.println("No Priors: Prediction=" + npPredAccuracy + ", Features=" + npFeatAccuracy + ", Cluster=" + npClustAccuracy + ", Representation=" + npRepAccuracy);

                                    out.get(curr).print(j + "\t" +repAccuracy + "\t" + npRepAccuracy + "\t" +  predAccuracy + "\t" + npPredAccuracy + "\t" + clustAccuracy + "\t" + npClustAccuracy + "\t" + friendlyClust + "\t" + npFriendly + "\t");
                                    out.get(curr).println(p.getLastRadius()[0] + "\t" + p.getLastRadiusWP() + "\t" + noPrior.getLastRadius()[0]);


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
                                    System.out.println("True clusters: " + clusters);



                                    System.out.println("Run #" + j + ": Prediction=" + predAccuracy + ", Features=" + featAccuracy + ", Cluster=" + clustAccuracy + ", Representation=" + repAccuracy);

                                    double[] radii = p.getTestedRadii();


                                    double[] results;
                                    double[]npParam;
                                    //for every radii and threshold combination, give the three accuracies

                                    if (!parallel) {
                                        npParam = getBestNPParams(radii, numComponents, noPrior, clusters, g, toRun, test); //Feature Select, Cluster Acc, Prediction Acc
                                        results = getBestParams(radii, radii, numComponents, p, dFile,clusters, g, toRun, test);
                                    } else {
                                        npParam = getBestNPParamsPar(radii, numComponents, noPrior, clusters, g, toRun, test);
                                        results = getBestParamsPar(radii, radii, numComponents, p,dFile, clusters, g, toRun, test);
                                    }

                                    out.get(curr).println(j + "\t" + p.getLastRadius()[0] + "\t" + p.getLastRadius()[1] + "\t" + noPrior.getLastRadius()[0] + "\t" +
                                            results[0] + "\t" + results[1] + "\t" + npParam[0] + "\t"  + friendlyClust + "\t" +  npFriendly  + "\t" +
                                            results[2] + "\t" + npParam[1]);


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


                                    out.get(curr).print(j + "\t" + p.getLastRadius()[0] +  "\t" + p.getLastRadiusWP() +  "\t" + noPrior.getLastRadius()[0]);
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
                                        PrefDivComparator pdc = new PrefDivComparator(p,"Target",g,clusters);

                                        try {
                                            predAccuracy = pdc.getPredictionAccuracy(causes, causesTest);
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

                                        pdc = new PrefDivComparator(noPrior,"Target",g,clusters);

                                        try {
                                            predAccuracy = pdc.getPredictionAccuracy(causes, causesTest);
                                        }catch(Exception e)
                                        {
                                            predAccuracy = Double.NaN;
                                            System.out.println(causes);
                                        }
                                        System.out.println("Done");

                                        out.get(curr).print(predAccuracy + "\t");




                                    }

                                    /***Code to get optimal predictions***/
                                    List<Node> causes = g.getAdjacentNodes(g.getNode("Target"));
                                    ArrayList<Gene> causalGenes = new ArrayList<Gene>();

                                    for(Node n:causes)
                                    {
                                        Gene temp = new Gene(0);
                                        temp.symbol = n.getName();
                                        causalGenes.add(temp);
                                    }

                                    DataSet train = RunPrefDiv.summarizeData(toRun, causalGenes, lastCluster, RunPrefDiv.ClusterType.NONE);
                                    train = getCausalFeatures(train,"Target");

                                    DataSet summarizedTest = RunPrefDiv.summarizeData(test, causalGenes, lastCluster, RunPrefDiv.ClusterType.NONE);
                                    DataSet causesTest = subsetByCols(train,summarizedTest);

                                    PiPrefDiv4 ignore = new PiPrefDiv4(d,"Target",10);

                                    PrefDivComparator pdc = new PrefDivComparator(ignore,"Target",g,clusters);


                                    out.get(curr).println(pdc.getPredictionAccuracy(train,causesTest));
                                }
                                /***Summarization Method Tests***/
                                else if (x == 3) {
                                    out.get(curr).print(j + "\t");
                                    Map<Gene, List<Gene>> lastCluster = p.getLastCluster();
                                    Map<Gene,List<Gene>> npCluster = noPrior.getLastCluster();

                                    for (int i = 0; i < allTypes.length; i++) {

                                        DataSet summarized = RunPrefDiv.summarizeData(toRun, selected, lastCluster, allTypes[i]);
                                        DataSet summarizedTest = RunPrefDiv.summarizeData(test, selected, lastCluster, allTypes[i]);

                                        PrefDivComparator pdc = new PrefDivComparator(p,"Target",g,clusters);
                                        try {
                                            predAccuracy = pdc.getPredictionAccuracy(summarized, summarizedTest);
                                        }
                                        catch(Exception e)
                                        {
                                            System.out.println(summarized);
                                            predAccuracy = Double.NaN;
                                        }

                                        out.get(curr).print(predAccuracy + "\t");

                                        summarized = RunPrefDiv.summarizeData(toRun, npGenes, npCluster, allTypes[i]);
                                        summarizedTest = RunPrefDiv.summarizeData(test, npGenes, npCluster, allTypes[i]);

                                        pdc = new PrefDivComparator(noPrior,"Target",g,clusters);
                                        try {
                                            predAccuracy = pdc.getPredictionAccuracy(summarized, summarizedTest);
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
                                /***Experiment changing the selected # of clusters***/
                                else if(x==6)
                                {
                                    for(int idx = 0; idx < numFeatsTotal;idx++)
                                    {
                                        /***Run with priors***/

                                        PiPrefDiv4 wp = new PiPrefDiv4(toRun, "Target", numFeats[idx], numParams);
                                        wp.setSubsamples(subs);
                                        wp.setParallel(false);
                                        wp.setPartialCorrs(partialCorr);
                                        wp.selectGenes(boot,numSamples,dFile);

                                        PrefDivComparator comp2 = new PrefDivComparator(wp,"Target",g,clusters);

                                        out.get(curr).print(j + "\tT\t" + numFeats[idx] + "\t" + comp2.getNumClusters()/(double)numComponents + "\t" + comp2.getTrueClusters(numComponents)/(double)minTargetParents);
                                        out.get(curr).println("\t" + comp2.getFriendlyCluster(numComponents));


                                        if(ii==0 && jj==0) {
                                            /***Run without priors if we haven't played with the number of priors/amount reliable yet***/

                                            PiPrefDiv4 np = new PiPrefDiv4(toRun, "Target", numFeats[idx], numParams);
                                            np.setSubsamples(subs);
                                            np.setParallel(false);
                                            np.setPartialCorrs(partialCorr);
                                            np.selectGenes(boot, numSamples);
                                            PrefDivComparator npComp2 = new PrefDivComparator(np, "Target", g, clusters);

                                            out.get(curr).print(j + "\tF\t" + numFeats[idx] + "\t" + npComp2.getNumClusters() / (double) numComponents + "\t" + npComp2.getTrueClusters(numComponents) / (double) minTargetParents);
                                            out.get(curr).println("\t" + npComp2.getFriendlyCluster(numComponents));
                                        }


                                    }
                                    //Run, Are Priors Included?, # of clusters selected, % of all clusters picked up, % of true clusters picked up, friendly cluster accuracy


                                }
                                /***Experiment changing the wp Cutoff***/
                                else if(x==7)
                                {
                                    for(int idx = 0; idx < numWp;idx++)
                                    {
                                        /***Run with priors***/

                                        PiPrefDiv4 wp = new PiPrefDiv4(toRun, "Target", minTargetParents, numParams);
                                        wp.setSubsamples(subs);
                                        wp.setParallel(false);
                                        wp.setPartialCorrs(partialCorr);
                                        wp.setWPCutoff(wpCutoff[idx]);
                                        wp.selectGenes(boot,numSamples,dFile);


                                        PrefDivComparator comp2 = new PrefDivComparator(wp,"Target",g,clusters);


                                        //Run, Are priors included?, Parameter, Representation Acc, Cluster Acc, Friendly Cluster Acc, Prediction Acc

                                        out.get(curr).print(j + "\tT\t" + wpCutoff[idx] + "\t" + comp2.getFeatAccuracy("ACC")+ "\t" + comp2.getClusterAccuracy(numComponents));
                                        out.get(curr).println("\t" + comp2.getFriendlyCluster(numComponents) + "\t" + comp2.getPredictionAccuracy(toRun,test));



                                    }
                                    out.get(curr).print(j + "\tF\t" + wpLow + "\t" + npCompare.getFeatAccuracy("ACC") + "\t" + npCompare.getClusterAccuracy(numComponents));
                                    out.get(curr).println("\t" + npCompare.getFriendlyCluster(numComponents) + "\t" + npCompare.getPredictionAccuracy(toRun,test));
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


    /***
     *
     * Standardizes both datasets based on the mean and sd of the training set (prevent information leak to testing set)
     * @param train Trainind Dataset
     * @param test Testing Dataset
     */
    public static void standardizeData(DataSet train, DataSet test)
    {
        TetradMatrix trainData = train.getDoubleData();
        for(int i = 0; i < train.getNumColumns();i++)
        {
            double [] col = trainData.getColumn(i).toArray();
            double mean = StatUtils.mean(col);
            double sd = StatUtils.sd(col);
            for(int j = 0; j < train.getNumRows();j++)
            {
                train.setDouble(j,i,(train.getDouble(j,i)-mean)/sd);
            }

            for(int j = 0; j < test.getNumRows();j++)
            {
                test.setDouble(j,i,(test.getDouble(j,i)-mean)/sd);
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


    /*** Run graphical modeling algorithm and get the features related to the target**/
    public static DataSet getCausalFeatures(DataSet input,String target)
    {
        try{
            input = runPriors.addDummy(input);
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
     * @param dPriors Which relationships do we have prior information for?
     * @param clusters True clusters in the data generating graph
     * @param g Data generating graph
     * @param toRun Training Data
     * @param test Testing Data
     * @return 1-D Array with two entries: Best Radii NP, Best Radii -> Both are with respect to cluster accuracy
     */
    public static double[] getBestParamsPar(final double [] radiiNP, final double[] radii, final int numComponents,final PiPrefDiv4 p, final String [] dPriors, final Map<String,List<Integer>>clusters, final Graph g, final DataSet toRun,final DataSet test)
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


            private synchronized void addAccuracy(PrefDivComparator comp, int i,int j)
            {
                result[i][j] = comp.getClusterAccuracy(numComponents);

            }

            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        //RadiiWP is the columns, RadiiNP is the rows
                        int j = s % radii.length; //Column number (radiiWP)
                        int i = s / radii.length; //Row Number (radiiNP)
                        try {
                            p.selectGenes(radiiNP[i],radii[j],dPriors);
                            PrefDivComparator comp = new PrefDivComparator(p,"Target",g,clusters);
                            addAccuracy(comp,i,j);


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
     * @param clusters The true clusters in the data generating graph
     * @param g The data generating graph
     * @param toRun Training dataset
     * @param test Testing dataset
     * @return The radii corresponding to the best cluster accuracy
     */
    private static double [] getBestNPParamsPar(final double [] radii,final int numComponents,final PiPrefDiv4 p, final  Map<String,List<Integer>>clusters, final Graph g, final DataSet toRun,final DataSet test)
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


            private synchronized void addAccuracy(PrefDivComparator pc, int i, String type)
            {
                result[i][2] = pc.getClusterAccuracy(numComponents);


                result[i][0] = pc.getFeatAccuracy(type);
                result[i][1] = pc.getPredictionAccuracy(toRun,test);
            }

            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        try {
                            List<Gene> selected = p.selectGenes(radii[s]);
                            PrefDivComparator pdc = new PrefDivComparator(p,"Target",g,clusters);
                            String type = "F1";
                            if(selected.size()==g.getAdjacentNodes(g.getNode("Target")).size())
                                type = "ACC";
                            addAccuracy(pdc,s,type);


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
    private static double[] getBestNPParams(double [] radii, int numComponents,PiPrefDiv4 p, Map<String,List<Integer>> clusters,Graph g, DataSet toRun,DataSet test)
    {
        double[][] result = new double[radii.length][3];
                for (int i = 0; i < radii.length; i++)
                {
                        List<Gene> selected = p.selectGenes(radii[i]);

                        PrefDivComparator pdc = new PrefDivComparator(p,"Target",g,clusters);

                        result[i][2] = pdc.getFriendlyCluster(numComponents);

                        String type = "F1";
                        if (selected.size() == g.getAdjacentNodes(g.getNode("Target")).size())
                            type = "ACC";
                        result[i][0] = pdc.getFeatAccuracy(type);
                        result[i][1] = pdc.getPredictionAccuracy(toRun, test);

              }
        return bestAcc(result,radii);
    }

    /*** Serial implementation of the above parallel method***/
    private static double[] getBestParams(double [] radiiNP, double [] radii, int numComponents,PiPrefDiv4 p, String [] dFile,Map<String,List<Integer>> clusters,Graph g, DataSet toRun,DataSet test)
    {

            double[] result = new double[3];

            double [][]allResults = new double[radii.length][radiiNP.length];

            //Find optimal parameters for no priors, and then loop again for with priors
            for (int i = 0; i < radii.length; i++) {


                for (int j = 0; j < radiiNP.length; j++) {

                    p.selectGenes(radiiNP[j], radii[i], dFile);

                    PrefDivComparator pdc = new PrefDivComparator(p, "Target", g, clusters);
                    double acc = pdc.getFriendlyCluster(numComponents);
                    allResults[i][j] = acc;


                    if (acc > result[2]) {
                        result[2] = acc;
                        result[1] = radii[i];
                        result[0] = radiiNP[j];
                    }
                }
            }
        try {
            PrintStream out =new PrintStream("Parameter_Search.txt");
            for(int j = 0; j < radiiNP.length;j++)
              out.print(radiiNP[j] + "\t");
            out.println();

            for(int i = 0; i < radii.length;i++)
            {
                out.print(radii[i] + "\t");
                for(int j =0 ; j < radiiNP.length;j++)
                {
                    out.print(allResults[i][j] + "\t");

                }
                out.println();
            }

            return result;
        }
        catch(Exception e)
        {
            e.printStackTrace();
            System.exit(-1);
        }
        return null;
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
