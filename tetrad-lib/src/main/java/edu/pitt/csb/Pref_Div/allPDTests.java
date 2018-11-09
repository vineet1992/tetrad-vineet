package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * The Purpose of this class is to produce the first figure for the experiments section of the Pref-Div Paper
 * How well can we evaluate the prior knowledge we are getting from both intensity and dissimilarity files?
 */
public class allPDTests {

    public static void main(String [] args) {


        int experiment = 0; //0 -> prior evaluation, 1 -> Accuracy of chosen parameters vs Optimal, 2-> Comparing different summarization methods for clustered genes
        //3 -> Evaluating how good cross validation is to select alpha parameters, 4-> Sensitivity analysis on number of folds, number of subsamples, number of parameters, ss, numvars, etc.


        int numRuns = 20;
        int numGenes = 200;
        int sampleSize = 200;
        double amountPrior = 0.6;//percentage of edges to have prior knowledge for
        boolean boot = false; //Should we use bootstrap samples for PiPrefDiv
        boolean loocv = false; //Should we use leave-one-out CV for PiPrefDiv
        boolean useCausalGraph = false; //Are we selecting genes to then use a causal graph to find connections or direct selection?
        int numSamples = 20; //Number of bootstrap/sub-sampled samples to use
        int numParams = 10;//Number of parameters to sweep over
        boolean noiseRandom = true;//Is the reliability range of reliable and unreliable priors set randomly?



        int numPriors = 10; //Number of prior knowledge sources
        int numReliable = 5; //Number of reliable sources
        int numComponents = 10; //How many components do we have for cluster simulation?
        int minTargetParents = 5; //How many true parents of the target are there?
        boolean amountRandom = false; //Should the priors have a random amount of prior knowledge?
        boolean targetContinuous = true; //Is the target variable continuous?
        boolean evenDistribution = true; //Is the distribution of nodes in each cluster even?
        int numCategories = 4; //number of categories for discrete variables

        String [] types = new String[]{"FS","Prediction","Cluster"}; //Help to create file header for experiment 1
        String [] starts = new String[]{"Best_Radius","Best_Threshold","Best_Accuracy","Pred_Accuracy"}; //More help


        File rFile = new File("Results");
        if (!rFile.isDirectory())
            rFile.mkdir();

        try {
            String start = "Prior_Evaluation_";
            if(experiment==1)
            {
                start = "Parameter_Accuracy_";
            }
            PrintStream out = new PrintStream("Results/" + start + numGenes + "_" + sampleSize + "_" + numRuns + "_" + numPriors + "_" + amountPrior + "_" + amountRandom + "_" + boot + "_" + loocv + "_"  + evenDistribution + "_" + targetContinuous + "_" + useCausalGraph +  "_" + numReliable + "_" + numComponents + "_" + minTargetParents + "_" + numParams + ".txt");

            if(experiment==0)
            {
                out.println("Run\tAmount_Prior\tReliable?\tIntensity?\tPredicted_Weight\tActual_Reliability");
            }
            else if(experiment==1)
            {
                out.print("Run\tPredicted_Radius\tPredicted_Threshold");

                for(int i = 0; i < types.length;i++)
                {
                    for(int j = 0; j < starts.length;j++)
                    {
                        out.print(starts[j] + "_" + types[i] + "\t");
                    }
                }
                out.println();
            }
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
                String priorIntensity = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_"  + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" +  noiseRandom + "_" +  j + "_intensity.txt";
                String[] dFile = new String[numPriors];
                for (int k = 0; k < numPriors; k++) {
                    dFile[k] = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_"  + noiseRandom + "_" + k + "_" + j +  "_dissimilarity.txt";
                }

                /***LOAD RELIABILITY SCORES***/
                File rd = new File("Reliabilities_D_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_"  + j + ".txt");
                File ri = new File("Reliabilities_I_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");
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
                DataSet toRun = d;
                DataSet test = d;
                if(experiment==1)
                {
                    int trainLength = (int)(toRun.getNumRows()*0.9);
                    int testLength = toRun.getNumRows()-trainLength;
                    int [] trainInds = new int[trainLength];
                    for(int i = 0; i < trainLength;i++)
                    {
                        trainInds[i] = i;
                    }
                    toRun = d.subsetRows(trainInds);
                    test.removeCols(trainInds);
                }

                PiPrefDiv p = new PiPrefDiv(toRun, "Target", minTargetParents,numParams);
                p.setSubsamples(subs);
                p.setVerbose();
                File f = new File("Results/Detailed_Score_Evaluation/");
                if(experiment==1)
                {
                    if(!f.exists())
                        f.mkdir();
                    p.setOutputScores(new PrintStream(f.getName() + "Parameter_Scores_" + j + ".txt"));
                }


                if(experiment==0) //Evaluate prior knowledge
                {
                    double[][] weights = p.evaluatePriors(boot, numSamples, priorIntensity, dFile, useCausalGraph);

                    /***PRINT RESULTS TO FILE***/
                    for (int i = 0; i < numPriors; i++) {
                        out.println(j + "\t" + getAmountPrior(priorIntensity, i) / ((double) g.getNumNodes() - 1) + "\t" + (i < numReliable) + "\tT\t" + weights[0][i] + "\t" + reliInt[i]);
                    }
                    for (int i = 0; i < numPriors; i++) {
                        out.println(j + "\t" + (double) (getAmountPrior(dFile[i])) / (g.getNumNodes() * (g.getNumNodes() - 1) / 2) + "\t" + (i < numReliable) + "\tF\t" + weights[1][i] + "\t" + reliDis[i]);
                    }
                }
                else if(experiment==1) //Evaluate Parameter accuracy
                {
                    ArrayList<Gene> selected = p.selectGenes(boot,numSamples,priorIntensity,dFile,useCausalGraph);
                    Map<Gene,List<Gene>> lastCluster = p.getLastCluster();
                    double clustAccuracy = getClusterAccuracy(lastCluster,clusters,numComponents);
                    double featAccuracy = -1;
                    if(minTargetParents==g.getAdjacentNodes(g.getNode("Target")).size())
                    {
                        //Correct amount of selections so use accuracy
                        featAccuracy = getFeatAccuracy(selected,g,"Target","ACC");

                    }
                    else
                    {
                        //Incorrect amount of selected variables so use F1
                        featAccuracy = getFeatAccuracy(selected,g,"Target","F1");

                    }
                    double predAccuracy = getAccuracy(toRun,test,selected);

                    double [] radii = p.getTestedRadii();
                    double [] thresholds = p.getTestedThresholds();

                    double [] intWeights = p.getLastIntensityWeights();
                    double [] simWeights = p.getLastSimilarityWeights();

                    ArrayList<Gene> meanGenes = PiPrefDiv.loadGenes(priorIntensity,intWeights);
                    float [] meanDis = PiPrefDiv.loadTheoryFiles(dFile,simWeights,numGenes);
                    //Use the weights from PiPrefDiv, we are interested in determining how good param selection is given that weights are good
                    //Because experiment 1 demonstrated that weights are good

                    RunPrefDiv rpd = new RunPrefDiv(meanDis,meanGenes,toRun,"Target",loocv);

                    rpd.setTopK(minTargetParents);
                    rpd.setAccuracy(0);
                    rpd.setNumAlphas(numParams);
                    rpd.setNS(subs.length);
                    rpd.setNumFolds(subs.length);
                    rpd.setCausalGraph(useCausalGraph);
                    rpd.useStabilitySelection();
                    rpd.clusterByCorrs();


                    //for every radii and threshold combination, give the three accuracies
                     double[][][] results = getBestParams(radii,thresholds,rpd,clusters,g,toRun,test); //Feature Select, Cluster Acc, Prediciton Acc



                    out.print(j + "\t" + p.getLastRadius() + "\t" + p.getLastThreshold());

                    //print out results to file
                    printResults(results, out,radii,thresholds,new double[]{featAccuracy,predAccuracy,clustAccuracy});

                    printDetailedScores(results,new PrintStream(f.getName() + "Detailed_Results_" + j));



                    //Roll through all tested radii,threshold combos and RunPrefDiv directly to get accuracies
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

    public static void printResults(double[][][]results, PrintStream out, double [] radii, double[] thresholds, double [] predAcc)
    {
        String [] types = new String[]{"FS","Prediction","Cluster"}; //Help to create file header for experiment 1
        String [] starts = new String[]{"Best_Radius","Best_Threshold","Best_Accuracy","Pred_Accuracy"}; //More help

        /**Feature Selection Results Printing**/
        double[] bestRadii = new double[3];
        double []maxValue = new double[3];
        double [] bestThreshold = new double[3];
        for(int i = 0; i < results.length;i++)
        {
            for(int j = 0; j < results[i].length;j++)
            {
                for(int k = 0; k < 3;k++)
                {
                    if(results[i][j][k]>maxValue[k])
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

    public static double[][][] getBestParams(double [] radii, double [] threshold, RunPrefDiv rpd, Map<String,Integer> clusters,Graph g, DataSet toRun,DataSet test)
    {
        double[][][] result = new double[radii.length][threshold.length][3];
        for(int i = 0; i < radii.length;i++)
        {
            for(int j =0 ; j < threshold.length;j++)
            {
                rpd.setRadius(radii[i]);
                rpd.setThreshold(threshold[j]);
                List<Gene> selected = rpd.runPD();
                Map<Gene,List<Gene>> lastCluster = rpd.getClusters();
               result[i][j][2] = getClusterAccuracy(lastCluster,clusters,selected.size());

                String type = "F1";
                if(selected.size()==g.getAdjacentNodes(g.getNode("Target")).size())
                    type = "ACC";
               result[i][j][0] = getFeatAccuracy(selected,g,"Target",type);
                result[i][j][1] = getAccuracy(toRun,test,selected);



            }
        }
        return result;
    }


    //TODO Incorporate discrete target variables
    public static double getAccuracy(DataSet train, DataSet test, List<Gene> selected)
    {
        RegressionDataset rd = new RegressionDataset(train);
        List<Node> regressors = new ArrayList<Node>();
        for(int i = 0; i < selected.size();i++)
            regressors.add(train.getVariable(selected.get(i).symbol));
        RegressionResult res = rd.regress(train.getVariable("Target"),regressors);
        int [] cols = new int[selected.size()];

        for(int i = 0; i < selected.size();i++)
        {
            cols[i] = test.getColumn(test.getVariable(selected.get(i).symbol));
        }
        double[]pred = new double[test.getNumRows()];
        double [] actual = new double[test.getNumRows()];
        for(int i = 0; i < test.getNumRows();i++)
        {
            double [] x = new double[cols.length];
            for(int j = 0; j < cols.length;j++) {
                x[j] = test.getDouble(i, j);

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
        return Math.sqrt(err);
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
