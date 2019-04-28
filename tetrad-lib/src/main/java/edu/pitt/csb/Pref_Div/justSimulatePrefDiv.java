package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.calculator.expression.Expression;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.GeneralizedSemIm;
import edu.cmu.tetrad.sem.GeneralizedSemPm;
import edu.cmu.tetrad.util.StatUtils;
import edu.pitt.csb.mgm.MixedUtils;
import nu.xom.Nodes;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Created by vinee_000 on 3/2/2018.
 */
public class justSimulatePrefDiv {


    public static double[] reliableParams = new double[]{10,2,2,10}; //True Low, True High, False Low, False High

    public static double[] unreliableParams = new double[]{4,4,4,4}; //True Low, True High, False Low, False High

    /***For when random priors are used***/
    private static double alphaMean = 10;
    private static double betaMean = 2;
    private static double unreliAlphaMean = 4;
    private static double unreliBetaMean = 4;




    public static double amountLimit = 1.0; //Specifies the total percentage of prior information that can be given by a single source

    public static boolean perfectPriors = false; // Priors give perfect information only about true relationships

    public static boolean simulateTogether = true; //Simulate both intensity and dissimilarity priors in a single file

    public static boolean latentControl = true; //Includes a master regulator for each cluster


    public static boolean realPathways = false; //Pathways are more realistic (not tightly clustered, distinct groups)
    public static boolean priorsToClusters = true; //Only important if real pathway simulation -> Is "true" info cluster membership? Or is it edge existence?

    private static boolean unreliableRandom = true; //FOR PRIOR KNOLWEDGE EVALUATION ONLY


    public static Random rand = new Random();

    public static int numGenes = 500;

    public static void main(String[] args) {
        double[] ap = new double[]{0.5};
        int[] nr = new int[]{1,3,5};
        int[] ss = new int[]{50,200};

        for (int xx = 0; xx < ap.length; xx++) {
            for (int y = 0; y < nr.length; y++) {

                for (int z = 0; z < ss.length; z++) {


                    int numRuns = 15;
                    int sampleSize = ss[z];
                    double amountPrior = ap[xx];//percentage of edges to have prior knowledge for
                    boolean boot = false; //Should we use bootstrap samples for PiPrefDiv
                    boolean loocv = false; //Should we use leave-one-out CV for PiPrefDiv
                    int numSamples = 20; //Number of bootstrap/sub-sampled samples to use

                    int avgPathways = 1; //Average number of pathways a gene belongs to (only for realPathways simulation)
                    int avgCauses = 1; //Average number of genes in a relevant pathway causally related to the target

                    double varLow = 0.01; //Low value for variance of variables
                    double varHigh = 2;  //High value for variance of variables

                    int numPriors = 5; //Number of prior knowledge sources
                    int numReliable = nr[y]; //Number of reliable sources
                    int numComponents = 50; //How many components do we have for cluster simulation?
                    int minTargetParents = 25; //How many true parents of the target are there?
                    boolean noiseRandom = true; //Should the priors have random amounts of reliability?


                    boolean amountRandom = true; //Should the priors have a random amount of prior knowledge?
                    boolean targetContinuous = true; //Is the target variable continuous?
                    boolean evenDistribution = true; //Is the distribution of nodes in each cluster even?
                    int numCategories = 4; //number of categories for discrete variables


                    int index = 0;
                    while (index < args.length) {
                        if (args[index].equals("-runs")) {
                            numRuns = Integer.parseInt(args[index + 1]);
                            index += 2;
                        } else if (args[index].equals("-ng")) {
                            numGenes = Integer.parseInt(args[index + 1]);
                            index += 2;
                        } else if (args[index].equals("-ss")) {
                            sampleSize = Integer.parseInt(args[index + 1]);
                            index += 2;
                        } else if (args[index].equals("-numParents")) {
                            minTargetParents = Integer.parseInt(args[index + 1]);
                            index += 2;
                        } else if (args[index].equals("-nc")) {
                            numComponents = Integer.parseInt(args[index + 1]);
                            index += 2;
                        } else if (args[index].equals("-categorical")) {
                            targetContinuous = false;
                            index++;
                        } else if (args[index].equals("-even")) {
                            evenDistribution = true;
                            index++;
                        } else if (args[index].equals("-numPriors")) {
                            numPriors = Integer.parseInt(args[index + 1]);
                            index += 2;

                        } else if (args[index].equals("-numReliable")) {
                            numReliable = Integer.parseInt(args[index + 1]);
                            index += 2;
                        } else if (args[index].equals("-amountPrior")) {
                            amountPrior = Double.parseDouble(args[index + 1]);
                            index += 2;
                        } else if (args[index].equals("-numSamples")) {
                            numSamples = Integer.parseInt(args[index + 1]);
                            index += 2;
                        } else if (args[index].equals("-amountRandom")) {
                            amountRandom = true;
                            index++;
                        } else if (args[index].equals("-noiseRandom")) {
                            noiseRandom = true;
                            index++;
                        } else if (args[index].equals("-boot")) {
                            boot = true;
                            index++;
                        } else if (args[index].equals("-loocv")) {
                            loocv = true;
                            index++;
                        } else if (args[index].equals("-even")) {
                            evenDistribution = true;
                            index++;
                        } else {
                            System.err.println("Unidentified argument: " + args[index]);
                            System.exit(-1);
                        }
                    }


                    boolean saveData = true;

                    //out variables refer to outputs
                    //correctness column to the stability file for the first boolean

                    //If either of outAccuracy or outStabAcc are true, then

                    //Create Directories to save results
                    File rFile = new File("Results");
                    if (!rFile.isDirectory())
                        rFile.mkdir();
                    File cFile = new File("Clusters");
                    File gFile = new File("Graphs");
                    File dFile = new File("Data");
                    File pFile = new File("Priors");
                    File subFile = new File("Subsamples");
                    if (saveData) {
                        if (!cFile.isDirectory())
                            cFile.mkdir();
                        if (!gFile.isDirectory())
                            gFile.mkdir();
                        if (!dFile.isDirectory())
                            dFile.mkdir();
                        if (!pFile.isDirectory())
                            pFile.mkdir();
                        if (!subFile.isDirectory())
                            subFile.mkdir();
                    }

                    for (int j = 0; j < numRuns; j++) {

                        System.out.println("Simulate data for run number: " + j);
                        HashMap<String, List<Integer>> clusters = new HashMap<String, List<Integer>>();

                        Graph g = null;
                        DataSet d = null;
                        File graphFile = null;


                        /**GENERATE GRAPH**/

                        System.out.print("Generating Graph...");

                        graphFile = new File("Graphs/Graph_" + numGenes + "_" + minTargetParents + "_" + numComponents + "_" + evenDistribution + "_" + j + ".txt");
                        if (graphFile.exists()) {


                            g = GraphUtils.loadGraphTxt(graphFile);

                            for (Node n : g.getNodes()) {
                                if (n.getName().startsWith("L"))
                                    n.setNodeType(NodeType.LATENT);
                            }
                        }
                        File dataFile = null;
                        dataFile = new File("Data/Dataset_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + j + ".txt");
                        try {
                            if (dataFile.exists())
                                d = MixedUtils.loadDataSet2(dataFile.getAbsolutePath());
                        } catch (Exception e) {
                            System.err.println("Unable to load dataset");
                            System.exit(-1);
                        }


                        MixedLeeHastieSimulation m = new MixedLeeHastieSimulation();
                        m.setRealPathways(realPathways);

                        m.setLatentControl(latentControl);

                        Parameters p = new Parameters();

                        if (realPathways)
                            p.setValue("numMeasures", numGenes);
                        else
                            p.setValue("numMeasures", numGenes + 1);
                        p.setValue("sampleSize", sampleSize);
                        if (targetContinuous)
                            p.setValue("percentDiscreteForMixedSimulation", 0);
                        else
                            p.setValue("percentDiscreteForMixedSimulation", 100 / (double) numGenes);
                        p.setValue("numCategories", numCategories);

                        p.setValue("numConnectedComponents", minTargetParents);
                        p.setValue("numComponents", numComponents);
                        if (targetContinuous)
                            p.setValue("targetContinuous", 1);
                        else
                            p.setValue("targetContinuous", 0);

                        p.setValue("varLow",varLow);
                        p.setValue("varHigh",varHigh);

                        if (g != null) {
                            m.setTrueGraph(g);
                        }
                        if (realPathways) {
                            p.setValue("avgPathways", avgPathways);
                            p.setValue("avgCauses", avgCauses);
                        }
                        clusters = m.simulateClustered(p, evenDistribution);

                        System.out.println("Done");
                        System.out.println("Graph: " + m.getTrueGraph());

                        System.out.println("Data: " + m.getDataSet(0));

                        System.out.println("Clusters: " + clusters);

                        File clusterFile = null;
                        clusterFile = new File("Clusters/Cluster_" + numGenes + "_" + minTargetParents + "_" + numComponents + "_" + evenDistribution + "_" + j + ".txt");
                        if (clusterFile.exists()) {
                            clusters = new HashMap<String, List<Integer>>();
                            try {
                                BufferedReader bTemp = new BufferedReader(new FileReader(clusterFile));
                                while (bTemp.ready()) {
                                    String[] line = bTemp.readLine().split("\t");
                                    List<Integer> temp = new ArrayList<Integer>();
                                    for (int x = 1; x < line.length; x++)
                                        temp.add(Integer.parseInt(line[x]));
                                    clusters.put(line[0], temp);
                                }
                            } catch (Exception e) {
                                System.err.println("Couldn't load clusters from file");
                                e.printStackTrace();
                                System.exit(-1);
                            }

                        }


                        String target = "Target";
                        if (d == null) {
                            d = m.getDataSet(0);

                        } else
                            m.setDataSet(d, 0);
                        g = m.getTrueGraph();


                        if (!simulateTogether) {
                            /***Consistent ordering is whatever is returned by g.getNodes() (ignoring the target variable)**/
                            File priorFile = null;
                            priorFile = new File("Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + "_intensity.txt");
                            File iFile = new File("Reliabilities_I_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");
                            if (!priorFile.exists()) {
                                if (perfectPriors)
                                    generatePerfectIntensity(priorFile, iFile, g, target, numPriors, amountPrior, amountRandom);
                                else
                                    generateIntensity(priorFile, iFile, g, target, numPriors, numReliable, amountPrior, amountRandom, noiseRandom);

                            }
                            System.out.println("Done");
                        }


                        /**Use this if prior file exists Make sure intensity is at the end**/
                        String[] disFiles = new String[numPriors];
                        try {
                            File temp = new File("Reliabilities_D_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");
                            if (simulateTogether)
                                temp = new File("Reliabilities_D_All_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");

                            if (!temp.exists()) {
                                PrintStream out2;
                                if (simulateTogether)
                                    out2 = new PrintStream("Reliabilities_D_All_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");
                                else
                                    out2 = new PrintStream("Reliabilities_D_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");

                                /***First decide the probability for each edge to be included in the information source***/
                                float[] prob;
                                if (simulateTogether)
                                    prob = generateDisAmounts(((numGenes + 1) * numGenes) / 2);
                                else
                                    prob = generateDisAmounts((numGenes * (numGenes - 1)) / 2);

                                for (int k = 0; k < numPriors; k++) {

                                    System.out.print("Generating Dissimilarity File #" + k + "...");


                                    if (simulateTogether)
                                        disFiles[k] = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + k + "_" + j + ".txt";
                                    else
                                        disFiles[k] = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + k + "_" + j + "_dissimilarity.txt";
                                    File f = new File(disFiles[k]);
                                    if (!f.exists()) {
                                        double reli = -1;

                                        /***Correct relationship means do they belong to the same cluster?***/
                                        if (priorsToClusters || latentControl) {
                                            reli = generateDissimilarityClustered(f, g, d, clusters, target, numReliable, amountPrior, amountRandom, noiseRandom, k, prob, numComponents);
                                        }
                                        /****Correct relationship means are they connected in the graph?***/
                                        else {
                                            if (perfectPriors)
                                                reli = generatePerfectDis(f, g, d, target, amountPrior, amountRandom, k, prob, simulateTogether);
                                            else
                                                reli = generateDissimilarity(f, g, d, target, numReliable, amountPrior, amountRandom, noiseRandom, k, prob, simulateTogether);
                                        }
                                        out2.println(reli);
                                        out2.flush();
                                        System.out.println("Done");

                                    }

                                    //End Generation of Prior Knowledge
                                }

                                out2.close();

                            }

                        } catch (Exception e) {
                            e.printStackTrace();
                            System.exit(-1);
                        }


                        System.out.print("Generating Subsamples...");
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
                            subs = PiPrefDiv4.genSubsamples(boot, numSamples, d, loocv);
                        }


                        System.out.print("Generating training subsamples...");
                        File subsFile2 = new File("Subsamples/Subsample_Train_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + boot + "_" + numSamples + "_" + loocv + "_" + j + ".txt");
                        int[][] subs2 = null;
                        if (subsFile2.exists()) {
                            try {
                                BufferedReader b = new BufferedReader(new FileReader(subsFile2));
                                int lineCount = 0;
                                while (b.ready()) {
                                    b.readLine();
                                    lineCount++;
                                }
                                subs2 = new int[lineCount][];
                                b = new BufferedReader(new FileReader(subsFile2));
                                lineCount = 0;
                                while (b.ready()) {
                                    String[] line = b.readLine().split("\t");
                                    subs2[lineCount] = new int[line.length];
                                    for (int x = 0; x < line.length; x++) {
                                        subs2[lineCount][x] = Integer.parseInt(line[x]);
                                    }
                                    lineCount++;
                                }
                                b.close();
                            } catch (Exception e) {
                                System.err.println("Unable to load training subsamples");
                                e.printStackTrace();
                                System.exit(-1);
                            }
                        } else {
                            int trainLength = (int) (d.getNumRows() * 0.9);
                            int[] trainInds = new int[trainLength];
                            for (int i = 0; i < trainLength; i++) {
                                trainInds[i] = i;
                            }
                            DataSet toRun = d.subsetRows(trainInds);
                            subs2 = PiPrefDiv4.genSubsamples(boot, numSamples, toRun, loocv);
                        }

                        System.out.println("Done");
/*
                RunPrefDiv r = new RunPrefDiv(dissimilarity, genes, d, target, false);

                if (subs != null)
                    r.setSubs(subs);
                r.useStabilitySelection();
                r.setStabilityOutput(outStabAcc);
                r.setTargets(g.getAdjacentNodes(g.getNode(target)));
                r.setRun(j);
                r.usePartialCorrelation(false);
                r.setCausalGraph(useCausalGraph);
                r.setB(B);*/

                        try {


                            PrintStream temp = new PrintStream(graphFile);
                            temp.println(g);
                            temp.flush();
                            temp.close();


                            temp = new PrintStream(clusterFile);
                            for (String s : clusters.keySet()) {
                                List<Integer> paths = clusters.get(s);

                                temp.print(s);
                                for (int x = 0; x < paths.size(); x++)
                                    temp.print("\t" + paths.get(x));
                                temp.println();
                                temp.flush();
                            }
                            temp.close();


                            temp = new PrintStream(dataFile);
                            temp.println(d);
                            temp.flush();
                            temp.close();

                            /***PRINT SUBSAMPLES OUT***/
                            temp = new PrintStream(subsFile);
                            for (int i = 0; i < subs.length; i++) {
                                for (int x = 0; x < subs[i].length; x++) {
                                    if (x == subs[i].length - 1)
                                        temp.println(subs[i][x]);
                                    else
                                        temp.print(subs[i][x] + "\t");
                                }
                            }
                            temp.flush();
                            temp.close();

                            temp = new PrintStream(subsFile2);
                            for (int i = 0; i < subs2.length; i++) {
                                for (int x = 0; x < subs2[i].length; x++) {
                                    if (x == subs2[i].length - 1)
                                        temp.println(subs2[i][x]);
                                    else
                                        temp.print(subs2[i][x] + "\t");
                                }
                            }
                            temp.flush();
                            temp.close();


                        } catch (Exception e) {
                            System.err.println("Couldn't write information to files " + j);
                            e.printStackTrace();
                            System.exit(-1);
                        }
                    }
                }
            }
        }
    }


    /****
     *
     * @param prior
     * @param intense
     * @param g
     * @param target
     * @param numPriors
     * @param amountPrior
     * @param amountRandom
     * @return
     */
    public static void generatePerfectIntensity(File prior, File intense, Graph g, String target, int numPriors, double amountPrior, boolean amountRandom) {
        List<Node> neighbors = g.getAdjacentNodes(g.getNode(target));
        List<Node> allNodes = g.getNodes();
        double[] reliability = new double[numPriors];
        int[] counts = new int[numPriors];
        double[] amount = new double[numPriors];
        for (int i = 0; i < numPriors; i++) {
            amount[i] = amountPrior;
            if (amountRandom)
                amount[i] = rand.nextDouble() * amountLimit;
        }

        try {
            PrintStream out = new PrintStream(prior);
            out.print("Gene\t");
            for (int i = 0; i < numPriors; i++) {
                out.print("Prior_" + i + "\t");
            }
            out.println();
            PrintStream out2 = new PrintStream(intense);

            for (int i = 0; i < allNodes.size(); i++) {
                Node curr = allNodes.get(i);
                if (curr.getName().equals(target))
                    continue;
                out.print(allNodes.get(i).getName() + "\t");

                /***Probability for all priors to include information about this relationship***/
                double prob = rand.nextDouble();
                for (int j = 0; j < numPriors; j++) {
                    /***Only give information if probability is less than amount and this is a true relationship***/

                    if (neighbors.contains(curr) && prob < amount[j]) {
                        double val = 1.0;
                        reliability[j] += val;

                        counts[j]++;
                        out.print(val + "\t");
                    } else {
                        out.print("-1\t");
                    }

                }
                out.println();
            }
            for (int i = 0; i < numPriors; i++) {
                reliability[i] /= counts[i];
                out2.println(reliability[i]);
            }
            out2.flush();
            out2.close();
            out.flush();
            out.close();

        } catch (Exception e) {
            System.err.println("Couldn't create intensity file");
            System.exit(-1);
        }

    }


    /***
     *
     * @param prior File to write the prior knowledge out to
     * @param g True simulated graph
     * @param target String representation of the target variable name
     * @param amountPrior What proportion of prior information should the source provide?
     * @param amountRandom Should the proportion be randomly chosen?
     * @param k Which index of dissimilarity prior are we on?
     * @param prob Probability of edge appearence
     * @param simulateTogether Should both intensity and similarity be simulated together
     * @return The reliability of this prior knowledge source
     */
    public static double generatePerfectDis(File prior, Graph g, DataSet data, String target, double amountPrior, boolean amountRandom, int k, float[] prob, boolean simulateTogether) {
        List<Node> allNodes = data.getVariables();


        double[][] dissim = new double[numGenes][numGenes];
        if (simulateTogether)
            dissim = new double[numGenes + 1][numGenes + 1];
        int row = 0;
        double reliability = 0;
        int counts = 0;
        double amount = amountPrior;
        if (amountRandom) {
            amount = rand.nextDouble() * amountLimit;
        }

        try {
            int index = 0;

            for (int i = 0; i < allNodes.size(); i++) {
                Node curr = g.getNode(data.getVariable(i).getName());

                if (curr.getName().startsWith("L") || (curr.getName().equals(target) && !simulateTogether))
                    continue;
                int col = 0;
                for (int j = 0; j < allNodes.size(); j++) {


                    Node curr2 = g.getNode(data.getVariable(j).getName());
                    if (curr2.getName().startsWith("L") || (curr2.getName().equals(target) && !simulateTogether))
                        continue;


                    //Currently reliable priors give information about true neighbors in the range [0.6,1] and false neighbors in the range [0,0.5]
                    if (i == j) {
                        dissim[row][col] = -1;
                    } else if (i > j) {
                        dissim[row][col] = dissim[col][row];
                    } else if (g.getEdge(curr, curr2) != null && prob[index] < amount) {
                        double val = 1.0;
                        dissim[row][col] = val;
                        reliability += val;
                        counts++;
                        index++;
                    } else {
                        dissim[row][col] = -1;
                        dissim[col][row] = -1;
                        index++;
                    }

                    col++;
                }
                row++;
            }
            reliability /= counts;


            PrintStream out = new PrintStream(prior);
            DecimalFormat df = new DecimalFormat("#.#####");
            for (int i = 0; i < dissim.length; i++) {
                for (int j = 0; j < dissim[i].length; j++) {
                    if (dissim[i][j] == -1) {
                        if (j == dissim[i].length - 1)
                            out.print("-1");
                        else
                            out.print("-1" + "\t");
                    } else {
                        if (j == dissim[i].length - 1)
                            out.print(df.format(dissim[i][j]));
                        else
                            out.print(df.format(dissim[i][j]) + "\t");
                    }
                }
                out.println();
            }

            out.flush();
            out.close();
            return reliability;

        } catch (Exception e) {
            System.err.println("Couldn't create dissimilarity file number " + k);
            e.printStackTrace();
            System.exit(-1);
        }

        return -1;
    }


    /***Same dissimilarity generation script, but here we use clusters instead of the graph to determine "true" adjacencies***/
    public static double generateDissimilarityClustered(File prior, Graph trueGraph, DataSet data, Map<String, List<Integer>> clusters, String target, int numReliable, double amountPrior, boolean amountRandom, boolean noiseRandom, int k, float[] prob, int numComponents) {
        /****Generate graph based on the clusters and call the usual method***/
        List<Node> nodes = new ArrayList<Node>();
        for (String s : clusters.keySet()) {
            nodes.add(new GraphNode(s));
        }
        Graph temp = new EdgeListGraphSingleConnections(nodes);

        List<List<Node>> gConnections = new ArrayList<List<Node>>();
        for (int i = 0; i < numComponents; i++) {
            gConnections.add(new ArrayList<Node>());
        }
        /***Get all nodes in cluster one and connect them***/
        for (int i = 0; i < nodes.size(); i++) {
            Node curr = nodes.get(i);
            /***Find all clusters that curr belongs to***/
            List<Integer> components = clusters.get(curr.getName());

            /***Add the node to these cluster groups***/
            for (int j = 0; j < components.size(); j++) {
                gConnections.get(components.get(j)).add(curr);
            }
        }

        /***For each component, add an edge between all of the pairs of nodes in the component***/
        for (int i = 0; i < gConnections.size(); i++) {
            List<Node> curr = gConnections.get(i);

            for (int j = 0; j < curr.size(); j++) {
                Node one = temp.getNode(curr.get(j).getName());
                for (int x = j + 1; x < curr.size(); x++) {
                    Node two = temp.getNode(curr.get(x).getName());
                    if (temp.getEdge(one, two) == null)
                        temp.addUndirectedEdge(one, two);

                }
            }
        }

        /***Add target variable and add all neighbors from the original graph***/
        temp.addNode(new GraphNode(target));

        for (Node n : trueGraph.getAdjacentNodes(trueGraph.getNode(target))) {
            temp.addUndirectedEdge(temp.getNode(n.getName()), temp.getNode(target));
        }


        /***Call approriate prior generation script***/
        if (perfectPriors)
            return generatePerfectDis(prior, temp, data, target, amountPrior, amountRandom, k, prob, simulateTogether);
        else
            return generateDissimilarity(prior, temp, data, target, numReliable, amountPrior, amountRandom, noiseRandom, k, prob, simulateTogether);

    }

    /****
     *
     * @param prior File Prior knowledge file
     * @param g True graph
     * @param target Target variable (string representation)
     * @param numReliable Number of reliable priors
     * @param amountPrior Proportion of reliable priors
     * @param amountRandom Should the amount this prior gives be random?
     * @param noiseRandom Should the noise for this prior be random?
     * @param k The index of this prior
     * @param prob Probability for each edge to have prior information
     * @param simulateTogether Should intensity and dissimilarity be generated together
     * @return Reliability measure for this prior
     *
     */
    public static double generateDissimilarity(File prior, Graph g, DataSet data, String target, int numReliable, double amountPrior, boolean amountRandom, boolean noiseRandom, int k, float[] prob, boolean simulateTogether) {




        List<Node> allNodes = data.getVariables();
        double[][] dissim = new double[numGenes][numGenes];

        if (simulateTogether)
            dissim = new double[numGenes + 1][numGenes + 1];
        int row = 0;
        double reliability = 0;
        int counts = 0;
        double amount = amountPrior;
        if (amountRandom) {
            amount = rand.nextDouble() * amountLimit;
        }
        if (noiseRandom) {
            NormalDistribution nd = new NormalDistribution(alphaMean,1);
            reliableParams[0] = nd.sample();
            reliableParams[3] = nd.sample();


            nd = new NormalDistribution(betaMean,1);
            double x = nd.sample();
            while(x <=0)
                x = nd.sample();
            reliableParams[1] = x;
            x = nd.sample();
            while(x <=0)
                x= nd.sample();
            reliableParams[2] = x;
        }

        if(unreliableRandom)
        {
            NormalDistribution nd = new NormalDistribution (unreliAlphaMean,1);
            unreliableParams[0] = nd.sample();
            unreliableParams[3] = nd.sample();
            nd = new NormalDistribution(unreliBetaMean,1);
            unreliableParams[1] = nd.sample();
            unreliableParams[2] = nd.sample();
        }


        double [] params = reliableParams;
        if(k >= numReliable)
        {
            params = unreliableParams;
        }

        BetaDistribution trueDist = new BetaDistribution(params[0],params[1]);
        BetaDistribution falseDist = new BetaDistribution(params[2],params[3]);



        try {
            int index = 0;

            for (int i = 0; i < allNodes.size(); i++) {
                /***Must follow the order of the dataset***/
                Node curr = g.getNode(data.getVariable(i).getName());


                /***Exclude latent variables and target variable, if simulated separately***/
                if (curr.getName().startsWith("L") || (curr.getName().equals(target) && !simulateTogether))
                    continue;


                int col = 0;
                for (int j = 0; j < allNodes.size(); j++) {

                    /***Same as above***/
                    Node curr2 = g.getNode(data.getVariable(j).getName());


                    if (curr2.getName().startsWith("L") || (curr2.getName().equals(target) && !simulateTogether))
                        continue;
                    //Currently reliable priors give information about true neighbors in the range [0.6,1] and false neighbors in the range [0,0.5]
                    if (i == j) {
                        dissim[row][col] = -1;
                    } else if (i > j) {
                        dissim[row][col] = dissim[col][row];
                    } else {
                        if (g.getEdge(curr, curr2) != null) {
                            if (prob[index] < amount) {
                                double val = trueDist.sample();

                                dissim[row][col] = val;
                                reliability += val;
                                counts++;

                            } else {
                                dissim[row][col] = -1;
                                dissim[col][row] = -1;
                            }
                        } else {
                            if (prob[index] < amount) {
                                double val = falseDist.sample();

                                dissim[row][col] = val;
                                reliability -= val;
                                counts++;
                            } else {
                                dissim[row][col] = -1;
                            }
                        }
                        index++;
                    }
                    col++;
                }
                row++;
            }
            reliability /= counts;


            PrintStream out = new PrintStream(prior);
            DecimalFormat df = new DecimalFormat("#.#####");
            for (int i = 0; i < dissim.length; i++) {
                for (int j = 0; j < dissim[i].length; j++) {
                    if (dissim[i][j] == -1) {
                        if (j == dissim[i].length - 1)
                            out.print("-1");
                        else
                            out.print("-1" + "\t");
                    } else {
                        if (j == dissim[i].length - 1)
                            out.print(df.format(dissim[i][j]));
                        else
                            out.print(df.format(dissim[i][j]) + "\t");
                    }
                }
                out.println();
            }

            out.flush();
            out.close();
            return reliability;

        } catch (Exception e) {
            System.err.println("Couldn't create dissimilarity file number " + k);
            e.printStackTrace();
            System.exit(-1);
        }

        return -1;
    }


    /***
     *
     * @param prior The file to write the intensity information to
     * @param intense The file to write the reliability information to
     * @param g The true data generating graph
     * @param target The target variable name
     * @param numPriors The total number of priors giving intensity information
     * @param numReliable The total number of reliable prior information sources
     * @param amountPrior The amount of prior information each source is giving (0 to 1)
     * @param amountRandom Whether or not the amount of prior information should be randomly generated
     * @param noiseRandom Whether or not the noise for each prior is randomly determined
     */
    public static void generateIntensity(File prior, File intense, Graph g, String target, int numPriors, int numReliable, double amountPrior, boolean amountRandom, boolean noiseRandom) {
        List<Node> neighbors = g.getAdjacentNodes(g.getNode(target));
        List<Node> allNodes = g.getNodes();
        double[] reliability = new double[numPriors];
        int[] counts = new int[numPriors];
        double[] amount = new double[numPriors];
        for (int i = 0; i < numPriors; i++) {
            amount[i] = amountPrior;
            if (amountRandom)
                amount[i] = rand.nextDouble() * amountLimit;
        }

        try {
            PrintStream out = new PrintStream(prior);
            out.print("Gene\t");
            for (int i = 0; i < numPriors; i++) {
                out.print("Prior_" + i + "\t");
            }
            out.println();
            PrintStream out2 = new PrintStream(intense);

            for (int i = 0; i < allNodes.size(); i++) {
                Node curr = allNodes.get(i);
                if (curr.getName().equals(target))
                    continue;
                out.print(allNodes.get(i).getName() + "\t");

                /***Probability for all priors to include information about this relationship***/
                double prob = rand.nextDouble();
                for (int j = 0; j < numPriors; j++) {

                    if (noiseRandom) {
                        reliableParams[0] = rand.nextDouble() * 0.8 + 0.2;
                        reliableParams[3] = rand.nextDouble() * 0.8;
                    }
                    /**Currently reliable priors give information about true neighbors in the range [0.6,1] and false neighbors in the range [0,0.5]**/
                    if (j < numReliable) {
                        if (neighbors.contains(curr)) {
                            if (prob < amount[j]) {
                                double val = rand.nextDouble() * (reliableParams[1] - reliableParams[0]) + reliableParams[0];
                                out.print(val + "\t");
                                reliability[j] += val;
                                counts[j]++;
                            } else {
                                out.print("-1\t");
                            }
                        } else {
                            if (prob < amount[j]) {
                                double val = rand.nextDouble() * (reliableParams[3] - reliableParams[2]) + reliableParams[2];
                                out.print(val + "\t");
                                reliability[j] -= val;
                                counts[j]++;
                            } else {
                                out.print("-1\t");
                            }
                        }
                    }
                    /***This prior is unreliable, needs to give bad information***/
                    else {
                        if (prob < amount[j]) {
                            double val = -1;
                            if (neighbors.contains(curr)) {
                                val = rand.nextDouble() * (unreliableParams[1] - unreliableParams[0]) + unreliableParams[0];
                                reliability[j] += val;
                            } else {
                                val = rand.nextDouble() * (unreliableParams[3] - unreliableParams[2]) + unreliableParams[2];
                                reliability[j] -= val;
                            }
                            out.print(val + "\t");

                            counts[j]++;
                        } else {
                            out.print("-1\t");
                        }
                    }

                }
                out.println();
            }
            for (int i = 0; i < numPriors; i++) {
                reliability[i] /= counts[i];
                out2.println(reliability[i]);
            }
            out2.flush();
            out2.close();
            out.flush();
            out.close();

        } catch (Exception e) {
            System.err.println("Couldn't create intensity file");
            System.exit(-1);
        }

    }


    /**
     * @param N Number of entries in the array
     * @return A float [] with random values between 0 and 1 as elements
     */
    public static float[] generateDisAmounts(int N) {
        float[] temp = new float[N];
        for (int i = 0; i < temp.length; i++) {
            temp[i] = rand.nextFloat();
        }
        return temp;
    }
}