package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.calculator.expression.Expression;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.EdgeListGraphSingleConnections;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.sem.GeneralizedSemIm;
import edu.cmu.tetrad.sem.GeneralizedSemPm;
import edu.cmu.tetrad.util.StatUtils;
import edu.pitt.csb.mgm.MixedUtils;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by vinee_000 on 3/2/2018.
 */
public class simulatePrefDivData {

    public static void main(String [] args) {
        int numRuns = 20;
        int numGenes = 200;
        int sampleSize = 200;

        int ns = 20; //number of subsamples
        int nf = 5;//number of folds for cross-validation
        int minTargetParents = 25; //How many true parents of the target are there?
        //double[] noises = {-3,-2,-1,-0.5,-0.25,0};
        double[] noises = {0.001,0.1,0.25,0.5,0.75,1.0};
        //int [] noises = {10,20,50,75,100};
         //   double [] noises = {1.0};
        for (int nn = 0; nn < noises.length; nn++) {
            double intensityNoise = noises[nn];
            double dissimilarityNoise = noises[nn];
            int B = 50;
            boolean stabilitySelection = true; //Should we use the stability selection version or the StARS version
            boolean targetContinuous = true; //Is the target variable continuous?
            boolean clusterSimulation = false; //Are we simulating the data in clusters?
            boolean evenDistribution = false; //Is the distribution of nodes in each cluster even?
            boolean clusterStability = false; //Do we use the stability of the clusters to determine the choice of alpha?
            boolean useCausalGraph = false;//Do we use the causal graph to select regression features? (Only for CPSS version)
            int numComponents = 10; //How many components do we have for cluster simulation?
            double stabilityThreshold = 0.1; //How much instability do we tolerate to choose alpha?
            double percentMissing = 0;
            int numCategories = 4; //number of categories for discrete variables

            double radius = -1; //Radius of similarity for Pref-Div
            double accuracy = 0.5; //A parameter in Pref-Div
            int topK = minTargetParents; //How many genes are we choosing with our method?

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
                } else if (args[index].equals("-ns")) {
                    ns = Integer.parseInt(args[index + 1]);
                    index += 2;
                } else if (args[index].equals("-numParents")) {
                    minTargetParents = Integer.parseInt(args[index + 1]);
                    index += 2;
                } else if (args[index].equals("-nc")) {
                    numComponents = Integer.parseInt(args[index + 1]);
                    index += 2;
                }else if(args[index].equals("-stabSelect"))
                {
                    stabilitySelection = true;
                    index++;
                }else if (args[index].equals("-stabThresh")) {
                    stabilityThreshold = Double.parseDouble(args[index + 1]);
                    index += 2;
                } else if (args[index].equals("-noise")) {
                    intensityNoise = Double.parseDouble(args[index + 1]);
                    dissimilarityNoise = intensityNoise;
                    index += 2;
                } else if (args[index].equals("-categorical")) {
                    targetContinuous = false;
                    index++;
                } else if (args[index].equals("-clusterSim")) {
                    clusterSimulation = true;
                    index++;
                } else if (args[index].equals("-even")) {
                    evenDistribution = true;
                    index++;
                } else if (args[index].equals("-clusterStab")) {
                    clusterStability = true;
                    index++;
                } else if(args[index].equals("-useCG")){
                        useCausalGraph = true;
                        index++;
                } else {
                    System.err.println("Unidentified argument: " + args[index]);
                    System.exit(-1);
                }
            }


            String runName = "";
            if (clusterSimulation)
                runName = numGenes + "_" + sampleSize + "_" + ns + "_" + minTargetParents + "_" + intensityNoise + "_" + targetContinuous + "_" + stabilityThreshold + "_" + numComponents + "_" + radius + "_" + accuracy + "_" + B + "_" + evenDistribution + "_" + clusterStability + "_" + stabilitySelection + "_" + useCausalGraph; //Name to append at the beginning of all output files
            else
                runName = numGenes + "_" + sampleSize + "_" + ns + "_" + minTargetParents + "_" + intensityNoise + "_" + targetContinuous + "_" + stabilityThreshold + "_" + radius + "_" + accuracy + "_" + B + "_" + clusterStability; //Name to append at the beginning of all output files
            boolean saveData = true;

            //out variables refer to outputs
            boolean outStability = true;//Stability vs alpha data file for each run
            boolean outAccuracy = true;//Percentage of correct individual features for each run, this will be its own file
            boolean outStabAcc = true;//Relationship between stability and accuracy of individual features, basically just add a
            //correctness column to the stability file for the first boolean

            //Not sure about this boolean, need to think about it TODO
            boolean mapFeatureChange = false; //How do stable features change with increasing theory use?
            //Produce a file with # of times each feature appeared for each alpha valu


            //If either of outAccuracy or outStabAcc are true, then


            GeneralizedSemIm im = null;

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


            ArrayList<String> accuracyResults = new ArrayList<String>();
            ArrayList<String> stabilityResults = new ArrayList<String>();
            File resultAccFile = new File("Results/Result_" + runName + "_accuracy" + ".txt");
            if (resultAccFile.exists()) {
                try {
                    BufferedReader r = new BufferedReader(new FileReader(resultAccFile));
                    r.readLine();
                    while (r.ready()) {
                        accuracyResults.add(r.readLine());
                    }
                } catch (Exception e) {
                    System.err.println("Error reading the accuracy results from file");
                    e.printStackTrace();
                    System.exit(-1);
                }
            }
            PrintStream out = null;
            PrintStream stabPS = null;
            try {
                out = new PrintStream("Results/Result_" + runName + "_accuracy" + ".txt");
                if (clusterSimulation)
                    out.println("Run\tPredicted_Features\tActual_Features\tPrecision\tRecall\tCluster_Precision\tCluster_Recall\tCluster_Accuracy");
                else
                    out.println("Run\tPredicted_Features\tActual_Features\tPrecision\tRecall");

            } catch (Exception e) {
                System.err.println("Unable to create accuracy file");
                System.exit(-1);
            }
            File resultStabFile = new File("Results/Result_" + runName + "_stabilities.txt");
            try {
                if (resultStabFile.exists()) {
                    BufferedReader r = new BufferedReader(new FileReader(resultStabFile));
                    r.readLine();
                    while (r.ready()) {
                        stabilityResults.add(r.readLine());
                    }
                    r.close();
                }
                stabPS = new PrintStream("Results/Result_" + runName + "_stabilities.txt");
                if(stabilitySelection)
                {
                    if (outStabAcc)
                        stabPS.println("Run\tAlpha\tCV Accuracy\tCV Stability\tStability\tPrecision\tRecall\tPredicted\tActual\tClusters");
                    else
                        stabPS.println("Run\tAlpha\tCV Accuracy\tCV Stability\tStability\tClusters");
                }
                else {
                    if (outStabAcc)
                        stabPS.println("Run\tAlpha\tStability\tPrecision\tRecall\tPredicted\tActual\tClusters");
                    else
                        stabPS.println("Run\tAlpha\tStability\tClusters");
                }
            } catch (Exception e) {
                System.err.println("Unable to create stability file");
                System.exit(-1);

            }
            for (int j = 0; j < numRuns; j++) {

                HashMap<String, Integer> clusters = new HashMap<String, Integer>();

                Graph g = null;
                DataSet d = null;
                File graphFile = null;


                if (clusterSimulation)
                    graphFile = new File("Graphs/Graph_" + numGenes + "_" + minTargetParents + "_" + numComponents + "_" + evenDistribution + "_" + j + ".txt");
                else
                    graphFile = new File("Graphs/Graph_" + numGenes + "_" + minTargetParents + "_" + j + ".txt");
                if (graphFile.exists())
                    g = GraphUtils.loadGraphTxt(graphFile);
                File dataFile = null;

                if (clusterSimulation)
                    dataFile = new File("Data/Dataset_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + j + ".txt");
                else
                    dataFile = new File("Data/Dataset_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + j + ".txt");
                try {
                    if (dataFile.exists())
                        d = MixedUtils.loadDataSet2(dataFile.getAbsolutePath());
                } catch (Exception e) {
                    System.err.println("Unable to load dataset");
                    System.exit(-1);
                }


                MixedLeeHastieSimulation m = new MixedLeeHastieSimulation();
                Parameters p = new Parameters();

                p.setValue("numMeasures", numGenes + 1);
                p.setValue("sampleSize", sampleSize);
                if (targetContinuous)
                    p.setValue("percentDiscreteForMixedSimulation", 0);
                else
                    p.setValue("percentDiscreteForMixedSimulation", 100 / (double) numGenes);
                p.setValue("numCategories", numCategories);
                NormalDistribution n = new NormalDistribution(numGenes * 2, numGenes / 3);
                p.setValue("numEdges", (int) n.sample());
                p.setValue("numConnectedComponents", minTargetParents);
                p.setValue("numComponents", numComponents);
                if (targetContinuous)
                    p.setValue("targetContinuous", 1);
                else
                    p.setValue("targetContinuous", 0);

                if (g != null) {
                    m.setTrueGraph(g);
                }

                if (clusterSimulation) {
                    clusters = m.simulateClustered(p, evenDistribution);
                    im = m.getIm();
                } else
                    m.simulate(p);


                File clusterFile = null;
                if (clusterSimulation)
                    clusterFile = new File("Clusters/Cluster_" + numGenes + "_" + minTargetParents + "_" + numComponents + "_" + evenDistribution + "_" + j + ".txt");
                if (clusterSimulation && clusterFile.exists()) {
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

                }


                String target = "Target";
                if (d == null) {
                    if (!clusterSimulation) {
                        A:
                        while (true) //Repeat until we have one target variable with a sufficient number of parents
                        {
                            g = m.getTrueGraph();
                            d = m.getDataSet(0);
                            for (int i = 0; i < d.getNumColumns(); i++) {
                                if (!targetContinuous && d.getVariable(i) instanceof DiscreteVariable) {
                                    if (g.getParents(g.getNode(d.getVariable(i).getName())).size() >= minTargetParents && g.getChildren(g.getNode(d.getVariable(i).getName())).size() == 0) {
                                        EdgeListGraphSingleConnections gTemp = (EdgeListGraphSingleConnections) g;
                                        gTemp.renameNode(d.getVariable(i).getName(), "Target");
                                        d.getVariable(i).setName("Target");
                                        target = d.getVariable(i).getName();
                                        break A;
                                    }
                                } else if (g.getParents(g.getNode(d.getVariable(i).getName())).size() >= minTargetParents && g.getChildren(g.getNode(d.getVariable(i).getName())).size() == 0) {
                                    EdgeListGraphSingleConnections gTemp = (EdgeListGraphSingleConnections) g;
                                    gTemp.renameNode(d.getVariable(i).getName(), "Target");
                                    d.getVariable(i).setName("Target");
                                    target = d.getVariable(i).getName();
                                    break A;
                                }
                            }
                            m = new MixedLeeHastieSimulation();
                            m.simulate(p);
                        }
                    }
                    d = m.getDataSet(0);
                } else
                    m.setDataSet(d, 0);
                NormalDistribution n_intense = new NormalDistribution(0, intensityNoise);
                NormalDistribution n_dissim = new NormalDistribution(0, dissimilarityNoise);
                g = m.getTrueGraph();
                double[][] data = d.getDoubleData().transpose().toArray();
                int targetColumn = d.getColumn(d.getVariable(target));
                ArrayList<Gene> genes = new ArrayList<Gene>();
                int count = 0;

                //Use this if prior file exists MAke sure intensity is at the end
                File priorFile = null;
                if (clusterSimulation)
                    priorFile = new File("Priors/Prior_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + intensityNoise + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + j + "_intensity.txt");
                else
                    priorFile = new File("Priors/Prior_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + intensityNoise + "_" + targetContinuous + "_" + j + "_intensity.txt");
                ArrayList<Double> values = new ArrayList<Double>();
                if (priorFile.exists()) {
                    try {
                        BufferedReader b = new BufferedReader(new FileReader(priorFile));
                        while (b.ready())
                            values.add(Double.parseDouble(b.readLine()));
                        b.close();
                    } catch (Exception e) {
                        System.err.println("Couldn't load prior file " + j);
                        e.printStackTrace();
                        System.exit(-1);
                    }
                }


                //Use this if we need to generate the prior file

                HashMap<String, Double> causalEffects = new HashMap<>();
                if (clusterSimulation) {
                    GeneralizedSemPm pm = im.getSemPm();
                    String exp = pm.getNodeExpressionString(g.getNode("Target"));
                    String[] stuff = exp.split("\\+");
                    for (int x = 0; x < stuff.length - 1; x++) {
                        String[] line = stuff[x].split("\\*");
                        causalEffects.put(line[1].trim(), Math.abs(im.getParameterValue(line[0].trim())));
                    }
                }
                //Plus some error
                count = 0;
                for (int i = 0; i < d.getNumColumns(); i++) {
                    if (!d.getVariable(i).getName().equals(target)) {
                        Gene g2 = new Gene(count);
                        g2.symbol = d.getVariable(i).getName();
                        if (targetContinuous) {
                            double corr = Math.abs(StatUtils.correlation(data[i], data[targetColumn]));
                            if (values.isEmpty() && clusterSimulation) {
                                if (causalEffects.get(g2.symbol) == null)
                                    g2.theoryIntensity = n_intense.sample();
                                else
                                    g2.theoryIntensity = causalEffects.get(g2.symbol) + n_intense.sample();
                            } else if (values.isEmpty())
                                g2.theoryIntensity = corr + n_intense.sample();
                            else
                                g2.theoryIntensity = values.get(count);
                        } else {
                            double corr = Functions.mixedMI(data[i], data[targetColumn], numCategories);
                            if (values.isEmpty())
                                g2.theoryIntensity = corr + n_intense.sample();
                            else
                                g2.theoryIntensity = values.get(count);
                        }
                        genes.add(g2);
                        count++;
                    }
                }

                float[] allCorrs = new float[genes.size()];
                for (int i = 0; i < genes.size(); i++) {
                    allCorrs[i] = (float) genes.get(i).theoryIntensity;
                }
                allCorrs = Functions.NPN(allCorrs, true);
                for (int i = 0; i < genes.size(); i++) {
                    genes.get(i).theoryIntensity = allCorrs[i];
                }
                System.out.println("Correlations: " + Arrays.toString(allCorrs));

                File priorFileDis = null;
                if (clusterSimulation)
                    priorFileDis = new File("Priors/Prior_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + intensityNoise + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + j + "_dissimilarity.txt");
                else
                    priorFileDis = new File("Priors/Prior_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + intensityNoise + "_" + targetContinuous + "_" + j + "_dissimilarity.txt");
                BufferedReader b = null;
                try {
                    if (priorFileDis.exists())
                        b = new BufferedReader(new FileReader(priorFileDis));
                } catch (Exception e) {
                    System.err.println("Unable to load prior dissimilarity file, but it exists " + j);
                    e.printStackTrace();
                    System.exit(-1);
                }

                causalEffects = new HashMap<>();
                if (clusterSimulation) {
                    GeneralizedSemPm pm = im.getSemPm();
                    List<Node> nodes = g.getNodes();
                    for (Node no : nodes) {
                        String exp = pm.getNodeExpressionString(no);
                        if (exp.contains("+")) {
                            String[] stuff = exp.split("\\+");
                            for (int x = 0; x < stuff.length - 1; x++) {
                                String[] line = stuff[x].split("\\*");
                                causalEffects.put(line[1].trim() + "," + no.getName(), Math.abs(im.getParameterValue(line[0].trim())));
                                causalEffects.put(no.getName() + "," + line[1].trim(), Math.abs(im.getParameterValue(line[0].trim())));
                            }
                        }
                    }
                }
                float[] dissimilarity = new float[genes.size() * (genes.size() - 1) / 2];
                index = 0;
                for (int i = 0; i < genes.size(); i++) {
                    double[] one = data[d.getColumn(d.getVariable(genes.get(i).symbol))];
                    for (int k = i + 1; k < genes.size(); k++) {
                        double[] two = data[d.getColumn(d.getVariable(genes.get(k).symbol))];
                        double corr = Math.abs(StatUtils.correlation(one, two) + n_dissim.sample());
                        if (b == null) {
                            if (clusterSimulation) {
                                if (causalEffects.get(genes.get(i).symbol + "," + genes.get(k).symbol) == null)
                                    dissimilarity[index] = (float) corr;
                                else
                                    dissimilarity[index] = (float) (causalEffects.get(genes.get(i).symbol + "," + genes.get(k).symbol).doubleValue() + n_dissim.sample());
                            } else
                                dissimilarity[index] = (float) corr;
                        } else {
                            try {
                                dissimilarity[index] = Float.parseFloat(b.readLine());
                            } catch (Exception e) {
                                System.err.println("Couldn't load theory dissimilarity from line: " + index);
                                System.exit(-1);
                            }
                        }
                        index++;
                    }
                }
                try {
                    if (b != null)
                        b.close();
                } catch (Exception e) {
                    System.err.println("Couldn't close reader for theory dissimilarity");
                    System.exit(-1);
                }

                //DISSIMILARITY SHOULD BE A MEASURE OF HOW NOT SIMILAR THEY ARE, so negate the values here to flip the normal distribution
                dissimilarity = Functions.NPN(dissimilarity, true);
                for (int i = 0; i < dissimilarity.length; i++) {
                    dissimilarity[i] = -1 * dissimilarity[i];
                }
                File subsFile = null;
                if (clusterSimulation)
                    subsFile = new File("Subsamples/Subsample_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + stabilitySelection + "_" + nf + "_" +  j +  ".txt");
                else
                    subsFile = new File("Subsamples/Subsample_" + numGenes + "_" + sampleSize + "_" + minTargetParents + "_" + targetContinuous + "_" + stabilitySelection + "_" + nf + "_" + j + ".txt");
                int[][] subs = null;
                if (subsFile.exists()) {
                    try {
                        b = new BufferedReader(new FileReader(subsFile));
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
                }

                RunPrefDiv r = new RunPrefDiv(dissimilarity, genes, d, target, false);

                if (subs != null)
                    r.setSubs(subs);
                r.useStabilitySelection(true);
                r.setClusterStability(clusterStability);
                r.setRadius(radius);
                r.setAccuracy(accuracy);
                r.setNS(ns);
                r.setTopK(topK);
                r.setNumFolds(nf);
                r.setThreshold(stabilityThreshold);
                r.setStabilityOutput(outStabAcc);
                r.setTargets(g.getAdjacentNodes(g.getNode(target)));
                r.setRun(j);
                r.usePartialCorrelation(false);
                r.setCausalGraph(useCausalGraph);
                r.setB(B);

                if (accuracyResults.size() <= j) {
                    ArrayList<Gene> result = r.runPD();
                    System.out.println("Top Genes: " + result);
                    System.out.println("Clusters: " + r.getClusters());
                    if (clusterSimulation)
                        accuracyResults.add(j + "\t" + result + "\t" + g.getAdjacentNodes(g.getNode(target)) + "\t" + getAccuracy(result, target, g)[1] + "\t" + getAccuracy(result, target, g)[2] + "\t" + simpleClusterAccuracy(result, target, g, clusters)[0] + "\t" + simpleClusterAccuracy(result, target, g, clusters)[1] + "\t" + averageClusterSim(result, target, g, clusters, r.getClusters()));
                    else
                        accuracyResults.add(j + "\t" + result + "\t" + g.getAdjacentNodes(g.getNode(target)) + "\t" + getAccuracy(result, target, g)[1] + "\t" + getAccuracy(result, target, g)[2]);

                    stabilityResults.addAll(r.getStabilityOutput());
                    if (subs == null)
                        subs = r.getSubs();
                }


                try {
                    PrintStream temp = new PrintStream(graphFile);
                    temp.println(g);
                    temp.flush();
                    temp.close();
                    if (clusterSimulation) {
                        temp = new PrintStream(clusterFile);
                        for (String s : clusters.keySet()) {
                            temp.println(s + "\t" + clusters.get(s));
                            temp.flush();
                        }
                        temp.close();
                    }
                    temp = new PrintStream(dataFile);
                    temp.println(d);
                    temp.flush();
                    temp.close();
                    temp = new PrintStream(priorFile);
                    count = 0;
                    for (int i = 0; i < d.getNumColumns(); i++) {
                        if (d.getVariable(i).getName().equals("Target"))
                            continue;
                        temp.println(allCorrs[count]);
                        count++;
                    }
                    temp.flush();
                    temp.close();
                    temp = new PrintStream(priorFileDis);
                    for (int i = 0; i < dissimilarity.length; i++) {
                        temp.println(dissimilarity[i]);
                    }
                    temp.flush();
                    temp.close();
                    temp = new PrintStream(subsFile);
                    for (int i = 0; i < subs.length; i++) {
                        for (int k = 0; k < subs[i].length; k++) {
                            if (k == subs[i].length - 1)
                                temp.println(subs[i][k]);
                            else
                                temp.print(subs[i][k] + "\t");
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

            for (int i = 0; i < accuracyResults.size(); i++) {
                out.println(accuracyResults.get(i));
            }
            out.flush();
            out.close();


            for (int i = 0; i < stabilityResults.size(); i++) {
                stabPS.println(stabilityResults.get(i));
            }
            stabPS.flush();
            stabPS.close();

            try {
                out.close();
            } catch (Exception e) {
                System.err.println("Error closing output accuracy file");
                System.exit(-1);
            }


            try {
                stabPS.close();
            } catch (Exception e) {
                System.err.println("Error closing output stability file");
                System.exit(-1);
            }

        }
    }



    public static double averageClusterSim(ArrayList<Gene> result, String target, Graph g, HashMap<String,Integer> clusters,HashMap<Gene,List<Gene>> predClusts)
    {
        //average tanimoto set similarity across clusters that were correctly identified (true connected components)
        //if same connected component was identified twice, merge these clusters
        List<Node> trueParents = g.getAdjacentNodes(g.getNode(target));

        double tss = 0;
        int count = 0;

        //Cluster # -> All genes predicted to be in that cluster
        HashMap<Integer,Set<String>> finalClusters = new HashMap<Integer,Set<String>>();
        for(Gene gene: predClusts.keySet())
        {
            int x = clusters.get(gene.symbol);
            List<Gene> temp = predClusts.get(gene);
            List<String> temp2 = new ArrayList<String>();
            temp2.add(gene.symbol);
            //if temp is null then this is an empty cluster, so just add 'gene' to the final set
            if(temp!=null) {
                for (Gene y : temp)
                    temp2.add(y.symbol);
            }
            if(finalClusters.get(x)==null)
                finalClusters.put(x,new HashSet(temp2));
            else
            {
                Set<String> set = finalClusters.get(x);
                set.addAll(temp2);
                finalClusters.put(x,set);
            }

        }
        HashMap<Integer,Set<String>> trueClusters = new HashMap<>();
        for(String s: clusters.keySet())
        {
            if(trueClusters.get(clusters.get(s))==null)
            {
                Set<String> temp = new HashSet<String>();
                temp.add(s);
                trueClusters.put(clusters.get(s),temp);
            }
            else
            {
                Set<String> temp = trueClusters.get(clusters.get(s));
                temp.add(s);
                trueClusters.put(clusters.get(s),temp);
            }
        }

        for(int i = 0; i < trueParents.size();i++)
        {
            if(finalClusters.get(i)!=null)
            {
                count++;
                tss+= (CollectionUtils.intersection(finalClusters.get(i),trueClusters.get(i)).size()/(double)CollectionUtils.union(finalClusters.get(i),trueClusters.get(i)).size());
            }
        }

        return tss/count;
    }
    //Returns Precision and recall for how good identified 10 cluster representatives are vs true 10 cluster representatives
    //You only get a single TP even if you pick multiple representatives from the same correct cluster, whereas you get multiple fp if you pick more than one
    //representative from he same wrong cluster
    //E.G. Selected Vars are from clusters (1,3,6) True vars are from clusters (8,13,6) Prec = 1/3, Rec = 1/3
    public static double [] simpleClusterAccuracy(ArrayList<Gene> result,String target, Graph g,HashMap<String,Integer> clusters)
    {
        List<Node> trueParents = g.getAdjacentNodes(g.getNode(target));
        double tp = 0;
        double fp = 0;
        double fn = 0;
        ArrayList<Integer> estimated = new ArrayList<Integer>();
        for(int i = 0; i < result.size();i++)
        {
            estimated.add(clusters.get(result.get(i).symbol));
        }
        HashSet<Integer> truth = new HashSet<Integer>();
        for(int i = 0; i < trueParents.size();i++)
        {
            truth.add(clusters.get(trueParents.get(i).getName()));
        }


        for(int i = 0; i < truth.size();i++)
        {
            if(truth.contains(i) && estimated.contains(i)) {
                tp++;
            }
            else
                fn++;
        }

        for(int i = 0; i < estimated.size();i++)
        {
            if(i >= truth.size())
                fp++;
        }
        double [] output = {tp/(tp+fp),tp/(tp+fn)};
        return output;

    }
    public static double[] getAccuracy(ArrayList<Gene> result,String target,Graph g)
    {
        double tp = 0;
        double fp = 0;
        double fn = 0;
        double [] output = new double[3]; //Accuracy, Precision,Recall
        List<Node> truth = g.getAdjacentNodes(g.getNode(target));
        for(int i = 0; i < result.size();i++)
        {
            boolean found = false;
            for(Node n: truth)
                if(result.get(i).symbol.equals(n.getName()))
                    found = true;
            if(found)
                tp++;
            else
                fp++;
        }
        for(int i = 0; i < truth.size();i++)
        {
            boolean found = false;
            for(int j = 0; j < result.size();j++)
            {
                if(result.get(j).symbol.equals(truth.get(i).getName()))
                    found = true;
            }
            if(!found)
                fn++;
        }
        output[1] = tp/(tp+fp);
        output[2] = tp/(tp+fn);
        output[0] = (output[1] + output[2])/2; //If the number of elements are equal, then precision and recall should be the same
        return output;
    }
}
