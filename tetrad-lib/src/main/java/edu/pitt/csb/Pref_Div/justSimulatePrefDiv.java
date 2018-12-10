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
import java.text.DecimalFormat;
import java.util.*;

/**
 * Created by vinee_000 on 3/2/2018.
 */
public class justSimulatePrefDiv {


    public static double [] reliableParams = new double[]{0.7,1.0,0.0,0.3}; //True Low, True High, False Low, False High
    public static double [] unreliableParams = new double[]{0.0,1.0,0.0,1.0}; //True Low, True High, False Low, False High
    public static double amountLimit = 0.5; //Specifies the total percentage of prior information that can be given by a single source
    public static Random rand = new Random();

    public static void main(String [] args) {
        int numRuns = 20;
        int numGenes = 300;
        int sampleSize = 1000;
        double amountPrior = 0.3;//percentage of edges to have prior knowledge for
        boolean boot = false; //Should we use bootstrap samples for PiPrefDiv
        boolean loocv = false; //Should we use leave-one-out CV for PiPrefDiv
        int numSamples =  20; //Number of bootstrap/sub-sampled samples to use



        int numPriors = 10; //Number of prior knowledge sources
        int numReliable = 5; //Number of reliable sources
        int numComponents = 20; //How many components do we have for cluster simulation?
        int minTargetParents = 10; //How many true parents of the target are there?
        boolean noiseRandom = true; //Should the priors have random amounts of reliability?


        boolean amountRandom = false; //Should the priors have a random amount of prior knowledge?
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
            HashMap<String, Integer> clusters = new HashMap<String, Integer>();

            Graph g = null;
            DataSet d = null;
            File graphFile = null;


            /**GENERATE GRAPH**/

            System.out.print("Generating Graph...");

            graphFile = new File("Graphs/Graph_" + numGenes + "_" + minTargetParents + "_" + numComponents + "_" + evenDistribution + "_" + j + ".txt");
            if (graphFile.exists())
                g = GraphUtils.loadGraphTxt(graphFile);
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
            Parameters p = new Parameters();

            p.setValue("numMeasures", numGenes + 1);
            p.setValue("sampleSize", sampleSize);
            if (targetContinuous)
                p.setValue("percentDiscreteForMixedSimulation", 0);
            else
                p.setValue("percentDiscreteForMixedSimulation", 100 / (double) numGenes);
            p.setValue("numCategories", numCategories);
            NormalDistribution n = new NormalDistribution(numGenes*2 , numGenes / 3);
            int edges = (int)n.sample();
            int nodesPerCluster = numGenes/numComponents;
            int limit = nodesPerCluster*(nodesPerCluster-1)/2;
            while(edges> limit)
                edges = (int)n.sample();
            System.out.println(edges);
            p.setValue("numEdges", edges);
            p.setValue("numConnectedComponents", minTargetParents);
            p.setValue("numComponents", numComponents);
            if (targetContinuous)
                p.setValue("targetContinuous", 1);
            else
                p.setValue("targetContinuous", 0);

            if (g != null) {
                m.setTrueGraph(g);
            }
            clusters = m.simulateClustered(p, evenDistribution);

            System.out.println("Done");
            System.out.println("Graph: " + m.getTrueGraph());

            System.out.println("Data: " + m.getDataSet(0));

            System.out.println("Clusters: " + clusters);

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

            }


            String target = "Target";
            if (d == null) {
                d = m.getDataSet(0);
            } else
                m.setDataSet(d, 0);
            g = m.getTrueGraph();


            System.out.print("Generating Intensity File...");
            //Consistent ordering is whatever is returned by g.getNodes() (ignoring the target variable)
            File priorFile = null;
            priorFile = new File("Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_"  + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" +  j + "_intensity.txt");
            File iFile = new File("Reliabilities_I_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");
            if (!priorFile.exists()) {
                generateIntensity(priorFile,iFile, g, target, numPriors, numReliable, amountPrior,amountRandom,noiseRandom);
            }
            System.out.println("Done");
            //Use this if prior file exists Make sure intensity is at the end
            String[] disFiles = new String[numPriors];
            try {
                File temp = new File("Reliabilities_D_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");
                if(!temp.exists()) {
                    PrintStream out2 = new PrintStream("Reliabilities_D_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + j + ".txt");
                    for (int k = 0; k < numPriors; k++) {

                        System.out.print("Generating Dissimilarity File #" + k + "...");
                        disFiles[k] = "Priors/Prior_" + numGenes + "_" + minTargetParents + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_" + targetContinuous + "_" + numComponents + "_" + evenDistribution + "_" + amountRandom + "_" + noiseRandom + "_" + k + "_" + j + "_dissimilarity.txt";
                        File f = new File(disFiles[k]);
                        if (!f.exists()) {
                            double reli = generateDissimilarity(f, g, target, numReliable, amountPrior, amountRandom, noiseRandom, k);
                            out2.println(reli);
                            out2.flush();
                            System.out.println("Done");

                        }

                        //End Generation of Prior Knowledge
                    }

                    out2.close();

                }

            }
            catch(Exception e) {
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
                }
                else
                {
                    subs = PiPrefDiv.genSubsamples(boot,numSamples,d,loocv);
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
            }
            else
            {
                int trainLength = (int)(d.getNumRows()*0.9);
                int [] trainInds = new int[trainLength];
                for(int i = 0; i < trainLength;i++)
                {
                    trainInds[i] = i;
                }
                DataSet toRun = d.subsetRows(trainInds);
                subs2 = PiPrefDiv.genSubsamples(boot,numSamples,toRun,loocv);
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
                        temp.println(s + "\t" + clusters.get(s));
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



    public static double generateDissimilarity(File prior, Graph g, String target, int numReliable,double amountPrior, boolean amountRandom, boolean noiseRandom, int k)
    {
        List<Node> allNodes = g.getNodes();
        double[][] dissim = new double[allNodes.size()-1][allNodes.size()-1];
        int row = 0;
        double reliability = 0;
        int counts = 0;
        double amount = amountPrior;
                if(amountRandom)
                {
                    amount = rand.nextDouble()*amountLimit;
                }
        if(noiseRandom)
        {
            reliableParams[0] = rand.nextDouble()*0.8+0.2;
            reliableParams[3] = rand.nextDouble()*0.8;
        }
        try {

            for(int i = 0; i < allNodes.size();i++)
            {
                Node curr = allNodes.get(i);
                if(curr.getName().equals(target))
                    continue;
                int col = 0;
                for(int j = 0; j < allNodes.size();j++) {
                    Node curr2 = allNodes.get(j);
                    if (curr2.getName().equals(target))
                        continue;
                    //Currently reliable priors give information about true neighbors in the range [0.6,1] and false neighbors in the range [0,0.5]
                    if (i == j) {
                        dissim[row][col] = -1;
                    } else if (i > j)
                    {
                        dissim[row][col] = dissim[col][row];
                    }
                    else if(k < numReliable)
                    {
                        if(g.getEdge(curr,curr2)!=null)
                        {
                            if(rand.nextDouble()<amount)
                            {
                                double val = rand.nextDouble()*(reliableParams[1]-reliableParams[0])+reliableParams[0];
                                dissim[row][col] = val;
                                reliability+=val;
                                counts++;
                            }
                            else
                            {
                                dissim[row][col] = -1;
                                dissim[col][row] = -1;
                            }
                        }
                        else
                        {
                            if(rand.nextDouble()<amount)
                            {
                                double val = rand.nextDouble()*(reliableParams[3]-reliableParams[2]) + reliableParams[2];
                                dissim[row][col] = val;
                                reliability-= val;
                                counts++;
                            }
                            else
                            {
                                dissim[row][col] = -1;
                            }
                        }
                    }
                    //This prior is unreliable, needs to give bad information
                    else
                    {
                        if(rand.nextDouble()<amount)
                        {
                            double val = -1;
                            if(g.getEdge(curr,curr2)!=null)
                            {
                                val = rand.nextDouble()*(unreliableParams[1]-unreliableParams[0]) + unreliableParams[0];
                                reliability+=val;

                            }
                            else
                            {
                                val = rand.nextDouble()*(unreliableParams[3]-unreliableParams[2]) + unreliableParams[2];
                                reliability-=val;


                            }
                            dissim[row][col] = val;
                            counts++;
                        }
                        else
                        {
                            dissim[row][col] = -1;
                        }
                    }
                    col++;
                }
                row++;
            }
            reliability /= counts;


            PrintStream out = new PrintStream(prior);
            DecimalFormat df = new DecimalFormat("#.#####");
            for(int i = 0; i < dissim.length;i++)
            {
                for(int j = 0; j < dissim[i].length;j++)
                {
                    if(dissim[i][j]==-1)
                    {
                        if (j == dissim[i].length - 1)
                            out.print("-1");
                        else
                            out.print("-1" + "\t");
                    }
                    else {
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

        }
        catch (Exception e)
        {
            System.err.println("Couldn't create dissimilarity file number " + k);
            e.printStackTrace();
            System.exit(-1);
        }

        return -1;
    }

    public static void generateIntensity(File prior,File intense, Graph g, String target, int numPriors, int numReliable,double amountPrior, boolean amountRandom,boolean noiseRandom)
    {
        List<Node> neighbors = g.getAdjacentNodes(g.getNode(target));
        List<Node> allNodes = g.getNodes();
        double [] reliability = new double[numPriors];
        int [] counts = new int[numPriors];
        double [] amount = new double[numPriors];
        for(int i = 0; i < numPriors;i++)
        {
            amount[i] = amountPrior;
            if(amountRandom)
                amount[i] = rand.nextDouble()*amountLimit;
        }

        try {
            PrintStream out = new PrintStream(prior);
            out.print("Gene\t");
            for(int i = 0; i < numPriors;i++)
            {
                out.print("Prior_" + i + "\t");
            }
            out.println();
            PrintStream out2 = new PrintStream(intense);

            for(int i = 0; i < allNodes.size();i++)
            {
                Node curr = allNodes.get(i);
                if(curr.getName().equals(target))
                    continue;
                out.print(allNodes.get(i).getName()+"\t");
               for(int j = 0; j < numPriors;j++)
               {

                   if(noiseRandom)
                   {
                       reliableParams[0] = rand.nextDouble()*0.8+0.2;
                       reliableParams[3] = rand.nextDouble()*0.8;
                   }
                   //Currently reliable priors give information about true neighbors in the range [0.6,1] and false neighbors in the range [0,0.5]
                   if(j < numReliable)
                   {
                       if(neighbors.contains(curr))
                       {
                           if(rand.nextDouble()<amount[j])
                           {
                               double val = rand.nextDouble()*(reliableParams[1]-reliableParams[0])+reliableParams[0];
                               out.print(val + "\t");
                               reliability[j]+=val;
                               counts[j]++;
                           }
                           else
                           {
                               out.print("-1\t");
                           }
                       }
                       else
                       {
                           if(rand.nextDouble()<amount[j])
                           {
                               double val = rand.nextDouble()*(reliableParams[3]-reliableParams[2]) + reliableParams[2];
                               out.print(val + "\t");
                               reliability[j]-= val;
                               counts[j]++;
                           }
                           else
                           {
                               out.print("-1\t");
                           }
                       }
                   }
                   //This prior is unreliable, needs to give bad information
                   else
                   {
                       if(rand.nextDouble()<amount[j])
                       {
                           double val = -1;
                           if(neighbors.contains(curr))
                           {
                               val = rand.nextDouble()*(unreliableParams[1]-unreliableParams[0]) + unreliableParams[0];
                               reliability[j]+=val;
                           }
                           else
                           {
                               val = rand.nextDouble()*(unreliableParams[3]-unreliableParams[2]) + unreliableParams[2];
                               reliability[j]-=val;
                           }
                           out.print(val + "\t");

                           counts[j]++;
                       }
                       else
                       {
                           out.print("-1\t");
                       }
                   }

               }
                out.println();
            }
            for(int i = 0; i < numPriors;i++)
            {
                reliability[i] /= counts[i];
                out2.println(reliability[i]);
            }
            out2.flush();
            out2.close();
            out.flush();
            out.close();

        }
        catch (Exception e)
        {
            System.err.println("Couldn't create intensity file");
            System.exit(-1);
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
