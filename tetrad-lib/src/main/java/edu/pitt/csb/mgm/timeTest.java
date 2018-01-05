package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.search.DagToPag;
import edu.cmu.tetrad.search.Fci;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

/**
 * Created by vinee_000 on 10/24/2016.
 */
public class timeTest {

    public static void main(String[] args) throws Exception {
        int numGraphs = 5;
        boolean inds = false;
        //String path = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/Sedgewick_CausalMGM_data/Sedgewick_CausalMGM_data/simulated_data/";
        //String path2 = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/bigData";
        //String path3 = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/sampleSizeOutput";
        String path = "";

        long fci_time[][][] = new long[20][5][numGraphs];

        long mgm_time[][][] = new long[20][5][numGraphs];


        long max_time[][][] = new long[20][5][numGraphs];
        long mfm_time[][][] = new long[20][5][numGraphs];

        double lab = 0.25;
        double lab2 = 0.55;
        double[] alp = {.001, .005, .01, .03, .05};
        ArrayList[] toRemove = new ArrayList[numGraphs+3];
        PrintStream mgmTime = new PrintStream("mgm_time.txt");
        mgmTime.println("Graph\tTime Lambda1\tTime Lambda2");
        for (int j = 0; j < 5; j++) {
            PrintStream numTests = new PrintStream("tests_" + j + "_inds.txt");
            double alpha = alp[j];
            if(!inds)
                numTests = new PrintStream("tests_" + j + ".txt");
            numTests.println("Graph\tMGM-FCI (0.25)\tMGM-FCI (0.55)\tFCI\tMGM-FCI-MAX (0.25)\tMGM-FCI-MAX (0.55)");
            for (int fe = 3; fe < numGraphs+3; fe++) {

                System.out.println(fe);
                String filename = "SF_" + fe + "_data";
                String filename2 = "SF_" + fe + "_graph.txt";
                String extension = ".txt";
                Graph trueGraph = GraphUtils.loadGraphTxt(new File(filename2));
                DataSet data = MixedUtils.loadDataSet2(filename + extension);
                boolean done = false;
                double[] lambda = {lab, lab, lab};
                double[] lambda2 = {lab2, lab2, lab2};
                DataSet removeStuff = data.copy();
                /*while (!done) {
                    removeStuff = data.copy();

                    removeStuff.permuteRows();
                    int[] removal = new int[950];
                    for (int las = 0; las < 950; las++)
                        removal[las] = las;
                    removeStuff.removeRows(removal);
                    try {
                        MGM m = new MGM(removeStuff, lambda);
                    } catch (Exception e) {
                        continue;
                    }
                    done = true;
                }*/
                int nn = 10;
                int numNodesToRemove = nn;


                int hubLimit = trueGraph.getNumNodes() / 15;
                int totalNodes = trueGraph.getNumNodes();
                int numNodesRemoved = 0;
                Random rand = new Random();
                int count = 0;
                ArrayList<Integer> usedInts = new ArrayList<Integer>();
                ArrayList<Integer> removeInts = new ArrayList<Integer>();
                boolean reset = false;
                A:
                while (numNodesToRemove > 0) {
                    count++;
                    if (count == 800) {
                        count = 0;
                        if (hubLimit > 1) {
                            hubLimit = hubLimit - 1;
                            usedInts.clear();
                            usedInts = (ArrayList<Integer>) removeInts.clone();
                            continue;
                        } else {
                            System.out.println("Stuck on file " + fe + "_" + j);
                            reset = true;
                            break A;
                        }
                    }

                    int nod = rand.nextInt(totalNodes);
                    if (usedInts.contains(nod))
                        continue;
                    else
                        usedInts.add(nod);
                    if (trueGraph.getOutdegree(trueGraph.getNode("X" + (nod + 1))) > hubLimit) {
                        //      System.out.println("Removing Node " + (nod + 1));
                        // trueGraph.getNode("X" + (nod + 1)).setNodeType(NodeType.LATENT);
                        numNodesToRemove--;
                        numNodesRemoved++;
                        //data.removeColumn(nod);
                        removeInts.add(nod);
                        totalNodes--; //not sure if i need this

                    }

                }
                if (reset)
                    System.exit(0);
                //  continue;
                Collections.sort(removeInts);
                Collections.reverse(removeInts);
                if (toRemove[fe] == null)
                    toRemove[fe] = removeInts;
                System.out.println("Removing Ints: " + toRemove[fe] + " for graph " + fe);
                ArrayList<Node> latents = new ArrayList<Node>();
                //System.out.println(removeInts);
                removeInts = toRemove[fe];
                for (int curr = 0; curr < removeInts.size(); curr++) {
                    latents.add(trueGraph.getNode("X" + (removeInts.get(curr) + 1)));
                    trueGraph.getNode("X" + (removeInts.get(curr) + 1)).setNodeType(NodeType.LATENT);
                    //removeStuff.removeColumn(removeInts.get(curr));
                    removeStuff.removeColumn(removeStuff.getVariable("X" + (removeInts.get(curr) + 1)));
                }
                usedInts.clear();
                // removeInts.clear();
                final DagToPag dagToPag = new DagToPag(trueGraph);
                dagToPag.setCompleteRuleSetUsed(false);
                Graph p = dagToPag.convert();
                double tolerance = 1e-7; //convergeance tolerance
                int iterLimit = 1000; //iteration limit
                MGM model = new MGM(removeStuff, lambda);
                MGM model2 = new MGM(removeStuff, lambda2);
                long time1 = System.nanoTime();
                model.learnEdges(iterLimit);
                time1 = System.nanoTime()-time1;
                long time2 = System.nanoTime();
                model2.learnEdges(iterLimit);
                time2 = System.nanoTime()-time2;
                if(j==0) {
                    mgmTime.println(fe + "\t" + time1/Math.pow(10,9) + "\t" + time2/Math.pow(10,9));
                    mgmTime.flush();
                }
                    Graph mgmGraph = model.graphFromMGM();
                Graph mgmGraph2 = model2.graphFromMGM();

                IndTestMultinomialAJ indTest = new IndTestMultinomialAJ(removeStuff, alpha);
                Fci f = new Fci(indTest);
                f.setInitialGraph(mgmGraph);
                long fciTime = System.nanoTime();
                Graph gFinal = f.search();
                fciTime = System.nanoTime()-fciTime;
                if(inds)
                    numTests.print(fe + "\t" + indTest.reset() + "'t");
                else
                    numTests.print(fe + "\t" + (fciTime+time1)/Math.pow(10,9) + "\t");

                f = new Fci(indTest);
                f.setInitialGraph(mgmGraph2);
                fciTime = System.nanoTime();
                f.search();
                fciTime = System.nanoTime()-fciTime;
                if(inds)
                    numTests.print(indTest.reset() + "\t");
                else
                    numTests.print((fciTime+time2)/Math.pow(10,9) + "\t");
                Fci f2 = new Fci(indTest);


                FciMaxP f3 = new FciMaxP(indTest);


                Graph gFinal2 = gFinal;
                fciTime = System.nanoTime();
                gFinal2 = f2.search();
                fciTime = System.nanoTime()-fciTime;
                if(inds)
                    numTests.print(indTest.reset()+"\t");
                else
                    numTests.print(fciTime/Math.pow(10,9) + "\t");


                Graph gFinal3 = gFinal;

                //gFinal3 = f3.search();
                //Need to get statistics of this run here
                // PrintWriter q = new PrintWriter(new BufferedWriter(new FileWriter(path + filename + "_output" + extension)));
                //out.println("FCI: " + fciTime + "\n" + gFinal2);

                //System.out.println(p);
                //out.println("FCI-Max: " + maxTime + "\n" + gFinal3);
                // PerformanceTests.graphComparison(gFinal3, p, out);
                f3 = new FciMaxP(indTest);
                f3.setInitialGraph(mgmGraph);
                long maxTime = System.nanoTime();
                Graph gFinal4 = f3.search();
                maxTime = System.nanoTime()-maxTime;
                if(inds)
                    numTests.print(indTest.reset() + "\t");
                else
                    numTests.print((maxTime+time1)/Math.pow(10,9) + "\t");
                f3 = new FciMaxP(indTest);
                f3.setInitialGraph(mgmGraph2);
                maxTime = System.nanoTime();
                f3.search();
                maxTime = System.nanoTime()-maxTime;
                if(inds)
                    numTests.print(indTest.reset()+"\t");
                else
                    numTests.println((maxTime+time2)/Math.pow(10,9));
                numTests.flush();
            }
            numTests.flush();
            numTests.close();
            mgmTime.flush();
            mgmTime.close();
        }

    }
}
