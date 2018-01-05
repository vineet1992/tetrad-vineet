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
 * Created by vinee_000 on 12/10/2016.
 */
public class mfmTime {

    public static void main(String[] args) throws Exception {
        int numGraphs = 10;
        //String path = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/Sedgewick_CausalMGM_data/Sedgewick_CausalMGM_data/simulated_data/";
        //String path2 = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/bigData";
        //String path3 = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/sampleSizeOutput";
        String path = "";
        PrintStream output = new PrintStream("time_alpha_fci.txt");
        output.println("Algorithm\tGraph\tLambda\tAlpha\tTime\tTests");
        double[] alp = {.001, .01, .03, .05};
        ArrayList[] toRemove = new ArrayList[numGraphs];
        for (int j = 3; j < 4; j++) {
            double alpha = alp[j];
            for (int fe = 7; fe < numGraphs; fe++) {


                String filename = "SF_" + fe + "_data_p";
                String filename2 = "SF_" + fe + "_graph_p.txt";
                String extension = ".txt";
                Graph trueGraph = GraphUtils.loadGraphTxt(new File(filename2));
                DataSet data = MixedUtils.loadDataSet2(filename + extension);
                double[] lambda = {.35, .35, .35};
                double[] lambda2 = {.55, .55, .55};
                DataSet removeStuff = data.copy();
               /*
                while(!done) {
                    removeStuff = data.copy();

                    removeStuff.permuteRows();
                    int[] removal = new int[950];
                    for (int las = 0; las < 950; las++)
                        removal[las] = las;
                    removeStuff.removeRows(removal);
                    try{
                        MGM m = new MGM(removeStuff,lambda);
                    }
                    catch(Exception e)
                    {
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
                //trueGraph = GraphUtils.loadGraphTxt(new File(filename2));
                //System.out.println("True Graph " + p);
                // System.out.println("Data: " + data);
                //Need to ask about this when you don't remove any variables
                //  double[] lambda = {lab, lab, lab};
                //double tolerance = 1e-7; //convergeance tolerance
                int iterLimit = 1000; //iteration limit
                MGM model = new MGM(removeStuff, lambda);
                long mgmTime = 0;
                mgmTime = System.nanoTime();
                model.learnEdges(iterLimit);
                mgmTime = System.nanoTime() - mgmTime;
                MGM model2 = new MGM(removeStuff, lambda2);
                long mgmTime2 = System.nanoTime();
                model2.learnEdges(iterLimit);
                mgmTime2 = System.nanoTime() - mgmTime2;
                Graph mgmGraph = model.graphFromMGM();
                Graph mgmGraph2 = model2.graphFromMGM();

                IndTestMultinomialAJ indTest = new IndTestMultinomialAJ(removeStuff, alpha);


                FciMaxP f3 = new FciMaxP(indTest);
                f3.setInitialGraph(mgmGraph);
                long mgmmax = System.nanoTime();
                //Graph gFinal4 = f3.search();
                mgmmax = System.nanoTime() - mgmmax;
                int mfmTests = indTest.timesCalled;
                indTest = new IndTestMultinomialAJ(removeStuff,alpha);
                f3 = new FciMaxP(indTest);
                f3.setInitialGraph(mgmGraph2);
                long mgmmax2 = System.nanoTime();
                //f3.search();
                mgmmax2 = System.nanoTime() - mgmmax2;
                int mfmTests2 = indTest.timesCalled;
                indTest = new IndTestMultinomialAJ(removeStuff,alpha);
                Fci f = new Fci(indTest);
                long fciTime = System.nanoTime();
                f.search();
                fciTime = System.nanoTime()-fciTime;
                int fciTests = indTest.timesCalled;
                indTest = new IndTestMultinomialAJ(removeStuff,alpha);
                f = new Fci(indTest);
                f.setInitialGraph(mgmGraph);
                long mfciTime = System.nanoTime();
                //f.search();
                mfciTime = System.nanoTime()-mfciTime;
                int mgmfciTests = indTest.timesCalled;
                indTest = new IndTestMultinomialAJ(removeStuff,alpha);
                f = new Fci(indTest);
                f.setInitialGraph(mgmGraph2);
                long mfciTime2 = System.nanoTime();
                //f.search();
                mfciTime2 = System.nanoTime()-mfciTime2;
                int mgmfciTests2 = indTest.timesCalled;
              output.println("FCI\t" + fe + "\t" + lambda[0] + "\t" + alpha + "\t" + fciTime + "\t" + fciTests);
                 output.println("MGM-FCI\t" + fe + "\t" + lambda[0] + "\t" + alpha + "\t" + (mgmTime+mfciTime) + "\t" + mgmfciTests);
                output.println("MGM-FCI\t" + fe + "\t" + lambda2[0] + "\t" + alpha + "\t" + (mgmTime2 + mfciTime2) + "\t" + mgmfciTests2);
                output.println("MFM\t" + fe + "\t" + lambda[0] + "\t" + alpha + "\t" + (mgmTime + mgmmax) + "\t" + mfmTests);
                output.println("MFM\t" + fe + "\t" + lambda2[0] + "\t" + alpha + "\t" + (mgmTime2+mgmmax2) + "\t" + mfmTests2);

                System.out.println("Just finished " + alpha + " on graph " + fe);
                output.flush();
            }
        }
output.close();
    }
}
