package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

/**
 * Created by vinee_000 on 9/14/2016.
 */
public class testParallelMax {

    public static void main(String [] args) throws Exception{
        int parallelism = Integer.parseInt(args[0]);
        int start = Integer.parseInt(args[1]);
        int numGraphs = Integer.parseInt(args[2]);
        PrintStream out = new PrintStream("time_mfm_" + parallelism + "_" + start + "_" + (start+numGraphs-1) + ".txt");
        double lchoice = 0.35;
        double alpha = 0.05;
        double [] factors = {0.25,0.5,1,2,3};
        for(int j = 0; j < factors.length;j++) {
            double factor = factors[j];
            for (int fe = start; fe < numGraphs; fe++) {
                String filename = "SF_" + fe + "_data_p";
                String filename2 = "SF_" + fe + "_graph_p.txt";
                String extension = ".txt";
                Graph trueGraph = GraphUtils.loadGraphTxt(new File(filename2));
                DataSet data = MixedUtils.loadDataSet2(filename + extension);

                DataSet removeStuff = data.copy();
                int nn = 10;
                if(trueGraph.getNodes().size() < 100)
                    nn = 4;

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
                ArrayList<Node> latents = new ArrayList<Node>();
                //System.out.println(removeInts);
                for (int curr = 0; curr < removeInts.size(); curr++) {
                    latents.add(trueGraph.getNode("X" + (removeInts.get(curr) + 1)));
                    trueGraph.getNode("X" + (removeInts.get(curr) + 1)).setNodeType(NodeType.LATENT);
                    //removeStuff.removeColumn(removeInts.get(curr));
                    removeStuff.removeColumn(removeStuff.getVariable("X" + (removeInts.get(curr) + 1)));
                }
                usedInts.clear();

                double[] lambda = {lchoice, lchoice, lchoice};
                double tolerance = 1e-7; //convergeance tolerance
                int iterLimit = 1000; //iteration limit
                long mgmTime = System.nanoTime();
                MGM m = new MGM(removeStuff, lambda);
                m.learnEdges(iterLimit);
                mgmTime = System.nanoTime() - mgmTime;
                IndependenceTest i = new IndTestMultinomialAJ(removeStuff, alpha);
                long fciTime = System.nanoTime();
                FciMaxP f = new FciMaxP(i, factor,parallelism);
                f.setInitialGraph(m.graphFromMGM());
                f.search();
                fciTime = System.nanoTime() - fciTime;


                out.println("Graph\tFactor\tMGM\tFCI-MAX\tTotal");
                out.println(fe + "\t" + factor + "\t" + mgmTime + "\t" + fciTime + "\t" + (mgmTime + fciTime));

                out.flush();
            }
        }
        out.close();

    }
}
