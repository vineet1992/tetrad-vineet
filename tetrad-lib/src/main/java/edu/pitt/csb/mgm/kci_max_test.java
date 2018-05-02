package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.search.DagToPag;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.pitt.csb.KCI.KCI;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

/**
 * Created by vinee_000 on 11/29/2016.
 */
public class kci_max_test {
    public static void main(String [] args) throws Exception{
        double [] l = {0.25,0.35,0.55};
        double [] alp = {.001,.01,.05};
        for(int i = 0; i < 3; i++)
        {
        for(int j = 0; j < 3; j ++) {
            double alpha = alp[j];
            int numGraphs = 3;
            ArrayList[] toRemove = new ArrayList[numGraphs];
            double lab = l[i];
            for (int fe = 0; fe < numGraphs; fe++) {
                String filename = "data." + (fe+1);
                String filename2 = "graph." + (fe+1) +".txt";
                String extension = ".txt";
                Graph trueGraph = GraphUtils.loadGraphTxt(new File(filename2));
                DataSet data = MixedUtils.loadDataSet2(filename + extension);
                boolean done = false;
                double[] lambda = {lab, lab, lab};
                DataSet removeStuff = data.copy();

                while(!done) {
                    removeStuff = data.copy();

                    removeStuff.permuteRows();
                    int[] removal = new int[4000];
                    for (int las = 0; las < 4000; las++)
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
                }
                int nn = 5;
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
                            System.out.println("Stuck on file " + fe + "_");
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
                Graph findLatents = GraphUtils.loadGraphTxt(new File("mult_" + i + "_" + j + "_" + fe + "_cg.txt"));
               /* for(int curr = 48; curr >=1; curr--)
                {
                    if(findLatents.getNode("X" + curr)==null)
                        trueGraph.getNode("X" + curr).setNodeType(NodeType.LATENT);
                    removeStuff.removeColumn(removeStuff.getVariable("X" + curr));
                }*/
                //System.out.println(removeInts);
                removeInts = toRemove[fe];
                for (int curr = 0; curr < removeInts.size(); curr++) {
                    latents.add(trueGraph.getNode("X" + (removeInts.get(curr) + 1)));
                    trueGraph.getNode("X" + (removeInts.get(curr) + 1)).setNodeType(NodeType.LATENT);
                    //removeStuff.removeColumn(removeInts.get(curr));
                    removeStuff.removeColumn(removeStuff.getVariable("X" + (removeInts.get(curr) + 1)));
                }
                usedInts.clear();
                 removeInts.clear();
                final DagToPag dagToPag = new DagToPag(trueGraph);
                dagToPag.setCompleteRuleSetUsed(false);
                Graph p = dagToPag.convert();
                trueGraph = GraphUtils.loadGraphTxt(new File(filename2));
                //System.out.println("True Graph " + p);
                // System.out.println("Data: " + data);
                //Need to ask about this when you don't remove any variables
                //  double[] lambda = {lab, lab, lab};
                double tolerance = 1e-7; //convergeance tolerance
                int iterLimit = 1000; //iteration limit
                MGM model = new MGM(removeStuff, lambda);
                long mgmTime = 0;
                mgmTime = System.currentTimeMillis();
                model.learnEdges(iterLimit);

                Graph mgmGraph = model.graphFromMGM();
                IndependenceTest kci = new KCI(removeStuff, alpha, alpha,trueGraph,l[i]);

                IndependenceTest mult = new IndTestMultinomialAJ(removeStuff, alpha);
                FciMaxP f = new FciMaxP(kci);
                f.setInitialGraph(mgmGraph);
                FciMaxP f2 = new FciMaxP(mult);
                f2.setInitialGraph(mgmGraph);

                Graph one = f.search();
                Graph two = f2.search();
             PrintStream out = new PrintStream("kci_" + i + "_" + j + "_" + fe + "_cg.txt");
                out.println(one);
                out.flush();
                out.close();
                PrintStream out2 = new PrintStream("mult_" + i + "_" + j + "_" + fe + "_cg.txt");
                out2.println(two);
                out2.flush();
                out2.close();
               PrintStream graphWriter = new PrintStream("true_" + fe + "_cg.txt");
                graphWriter.println(trueGraph);
               graphWriter.flush();
               graphWriter.close();
                graphWriter = new PrintStream("PAG_" + i + "_" + fe + "_cg.txt");
                graphWriter.println(p);
                graphWriter.flush();
                graphWriter.close();
            }
        }
        }
    }
}
