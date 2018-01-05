package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.performance.PerformanceTests;
import edu.cmu.tetrad.search.DagToPag;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

/**
 * Created by vinee_000 on 10/22/2016.
 */
public class fciMaxTests {
    public static void main(String [] args) throws Exception
    {
        for(int i = 0; i < 10 ; i ++) {
            PrintStream write = new PrintStream("out_" + i + "_graph.txt");
            Graph trueGraph = GraphUtils.loadGraphTxt(new File("SF_" + i + "_graph.txt"));
            DataSet d = MixedUtils.loadDataSet2("SF_" + i + "_data_NonMono.txt");
            DataSet removeStuff = d.copy();

            int numNodesToRemove = 4;
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
                B:for(int nod = 0; nod < totalNodes;nod++) {
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
                        break B;
                    }
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
            DagToPag dtp = new DagToPag(trueGraph);
           Graph pag =  dtp.convert();
            usedInts.clear();

            IndependenceTest ii = new IndTestMultinomialAJ(removeStuff,.05);
            FciMaxP p = new FciMaxP(ii);
            Graph out = p.search();
            write.println(out);
            System.out.println(i);
            System.out.println(PerformanceTests.endpointMisclassification(out.getNodes(),out,pag));

        }
    }
}
