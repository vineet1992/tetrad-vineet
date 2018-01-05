package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;

import java.io.PrintStream;
import java.util.List;

/**
 * Created by vinee_000 on 6/16/2016.
 */
public class testLOOCV {
    public static void main(String [] args) throws Exception{
        DataSet dat = MixedUtils.loadDataSet2("final_data_600.txt");
        double [] lambda = {.1,.1,.1};
        double alpha = .5;
       /* MGM m = new MGM(dat,lambda);
        m.learnEdges(400);
        Graph g = m.graphFromMGM();
        List<Node> neighbors = g.getAdjacentNodes(dat.getVariable(0));
        if (neighbors.size() == 0) {
            System.out.println("No neighbors of response");
        }
        PrintStream out = new PrintStream("neighbors_mgm.txt");
        for (Node nod : neighbors) {
            out.println(dat.getColumn(nod));
        }
        out.flush();
        out.close();*/
        LOOCV l = new LOOCV(dat,"MGM",alpha,lambda,0);
        l.run();
        System.out.println(l.stability);
        System.out.println(l.geneRanks);
    }
}
