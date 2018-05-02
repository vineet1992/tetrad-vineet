package edu.pitt.csb.KCI;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.File;
import java.util.ArrayList;

/**
 * Created by vinee_000 on 6/2/2016.
 */
public class KCI_Ind_Debug_2 {
    public static void main(String [] args) throws Exception
    {
        int i = 2;
        Graph g = GraphUtils.loadGraphTxt(new File("C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/Sedgewick_CausalMGM_data/Sedgewick_CausalMGM_data/simulated_data/dataAndNetworks/DAG_" + i + "_graph.txt"));
        DataSet ds1 = MixedUtils.loadDataSet2("C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/Sedgewick_CausalMGM_data/Sedgewick_CausalMGM_data/simulated_data/updated_data/DAG_" + i + "_data_java.txt");
        DataSet ds2 = MixedUtils.loadDataSet2("C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2015/Mixed-Latent/Sedgewick_CausalMGM_data/Sedgewick_CausalMGM_data/simulated_data/updated_data/DAG_" + i + "_data_java.txt");
        int[] selectedRows = new int[4500];
        for (int j = 0; j < 4500; j++)
            selectedRows[j] = j;
        ds1.removeRows(selectedRows);
        ds2.removeRows(selectedRows);
        // System.out.println(ds.getNumRows());
        double alpha2 = 0.05;
        double alpha = 0.05;

        KCI_Ind k3 = new KCI_Ind(ds1,ds2,alpha,alpha2,g,Integer.MAX_VALUE);

        IndTestMultinomialAJ k4 = new IndTestMultinomialAJ(ds1,.05);
        ArrayList<Node> z = new ArrayList<Node>();
        z.add(ds2.getVariable(19));
       z.add(ds2.getVariable(24));
        System.out.println(ds2.getVariable(19));
        //z.clear();
        k3.isIndependent(ds2.getVariable(10),ds2.getVariable(29),z);
        System.out.println(k3.getPValue());
        k4.isIndependent(ds2.getVariable(10),ds2.getVariable(29),z);
        System.out.println(k4.getPValue());

    }
}
