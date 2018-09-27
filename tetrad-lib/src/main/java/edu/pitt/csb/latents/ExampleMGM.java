package edu.pitt.csb.latents;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.*;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.List;

/**
 * Created by vinee_000 on 7/16/2018.
 */
public class ExampleMGM {
    public static void main(String[]args) throws Exception
    {
        //First step is loading in the dataset

        //Replace file.txt with the filename of your data
        DataSet d = MixedUtils.loadDataSet2("file.txt");
        //Need to setup parameters for MGM, these tell the algorithm how many/few edges to include
        //A larger number means that fewer edges are included, a lower number means more edges are included
        double [] lambda = {0.2,0.2,0.2};

        //Create an MGM object with lambda parameters
        MGM m = new MGM(d,lambda);

        //Run MGM to get an output graph
        m.learnEdges(1000);

        //Extract the graph from MGM
        Graph g = m.graphFromMGM();

        //Create an independence test object for FCI-MAX to use
        IndependenceTest ind = new IndTestMultinomialAJ(d,0.05,true);
        //Run FCI-MAX to get edge directions
        FciMaxP f = new FciMaxP(ind);
        Graph out = f.search();

        //Save the graph created by FCI-MAX
        GraphUtils.saveGraph(out,new File("output.txt"),false);


        PrintStream output = new PrintStream("latents_tumor.txt");


        LatentPrediction lp = new LatentPrediction(d,4,0.01,15);
        List<LatentPrediction.Pair> latents = lp.runRegularAlgorithm("MGM-FCI-MAX",false);
        for(LatentPrediction.Pair lts: latents)
        {
            output.println(lts);
            output.flush();
        }
        output.close();

    }
}
