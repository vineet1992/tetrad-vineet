package edu.pitt.csb.latents;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.pitt.csb.latents.LatentPrediction;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import java.io.File;
import java.io.PrintStream;
import java.util.List;

public class MGM_project{
    public static void main(String[]args) throws Exception
    {
        //First step is loading in the data set

        //Replace file.txt with the filename of your data

        DataSet d = MixedUtils.loadDataSet2("finaltumordata_cleaned.txt");
        d.removeColumn(d.getVariable("patient.gender"));

        //Need to setup parameters for MGM, these tell the algorithm how many/few edges to include
        //A larger number means that fewer edges are included, a lower number means more edges are included
        double [] lambda = {0.2,0.2,0.2};

        PrintStream output = new PrintStream("latents_tumor.txt");


        LatentPrediction lp = new LatentPrediction(d,4,0.01,10);
        List<LatentPrediction.Pair> latents = lp.runRegularAlgorithm("MGM-FCI-MAX",false);
        for(LatentPrediction.Pair lts: latents)
        {
            output.println(lts);
            output.flush();
        }
        output.close();

    }
}

