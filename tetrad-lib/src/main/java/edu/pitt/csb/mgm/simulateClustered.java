package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.score.SemBicScore;
import edu.cmu.tetrad.algcomparison.simulation.ContinuousLinearGaussianSemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.performance.PerformanceTests;
import edu.cmu.tetrad.search.Fgs2;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.search.PcStable;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 12/12/2016.
 */
public class simulateClustered {
    public static void main(String [] args) throws Exception {
        int maxVars = 75;
        int numRuns = 20;
        PrintStream out = new PrintStream("data_dense.txt");
        out.println("Vars\tAPrec\tARec\tAHPrec\tAHRec");
        PrintStream out2 = new PrintStream("matrices_dense.txt");
        PrintStream out3 = new PrintStream("data_dense_fgs.txt");
        out3.println("Vars\tAPrec\tARec\tAHPrec\tAHRec");
        PrintStream out4 = new PrintStream("matrices_dense_fgs.txt");

for(int i = 0; i < maxVars; i +=10) {
//Parameters, how correlated are the variables? 0.01 to 0.05
    //How many deterministic variables are there?
    double [][] avgCounts = new double [4][4];
    double [][] avgCounts2 = new double[4][4];
    for(int j = 0; j < numRuns;j++) {


        Parameters p = new Parameters();
        p.setValue("numMeasures",maxVars);
        ContinuousLinearGaussianSemSimulation c = new ContinuousLinearGaussianSemSimulation();
        c.setruns(i);
        c.simulate(p);

// Get average endpoint misclassification matrices for increasing number of determinisitic variables



        PcStable p2 = new PcStable(new IndTestFisherZ(c.getDataSet(0), .05));
        edu.cmu.tetrad.search.Fgs2 f = new Fgs2(new edu.cmu.tetrad.search.SemBicScore2(new CovarianceMatrix(c.getDataSet(0))));
        Graph est = p2.search();
        Graph est2 = f.search();
        int[][] counts = PerformanceTests.endpointMisclassification2(p2.getNodes(), est, c.getTrueGraph());
        System.out.println("(" + i + "," + j + ")");
        System.out.println(PerformanceTests.convert(counts));
        for(int ii = 0; ii < 4;ii++) {
            for (int jj = 0; jj < 4; jj++)
            {
                avgCounts[ii][jj] = avgCounts[ii][jj] + counts[ii][jj]/(double)numRuns;
            }

        }

        int[][] counts2 = PerformanceTests.endpointMisclassification2(c.getDataSet(0).getVariables(), est2, c.getTrueGraph());
        System.out.println("(" + i + "," + j + ")");
        System.out.println(PerformanceTests.convert(counts2));
        for(int ii = 0; ii < 4;ii++) {
            for (int jj = 0; jj < 4; jj++)
            {
                avgCounts2[ii][jj] = avgCounts2[ii][jj] + counts2[ii][jj]/(double)numRuns;
            }

        }

    }
    out2.println(i + " Correlated Variables");
    out2.println(PerformanceTests.convert(avgCounts));
    out2.println();
    out2.flush();


    //Need to save variables of interest before they are reset
out.println(i + "\t" + PerformanceTests.adjPrec(avgCounts)+ "\t" + PerformanceTests.adjRec(avgCounts) + "\t" + PerformanceTests.AHPrec(avgCounts) + "\t" + PerformanceTests.AHRec(avgCounts));
    out.flush();

    out4.println(i + " Correlated Variables");
    out4.println(PerformanceTests.convert(avgCounts2));
    out4.println();
    out4.flush();


    //Need to save variables of interest before they are reset
    out3.println(i + "\t" + PerformanceTests.adjPrec(avgCounts2)+ "\t" + PerformanceTests.adjRec(avgCounts2) + "\t" + PerformanceTests.AHPrec(avgCounts2) + "\t" + PerformanceTests.AHRec(avgCounts2));
    out3.flush();
}
out.close();
        out2.close();
        out3.close();
        out4.close();



    }

}
