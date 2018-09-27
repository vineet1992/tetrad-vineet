package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.algcomparison.statistic.utilities.AdjacencyConfusion;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * Created by vinee_000 on 9/1/2017.
 */
public class testStepsPar {
    public static void main(String [] args) throws Exception
    {
        int numVariables = 50;
        int sampleSize = 30;
        int numCategories = 2;
        int ns = 20;
        MixedLeeHastieSimulation c = new MixedLeeHastieSimulation();
        Parameters p = new Parameters();
        p.setValue("numMeasures", numVariables);
        NormalDistribution n = new NormalDistribution(2*numVariables,numVariables*.75);
        p.setValue("numEdges", (int)n.sample());
        p.setValue("sampleSize", sampleSize);
        p.setValue("numCategories",numCategories);
        double low = .1;
        double high = .8;
        double[] initLambdas = new double[40];
        for (int i = 0; i < 40; i++) {
            initLambdas[i] = i * (high - low) / 40 + low;
        }
        int numRuns = 10;
        PrintStream out = new PrintStream("steps_log.txt");
        out.println("Run\tTime_Non_Par\tTime_Par\tPrec_Non_Par\tPrec_Par\tRec_Non_Par\tRec_Par\tLambdas_Non_Par\tLambdas_Par");
        for(int i = 0; i < numRuns;i++) {
            c.simulate(p);
            DataSet curr = c.getDataSet(0);
            Graph truth = c.getTrueGraph();

            STEPS s = new STEPS(curr,initLambdas,.1,ns,false);
            long time = System.nanoTime();
            Graph g = s.runSteps();
            double [] lbNon = s.lastLambda;
            time = System.nanoTime()-time;
            long time2 = System.nanoTime();
            Graph g2 = s.runStepsPar();
            double [] lb = s.lastLambda;
            time2 = System.nanoTime()-time2;
            out.print(i + "\t" + time + "\t" + time2 + "\t");
                    AdjacencyConfusion adjConfusion = new AdjacencyConfusion(truth, g);
            int adjTp = adjConfusion.getAdjTp();
            int adjFp = adjConfusion.getAdjFp();
            int adjFn = adjConfusion.getAdjFn();
            int adjTn = adjConfusion.getAdjTn();
            double precNonPar = adjTp / (double) (adjTp + adjFp);
            double recNonPar = adjTp/ (double) (adjTp + adjFn);

            adjConfusion = new AdjacencyConfusion(truth, g2);
            adjTp = adjConfusion.getAdjTp();
            adjFp = adjConfusion.getAdjFp();
            adjFn = adjConfusion.getAdjFn();
            adjTn = adjConfusion.getAdjTn();
            double precPar = adjTp / (double) (adjTp + adjFp);
            double recPar = adjTp/ (double) (adjTp + adjFn);
            out.println(precNonPar + "\t" + precPar + "\t" + recNonPar + "\t" + recPar + "\t" + Arrays.toString(lbNon) + "\t" + Arrays.toString(lb));
            out.flush();


        }
        out.close();
    }
}
