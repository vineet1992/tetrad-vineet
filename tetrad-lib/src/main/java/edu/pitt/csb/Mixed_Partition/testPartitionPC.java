package edu.pitt.csb.Mixed_Partition;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;
import edu.cmu.tetrad.search.PcStableExp;
import edu.cmu.tetrad.util.StatUtils;
import edu.pitt.csb.CompareGraphs;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;

/**
 * Created by vinee_000 on 9/8/2017.
 */
public class testPartitionPC {
    public static void main(String [] args)throws Exception
    {
        Parameters p = new Parameters();
        p.setValue("numMeasures",300);
    int numRuns = 10;
        p.put("numEdges", 2 * p.getInt("numMeasures"));
        p.put("sampleSize", 200);
        p.put("alpha", .01);
        p.put("numCategories",4);
        MixedLeeHastieSimulation m = new MixedLeeHastieSimulation();
        double [][] times = new double[2][numRuns];
        double [][] scores = new double[8][numRuns];

        for(int i = 0; i < numRuns;i++)
        {
            m.simulate(p);
            DataSet d  = m.getDataSet(0);
            Graph g = m.getTrueGraph();
            IndependenceTest ii = new IndTestMultinomialAJ(d,.01);
            PcStableExp p2 = new PcStableExp(ii,true);
            long time = System.nanoTime();
            Graph est = p2.search();
            double ftime = ((System.nanoTime()-time)/Math.pow(10,9));
            times[0][i] = ftime;
            PcStable p3 = new PcStable(ii);
            CompareGraphs c = new CompareGraphs(g,est);
            System.out.println("Partitioned:");
            System.out.println(c.getPrec() + "," + c.getRec());
            System.out.println(c.getArrowPrec() + "," + c.getArrowRec());

            long time2 = System.nanoTime();
            Graph est2 = p3.search();
            double gtime = ((System.nanoTime()-time2)/Math.pow(10,9));
            times[1][i] = gtime;
            c = new CompareGraphs(g,est2);
            System.out.println("Normal:");
            System.out.println(c.getPrec() + "," + c.getRec());
            System.out.println(c.getArrowPrec() + "," + c.getArrowRec());
            System.out.println("Partitioned: " + ftime + "," + "Default: " + gtime);
            System.out.println("\n\n\n");


        }
        System.out.println("Average Time Paritioned: " + StatUtils.mean(times[0]) + " (" + StatUtils.sd(times[0]) + " )");
        System.out.println("Average Time Default: " + StatUtils.mean(times[1]) + " (" + StatUtils.sd(times[1]) + " )");
    }
}
