package edu.pitt.csb.Mixed_Partition;

/**
 * Created by vinee_000 on 9/5/2017.
 */
public class testMixedData {
    public static void main(String [] args) throws Exception {
     /*       int numRuns = 10;
            int cutoffMax = 4;
            Graph[] trueGraphs = new Graph[numRuns];
            DataSet[] sets = new DataSet[numRuns];
            for(int i = 0; i < numRuns;i++)
            {
                Parameters p = new Parameters();
                //p.setValue("");
              //ContinuousLinearGaussianSemSimulation c = new ContinuousLinearGaussianSemSimulation();
                c.simulate(p);
                DataSet d =  c.getDataSet(0);
                Graph g = c.getTrueGraph();


          /*  ContinuousNonlinearNongaussianSimulation c = new ContinuousNonlinearNongaussianSimulation();
            c.simulate(p);
            Graph truth = c.getTrueGraph();
            DataSet d = c.getDataSet(0);
            trueGraphs[i] = truth;
            sets[i] = d;
            }
            Comparison c = new Comparison();

//IGNORE EVERYTHING BELOW HERE

            for(int index = 0; index < cutoffMax;index++) {
                PrintStream out = new PrintStream("time_" + (index) + ".txt");
                PrintStream out2 = new PrintStream("accuracy_" + (index) + ".txt");

                out.println("KCI TIME\tZTIME");
                out2.println("Test\tAdj_Prec\tAdj_Rec\tAH_Prec\tAH_Rec");
                long[] kciTime = new long[numRuns];
                long[] zTime = new long[numRuns];
                double[][][] avgCounts = new double[numRuns][4][4];
                double[][][] avgCounts2 = new double[numRuns][4][4];
                for (int i = 0; i < numRuns; i++) {
                    DataSet d = sets[i];
                    Graph truth = trueGraphs[i];
                    KCI k = new KCI(d, .05, .05);
                    k.setCutoff(index-1);
                    CpcStable p2 = new CpcStable(k);
                    kciTime[i] = System.nanoTime();
                    Graph est = p2.search();
                    kciTime[i] = System.nanoTime() - kciTime[i];
                    zTime[i] = System.nanoTime();
                    p2 = new CpcStable(new IndTestMultinomialAJ(d, .05));
                    zTime[i] = System.nanoTime() - zTime[i];
                    Graph est2 = p2.search();
                    System.out.println(PerformanceTests.endpointMisclassification(est.getNodes(),est,truth));
                    System.out.println(PerformanceTests.endpointMisclassification(est2.getNodes(),est2,truth));

                    int[][] counts = PerformanceTests.endpointMisclassification2(est.getNodes(), est, truth);
                    int[][] counts2 = PerformanceTests.endpointMisclassification2(est2.getNodes(), est2, truth);
                    for (int ii = 0; ii < 4; ii++) {
                        for (int jj = 0; jj < 4; jj++) {
                            avgCounts[i][ii][jj] = counts[ii][jj];
                            avgCounts2[i][ii][jj] = counts2[ii][jj];
                        }

                    }
                }
                double[][] stats = new double[4][numRuns];
                double[][] stats2 = new double[4][numRuns];

                for (int i = 0; i < numRuns; i++) {
                    stats[0][i] = PerformanceTests.adjPrec(avgCounts[i]);
                    stats[1][i] = PerformanceTests.adjRec(avgCounts[i]);
                    stats[2][i] = PerformanceTests.AHPrec(avgCounts[i]);
                    stats[3][i] = PerformanceTests.AHRec(avgCounts[i]);

                    stats2[0][i] = PerformanceTests.adjPrec(avgCounts2[i]);
                    stats2[1][i] = PerformanceTests.adjRec(avgCounts2[i]);
                    stats2[2][i] = PerformanceTests.AHPrec(avgCounts2[i]);
                    stats2[3][i] = PerformanceTests.AHRec(avgCounts2[i]);

                }

                out2.print("KCI\t");
                for (int i = 0; i < 4; i++) {
                    out2.print(StatUtils.mean(stats[i]) + "(" + StatUtils.sd(stats[i]) + ")");
                    if (i == 3)
                        out2.println();
                    else
                        out2.print("\t");
                }

                out2.print("Simple\t");
                for (int i = 0; i < 4; i++) {
                    out2.print(StatUtils.mean(stats2[i]) + "(" + StatUtils.sd(stats2[i]) + ")");
                    if (i == 3)
                        out2.println();
                    else
                        out2.print("\t");
                }
                out2.flush();
                out2.close();
                for (int i = 0; i < numRuns; i++) {
                    out.println(kciTime[i] / Math.pow(10, 9) + "\t" + zTime[i] / Math.pow(10, 9));
                }
                out.flush();
                out.close();
            }*/
        }
}
