package edu.pitt.csb.KCI;

import edu.cmu.tetrad.algcomparison.graph.RandomGraph;
import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.algcomparison.simulation.SemSimulation;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndTestDSep;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

/**
 * Created by vinee_000 on 4/12/2018.
 */
public class KCI_LinGauss_Tests {



    public static void main(String [] args)throws Exception
    {
        String runName = "Conditional_Independencies_Small_Graph";

        PrintStream out = new PrintStream(runName + ".txt");
        out.println("Run\tVar_1\tVar_2\tConditioning_Set\tKCI_PVal\tAJ_PVal\tSampleSize\tTruth\tConditioning_Size");
        int numRuns = 5;
        int samplingSize = 500;
        int [] conditioningSize = {0,1,2,3,4};
        List<Node> z = new ArrayList<Node>();
        Random rand = new Random();
       // int [] sampleSizes = {100,500,1000};
        int [] sampleSizes= {100};
            for(int j = 0; j < numRuns;j++) {
                System.out.println("Run " + j);
                Parameters p = new Parameters();
                p.setValue("numMeasures", 50);
                p.setValue("sampleSize", 5000);
                p.setValue("numEdges", 75);
                MixedLeeHastieSimulation c = new MixedLeeHastieSimulation();
                c.simulate(p);
                DataSet temp = c.getDataSet(0).copy();
                temp.permuteRows();
                PrintStream graph = new PrintStream("Graph_" + j + ".txt");
                graph.println(c.getTrueGraph());
                graph.flush();
                graph.close();
                graph = new PrintStream("Data_" + j + ".txt");
                graph.println(temp);
                graph.flush();
                graph.close();
                for(int i = 0; i < sampleSizes.length;i++) {
                    System.out.println("Sample Size: " + sampleSizes[i]);
                    int [] rows = new int[sampleSizes[i]];
                    for(int y = 0; y < rows.length;y++)
                        rows[y] = y;
                    DataSet toUse = temp.subsetRows(rows);
                IndependenceTest k = new KCI(toUse, 0.05);
                IndependenceTest ii = new IndTestMultinomialAJ(toUse, 0.05);
                List<Node> allVars = c.getTrueGraph().getNodes();
                for(int x = 0 ; x < conditioningSize.length;x++) {



                    for(int count = 0; count < samplingSize;count++) {
                        int cs = conditioningSize[x];
                        ArrayList<Integer> allNodes = new ArrayList<Integer>();
                        for (int y = 0; y < cs + 2; y++) {
                            int ints = rand.nextInt(allVars.size());
                            while (allNodes.contains(ints))
                                ints = rand.nextInt(allVars.size());
                            allNodes.add(ints);
                        }
                        List<Node> conditioningSetK = new ArrayList<Node>();
                        List<Node> conditioningSetII = new ArrayList<Node>();

                        for (int y = 2; y < allNodes.size(); y++) {
                            conditioningSetK.add(k.getVariable(toUse.getVariable(allNodes.get(y)).getName()));
                            conditioningSetII.add(ii.getVariable(toUse.getVariable(allNodes.get(y)).getName()));
                            z.add(c.getTrueGraph().getNode(toUse.getVariable(allNodes.get(y)).getName()));
                        }

                        k.isIndependent(k.getVariable(toUse.getVariable(allNodes.get(0)).getName()), k.getVariable(toUse.getVariable(allNodes.get(1)).getName()), conditioningSetK);
                        ii.isIndependent(ii.getVariable(toUse.getVariable(allNodes.get(0)).getName()), ii.getVariable(toUse.getVariable(allNodes.get(1)).getName()), conditioningSetII);
                        IndependenceTest dsep = new IndTestDSep(c.getTrueGraph());
                        /*TetradMatrix cov = toUse.getCovarianceMatrix();
                        int [] cols = new int[conditioningSize[x]+2];
                        cols[0] = allNodes.get(0);
                        cols[1] = allNodes.get(1);
                        for(int y = 2; y < cols.length;y++)
                        {
                            cols[y] = toUse.getColumn(toUse.getVariable(z.get(y-2).getName()));
                        }
                        DataUtils.subMatrix
                        double pc = StatUtils.partialCorrelation(cov.getSelection(cols,cols));*/
                        out.println(j + "\t" + toUse.getVariable(allNodes.get(0)) + "\t" + toUse.getVariable(allNodes.get(1)) + "\t" + conditioningSetK + "\t" + k.getPValue() + "\t" + ii.getPValue() + "\t" + sampleSizes[i] + "\t" + dsep.isDependent(c.getTrueGraph().getNode(temp.getVariable(allNodes.get(0)).getName()), c.getTrueGraph().getNode(temp.getVariable(allNodes.get(1)).getName()), z) + "\t" + conditioningSize[x]);
                    }
                }
            }
            out.flush();
        }
        out.close();
    }
}
