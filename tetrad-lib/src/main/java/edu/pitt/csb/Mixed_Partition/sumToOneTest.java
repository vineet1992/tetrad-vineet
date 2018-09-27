package edu.pitt.csb.Mixed_Partition;

import cern.jet.random.NegativeBinomial;
import cern.jet.random.engine.RandomEngine;
import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.EdgeListGraphSingleConnections;
import edu.cmu.tetrad.graph.Endpoint;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.StatUtils;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;

import java.io.PrintStream;
import java.util.Random;

/**
 * Created by vinee_000 on 4/13/2018.
 */
public class sumToOneTest {

    public static void main(String [] args)throws Exception
    {
        Random rand = new Random();

        PrintStream out = new PrintStream("Sum_to_one_tests_500.txt");
        out.println("Run\tIs_Edge\tcorr_nb\tcorr_scaled\tcorr_NPN\tp_val_nb\tp_val_scaled\tp_val_NPN\tsample_size");
        int numRuns = 1000;
        int [] samples = {100,500,1000};
        int numVars = 500;
        for(int i = 0; i < numRuns;i++) {
            Parameters p = new Parameters();
            p.setValue("sampleSize", 1000);
            p.setValue("numMeasures",numVars);
            MixedLeeHastieSimulation c = new MixedLeeHastieSimulation();
            //TODO


            Graph g = new EdgeListGraphSingleConnections();
            g.addNode(new ContinuousVariable("X"));
            g.addNode(new ContinuousVariable("Y"));
            boolean edge = false;
            if(rand.nextBoolean()) {
                g.addEdge(new Edge(g.getNode("X"), g.getNode("Y"), Endpoint.TAIL, Endpoint.ARROW));
                edge = true;
            }
            for(int j = 2;j < numVars;j++)
            {
                g.addNode(new ContinuousVariable("Z" + j));
            }
            c.setTrueGraph(g);
            c.simulate(p);
            double prob = rand.nextDouble()*0.5 + 0.49;
            NegativeBinomial n = new NegativeBinomial(1,prob, RandomEngine.makeDefault());
            prob = rand.nextDouble()*0.5+0.49;
            NegativeBinomial n_error = new NegativeBinomial(1,prob,RandomEngine.makeDefault());
            double beta = 0;
            if(rand.nextBoolean())
            {
                beta = rand.nextDouble()+0.5;
            }
            else
            {
                beta = rand.nextDouble()-1.5;
            }
            NegativeBinomial [] others = new NegativeBinomial[numVars-2];
            for(int k = 0; k < others.length;k++)
            {
                prob = rand.nextDouble()*0.5+0.49;
                others[k] = new NegativeBinomial(1,prob,RandomEngine.makeDefault());
            }
            DataSet data = c.getDataSet(0);
            for(int j = 0; j < data.getNumRows();j++) {
                if (edge) {
                    data.setDouble(j, 0, n.nextDouble());
                    data.setDouble(j, 1, data.getDouble(j, 0) * beta + n_error.nextDouble());
                    while(data.getDouble(j,0)==0&&data.getDouble(j,1)==0)
                    {
                        data.setDouble(j,0,n.nextDouble());
                        data.setDouble(j,1,data.getDouble(j,0)*beta*n_error.nextDouble());
                    }
                    //   if(data.getDouble(j,1)<0)
                    //     data.setDouble(j,1,0);
                } else {
                    data.setDouble(j, 0, n.nextDouble());
                    data.setDouble(j, 1, n_error.nextDouble());
                    while(data.getDouble(j,0)==0 && data.getDouble(j,1)==0)
                    {
                        data.setDouble(j,0,n.nextDouble());
                        data.setDouble(j,1,n_error.nextDouble());
                    }
                }

                for(int k = 2; k < data.getNumColumns();k++)
                {
                    data.setDouble(j,k,others[k-2].nextDouble());
                }
            }
            DataSet d2 = data.copy();

            for(int j = 0; j < d2.getNumRows();j++)
            {
                double denom = 0;
                for(int k = 0; k < d2.getNumColumns();k++)
                {
                    denom+=d2.getDouble(j,k);
                }
                for(int k = 0; k < d2.getNumColumns();k++)
                {
                    d2.setDouble(j,k,d2.getDouble(j,k)/denom);
                }
            }

            DataSet d3 = DataUtils.getNonparanormalTransformed(d2);
            for(int j = 0; j < samples.length;j++)
            {
                int [] rows = new int[samples[j]];
                for(int k = 0; k < rows.length;k++)
                    rows[k]=k;
                DataSet data_mod = data.subsetRows(rows);
                DataSet d2_mod = d2.subsetRows(rows);
                DataSet d3_mod = d3.subsetRows(rows);
                IndependenceTest ii = new IndTestMultinomialAJ(data_mod,0.05);
                IndependenceTest ii2 = new IndTestMultinomialAJ(d2_mod,0.05);
                IndependenceTest ii3 = new IndTestMultinomialAJ(d3_mod,0.05);
                ii.isIndependent(ii.getVariable("X"),ii.getVariable("Y"));
                ii2.isIndependent(ii2.getVariable("X"),ii2.getVariable("Y"));
                ii3.isIndependent(ii3.getVariable("X"),ii3.getVariable("Y"));


                out.println(i + "\t" + edge + "\t" + data_mod.getCorrelationMatrix().get(data_mod.getColumn(data_mod.getVariable("X")),data_mod.getColumn(data_mod.getVariable("Y"))) + "\t" + d2_mod.getCorrelationMatrix().get(d2_mod.getColumn(d2_mod.getVariable("X")),d2_mod.getColumn(d2_mod.getVariable("Y"))) + "\t" + d3_mod.getCorrelationMatrix().get(d3_mod.getColumn(d3_mod.getVariable("X")),d3_mod.getColumn(d3_mod.getVariable("Y"))) + "\t" + ii.getPValue() + "\t" + ii2.getPValue() + "\t" + ii3.getPValue() + "\t" +samples[j]);
                out.flush();
            }


        }
        out.flush();
        out.close();
    }
}
