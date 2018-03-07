package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.StatUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by vinee_000 on 3/2/2018.
 */
public class simulatePrefDivData {

    public static void main(String [] args)
    {
        int numRuns = 1;
        int numGenes = 100;
        int sampleSize = 10000;
        int ns = 20;
        int minTargetParents = 5;
        double intensityNoise = 0.1;
        double dissimilarityNoise = 0.1;
        boolean targetContinuous = true;
        double percentMissing = 0;
        int numCategories = 4;
        MixedLeeHastieSimulation m = new MixedLeeHastieSimulation();
        Parameters p = new Parameters();

        p.setValue("numMeasures",numGenes+1);
        p.setValue("sampleSize",sampleSize);
        if(targetContinuous)
            p.setValue("percentDiscreteForMixedSimulation",0);
        else
            p.setValue("percentDiscreteForMixedSimulation",100/(double)numGenes);
        p.setValue("numCategories",numCategories);

        for(int j = 0; j < numRuns;j++) {
            NormalDistribution n = new NormalDistribution(numGenes * 2, numGenes / 2);
            p.setValue("numEdges", (int) n.sample());


            m.simulate(p);
            String target = "";
            A: while (true) //Repeat until we have one target variable with a sufficient number of parents
            {
                Graph g = m.getTrueGraph();
                DataSet d = m.getDataSet(0);
                for (int i = 0; i < d.getNumColumns(); i++) {
                    if (!targetContinuous && d.getVariable(i) instanceof DiscreteVariable) {
                        if (g.getParents(g.getNode(d.getVariable(i).getName())).size() >= minTargetParents && g.getChildren(g.getNode(d.getVariable(i).getName())).size()==0) {
                            target = d.getVariable(i).getName();
                            break A;
                        }
                    } else if (g.getParents(g.getNode(d.getVariable(i).getName())).size() >= minTargetParents && g.getChildren(g.getNode(d.getVariable(i).getName())).size()==0) {
                        target = d.getVariable(i).getName();
                        break A;
                    }
                }
                m.simulate(p);
            }
            DataSet d = m.getDataSet(0);
            NormalDistribution n_intense = new NormalDistribution(0,intensityNoise);
            NormalDistribution n_dissim = new NormalDistribution(0,dissimilarityNoise);
            Graph g = m.getTrueGraph();
            System.out.println("Data: " + d);
            System.out.println("Truth:" + g + "\n Target: " + target);
            double [][] data = d.getDoubleData().transpose().toArray();
            int targetColumn = d.getColumn(d.getVariable(target));
            ArrayList<Gene> genes = new ArrayList<Gene>();
            int count = 0;
            for(int i = 0; i < d.getNumColumns();i++)
            {
                if(!d.getVariable(i).getName().equals(target))
                {
                    Gene g2 = new Gene(count);
                    g2.symbol = d.getVariable(i).getName();
                    if(targetContinuous)
                    {
                        double corr = Math.abs(StatUtils.correlation(data[i],data[targetColumn]));
                        g2.theoryIntensity = corr + n_intense.sample();
                    }
                    else
                    {
                        double corr = Functions.mixedMI(data[i],data[targetColumn],3);
                        g2.theoryIntensity = corr+n_intense.sample();
                    }
                    genes.add(g2);
                    count++;
                }
            }

            float [] allCorrs = new float[genes.size()];
            for(int i = 0; i < genes.size();i++)
            {
                allCorrs[i] = (float)genes.get(i).theoryIntensity;
            }
            allCorrs = Functions.NPN(allCorrs,true);
            for(int i = 0; i < genes.size();i++)
            {
                genes.get(i).theoryIntensity = allCorrs[i];
            }
            System.out.println("Correlations: " + Arrays.toString(allCorrs));

            float [] dissimilarity = new float[genes.size()*(genes.size()-1)/2];
            int index = 0;
            for(int i = 0; i < genes.size();i++)
            {
                double [] one = data[d.getColumn(d.getVariable(genes.get(i).symbol))];
                for(int k = i+1; k < genes.size();k++)
                {
                    double [] two = data[d.getColumn(d.getVariable(genes.get(k).symbol))];
                    double corr = Math.abs(StatUtils.correlation(one,two));
                    dissimilarity[index] = (float) corr;
                    index++;
                }
            }
            dissimilarity = Functions.NPN(dissimilarity,true);

            RunPrefDiv r = new RunPrefDiv(dissimilarity,genes,d,target,false);
            r.setApproxCorrelations(true);
            r.setClusterStability(true);
            r.setRadius(0.25);
            r.setAccuracy(0);
            r.setNS(ns);
            r.setTopK(minTargetParents);
            r.setThreshold(0.01);
            ArrayList<Gene> result = r.runPD();
            System.out.println("Top Genes: " + result);
            System.out.println("Clusters: " + r.getClusters());
        }
    }
}
