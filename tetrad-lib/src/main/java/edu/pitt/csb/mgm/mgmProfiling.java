package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;

import java.io.PrintStream;
import java.util.HashMap;

/**
 * Created by vinee_000 on 7/20/2018.
 */
public class mgmProfiling {
    public static void main(String [] args)throws Exception
    {
        //int [] numVars = {50,100,200,500,1000,1500};
        int [] numVars = {100,200,500,1000,1500};
        int [] sampleSize = {100,500,1000};
        int [] density = {2,4,6};
        int numRuns = 5;
        double [] lowLambdas = {0.05,0.08,0.1,0.15,0.2};
        double [] highLambdas = {0.15,0.2,0.35,0.5,0.8};
        PrintStream out = new PrintStream("Runtimes.txt");
        out.println("Run\tNumVars\tNumSamples\tGraph_Density\tLambda\tTime");
        for(int i = 0; i < numVars.length;i++)
        {
            for(int j = 0; j < sampleSize.length;j++)
            {
                for(int d = 0; d < density.length;d++) {
                    for (int k = 0; k < numRuns; k++) {
                        Parameters p = new Parameters();
                        p.setValue("numMeasures", numVars[i]);
                        p.setValue("sampleSize", sampleSize[j]);
                        p.setValue("numEdges",density[d]*numVars[i]);
                        p.setValue("numCategories",4);
                        p.setValue("percentDiscreteForMixedSimulation",50);
                        MixedLeeHastieSimulation m = new MixedLeeHastieSimulation();
                        m.simulate(p);
                        while(!checkData(m.getDataSet(0),4))
                        {
                            m.simulate(p);
                        }


                        if(numVars[i]>=500)
                        {
                            for(int l = 0; l < highLambdas.length;l++)
                            {
                                System.out.println("Running MGM with " + numVars[i] + " variables, " + sampleSize[j] + " samples, " + density[d] + " graph density, for run #" + k + ", with lambda " + lowLambdas[l]);
                                //System.out.println(m.getTrueGraph());
                                long time = System.nanoTime();
                                MGM mgm = new MGM(m.getDataSet(0),new double[]{highLambdas[l],highLambdas[l],highLambdas[l]});
                                mgm.learnEdges(1000);
                                time = System.nanoTime()-time;
                                System.out.println("Runtime: " + time/Math.pow(10,9));
                                out.println(k + "\t" + numVars[i] + "\t" + sampleSize[j] + "\t" + density[d] + "\t" + highLambdas[l] + "\t" + time/Math.pow(10,9));
                            }
                        }
                        else
                        {
                            for(int l = 0; l < lowLambdas.length;l++)
                            {
                                System.out.println("Running MGM with " + numVars[i] + " variables, " + sampleSize[j] + " samples, " + density[d] + " graph density, for run #" + k + ", with lambda " + lowLambdas[l]);
                                //System.out.println(m.getTrueGraph());
                                long time = System.nanoTime();
                                MGM mgm = new MGM(m.getDataSet(0),new double[]{lowLambdas[l],lowLambdas[l],lowLambdas[l]});
                                mgm.learnEdges(1000);
                                time = System.nanoTime()-time;
                                System.out.println("Runtime: " + time/Math.pow(10,9));
                                out.println(k + "\t" + numVars[i] + "\t" + sampleSize[j] + "\t" + density[d] + "\t" + lowLambdas[l] + "\t" + time/Math.pow(10,9));
                            }
                        }
                    }
                }
            }
        }
    }
    public static boolean checkData(DataSet d, int bound)
    {
        for(int i = 0; i < d.getNumColumns();i++)
        {
            if(d.getVariable(i)instanceof DiscreteVariable)
            {
                HashMap<Integer,Integer> map = new HashMap<Integer,Integer>();
                for(int j = 0; j < d.getNumRows();j++)
                {
                    if(map.get(d.getInt(j,i))==null)
                        map.put(d.getInt(j,i),1);
                    else
                        map.put(d.getInt(j,i),map.get(d.getInt(j,i))+1);
                }
                for(int j = 0; j < 4;j++)
                {
                    if(map.get(j)==null || map.get(j)<=bound)
                        return false;
                }
            }
        }
        return true;
    }
}
