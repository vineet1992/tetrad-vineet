package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.util.StatUtils;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Random;

/**
 * Created by vinee_000 on 3/23/2018.
 */
public class testMGM {

    public static void main(String [] args) throws Exception
    {
        int numRuns = 2000;
        int numVars = 2500;
        int samples = 1500;
        PrintStream out = new PrintStream("Runtime.txt");
        out.println("Num_Variables\tGraph_Density\tSample_Size\tLambda\tRuntime\tTime_Per_Iteration\tNum_Iterations");
        Parameters p = new Parameters();
        Random rand = new Random();
        A:for(int i = 0; i < numRuns;i++) {


            int numMeasures = rand.nextInt(numVars) + 50;
            int numEdges = (int) (rand.nextDouble() * 6 * numMeasures);
            int sampleSize = rand.nextInt(samples) + 100;
            p.setValue("numMeasures", numMeasures);
            p.setValue("numEdges", numEdges);
            p.setValue("percentDiscreteForMixedSimulation", 50);
            p.setValue("numCategories", 4);
            p.setValue("sampleSize", sampleSize);
            try {
            MixedLeeHastieSimulation m = new MixedLeeHastieSimulation();
            m.simulate(p);

            DataSet d = m.getDataSet(0);

            int attempts = 5000;
            while (!checkForVariance(d)) {
                m = new MixedLeeHastieSimulation();
                m.simulate(p);
                attempts--;
                d = m.getDataSet(0);
                if(attempts==0)
                    continue A;

            }
            double lambda = rand.nextDouble()*0.85+0.05;
            System.out.println(numMeasures + "\t" + numEdges/(double)numMeasures + "\t" + sampleSize + "\t" + lambda);
            long time = -1;

                MGM m2 = new MGM(d, new double[]{lambda, lambda, lambda});
                m2.setTimeout(43200*(long)Math.pow(10,9));
                time = System.nanoTime();
                m2.learnEdges(1000);
                time = System.nanoTime()-time;
                System.out.println(m2.iterCount + "\t" + m2.timePerIter);
                out.println(numMeasures + "\t" + numEdges/(double)numMeasures + "\t" + sampleSize + "\t" + lambda + "\t" + time/Math.pow(10,9) + "\t" + m2.timePerIter/((double)m2.iterCount) + "\t" + m2.iterCount);
                out.flush();
            }
            catch(Exception e)
            {
                continue A;
            }

        }
        out.close();

    }
    public static boolean checkForVariance(DataSet d)
    {
        for(int i = 0; i < d.getNumColumns();i++)
        {
            if(d.getVariable(i)instanceof ContinuousVariable)
            {
                double [] temp = d.getDoubleData().getColumn(i).toArray();
                if(StatUtils.variance(temp)<0.0001)
                    return false;
            }else
            {
                HashMap<Integer,Integer> map = new HashMap<Integer,Integer>();
                for(int j = 0; j < d.getNumRows();j++)
                {
                    if(map.get(d.getInt(j,i))==null)
                    {
                        map.put(d.getInt(j,i),1);
                    }else
                    {
                        map.put(d.getInt(j,i),map.get(d.getInt(j,i))+1);
                    }
                }
                for(Integer ii:map.keySet())
                {
                    if(map.get(ii)<4)
                        return false;
                }
            }
        }
        return true;
    }
}
