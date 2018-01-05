package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.DataSet;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 9/11/2017.
 */
public class debugSteps {
    public static void main(String [] args)throws Exception
    {
        /*MixedLeeHastieSimulation m = new MixedLeeHastieSimulation();
        Parameters p = new Parameters();
        p.setValue("numMeasures",40);
        p.setValue("percentDiscreteForMixedSimulation",10);

        m.simulate(p);
        System.out.println(m.getDataSet(0));*/
        DataSet d = MixedUtils.loadDataSet2("C:/Users/vinee_000/Downloads/Lung-tetrad_hv-14vars.txt");
        double low = .1;
        double high = .8;
        double[] initLambdas = new double[40];
        for (int i = 0; i < 40; i++) {
            initLambdas[i] = i * (high - low) / 40 + low;
        }

      //  STEPS s = new STEPS(m.getDataSet(0),initLambdas,0.05,20);
        STEPS s = new STEPS(d,initLambdas,0.01,20);
        s.b = 150;
        double [][] arr1 = s.runStepsArrayPar();
        PrintStream o = new PrintStream("parallel_stab.txt");
        for(int i = 0; i < arr1.length;i++)
        {
            for(int j = 0; j < arr1[0].length;j++)
            {
                o.print(arr1[i][j]+"\t");
            }
            o.println();
        }
        o.flush();
        o.close();
        System.out.println("Finished parallel search");
        double [][] arr = s.runStepsArray();
        PrintStream out = new PrintStream("stab.txt");
        for(int i = 0; i < arr.length;i++)
        {
            for(int j = 0; j < arr[0].length;j++)
            {
                out.print(arr[i][j]+"\t");
            }
            out.println();
        }

        out.flush();
        out.close();

    }

}
