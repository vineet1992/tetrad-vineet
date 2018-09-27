package edu.pitt.csb.Survival;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.pitt.csb.mgm.MGM;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 5/31/2018.
 */
public class mgmTest {
    public static void main(String [] args)throws Exception
    {
        MixedLeeHastieSimulation m = new MixedLeeHastieSimulation();
        Parameters p = new Parameters();
        p.setValue("numMeasures",6);
        p.setValue("percentDiscreteForMixedSimulation",50);
       p.setValue("numEdges", 1);

      p.setValue("sampleSize", 7);
        p.setValue("numCategories",3);
        m.simulate(p);
        PrintStream p2 = new PrintStream("test_data.txt");
        p2.println(m.getDataSet(0));
        System.out.println(m.getDataSet(0));
        MGM mg = new MGM(m.getDataSet(0),new double[]{0.2,0.4,0.5});
        mg.learnEdges(1000);
        mg.graphFromMGM();
    }
}
