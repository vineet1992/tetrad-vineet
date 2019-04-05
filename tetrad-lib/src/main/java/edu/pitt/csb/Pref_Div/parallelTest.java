package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.data.DataSet;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

/**
 * Created by vinee_000 on 4/2/2019.
 */
public class parallelTest {


    /***TODO Change this script to use with prior parallelism to test***/
    public static void main(String [] args) throws Exception
    {
        int [] minTargetParents = new int[]{5,10,25,75,75};
        int [] numComponents = new int[]{10,25,50,300,300};
        int numRuns = 5;
        int [] numVariables = new int[]{50,200,500,3000,7000};
        PrintStream out = new PrintStream("Parallel_Results.txt");
        out.println("Run\tVariables\tTime\tParallel_Time\tErrors\tGene_Errors");
        int numParams = 20;
        int numPriors = 5;
        int numSamples = 20;
        int numReliable = 3;
        double amountPrior = 0.5;

        for(int i = 0; i < numVariables.length;i++)
        {
            for(int j = 0; j < numRuns;j++)
            {

                String[] dFile = new String[numPriors];
                for (int k = 0; k < numPriors; k++) {
                    dFile[k] = "Priors/Prior_" + numVariables[i] + "_" + minTargetParents[i] + "_" + numPriors + "_" + numReliable + "_" + amountPrior + "_true_" + numComponents[i] + "_true_false_true" + "_" + k + "_" + j + ".txt";
                  }
                File dataFile = new File("Data/Dataset_" + numVariables[i] + "_" + 200 + "_" + minTargetParents[i] + "_true_" + numComponents[i] + "_true_" + j + ".txt");

                out.print(j + "\t" + numVariables[i] + "\t");
                DataSet data = MixedUtils.loadDataSet2(dataFile.getAbsolutePath());
                int [][] samps = StabilityUtils.subSampleNoReplacement(data.getNumRows(),numSamples);

                System.out.println("**********SERIAL***********");
                PiPrefDiv4 noPrior = new PiPrefDiv4(data,"Target",minTargetParents[i],numParams);
                noPrior.setParallel(false);
                noPrior.setPartialCorrs(false);
                noPrior.setSubsamples(samps);
                long time = System.nanoTime();
                ArrayList<Gene> genes = noPrior.selectGenes(false,20,dFile);
                time = System.nanoTime()-time;
                out.print(time/Math.pow(10,9) + "\t");

                System.out.println("*********PARALLEL**********");

                PiPrefDiv4 noPrior2 = new PiPrefDiv4(data,"Target",minTargetParents[i],numParams);
                noPrior2.setParallel(true);
                noPrior2.setPartialCorrs(false);
                noPrior2.setSubsamples(samps);
                time = System.nanoTime();
                ArrayList<Gene> genes2 = noPrior2.selectGenes(false,20,dFile);
                time = System.nanoTime()-time;


                int mistakes = 0;
                for(int x = 0;x  < noPrior.constrictCorrs.length;x++)
                {
                    if(Math.abs(noPrior2.constrictCorrs[x]-noPrior.constrictCorrs[x]) > 0.00001)
                    {
                        mistakes++;
                    }
                }
                int geneMistakes = 0;
                for(int x = 0; x < genes.size();x++)
                {
                    if(!genes.get(x).symbol.equals(genes2.get(x).symbol))
                        geneMistakes++;
                }
                out.println(time/Math.pow(10,9) + "\t" + mistakes + "\t" + geneMistakes);
                out.flush();
            }
        }
        out.close();

    }
}
