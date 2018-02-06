package edu.pitt.csb.stability;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MGM_Priors;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.Arrays;

/**
 * Created by vinee_000 on 1/25/2018.
 * The purpose of this class is to generate cross validation datasets for a downstream ML Task, and to learn an MGM model on each of the training datasets for this task
 * Expected Inputs: Full Dataset, MGM Lambda parameters, directory to print Subsampled datasets to, Target Variable, (Optional) Edges with prior information, (Optional) Lambdas for edges with prior information,
 * Outputs: Learned neighbors of the target as a 2d array in a text file, All of the subsampled datasets in a separate directory
 */
public class CrossValidationSets {
    private DataSet data;
    private double [] lambda;
    private double [] lambdaWP;
    private boolean prior;
    private String directory;  //Input Directory files
    private String outDirectory; //Output directory files
    private String target;
    private boolean [][] havePrior;
    private int k; //K-Fold cross validation
    private String runName;


    public CrossValidationSets(DataSet d, double [] lambda, String dir,String dirOut, String t, int k, String runName)
    {
        prior = false;
        data = d;
        this.lambda = lambda;
        this.directory = dir;
        this.target = t;
        this.k = k;
        this.runName = runName;
        this.outDirectory = dirOut;
    }
    public CrossValidationSets(DataSet d, double [] lambda, double [] lambdaWP, String dir, String dirOut, String t, boolean [][] havePrior, int k, String runName)
    {
        this.data = d;
        this.lambda = lambda;
        this.lambdaWP = lambdaWP;
        directory = dir;
        target = t;
        this.havePrior = havePrior;
        prior = true;
        this.k = k;
        this.runName = runName;
        this.outDirectory = dirOut;
    }


    //Writes the generated features and generated subsampled sets directly to the output directory
    public void crossValidate()
    {
        try {
            File f = new File(directory);
            if (!f.exists()) {
                System.out.println("Creating directory: " + directory + " for datasets");
                f.mkdir();
            }
        }
        catch(Exception e)
        {
            System.out.println("Couldn't create subsample directory");
            return;
        }
        //Generate Cross-Validation Subsamples
        int [][] samps = StabilityUtils.generateSubsamples(k,data.getNumRows());


        //Write subsamples to file
        for(int i = 0; i < samps.length;i++)
        {

            Arrays.sort(samps[i]);
            try {
                File f = new File(directory + "/Cross_Validation_Test_" + i + ".txt");
                if(!f.exists()) {
                    System.out.println("Generating CV File " + i);
                    PrintStream out = new PrintStream(directory + "/Cross_Validation_Test_" + i + ".txt");
                    out.println(data.subsetRows(samps[i]));
                    out.flush();
                    out.close();
                    out = new PrintStream(directory + "/Cross_Validation_Train_" + i + ".txt");
                    DataSet temp = data.copy();

                    temp.removeRows(samps[i]);
                    out.println(temp);
                    out.flush();
                    out.close();

                }
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.out.println("Couldn't create CV dataset " + i);
                return;
            }
        }

        //For each subsample, run MGM with specified lambda parameter to generate Graph and identify neighbors of target
        PrintStream out;
        try{
            out = new PrintStream(outDirectory + "/Selected_Features_" + runName + ".txt");
        }
        catch(Exception e)
        {
            System.out.println("Couldn't create Selected Features file");
            return;
        }
        for(int i = 0; i < samps.length;i++)
        {
            System.out.println("Running Cross Validation: " + i);

            //Load training dataset from the file that was just created or an old file if it already existed
            DataSet train;
            try {
                train = MixedUtils.loadDataSet2(directory + "/Cross_Validation_Train_" + i + ".txt");   }
            catch(Exception e)
            {
                System.out.println("Couldn't load Training Data");
                return;
            }

            //When running MGM with prior information
            if(prior)
            {
                MGM_Priors m = new MGM_Priors(train,lambda,lambdaWP,havePrior);
                m.learnEdges(1000);
                Graph g = m.graphFromMGM();
                try{
                    for(Node n: g.getAdjacentNodes(g.getNode(target)))
                        out.print(n + "\t");
                    out.println();
                    out.flush();
                }
                catch(Exception e)
                {
                    System.out.println("Error writing to output file ");
                    return;
                }
            }
            //Otherwise
            else
            {
                MGM m = new MGM(train,lambda);
                m.learnEdges(1000);
                Graph g = m.graphFromMGM();
                try {
                    for (Node n : g.getAdjacentNodes(g.getNode(target))) {
                        out.print(n + "\t");
                    }
                    out.println();
                    out.flush();
                }
                catch(Exception e)
                {
                    System.out.println("Error Writing to output file for iteration: " + i);
                    return;
                }
            }
        }
        try {
            out.close();
        }
        catch(Exception e)
        {
            System.out.println("Error closing output file");
            return;
        }

    }

}
