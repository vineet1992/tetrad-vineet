package edu.pitt.csb.stability;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.Priors.mgmPriors;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.File;

import static edu.pitt.csb.Priors.realDataPriorTest.loadPAM50;

/**
 * Created by vinee_000 on 1/25/2018.
 */
public class quickCVTest {
    public static void main(String [] args) throws Exception
    {


        //No prior [0.573076923076923, 0.3551282051282051, 0.09358974358974359]

        //Irrelevant Prior [Lambdas no prior: [0.3808349769888231, 0.3702169625246548, 0.15785667324128863]
        //Lambdas with prior: [0.11538461538461538, 0.4870151216305062, 0.11538461538461538]

        //Relevant Prior:Lambdas no prior: [0.3730111768573307, 0.33948060486522025, 0.149474030243261]
        //Lambdas with prior: [0.09358974358974359, 0.518310322156476, 0.09358974358974359]

        if(args[0].equals("-np")) {

            DataSet data = MixedUtils.loadDataSet2("genes_with_clinical.txt");
            data = MixedUtils.completeCases(data);
            double[] lambda = {0.5948717948717949, 0.4423076923076923, 0.11538461538461538};
            CrossValidationSets cv = new CrossValidationSets(data, lambda, "Subsamples2", ".","Subtype",10, "No_Prior");
            cv.crossValidate();
        }
        else if(args[0].equals("-ip"))
        {
            DataSet data = MixedUtils.loadDataSet2("genes_with_clinical.txt");
            data = MixedUtils.completeCases(data);
            double [] lambdaNP = {0.3808349769888231,0.3702169625246548,0.15785667324128863};
            double [] lambdaWP = {0.11538461538461538,0.4870151216305062,0.11538461538461538};
            int numPriors = 5;
            TetradMatrix[] priors = new TetradMatrix[numPriors];
            for(int i = 0; i < numPriors;i++) {
                //priors[i] = new TetradMatrix(loadPrior(new File("prior_sources/Irr_Prior_" + i + ".txt"),data.getNumColumns()));
                priors[i] = new TetradMatrix(loadPAM50(new File("prior_sources/Irr_PAM50_" + i + ".txt"),data.getNumColumns()));
            }
            boolean [][] temp = new boolean[priors[0].rows()][priors[0].columns()];
            for(int i = 0; i < priors.length;i++)
            {

                TetradMatrix curr = priors[i];

                for(int j = 0; j < curr.rows();j++)
                {
                    for(int k = j +1; k < curr.columns();k++)
                    {
                        if(curr.get(j,k)!=0) {
                            temp[j][k] = true;
                        }
                    }
                }
            }
            CrossValidationSets cv = new CrossValidationSets(data,lambdaNP,lambdaWP,"Subsamples",".","Subtype",temp,10,"Irrelevant_Prior");
            cv.crossValidate();

        }
        else if(args[0].equals("-rp"))
        {
            DataSet data = MixedUtils.loadDataSet2("genes_with_clinical.txt");
            data = MixedUtils.completeCases(data);


            //Lambdas with prior: [0.09358974358974359, 0.518310322156476, 0.09358974358974359]

            double [] lambdaNP =  {0.3730111768573307, 0.33948060486522025, 0.149474030243261};
            double [] lambdaWP = {0.09358974358974359, 0.518310322156476, 0.09358974358974359};

            int numPriors = 5;
            TetradMatrix [] priors = new TetradMatrix[numPriors+1];
            for(int i = 0; i < numPriors;i++) {
                priors[i] = new TetradMatrix(loadPAM50(new File("prior_sources/Irr_PAM50_" + i + ".txt"),data.getNumColumns()));

            }
            priors[priors.length-1] = new TetradMatrix(loadPAM50(new File("prior_sources/Prior_PAM50.txt"),data.getNumColumns()));

            boolean [][] temp = new boolean[priors[0].rows()][priors[0].columns()];
            for(int i = 0; i < priors.length;i++)
            {

                TetradMatrix curr = priors[i];

                for(int j = 0; j < curr.rows();j++)
                {
                    for(int k = j +1; k < curr.columns();k++)
                    {
                        if(curr.get(j,k)!=0) {
                            temp[j][k] = true;
                        }
                    }
                }
            }
            CrossValidationSets cv = new CrossValidationSets(data,lambdaNP,lambdaWP,"Subsamples",".","Subtype",temp,10,"Relevant_Prior");
            cv.crossValidate();
        }
        else if(args[0].equals("-orp"))
        {
            DataSet data = MixedUtils.loadDataSet2("genes_with_clinical.txt");
            data = MixedUtils.completeCases(data);

//Lambdas no prior: [0.35345167652859955, 0.3299802761341223, 0.15394477317554242]
            //Lambdas with prior: [0.07179487179487179, 0.07179487179487179, 0.07179487179487179]
            double [] lambdaNP =  {0.35345167652859955, 0.3299802761341223, 0.15394477317554242};
            double [] lambdaWP = {0.07179487179487179, 0.07179487179487179, 0.07179487179487179};

            int numPriors = 1;
            TetradMatrix [] priors = new TetradMatrix[numPriors];
            priors[0] = new TetradMatrix(loadPAM50(new File("prior_sources/Prior_PAM50.txt"),data.getNumColumns()));

            boolean [][] temp = new boolean[priors[0].rows()][priors[0].columns()];
            for(int i = 0; i < priors.length;i++)
            {

                TetradMatrix curr = priors[i];

                for(int j = 0; j < curr.rows();j++)
                {
                    for(int k = j +1; k < curr.columns();k++)
                    {
                        if(curr.get(j,k)!=0) {
                            temp[j][k] = true;
                        }
                    }
                }
            }
            CrossValidationSets cv = new CrossValidationSets(data,lambdaNP,lambdaWP,"Subsamples",".","Subtype",temp,10,"Only_Relevant_Prior");
            cv.crossValidate();
        }

    }
}
