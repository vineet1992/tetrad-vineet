package edu.pitt.csb.Priors;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.MGM_Priors;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by vinee_000 on 12/6/2017.
 */
public class realDataPriorTest {

    public static void main(String [] args) throws Exception
    {

        //Three Groups, No Prior = run with the data alone (just to compute the PAM50 enrichment and prediction accuracy)
        //Irrelevant Prior = run with the Irrelevant PAM50 Genesets alone (also to compute PAM50 enrichment and prediction accuracy)
        //Relevant Prior = run with Irrelevant and actual PAM50 to compare weighting scores, and to compute Prediction accuracy of this network
        //Pathway Prior = run with just all of the pathways to compute weighting


        boolean runNoPrior = false;
        boolean runIrrelevantPrior = false;
        boolean runRelevantPrior = false;
        boolean runBoth = false; //Pathway Priors
        boolean tumors = false;
        boolean computeMetrics = false;
        int numLambda = 40;
        int numPriors = 9;
        double g = 0.01;

        double low = .05;
        double high = .9;



        int index = 0;
        while(index < args.length)
        {
            if(args[index].equals("-np"))
            {
                runNoPrior = true;
                index++;
            }
            else if(args[index].equals("-ip"))
            {
                runIrrelevantPrior = true;
                index++;
            }
            else if(args[index].equals("-rp"))
            {
                runRelevantPrior = true;
                index++;
            }
            else if(args[index].equals("-pp"))
            {
                runBoth = true;
                index++;
            }
            else if(args[index].equals("-t"))
            {
                tumors = true;
                index++;
            }
            else if(args[index].equals("-nl"))
            {
                numLambda = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-g"))
            {
                g = Double.parseDouble(args[index+1]);
                index+=2;
            }
        }

        double [] lambda = new double[numLambda];
        double inc = (high-low)/(numLambda-1);
        for(int i = 0; i < numLambda;i++)
        {
            lambda[i] = i*inc + low;
        }
        DataSet data = null;
        if(tumors)
            data = MixedUtils.loadDataSet2("genes_with_clinical.txt");
        else
            data = MixedUtils.loadDataSet2("genes_with_clinical_normals.txt");
        //System.out.println(data);
        data = MixedUtils.completeCases(data);
        //System.out.println(data);
        File f = new File("Subsamples");
        if(!tumors)
            f = new File("Normal_Subsamples");
        if(!f.exists())
            f.mkdir();
        int numSub = 20;

        int b = (int) Math.floor(10 * Math.sqrt(data.getNumRows()));
        if (b > data.getNumRows())
            b = data.getNumRows() / 2;
        int[][] samps = StabilityUtils.subSampleNoReplacement(data.getNumRows(), b, numSub);
        DataSet [] subsamples = new DataSet[numSub];
        for(int i = 0; i< numSub;i++)
        {
            //File temp = new File("Subsamples/Subsample_" + i + ".txt");
            File temp = new File("Subsamples/Subsample_" + i + ".txt");
            if(!tumors)
                temp = new File("Normal_Subsamples/Normal_Subsample_" + i + ".txt");
            if(!temp.exists())
            {
                DataSet tData = data.subsetRows(samps[i]);
                PrintStream out2 = new PrintStream(temp.getAbsolutePath());
                out2.println(tData);
                out2.flush();
                out2.close();
                subsamples[i] = tData;
            }
            else
            {
                subsamples[i] = MixedUtils.loadDataSet2(temp.getAbsolutePath());
            }
        }
        //Test no Prior situation
        if(runNoPrior) {
            STEPS s = new STEPS(data, lambda, g, subsamples);
            Graph graph = s.runStepsPar();
            double[][] stab = s.stabilities;
            PrintStream out;
            if(tumors)
                out = new PrintStream("Stabilities_No_Prior.txt");
            else
                out = new PrintStream("Stabilities_No_Prior_Normals.txt");
            for (int i = 0; i < data.getNumColumns(); i++) {
                if (i == data.getNumColumns() - 1)
                    out.println(data.getVariable(i).getName());
                else
                    out.print(data.getVariable(i).getName() + "\t");
            }
            for (int i = 0; i < stab.length; i++) {
                out.print(data.getVariable(i).getName() + "\t");
                for (int j = 0; j < stab[i].length; j++) {
                    if (j == stab[i].length - 1)
                        out.println(stab[i][j]);
                    else
                        out.print(stab[i][j] + "\t");
                }
            }
            out.flush();
            out.close();
            if(tumors)
            out = new PrintStream("Graph_No_Prior.txt");
            else
                out = new PrintStream("Graph_No_Prior_Normals.txt");
            out.println(graph);
            out.flush();
            out.close();
        }

        //Test Irrelevant Prior Situation

        if(runIrrelevantPrior)
        {

            numPriors = 3;
            TetradMatrix [] priors = new TetradMatrix[numPriors];
            for(int i = 0; i < numPriors;i++) {
                //priors[i] = new TetradMatrix(loadPrior(new File("prior_sources/Irr_Prior_" + i + ".txt"),data.getNumColumns()));
                priors[i] = new TetradMatrix(loadPrior(new File("prior_sources/Irr_PAM50_" + i + ".txt"),data.getNumColumns()));
                System.out.println(priors[i].rows() + "," + priors[i].columns());
            }
            System.out.println("Constructing lambdas...");
            mgmPriors m = new mgmPriors(subsamples.length,lambda,data,priors,subsamples);
           /* PrintStream lb = new PrintStream("Irrelevant_Lambdas.txt");
            for(int i = 0; i < m.getLambdas().length;i++)
            {
                lb.println(Arrays.toString(m.getLambdas()[i]));
            }
            lb.flush();
            lb.close();*/
            System.out.println("Running Priors...");
            Graph g2 = m.runPriors();
            double [][] stab = m.edgeScores;
            double [] weights = m.expertWeights;
            PrintStream w;
            if(tumors)
                w = new PrintStream("Weights_Irrelevant_Prior.txt");
            else
                w = new PrintStream("Weights_Irrelevant_Prior_Normals.txt");
            for(int i = 0; i < weights.length;i++)
            {
                w.println("Irr_PAM50_" + i  + "\t" + weights[i]);
            }
            w.flush();
            w.close();
            PrintStream out;
            if(tumors)
                out = new PrintStream("Stabilities_Irrelevant_Prior.txt");
            else
                out = new PrintStream("Stabilities_Irrelevant_Prior_Normals.txt");
            for (int i = 0; i < data.getNumColumns(); i++) {
                if (i == data.getNumColumns() - 1)
                    out.println(data.getVariable(i).getName());
                else
                    out.print(data.getVariable(i).getName() + "\t");
            }
            for (int i = 0; i < stab.length; i++) {
                out.print(data.getVariable(i).getName() + "\t");
                for (int j = 0; j < stab[i].length; j++) {
                    if (j == stab[i].length - 1)
                        out.println(stab[i][j]);
                    else
                        out.print(stab[i][j] + "\t");
                }
            }
            out.flush();
            out.close();
            if(tumors)
            out = new PrintStream("Graph_Irrelevant_Prior.txt");
            else
                out = new PrintStream("Graph_Irrelevant_Prior_Normals.txt");
            out.println(g2);
            out.flush();
            out.close();


        }




        //Test Relevant PAM 50 Situation
        if(runRelevantPrior)
        {

            numPriors = 3;
            TetradMatrix [] priors = new TetradMatrix[numPriors+1];
            for(int i = 0; i < numPriors;i++) {
                priors[i] = new TetradMatrix(loadPrior(new File("prior_sources/Irr_PAM50_" + i + ".txt"),data.getNumColumns()));

            }
            priors[priors.length-1] = new TetradMatrix(loadPrior(new File("prior_sources/Prior_PAM50.txt"),data.getNumColumns()));

            mgmPriors m = new mgmPriors(subsamples.length,lambda,data,priors,subsamples);
            PrintStream lb;
            if(tumors)
                lb = new PrintStream("Relevant_Lambdas.txt");
            else
                lb = new PrintStream("Relevant_Lambdas_Normals.txt");
            for(int i = 0; i < m.getLambdas().length;i++)
            {
                lb.println(Arrays.toString(m.getLambdas()[i]));
            }
            lb.flush();
            lb.close();
            System.out.println("Running Priors...");
            Graph g2 = m.runPriors();
            double [][] stab = m.edgeScores;
            double [] weights = m.expertWeights;
            PrintStream w;
            if(tumors)
                 w = new PrintStream("Weights_Relevant_Prior.txt");
            else
                w = new PrintStream("Weights_Relevant_Prior_Normals.txt");
            for(int i = 0; i < weights.length;i++)
            {
                if(i==weights.length-1)
                    w.println("PAM50\t" + weights[i]);
                else
                    w.println("Irr_PAM50_" + i + "\t" + weights[i]);
            }
            w.flush();
            w.close();
            PrintStream out;
            if(tumors)
                out = new PrintStream("Stabilities_Relevant_Prior.txt");
            else
                 out = new PrintStream("Stabilities_Relevant_Prior_Normals.txt");
            for (int i = 0; i < data.getNumColumns(); i++) {
                if (i == data.getNumColumns() - 1)
                    out.println(data.getVariable(i).getName());
                else
                    out.print(data.getVariable(i).getName() + "\t");
            }
            for (int i = 0; i < stab.length; i++) {
                out.print(data.getVariable(i).getName() + "\t");
                for (int j = 0; j < stab[i].length; j++) {
                    if (j == stab[i].length - 1)
                        out.println(stab[i][j]);
                    else
                        out.print(stab[i][j] + "\t");
                }
            }
            out.flush();
            out.close();
            if(tumors)
                 out = new PrintStream("Graph_Relevant_Prior.txt");
            else
                out = new PrintStream("Graph_Relevant_Prior_Normals.txt");
            out.println(g2);
            out.flush();
            out.close();


        }

        if(runBoth)
        {

            f = new File("prior_sources");
            File [] stuff = f.listFiles();
            ArrayList<File> files = new ArrayList<File>();
            for(int i = 0; i < stuff.length;i++)
            {
                if(!stuff[i].getName().contains("PAM50"))
                    files.add(stuff[i]);
            }
            TetradMatrix [] priors = new TetradMatrix[files.size()];
            for(int i = 0; i < files.size();i++) {
                priors[i] = new TetradMatrix(loadPrior(files.get(i),data.getNumColumns()));
            }
            mgmPriors m = new mgmPriors(subsamples.length,lambda,data,priors,subsamples);
            PrintStream lb;
            if(tumors)
                lb = new PrintStream("Pathway_Lambdas.txt");
            else
                lb = new PrintStream("Pathway_Lambdas_Normals.txt");
            for(int i = 0; i < m.getLambdas().length;i++)
            {
                lb.println(Arrays.toString(m.getLambdas()[i]));
            }
            lb.flush();
            lb.close();


            Graph g2 = m.runPriors();
            double [][] stab = m.edgeScores;
            double [] weights = m.expertWeights;
            PrintStream w;
            if(tumors)
                 w = new PrintStream("Weights_Pathway_Prior.txt");
            else
                 w = new PrintStream("Weights_Pathway_Prior_Normals.txt");
            for(int i = 0; i < weights.length;i++)
            {
                w.println(files.get(i).getName().replace(".txt","") + "\t" + weights[i]);
            }
            w.flush();
            w.close();
            PrintStream out;
            if(tumors)
                out = new PrintStream("Stabilities_Pathway_Prior.txt");
            else
                out = new PrintStream("Stabilities_Pathway_Prior_Normals.txt");
            for (int i = 0; i < data.getNumColumns(); i++) {
                if (i == data.getNumColumns() - 1)
                    out.println(data.getVariable(i).getName());
                else
                    out.print(data.getVariable(i).getName() + "\t");
            }
            for (int i = 0; i < stab.length; i++) {
                out.print(data.getVariable(i).getName() + "\t");
                for (int j = 0; j < stab[i].length; j++) {
                    if (j == stab[i].length - 1)
                        out.println(stab[i][j]);
                    else
                        out.print(stab[i][j] + "\t");
                }
            }
            out.flush();
            out.close();
            if(tumors)
                out = new PrintStream("Graph_Pathway_Prior.txt");
            else
                out = new PrintStream("Graph_Pathway_Prior_Normals.txt");
            out.println(g2);
            out.flush();
            out.close();


        }


        //Output metrics, precision and recall of PAM50 Genes (Check the neighborhood of the Subtype variable)
        //, Recovery of relevant pathway connections
    }
    public static double [][] loadPrior(File data, int numVariables) throws Exception
    {
        double [][] temp = new double[numVariables][numVariables];
        BufferedReader b = new BufferedReader(new FileReader(data));
        b.readLine();
        for(int i = 0; i < numVariables;i++)
        {
            String [] line = b.readLine().split("\t");
            for(int j = 0; j < numVariables;j++)
            {
                temp[i][j] = Double.parseDouble(line[j]);
            }
        }
        return temp;
    }
}
