package edu.pitt.csb.Priors;

import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;

/**
 * Created by vinee_000 on 8/31/2017.
 */
public class priorTestResults{
public static void main(String [] args) throws Exception
    {
       // String [] algs = {"mgm_priors","mgm_priors_split"};
        String [] algs = {"mgm_priors","STEPS","oracle","oracle_one","mgm_one_steps"};
        String [] type = {"CC","CD","DD","ALL"};
        int numSubsamples = 10;
        double [] priors = {0.1,0.3,0.6};
     //   int [] sampleSize = {500,1000,2000,3000};
        int [] sampleSize = {200};
        //double [] priors = {0.6};
        int numVars = 100;
        //int [] numExperts = {1,2,3,4,5};
        int [] numExperts = {5};
        int [] reliableExperts = {1,3,5};


        boolean priorInfo = false; //Compute information about each prior source
        boolean gatherWeights = false; //Compute information about weights computed by the SP MGM algorithm
        boolean unreliable = true; //Compute results for experiment with unreliable prior information sources



        //only matters for gathering weights
       int numRuns = 15;
        int rexp = 3; //number of reliable experts to use for analysis
            for(int s = 0; s < sampleSize.length;s++)
            {
                if(gatherWeights) {
                    //Priors_18_50_0.1_3_1
                    PrintStream out = new PrintStream("Prior Weights_" + numVars + "_" + sampleSize[s] + ".txt");
                    out.println("Graph\tPrior_Index\tAmount_Prior\tScore_All\tScore_CC\tScore_CD\tScore_DD\tPercent_Prior_All\tPercent_Prior_CC\tPercent_Prior_CD\tPercent_Prior_DD\tSP_MGM\tSP_MGM_F1_All\tSP_MGM_F1_CC\tSP_MGM_F1_CD\tSP_MGM_F1_DD\tSplit\tSplit_F1_All\tSplit_F1_CC\tSplit_F1_CD\tSplit_F1_DD");
                    for(int i = 0; i < numRuns;i++) { // index
                        Graph graph = GraphUtils.loadGraphTxt(new File("../Graphs/Graph_" + i + "_" + numVars + ".txt"));
                        DataSet d = MixedUtils.loadDataSet2("../Data/Data_" + i + "_" + numVars + "_" + sampleSize[s] + ".txt");

                            for(int k = 0; k < numExperts[0];k++) // which expert?
                            {

                                for(int m = 0; m < priors.length;m++) {
                                    out.print(i + "\t" + k + "\t");
                                    out.print(priors[m] + "\t");
                                    File f = new File("../Priors/Priors_" + i + "_" + numVars + "_" + priors[m] + "_" + rexp + "_" + k + ".txt");
                                    double [][] currPrior = priorTest.loadPrior(f, numVars);
                                    out.print(getScore(graph,currPrior,d,"All")[0] + "\t" + getScore(graph,currPrior,d,"CC")[0] + "\t" + getScore(graph,currPrior,d,"CD")[0] + "\t" + getScore(graph,currPrior,d,"DD")[0] + "\t" + getScore(graph,currPrior,d,"All")[1] + "\t" + getScore(graph,currPrior,d,"CC")[1] + "\t" + getScore(graph,currPrior,d,"CD")[1] + "\t" + getScore(graph,currPrior,d,"DD")[1] + "\t");

                                    for (int j = 0; j < algs.length; j++) { //which algorithm?
                                        BufferedReader b = new BufferedReader(new FileReader("../Weights/Weight_" + algs[j] + "_" + i + "_" + numVars + "_" + sampleSize[s] + "_" + priors[m] + "_" + rexp + "_" + numExperts[0] + ".txt" ));
                                        String [] line = b.readLine().replace("[","").replace("]","").split(",");
                                        out.print(line[k].trim() + "\t");
                                        b.close();

                                       b = new BufferedReader(new FileReader(algs[j] + "_" + priors[m] + "_" + rexp + "_" + numExperts[0] + "_" + numVars + "_" + sampleSize[s] + "_" + numSubsamples + ".txt"));
                                        b.readLine();
                                        for(int x = 0; x < i;x++)
                                        {
                                            b.readLine();
                                        }
                                        line = b.readLine().split("\t");
                                        out.print(line[4] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3]);
                                        if(j==0)
                                            out.print("\t");




                                }
                                out.println();
                                    out.flush();
                            }


                        }
                    }
                    out.flush();
                    out.close();
                }
                else if(!priorInfo && !unreliable) {
                    PrintStream out = new PrintStream("results_" + numVars + "_" + sampleSize[s] + ".txt");
                    out.println("Algorithm\tCC\tCC_STD\tCD\tCD_STD\tDD\tDD_STD\tALL\tALL_STD");
                    for (int j = 0; j < algs.length; j++) {
                        if (algs[j].equals("mgm_priors") || algs[j].equals("mgm_priors_split")) {
                            for(int e = 0; e < numExperts.length;e++) {
                                for (int i = 0; i < priors.length; i++) {
                                    out.print(algs[j] + "_" + priors[i] + "_" + numExperts[e] + "\t");
                                    for (int k = 0; k < type.length; k++) {
                                        BufferedReader b = new BufferedReader(new FileReader(algs[j] + "_" + priors[i] + "_" + numExperts[e] + "_" + numVars + "_" + sampleSize[s] + "_" + numSubsamples + ".txt"));
                                        b.readLine();
                                        ArrayList<Double> results = new ArrayList<Double>();
                                        while (b.ready()) {
                                            String[] line = b.readLine().split("\t");
                                            if (line[k + 1].equals("NaN"))
                                                continue;
                                            results.add(Double.parseDouble(line[k + 1]));
                                        }
                                        b.close();
                                        out.print(mean(results) + "\t");


                                        if (k == type.length - 1)
                                            out.println(std(results));
                                        else
                                            out.print(std(results) + "\t");
                                    }
                                }
                            }
                        } else {
                            int i = 0;
                            out.print(algs[j] + "\t");
                            for (int k = 0; k < type.length; k++) {
                                BufferedReader b = new BufferedReader(new FileReader(algs[j] + "_" + priors[i] + "_" + 5 + "_" + numVars + "_" + sampleSize[s] + "_" + numSubsamples + ".txt"));
                                b.readLine();
                                ArrayList<Double> results = new ArrayList<Double>();
                                while (b.ready()) {
                                    String[] line = b.readLine().split("\t");
                                    if (line[k + 1].equals("NaN"))
                                        continue;
                                    results.add(Double.parseDouble(line[k + 1]));
                                }
                                b.close();
                                out.print(mean(results) + "\t");
                                if (k == type.length - 1)
                                    out.println(std(results));
                                else
                                    out.print(std(results) + "\t");
                            }

                        }

                    }

                    out.flush();
                    out.close();
                }
                else if(unreliable)
                {
                    PrintStream out = new PrintStream("results_" + numVars + "_" + sampleSize[s]  + "_unreliable.txt");
                    out.println("Algorithm\tCC\tCC_STD\tCD\tCD_STD\tDD\tDD_STD\tALL\tALL_STD");
                    for (int j = 0; j < algs.length; j++) {
                        if (algs[j].equals("mgm_priors") || algs[j].equals("mgm_priors_split")) {
                            for(int e = 0; e < numExperts.length;e++) {
                                for (int i = 0; i < priors.length; i++) {
                                    for(int re = 0; re < reliableExperts.length;re++)
                                    {
                                        out.print(algs[j] + "_" + priors[i] + "_" + reliableExperts[re] + "\t");
                                        for (int k = 0; k < type.length; k++) {
                                            BufferedReader b;
                                            if(reliableExperts[re]==5)
                                                b = new BufferedReader(new FileReader(algs[j] + "_" + priors[i] + "_" + numExperts[e] + "_" + numVars + "_" + sampleSize[s] + "_10.txt"));
                                            else
                                                b = new BufferedReader(new FileReader(algs[j] + "_" + priors[i] + "_"+ reliableExperts[re] + "_" + numExperts[e] + "_" + numVars + "_" + sampleSize[s] + "_10.txt"));
                                            b.readLine();
                                            ArrayList<Double> results = new ArrayList<Double>();
                                            while (b.ready()) {
                                                String[] line = b.readLine().split("\t");
                                                if (line[k + 1].equals("NaN"))
                                                    continue;
                                                results.add(Double.parseDouble(line[k + 1]));
                                            }
                                            b.close();
                                            out.print(mean(results) + "\t");


                                            if (k == type.length - 1)
                                                out.println(std(results));
                                            else
                                                out.print(std(results) + "\t");
                                         }
                                    }
                                }
                            }
                        } else {
                            int i = 0;
                            out.print(algs[j] + "\t");
                            for (int k = 0; k < type.length; k++) {
                                BufferedReader b = new BufferedReader(new FileReader(algs[j] + "_" + priors[i] + "_5_" + numVars + "_" + sampleSize[s] + "_10.txt"));
                                b.readLine();
                                ArrayList<Double> results = new ArrayList<Double>();
                                while (b.ready()) {
                                    String[] line = b.readLine().split("\t");
                                    if (line[k + 1].equals("NaN"))
                                        continue;
                                    results.add(Double.parseDouble(line[k + 1]));
                                }
                                b.close();
                                out.print(mean(results) + "\t");
                                if (k == type.length - 1)
                                    out.println(std(results));
                                else
                                    out.print(std(results) + "\t");
                            }

                        }

                    }

                    out.flush();
                    out.close();
                }
                else
                {
                    int totalExperts = 5;
                    numRuns = 20;
                    PrintStream out = new PrintStream("Prior_Description_" + numVars + "_" + sampleSize[s] + ".txt");
                    PrintStream out2 = new PrintStream("Graph_Description_" + numVars + "_" + sampleSize[s] +  ".txt");
                    out.println("CC\tCD\tDD");
                    out2.println("Average_CC\tAverage_CD\tAverage_DD");
                    double [] averages = new double[3];
                    double [] gAverages = new double[3];
                    for(int k = 0; k < priors.length;k++) {
                        for (int i = 0; i < numRuns; i++) {

                            File f = new File("Graphs/Graph_" + i + "_" + numVars + ".txt");
                                Graph g = GraphUtils.loadGraphTxt(f);
                            for (int j = 0; j < totalExperts; j++) {

                                TetradMatrix currPrior;
                                f = new File("Priors/Priors_" + i + "_" + numVars + "_" + priors[k] + "_" + j + ".txt");
                                currPrior = new TetradMatrix(priorTest.loadPrior(f, numVars));

                                f = new File("Data/Data_" + i + "_" + numVars + "_" + sampleSize[s] + ".txt");
                                DataSet d = MixedUtils.loadDataSet2("Data/Data_" + i + "_" + numVars + "_" + sampleSize[s] + ".txt");
                                if(j==0 &&k ==0)
                                {
                                    for(Edge e: g.getEdges())
                                    {
                                        if(d.getVariable(e.getNode1().getName())instanceof ContinuousVariable && d.getVariable(e.getNode2().getName())instanceof ContinuousVariable)
                                            gAverages[0]++;
                                        else if(d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable)
                                            gAverages[2]++;
                                        else
                                            gAverages[1]++;
                                    }
                                }
                                for(int ii = 0; ii < currPrior.rows();ii++)
                                {
                                    for(int jj = ii + 1; jj < currPrior.columns();jj++)
                                    {
                                        if(currPrior.get(ii,jj)!=0)
                                        {
                                            if(d.getVariable(ii)instanceof ContinuousVariable && d.getVariable(jj)instanceof ContinuousVariable)
                                                averages[0]++;
                                            else if(d.getVariable(ii)instanceof DiscreteVariable && d.getVariable(jj) instanceof DiscreteVariable)
                                                averages[2]++;
                                            else
                                                averages[1]++;
                                        }
                                    }
                                }
                            }
                        }
                        int total = numRuns*totalExperts;
                        out.println(averages[0]/total + "\t" + averages[1]/total + "\t" + averages[2]/total);
                        out.println();
                    }
                    //Load in prior for each expert
                    //Load in DataSet
                    //Compute Distribution of prior across each edge type
                    //Output result in  "numVars_SampleSize_prior_description.txt"
                    //rows are CC, CD, DD
                    //columns are 10%, 30%, 60% prior info

                    out.flush();
                    out.close();
                    out2.println(gAverages[0]/20 + "\t" + gAverages[1]/20 + "\t" + gAverages[2]/20);
                    out2.flush();
                    out2.close();
                }
        }
    }
    public static double mean(ArrayList<Double> results)
    {
        double x = 0;
        for(Double d:results)
            x+= d.doubleValue();
        return x/results.size();
    }
    public static double[] getScore(Graph g, double [][] prior, DataSet d, String type)
    {
        double sum = 0;
        int count = 0;
        int totCount = 0;
        for(int i = 0; i < prior.length;i++)
        {
            for(int j = i+1; j < prior[i].length;j++)
            {
                if(prior[i][j]!=0)
                {
                    if(type.equals("All"))
                    {
                        count++;
                        sum+=prior[i][j];
                    }
                    else if(d.getVariable(i)instanceof ContinuousVariable && d.getVariable(j)instanceof ContinuousVariable && type.equals("CC"))
                    {
                        count++;
                        sum+=prior[i][j];
                    }
                    else if(d.getVariable(i)instanceof DiscreteVariable && d.getVariable(j)instanceof DiscreteVariable && type.equals("DD"))
                    {
                        count++;
                        sum+=prior[i][j];
                    }
                    else if((d.getVariable(i)instanceof DiscreteVariable != d.getVariable(j)instanceof DiscreteVariable) && type.equals("CD"))
                    {
                        count++;
                        sum+=prior[i][j];
                    }
                }
            }
        }
        for(Edge e:g.getEdges())
        {
            String i = e.getNode1().getName();
            String j = e.getNode2().getName();
            if(type.equals("All"))
            {
                totCount++;
            }
            else if(d.getVariable(i)instanceof ContinuousVariable && d.getVariable(j)instanceof ContinuousVariable && type.equals("CC"))
            {
                totCount++;
            }
            else if(d.getVariable(i)instanceof DiscreteVariable && d.getVariable(j)instanceof DiscreteVariable && type.equals("DD"))
            {
                totCount++;
            }
            else if((d.getVariable(i)instanceof DiscreteVariable != d.getVariable(j)instanceof DiscreteVariable) && type.equals("CD"))
            {
                totCount++;
            }
        }
        double [] temp = {sum/count,count/(double)totCount};
        return temp;
    }
    public static double std(ArrayList<Double>results)
    {
        double [] val = new double[results.size()];
        for(int i =0; i < val.length;i++)
            val[i] = results.get(i).doubleValue();

        return StatUtils.sd(val)/Math.sqrt(val.length);
    }
}
