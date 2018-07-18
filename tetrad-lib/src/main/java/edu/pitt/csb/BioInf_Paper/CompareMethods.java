package edu.pitt.csb.BioInf_Paper;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.*;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

/**
 * Created by vinee_000 on 7/10/2018.
 */
public class CompareMethods {
    public static void main(String [] args)throws Exception
    {
        boolean justGenerate = false;
        //int numRuns = 50;
        int numRuns = 50;
        int [] numVariables = {50,200};
        int [] numSamples = {500,100};
        int numCats = 3;
        int maxDegree = 10;
        int avgDegree = 2;
        String [] algs = {"MGMPCS","MGMCPCS"};
        //String [] algs = {"CPCS","CPSS","PCS","CG","MGMPCS","MGMCPCS"};
        //String [] algs = {"CPCS","CPSS","PCS","CG","MGMPCS","MGMCPCS","Copula"};
        //String [] algs = {"Copula"};
        //String [] algs = {"CPSS"};
        String [] types = {"CC","CD","DD","All"};
        double [] alphas = {0.001,0.01,0.05,0.1};
        //double [] lambdas = {0.2,0.28,0.4,0.57,0.8};
        double [] lambdas = {0.05,0.071,0.1,0.14,0.2};
        boolean reload = false;
        for(int i = 0; i < 1;i++)
        {
            int numMeasures = numVariables[i];
            int samples = numSamples[i];
            String directory = "";
            if(i==0)
            {
                //create directory for low dimensional dataset and set dir to this
                directory = "LD";
                File f = new File(directory);
                if(!f.exists() || !f.isDirectory())
                    f.mkdir();
            }
            else
            {
                directory = "HD";
                File f = new File(directory);
                if(!f.exists() || !f.isDirectory())
                    f.mkdir();
            }
            File f = new File(directory + "/Results");
            if(!f.exists())
                f.mkdir();
            f = new File(directory + "/Estimated_Graphs");
            if(!f.exists())
                f.mkdir();
            Parameters p = new Parameters();
            f = new File(directory + "/Graphs");
            if(!f.exists())
                f.mkdir();
            f = new File(directory + "/Data");
            if(!f.exists())
                f.mkdir();
            p.setValue("numMeasures",numMeasures);
            p.setValue("sampleSize",samples);
            p.setValue("numCategories",numCats);
            p.setValue("maxDegree",maxDegree);
            p.setValue("numEdges",numMeasures*2);

            MixedLeeHastieSimulation m = new MixedLeeHastieSimulation();
            PrintStream [] result = new PrintStream[algs.length];
            if(!justGenerate) {
                for (int x = 0; x < algs.length; x++) {
                    if(reload)
                        result[x] = new PrintStream(new FileOutputStream(directory + "/Results/"  + algs[x] + ".txt",true));
                    else {
                        result[x] = new PrintStream(directory + "/Results/" + algs[x] + ".txt");
                        result[x].println("Alpha\tLambda\tType\tRun\tAP\tAR\tAHP\tAHR\tSHD\tTime");
                    }
                    result[x].flush();
                }
            }
            RUN:for(int j = 0; j < numRuns;j++)
            {
                f = new File(directory + "/Graphs/DAG_" + j + ".txt");
                Graph trueGraph;
                DataSet data;
                Graph truePattern;
                if(f.exists())
                {
                    trueGraph = GraphUtils.loadGraphTxt(f);
                    truePattern = GraphUtils.loadGraphTxt(new File(directory + "/Graphs/Pattern_" + j + ".txt"));
                    data = MixedUtils.loadDataSet2(directory + "/Data/Data_" + j + ".txt");
                }
                else {
                    System.out.print("Simulating data for run " + j + " ...");
                    m.simulate(p);
                    trueGraph = m.getTrueGraph();
                    //DagToPattern dtp = new DagToPattern(trueGraph);
                    IndependenceTest itd = new IndTestDSep(trueGraph);
                    PcStable pcs = new PcStable(itd);
                    System.out.println("Converting to pattern");
                    truePattern = pcs.search();
                    //truePattern = dtp.convert();
                    data = m.getDataSet(0);
                    System.out.println("Done");
                    //Save everything
                    PrintStream graphOut = new PrintStream(directory + "/Graphs/DAG_" + j + ".txt");
                    graphOut.println(trueGraph);
                    graphOut.flush();
                    graphOut.close();
                    PrintStream patternOut = new PrintStream(directory + "/Graphs/Pattern_" + j + ".txt");
                    patternOut.println(truePattern);
                    patternOut.flush();
                    patternOut.close();
                    PrintStream dataOut = new PrintStream(directory + "/Data/Data_" + j + ".txt");
                    dataOut.println(data);
                    dataOut.flush();
                    dataOut.close();
                }
                if(justGenerate)
                    continue RUN;
                //Run each algorithm on the data and save the output to Estimated folder
                for(int k = 0; k < algs.length;k++)
                {
                    A:for(int a = 0; a < alphas.length;a++)
                    {
                        IndependenceTest ind = new IndTestMultinomialAJ(data,alphas[a],true);
                        for(int b = 0; b < lambdas.length;b++)
                        {
                            System.out.print("Running " + algs[k] + " with alpha " + alphas[a] + " and lambda " + lambdas[b] + " for run " + j + "...");
                            Graph est = null;
                            double [] lambda = {lambdas[b],lambdas[b],lambdas[b]};
                            long time = System.nanoTime();
                            File tempGraph = new File(directory + "/Estimated_Graphs/" + algs[k] + "_" + alphas[a] + "_" + lambdas[b] + "_" + j + ".txt");
                            if(tempGraph.exists() && reload)
                            {
                                est = GraphUtils.loadGraphTxt(tempGraph);
                            }
                            else if(algs[k].equals("CPCS"))
                            {
                                if(b!=0)
                                    continue A;
                                CpcStable cp = new CpcStable(ind);
                                est = cp.search();
                            }
                            else if(algs[k].equals("CPSS"))
                            {
                                if(b!=0)
                                    continue A;
                                double alp = 0.05;
                                double lb = 0.1;
                                if(i!=0)
                                    lb = 0.2;
                                //TODO CPC should be part of the estimation procedure I believe
                                CPSS cp = new CPSS(data,new double[]{lb,lb,lb},alphas[a]);
                                Graph init = cp.runCPSSPar();
                                IndependenceTest cpssTest = new IndTestMultinomialAJ(data,alp);
                                CpcStable cpc = new CpcStable(cpssTest);
                                if(init.getNumEdges()>0)
                                    cpc.setInitialGraph(init);
                                est = cpc.search();
                            }
                            else if(algs[k].equals("Copula"))
                            {
                                //REMEMBER COPULA CAN'T RUN WITH P > N
                                if(b!=0)
                                    continue A;
                                est = GraphUtils.loadGraphTxt(new File(directory + "/Estimated_Graphs/" + algs[k] + "_" + alphas[a] + "_" + lambdas[b] + "_" + j + ".txt"));
                            }
                            else if(algs[k].equals("PCS"))
                            {
                                if(b!=0)
                                    continue A;
                                PcStable cp = new PcStable(ind);
                                est = cp.search();
                            }
                            else if(algs[k].equals("CG"))
                            {
                                if(b!=0)
                                    continue A;
                                IndependenceTest cg = new IndTestConditionalGaussianLrt(data,alphas[a]);
                                PcStable cp = new PcStable(cg);
                                est = cp.search();
                            }
                            else if(algs[k].equals("MGMPCS"))
                            {
                                MGM mgm = new MGM(data,lambda);
                                mgm.learnEdges(1000);
                                Graph temp = mgm.graphFromMGM();
                                PcStable cp = new PcStable(ind);
                                cp.setInitialGraph(temp);
                                est = cp.search();
                            }
                            else if(algs[k].equals("MGMCPCS"))
                            {
                                MGM mgm = new MGM(data,lambda);
                                mgm.learnEdges(1000);
                                Graph temp = mgm.graphFromMGM();
                                CpcStable cp = new CpcStable(ind);
                                cp.setInitialGraph(temp);
                                est = cp.search();
                            }
                            time = System.nanoTime()-time;
                            if(!(tempGraph.exists() && reload)) {
                                PrintStream temp = new PrintStream(directory + "/Estimated_Graphs/" + algs[k] + "_" + alphas[a] + "_" + lambdas[b] + "_" + j + ".txt");
                                temp.println(est);
                                temp.flush();
                                temp.close();


                                int[][] stats = MixedUtils.allEdgeStatsBioInf(truePattern, est,data);


                                //Compute ap, ar, ahp, ahr, shd for each edge type, and time

                                //Print out to file each line (for each edge type)
                                for (int t = 0; t < types.length; t++) {
                                    result[k].println(alphas[a] + "\t" + lambdas[b] + "\t" + types[t] + "\t" + j + "\t" + (stats[t][0] / (double) (stats[t][0] + stats[t][1])) + "\t" + (stats[t][0] / (double) (stats[t][0] + stats[t][2])) + "\t" + (stats[t][3] / (double) (stats[t][3] + stats[t][4])) + "\t" + (stats[t][3] / (double) (stats[t][3] + stats[t][5])) + "\t" + stats[t][6] + "\t" + time / Math.pow(10, 9));
                                }
                                result[k].flush();
                            }
                            System.out.println("Done");
                        }
                    }

                }

            }

        }

    }
}
