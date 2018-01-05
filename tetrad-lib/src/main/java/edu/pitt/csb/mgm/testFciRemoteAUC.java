package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.performance.PerformanceTests;
import edu.cmu.tetrad.search.*;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

/**
 * Created by vinee_000 on 1/6/2016.
 */
public class testFciRemoteAUC {


    /**
     * Created by vinee_000 on 12/29/2015.
     */

    //Script to test algorithms for latent variable causal inference!

    public static void main(String[] args) throws Exception {

        //Make results folder if it doesn't exist
        boolean genLatents = false;
        boolean sf = false;
        boolean hard = true;
        ArrayList<String> algs = new ArrayList<String>();
        algs.add("FCI");
        algs.add("FCI-MAX");
        algs.add("MGM-FCI");
        algs.add("MGM-FCI-MAX");
        algs.add("CFCI");
        algs.add("MGM-CFCI");
        int numGraphs = 10;
        int numVars = 50;
        int sampleSize = 1000;
        int numLatents = 5;
        int [] iarr = {10,4,2,1};
        double [] alp = {.001,.01,.05,.1};
        PrintStream [] times = new PrintStream[algs.size()];

        if(numVars==500)
        {
            iarr = new int[3];
            iarr[0] = 10;
            iarr[1] = 6;
            iarr[2] = 4;
            alp = new double[4];
            alp[0] = .001;
            alp[1] = .01;
            alp[2] = .03;
            alp[3] = .05;
        }
        File dir = new File("Output_Graphs_" + numVars + "_Nodes_" + sampleSize + "_Samples" ); //directory for produced graphs

        if(sampleSize==1000)
            dir = new File("Output_Graphs_" + numVars + "_Nodes");
        if(sf)
            dir = new File(dir.getPath() + "_SF");
        if(hard)
            dir = new File(dir.getPath() + "_hard");
       if(!dir.exists())
            dir.mkdir();

        if(genLatents) {
            generatePAGS(numGraphs, numLatents, dir, sf,hard, numVars);
            System.exit(0);
        }


        for(int j = 0; j < times.length;j++)
        {
            if(sf)
                times[j] = new PrintStream(dir + "/" + algs.get(j) + "_" + numVars + "_Nodes_" + sampleSize + "_samples_SF_time.txt");
            else
                times[j] = new PrintStream(dir + "/" + algs.get(j) + "_" + numVars + "_Nodes_" + sampleSize + "_samples_time.txt");

            times[j].println("Lambda\tAlpha\tGraph\tTime\tInd_Tests");
            times[j].flush();
        }





        for(int ii = 0; ii < iarr.length;ii++)
        {
        int i = iarr[ii]; //10, 4, 2, 1

        double lab = 0.05 * (i + 1);

        ArrayList[] toRemove = new ArrayList [numGraphs];
        for (int j = 0; j < alp.length; j++) {

            double alpha = alp[j];
            /*PrintWriter output_CC = new PrintWriter(new BufferedWriter(new FileWriter("CC_output_" + i + "_" + j + ".txt")));
            PrintWriter output_CD = new PrintWriter(new BufferedWriter(new FileWriter("CD_output_" + i + "_" + j + ".txt")));
            PrintWriter output_DD = new PrintWriter(new BufferedWriter(new FileWriter("DD_output_" + i + "_" + j + ".txt")));

            output_CC.println("Algorithm\tLambda\tAlpha\tGraph_Num\tTP\tFP\tFN\tTime\tOrientation TP\tOrientation FP\tOrientation FN\tBidirected TP\tBidirected FP\t Bidirected FN");
            output_CC.println();
            output_CD.println("Algorithm\tLambda\tAlpha\tGraph_Num\tTP\tFP\tFN\tTime\tOrientation TP\tOrientation FP\tOrientation FN\tBidirected TP\tBidirected FP\tBidirected FN");
            output_CD.println();
            output_DD.println("Algorithm\tLambda\tAlpha\tGraph_Num\tTP\tFP\tFN\tTime\tOrientation TP\tOrientation FP\tOrientation FN\tBidirected TP\tBidirected FP\tBidirected FN");
            output_DD.println();*/
            for (int fe = 0; fe < numGraphs; fe++) {
                System.out.println(i + "," + fe);
                String filename = "Data_" + numVars + "_Nodes/DAG_" + fe + "_data";
                if(sf)
                    filename = "Data_" + numVars + "_Nodes_SF/SF_" + fe + "_data";
                String filename2 = "Graphs_" + numVars + "_Nodes/DAG_" + fe + ".txt";
                if(sf)
                    filename2 = "Graphs_" + numVars + "_Nodes_SF/SF_" + fe + ".txt";
                String extension = ".txt";
                Graph trueGraph = GraphUtils.loadGraphTxt(new File(filename2));
                DataSet data = MixedUtils.loadDataSet2(filename + extension);
                boolean done = false;
                double[] lambda = {lab, lab, lab};
                DataSet removeStuff = data.copy();

                A:
                while (!done) {
                    removeStuff = data.copy();
                    if (sampleSize >= removeStuff.getNumRows())
                        break A;
                    int[] removal = new int[removeStuff.getNumRows() - sampleSize];
                    for (int las = 0; las < removeStuff.getNumRows() - sampleSize; las++)
                        removal[las] = las;
                    removeStuff.removeRows(removal);
                    try {
                        MGM m = new MGM(removeStuff, lambda);
                    } catch (Exception e) {
                        continue;
                    }
                    done = true;
                }
                Graph p;
                Graph PAG = GraphUtils.loadGraphTxt(new File(dir + "/PAG_10_" + fe + ".txt"));
                for (int index = 1; index < numVars; index++) {
                    if (PAG.getNode("X" + index) == null) {
                        trueGraph.getNode("X" + index).setNodeType(NodeType.LATENT);
                        removeStuff.removeColumn(removeStuff.getVariable("X" + index));
                    }
                }
                final DagToPag dagToPag = new DagToPag(trueGraph);
                dagToPag.setCompleteRuleSetUsed(false);
                p = dagToPag.convert();
                double tolerance = 1e-7; //convergeance tolerance
                int iterLimit = 1000; //iteration limit
                MGM model = new MGM(removeStuff, lambda);
                long mgmTime = 0;
                mgmTime = System.nanoTime();
                model.learnEdges(iterLimit);

                Graph mgmGraph = model.graphFromMGM();
                mgmTime = System.nanoTime() - mgmTime;


                IndTestMultinomialAJ indTest = new IndTestMultinomialAJ(removeStuff.copy(), alpha);
                PrintStream graphWriter;

                if(algs.contains("MGM-FCI")) {
                     graphWriter = new PrintStream(dir + "/mgmfci_" + i + "_" + j + "_" + fe + ".txt");
                    Fci f = new Fci(indTest);
                    f.setInitialGraph(mgmGraph);
                    long mgmFciTime = System.nanoTime();
                    Graph gFinal = f.search();

                    times[algs.indexOf("MGM-FCI")].println(i + "\t" + j + "\t" + fe + "\t" + (System.nanoTime() - mgmFciTime + mgmTime) + "\t" + indTest.timesCalled);
                    times[algs.indexOf("MGM-FCI")].flush();
                    graphWriter.println(gFinal);
                    graphWriter.flush();
                    graphWriter.close();

                }

                if(algs.contains("FCI")&&i==10) {
                    indTest = new IndTestMultinomialAJ(removeStuff.copy(), alpha);

                    Fci f2 = new Fci(indTest);


                    long fciTime = System.nanoTime();
                    Graph gFinal2 = f2.search();
                    fciTime = System.nanoTime() - fciTime;
                    times[algs.indexOf("FCI")].println(i + "\t" + j + "\t" + fe + "\t" + fciTime + "\t" + indTest.timesCalled);
                    times[algs.indexOf("FCI")].flush();
                    graphWriter = new PrintStream(dir + "/fci_" + i + "_" + j + "_" + fe + ".txt");
                    graphWriter.println(gFinal2);
                    graphWriter.flush();
                    //double [][] stats_fci = MixedUtils.allEdgeStatsLatentNew(p, gFinal2,trueGraph,latents);

                    graphWriter.close();
                }

                if(algs.contains("FCI-MAX")&&i==10) {
                    indTest = new IndTestMultinomialAJ(removeStuff.copy(), alpha);
                    FciMaxP f3 = new FciMaxP(indTest);


                    long maxTime = System.nanoTime();

                    Graph gFinal3 = f3.search();
                    maxTime = System.nanoTime() - maxTime;
                    times[algs.indexOf("FCI-MAX")].println(i + "\t" + j + "\t" + fe + "\t" + (maxTime) + "\t" + indTest.timesCalled);
                    graphWriter = new PrintStream(dir + "/max_" + i + "_" + j + "_" + fe + ".txt");
                    graphWriter.println(gFinal3);
                    graphWriter.flush();
                    graphWriter.close();

                }
                if(algs.contains("MGM-FCI-MAX")) {
                    indTest = new IndTestMultinomialAJ(removeStuff.copy(), alpha);
                    FciMaxP f3 = new FciMaxP(indTest);
                    f3.setInitialGraph(mgmGraph);
                    graphWriter = new PrintStream(dir + "/mfm_" + i + "_" + j + "_" + fe + ".txt");
                    long mgmmax = System.nanoTime();
                    Graph gFinal4 = f3.search();
                    mgmmax = System.nanoTime() - mgmmax;
                    times[algs.indexOf("MGM-FCI-MAX")].println(i + "\t" + j + "\t" + fe + "\t" + (mgmmax+mgmTime) + "\t" + indTest.timesCalled);
                    graphWriter.println(gFinal4);
                    graphWriter.flush();
                    graphWriter.close();
                }
                //double [][] stats_mfm = MixedUtils.allEdgeStatsLatentNew(p, gFinal4,trueGraph,latents);
                if(algs.contains("CFCI")&&i==10) {
                    indTest = new IndTestMultinomialAJ(removeStuff.copy(), alpha);

                    Cfci f5 = new Cfci(indTest);
                    long time = System.nanoTime();
                    Graph gfinal5 = f5.search();
                    time = System.nanoTime() - time;
                    times[algs.indexOf("CFCI")].println(i + "\t" + j + "\t" + fe + "\t" + (time) + "\t" + indTest.timesCalled);
                    graphWriter = new PrintStream(dir + "/cfci_" + i + "_" + j + "_" + fe + ".txt");
                    graphWriter.println(gfinal5);
                    graphWriter.flush();
                    graphWriter.close();
                }

                if(algs.contains("MGM-CFCI")) {

                    indTest = new IndTestMultinomialAJ(removeStuff.copy(), alpha);

                    Cfci f6 = new Cfci(indTest);
                    f6.setInitialGraph(mgmGraph);
                    long time = System.nanoTime();
                    Graph gfinal6 = f6.search();
                    time = System.nanoTime()-time;
                    times[algs.indexOf("MGM-CFCI")].println(i + "\t" + j + "\t" + fe + "\t" + (time + mgmTime) + "\t" + indTest.timesCalled);
                    graphWriter = new PrintStream(dir + "/mgmcfci_" + i + "_" + j + "_" + fe + ".txt");
                    graphWriter.println(gfinal6);
                    graphWriter.flush();
                    graphWriter.close();

                }
               /* int TPU = 0;
                int FPU = 1;
                int FNU = 2;
                int OTP = 3;
                int OFP = 4;
                int OFN = 5;
                int LTP = 6;
                int LFP = 7;
                int LFN = 8;
                DecimalFormat df = new DecimalFormat("#.##");
                output_CC.print("MGM-FCI\t" + df.format(.05 * (i + 1)) + "\t" + df.format(alp[j]) + "\t" + fe + "\t" + stats_mgm[0][TPU] + "\t" + stats_mgm[0][FPU] + "\t" + stats_mgm[0][FNU]);
                output_CC.println("\t" + mgm_time[i][j][fe] + "\t" + stats_mgm[0][OTP] + "\t" + stats_mgm[0][OFP] + "\t" + stats_mgm[0][OFN] + "\t" + stats_mgm[0][LTP] + "\t" + stats_mgm[0][LFP] + "\t" + stats_mgm[0][LFN]);
                output_CD.print("MGM-FCI\t" + df.format(.05 * (i + 1)) + "\t" + df.format(alp[j]) + "\t" + fe + "\t" + stats_mgm[1][TPU] + "\t" + stats_mgm[1][FPU] + "\t" + stats_mgm[1][FNU]);
                output_CD.println("\t" + mgm_time[i][j][fe] + "\t" + stats_mgm[1][OTP] + "\t" + stats_mgm[1][OFP] + "\t" + stats_mgm[1][OFN] + "\t" + stats_mgm[1][LTP] + "\t" + stats_mgm[1][LFP] + "\t" + stats_mgm[1][LFN]);
                output_DD.print("MGM-FCI\t" + df.format(.05 * (i + 1)) + "\t" + df.format(alp[j]) + "\t" + fe + "\t" + stats_mgm[2][TPU] + "\t" + stats_mgm[2][FPU] + "\t" + stats_mgm[2][FNU]);
                output_DD.println("\t" + mgm_time[i][j][fe] + "\t" + stats_mgm[2][OTP] + "\t" + stats_mgm[2][OFP] + "\t" + stats_mgm[2][OFN] + "\t" + stats_mgm[2][LTP] + "\t" + stats_mgm[2][LFP] + "\t" + stats_mgm[2][LFN]);




                output_CC.print("FCI\t" + df.format(.05 * (i + 1)) + "\t" + df.format(alp[j]) + "\t" + fe + "\t" + stats_fci[0][TPU] + "\t" + stats_fci[0][FPU] + "\t" + stats_fci[0][FNU]);
                output_CC.println("\t" + fci_time[i][j][fe] + "\t" + stats_fci[0][OTP] + "\t" + stats_fci[0][OFP] + "\t" + stats_fci[0][OFN] + "\t" + stats_fci[0][LTP] + "\t" + stats_fci[0][LFP] + "\t" + stats_fci[0][LFN]);
                output_CD.print("FCI" + "\t" + df.format(.05 * (i + 1)) + "\t" + df.format(alp[j]) + "\t" + fe + "\t" + stats_fci[1][TPU] + "\t" + stats_fci[1][FPU] + "\t" + stats_fci[1][FNU]);
                output_CD.println("\t" + fci_time[i][j][fe] + "\t" + stats_fci[1][OTP] + "\t" + stats_fci[1][OFP] + "\t" + stats_fci[1][OFN] + "\t" + stats_fci[1][LTP] + "\t" + stats_fci[1][LFP] + "\t" + stats_fci[1][LFN]);
                output_DD.print("FCI\t" + df.format(.05 * (i + 1)) + "\t" + df.format(alp[j]) + "\t" + fe + "\t" + stats_fci[2][TPU] + "\t" + stats_fci[2][FPU] + "\t" + stats_fci[2][FNU]);
                output_DD.println("\t" + fci_time[i][j][fe] + "\t" + stats_fci[2][OTP] + "\t" + stats_fci[2][OFP] + "\t" + stats_fci[2][OFN] + "\t" + stats_fci[2][LTP] + "\t" + stats_fci[2][LFP] + "\t" + stats_fci[2][LFN]);



                output_CC.print("MFM\t" + df.format(.05 * (i + 1)) + "\t" + df.format(alp[j]) + "\t" + fe + "\t" + stats_mfm[0][TPU] + "\t" + stats_mfm[0][FPU] + "\t" + stats_mfm[0][FNU]);
                output_CC.println("\t" + mfm_time[i][j][fe] + "\t" + stats_mfm[0][OTP] + "\t" + stats_mfm[0][OFP] + "\t" + stats_mfm[0][OFN] + "\t" + stats_mfm[0][LTP] + "\t" + stats_mfm[0][LFP] + "\t" + stats_mfm[0][LFN]);
                output_CD.print("MFM\t" + df.format(.05 * (i + 1)) + "\t" + df.format(alp[j]) + "\t" + fe + "\t" + stats_mfm[1][TPU] + "\t" + stats_mfm[1][FPU] + "\t" + stats_mfm[1][FNU]);
                output_CD.println("\t" + mfm_time[i][j][fe] + "\t" + stats_mfm[1][OTP] + "\t" + stats_mfm[1][OFP] + "\t" + stats_mfm[1][OFN] + "\t" + stats_mfm[1][LTP] + "\t" + stats_mfm[1][LFP] + "\t" + stats_mfm[1][LFN]);
                output_DD.print("MFM\t" + df.format(.05 * (i + 1)) + "\t" + df.format(alp[j]) + "\t" + fe + "\t" + stats_mfm[2][TPU] + "\t" + stats_mfm[2][FPU] + "\t" + stats_mfm[2][FNU]);
                output_DD.println("\t" + mfm_time[i][j][fe] + "\t" + stats_mfm[2][OTP] + "\t" + stats_mfm[2][OFP] + "\t" + stats_mfm[2][OFN] + "\t" + stats_mfm[2][LTP] + "\t" + stats_mfm[2][LFP] + "\t" + stats_mfm[2][LFN]);
*/
                //need to edit the close settings, and play with the all edge stats method
            }
          /*  output_CC.flush();
            output_CC.close();
            output_CD.flush();
            output_CD.close();
            output_DD.flush();
            output_DD.close();*/
        }
        }
        for(int j = 0; j < times.length;j++)
        {
            times[j].flush();
            times[j].close();
        }
    }
    public static void generatePAGS(int numGraphs, int numLatents, File dir,boolean sf,boolean hard,int numVars) throws Exception
    {

        for (int fe = 0; fe < numGraphs; fe++) {
        String filename2 = "Graphs_" + numVars + "_Nodes/DAG_" + fe + ".txt";
        if(sf)
            filename2 = "Graphs_" + numVars + "_Nodes_SF/SF_" + fe + ".txt";
        String extension = ".txt";
            String filename = "Data_" + numVars + "_Nodes/DAG_" + fe + "_data";
            if(sf)
                filename = "Data_" + numVars + "_Nodes_SF/SF_" + fe + "_data";
            if(hard)
                filename += "_hard";
            Graph trueGraph = GraphUtils.loadGraphTxt(new File(filename2));
            DataSet data = MixedUtils.loadDataSet2(filename + extension);
            boolean done = false;

            Graph p;
            int numNodesToRemove = numLatents;

                int totalNodes = trueGraph.getNumNodes();
                int numNodesRemoved = 0;
                Random rand = new Random();
                int count = 0;
                ArrayList<Integer> usedInts = new ArrayList<Integer>();
                ArrayList<Integer> removeInts = new ArrayList<Integer>();
                boolean reset = false;
                A:
                while (numNodesToRemove > 0) {
                    count++;
                    if (count == 1000) {
                        count = 0;
                        System.out.println("Stuck on file " + fe );
                        reset = true;
                        break A;
                    }

                    int nod = rand.nextInt(totalNodes);
                    if (usedInts.contains(nod))
                        continue;
                    else
                        usedInts.add(nod);
                    if (trueGraph.getOutdegree(trueGraph.getNode("X" + (nod + 1))) > 1) {

                        numNodesToRemove--;
                        numNodesRemoved++;

                        removeInts.add(nod);
                        totalNodes--; //not sure if i need this

                    }

                }
                if (reset)
                    System.exit(0);
                //  continue;
                Collections.sort(removeInts);
                Collections.reverse(removeInts);

                ArrayList<Node> latents = new ArrayList<Node>();
                for (int curr = 0; curr < removeInts.size(); curr++) {
                    latents.add(trueGraph.getNode("X" + (removeInts.get(curr) + 1)));
                    trueGraph.getNode("X" + (removeInts.get(curr) + 1)).setNodeType(NodeType.LATENT);
                }
                usedInts.clear();

                final DagToPag dagToPag = new DagToPag(trueGraph);
                dagToPag.setCompleteRuleSetUsed(false);
                p = dagToPag.convert();
                PrintStream out = new PrintStream(dir + "/PAG_10_" + fe + ".txt");
                if(hard)
                    out = new PrintStream(dir + "/PAG_10_" + fe + "_hard.txt");
                out.println(p);
                out.flush();
                out.close();
            }
    }

}
