package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;

/**
 * Created by vinee_000 on 7/6/2016.
 */
public class mfmResults {
    public static void main(String[] args) throws Exception {
        String[] algs = {"mgmfci","mfm","fci","cfci","mgmcfci","max"};
        boolean sf = false;
        boolean includeBCCD = false;
        boolean hard = true;
        String graphFolder = "Graphs_50_Nodes";
        String dataFolder = "Data_50_Nodes";
        String outputFolder = "Output_Graphs_50_Nodes";

        if(hard)
        {
            outputFolder+="_hard";
        }
        if(sf)
        {
            graphFolder+="_SF";
            dataFolder+="_SF";
            outputFolder+="_SF";
        }
        int numVars = 50;
        int numGraphs = 10;
        int numLatents = 5;
        int[] iarr = {1, 2, 4};
        double[] alp = {.001, .01, .05, .1};
      //  int [] iarr = {4,6};
       // double [] alp = {.001,.01,.03,.05};
        PrintStream[] streams = new PrintStream[algs.length];
        for (int i = 0; i < algs.length; i++) {
            PrintStream temp = new PrintStream(outputFolder + "/" + algs[i] + "_results.txt");

            streams[i] = temp;
        }


        String header = "Alpha\tLambda\tGraph\tType\tAdjacency TP\tAdjacency FP\tAdjacency FN\tExperimental TP\tExperimental TN\tExperimental FP\tExperimental FN\tLatent TP\tLatent FP\tLatent FN";
        int TPU = 0;
        int FPU = 1;
        int FNU = 2;
        int ETP = 3;
        int EFP = 4;
        int EFN = 5;
        int ETN = 6;
        int ITP = 7;
        int IFP = 8;
        int IFN = 9;
        int ITN = 10;
        int LTP = 11;
        int LFP = 12;
        int LFN = 13;
        for (PrintStream p : streams)
            p.println(header);
        for (int i = 0; i < iarr.length; i++) {
            for (int j = 0; j < alp.length; j++) {
                for (int fe = 0; fe < numGraphs; fe++) {
                    int ind = iarr[i];
                    // Graph trueDAG = GraphUtils.loadGraphTxt(new File("DAG_" + fe + "_graph.txt"));
                    Graph trueDAG;
                    if(sf)
                        trueDAG = GraphUtils.loadGraphTxt(new File(graphFolder + "/SF_" + fe + ".txt"));
                    else
                        trueDAG = GraphUtils.loadGraphTxt(new File(graphFolder + "/DAG_" + fe + ".txt"));
                    Graph truePAG = GraphUtils.loadGraphTxt(new File(outputFolder + "/PAG_10_" + fe + ".txt"));
                    ArrayList<Node> latents = new ArrayList<Node>();
                    for (Node n : trueDAG.getNodes()) {
                        if (truePAG.getNode(n.getName()) == null)
                            latents.add(n);
                    }
                    Graph[] resultPAGS = new Graph[algs.length];
                    for (int al = 0; al < algs.length; al++) {
                        if (algs[al].contains("mfm") || algs[al].contains("mgm"))
                            resultPAGS[al] = GraphUtils.loadGraphTxt(new File(outputFolder + "/" + algs[al] + "_" + ind + "_" + j + "_" + fe + ".txt"));
                        else
                            resultPAGS[al] = GraphUtils.loadGraphTxt(new File(outputFolder + "/" + algs[al] + "_10_" + j + "_" + fe + ".txt"));
                    }


                    // DataSet data = MixedUtils.loadDataSet2("SF_" + fe + "_data_NonMono.txt");
                    DataSet data;
                    if(!sf)
                        data = MixedUtils.loadDataSet2(dataFolder + "/DAG_" + fe + "_data.txt");
                    else
                        data = MixedUtils.loadDataSet2(dataFolder + "/SF_" + fe + "_data.txt");
                    ArrayList<double[][]> stats = new ArrayList<double[][]>();
                    for (int al = 0; al < algs.length; al++) {
                        stats.add(MixedUtils.allEdgeStatsLatentNew(truePAG, resultPAGS[al], trueDAG, latents, data));
                    }

                    double lab = 0.05 * (ind + 1);
                    for (int type = 0; type < 3; type++) {
                        if (type == 0) {
                            for (PrintStream x : streams)
                                x.print(alp[j] + "\t" + lab + "\t" + fe + "\tCC\t");
                        } else if (type == 1) {
                            for (PrintStream x : streams)
                                x.print(alp[j] + "\t" + lab + "\t" + fe + "\tCD\t");
                        } else {
                            for (PrintStream x : streams)
                                x.print(alp[j] + "\t" + lab + "\t" + fe + "\tDD\t");
                        }

                        for (int k = 0; k < algs.length; k++) {
                            streams[k].println(stats.get(k)[type][TPU] + "\t" + stats.get(k)[type][FPU] + "\t" + stats.get(k)[type][FNU] + "\t" + stats.get(k)[type][ETP] + "\t" + stats.get(k)[type][ETN] + "\t" + stats.get(k)[type][EFP] + "\t" + stats.get(k)[type][EFN] + "\t" + stats.get(k)[type][LTP] + "\t" + stats.get(k)[type][LFP] + "\t" + stats.get(k)[type][LFN]);
                            streams[k].flush();
                        }

                    }

                }
            }
        }

        for (int k = 0; k < algs.length; k++)
            streams[k].close();


        if (includeBCCD) {
            int[][] latents = new int[numGraphs][numLatents];
            BufferedReader b = new BufferedReader(new FileReader(outputFolder + "/latents.txt"));
            int graph = 0;
            while (b.ready()) {
                String[] line = b.readLine().split("\t");
                for (int i = 0; i < line.length; i++) {
                    latents[graph][i] = Integer.parseInt(line[i]);
                }
                graph++;
            }
            double[] p_loc = {.25, .5, .75};
            double[] p_edge = {.25, .5, .75};
            header = "Loc_Probability\tEdge_Probability\tGraph Number\tType\tAdjacency TP\tAdjacency FP\t Adjacency FN\tOrientation TP\tOrientation FP\tOrientation FN";
            PrintStream out = new PrintStream(outputFolder + "/bccd_results.txt");
            out.println(header);



                for (int i = 0; i < p_loc.length; i++) {
                    for (int j = 0; j < p_edge.length; j++) {
                        for (int g = 0; g < numGraphs; g++) {
                        b = new BufferedReader(new FileReader(outputFolder + "/BCCD_" + g + "_" + (i + 1) + "_" + (j + 1) + ".txt"));
                        int[][] data = new int[numVars - numLatents][numVars - numLatents];
                        int row = 0;
                        while (b.ready()) {
                            String[] stuff = b.readLine().split("\t");
                            for (int col = 0; col < stuff.length; col++) {
                                data[row][col] = Integer.parseInt(stuff[col]);
                            }

                            row++;
                        }
                        Graph pag = convertToGraph(data, latents[g]);

                        Graph truePAG = GraphUtils.loadGraphTxt(new File(outputFolder + "/PAG_10_" + g + ".txt"));
                        Graph trueDAG;
                        if(sf)
                            trueDAG = GraphUtils.loadGraphTxt(new File(graphFolder + "/SF_" + g + ".txt"));
                        else
                            trueDAG = GraphUtils.loadGraphTxt(new File(graphFolder + "/DAG_" + g + ".txt"));
                        if(i==0&&j==0) {
                            System.out.println(pag);
                            System.out.println(truePAG);
                        }
                        ArrayList<Node> latentSet = new ArrayList<Node>();
                        for (Node n : trueDAG.getNodes()) {
                            if (truePAG.getNode(n.getName()) == null)
                                latentSet.add(n);
                        }
                        DataSet dataSet;
                        if(sf)
                        dataSet = MixedUtils.loadDataSet2(dataFolder + "/SF_" + g + "_data.txt");
                        else
                            dataSet = MixedUtils.loadDataSet2(dataFolder + "/DAG_" + g + "_data.txt");
                        double[][] stats = MixedUtils.allEdgeStatsLatentNew(truePAG, pag, trueDAG, latentSet, dataSet);
                        for (int type = 0; type < 3; type++) {
                            out.print(p_loc[i] + "\t" + p_edge[j] + "\t" + g + "\t");
                            if (type == 0) {
                                out.print("CC\t");
                            } else if (type == 1) {
                                out.print("CD\t");
                            } else {
                                out.print("DD\t");
                            }
                            out.println(stats[type][TPU] + "\t" + stats[type][FPU] + "\t" + stats[type][FNU] + "\t" + stats[type][ETP] + "\t" + stats[type][EFP] + "\t" + stats[type][EFN]);
                            out.flush();
                        }
                    }
                }
            }


        }
    }
    public static Graph convertToGraph(int[][]data,int[]latents) {
        Graph temp = new EdgeListGraph();
        ArrayList<GraphNode> nodes = new ArrayList<>();
        A: for(int i = 1; i < data.length + latents.length + 1;i++)
        {
            for(int j = 0; j < latents.length;j++) {
                if (i == latents[j])
                    continue A;
            }
            GraphNode tempNode = new GraphNode("X" + i);
            nodes.add(tempNode);
            temp.addNode(tempNode);
        }
        for(int i = 0; i < data.length;i++)
        {
            for(int j = i + 1; j < data.length;j++)
            {
                if(data[i][j]!=0)
                {
                    Node x = nodes.get(i);
                    Node y = nodes.get(j);
                    Endpoint one = convertToEndpoint(data[i][j]);
                    Endpoint two = convertToEndpoint(data[j][i]);
                    if(one==Endpoint.TAIL && (two==Endpoint.CIRCLE || two==Endpoint.TAIL))
                    {
                        one = Endpoint.CIRCLE;
                        two = Endpoint.CIRCLE;
                    }
                    if(two==Endpoint.TAIL && (one==Endpoint.CIRCLE || one == Endpoint.TAIL))
                    {
                        one = Endpoint.CIRCLE;
                        two = Endpoint.CIRCLE;
                    }
                    Edge e = new Edge(x,y,one,two);
                    temp.addEdge(e);
                }
            }
        }
        return temp;
    }
    public static Endpoint convertToEndpoint(int i )
    {
        if(i==0)
            return Endpoint.NULL;
        else if(i==1)
            return Endpoint.TAIL;
        else if(i==2)
            return Endpoint.ARROW;
        else
            return Endpoint.CIRCLE;
    }
}