package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by vinee_000 on 9/13/2016.
 */
public class bccd_results {

    public static void main(String [] args) throws Exception {
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

        int numGraphs = 10;
        int numNodes = 48;
        int numLatents = 3;

        int[][] latents = new int[numGraphs][numLatents];
        BufferedReader b = new BufferedReader(new FileReader("latents.txt"));
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
        String header = "Loc_Probability\tEdge_Probability\tGraph Number\tType\tAdjacency TP\tAdjacency FP\t Adjacency FN\tOrientation TP\tOrientation FP\tOrientation FN";
        PrintStream out = new PrintStream("bccd_results_11_18.txt");
        out.println(header);
        for (int g = 0; g < numGraphs; g++) {


            for (int i = 0; i < p_loc.length; i++) {
                for (int j = 0; j < p_edge.length; j++) {
                    b = new BufferedReader(new FileReader("BCCD_" + g + "_" + (i+1) + "_" + (j+1) + ".txt"));
                    int[][] data = new int[numNodes - numLatents][numNodes - numLatents];
                    int row = 0;
                    while (b.ready()) {
                        String[] stuff = b.readLine().split("\t");
                        for (int col = 0; col < stuff.length; col++) {
                            data[row][col] = Integer.parseInt(stuff[col]);
                        }

                        row++;
                    }
                    Graph pag = convertToGraph(data, latents[g]);
                    Graph truePAG = GraphUtils.loadGraphTxt(new File("PAG_1_" + g + ".txt"));
                    Graph trueDAG = GraphUtils.loadGraphTxt(new File("SF_" + g + "_graph.txt"));
                    ArrayList<Node> latentSet = new ArrayList<Node>();
                    for (Node n : trueDAG.getNodes()) {
                        if (truePAG.getNode(n.getName()) == null)
                            latentSet.add(n);
                    }
                    DataSet dataSet = MixedUtils.loadDataSet2("SF_" + g + "_data_NonMono.txt");
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
