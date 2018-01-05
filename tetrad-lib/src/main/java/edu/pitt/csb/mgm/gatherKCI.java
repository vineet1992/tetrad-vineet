package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

/**
 * Created by vinee_000 on 12/11/2016.
 */
public class gatherKCI {
    public static void main(String [] args) throws Exception {

        int numGraphs = 10;
        double [] iarr = {.25,.35,.55};
        double [] alp = {.001,.01,.05};
        PrintStream kci_results;
        PrintStream mult_results;

        kci_results = new PrintStream("kci_results_12_11.txt");
        mult_results = new PrintStream("mult_results_12_11.txt");
        String header = "Alpha\tLambda\tGraph\tType\tAdjacency TP\tAdjacency FP\tAdjacency FN\tExperimental TP\tExperimental TN\tExperimental FP\tExperimental FN\tLatent TP\tLatent FP\tLatent FN";
        int TPU = 0;
        int FPU = 1;
        int FNU = 2;
        int ETP = 3;
        int EFP = 4;
        int EFN = 5;
        int ETN = 6;
        int ITP =  7;
        int IFP = 8;
        int IFN = 9;
        int ITN = 10;
        int LTP = 11;
        int LFP = 12;
        int LFN = 13;
        kci_results.println(header);
        mult_results.println(header);

        for (int i = 0; i < iarr.length; i++) {
            for (int j = 0; j < alp.length; j++) {
                for (int fe = 0; fe < numGraphs; fe++) {
                    System.out.println("Here " + fe);
                    // Graph trueDAG = GraphUtils.loadGraphTxt(new File("DAG_" + fe + "_graph.txt"));
                    Graph trueDAG = GraphUtils.loadGraphTxt(new File("SF_" + fe + "_graph.txt"));
                    Graph truePAG = GraphUtils.loadGraphTxt(new File("PAG_" + i + "_" + fe + ".txt"));
                    ArrayList<Node> latents = new ArrayList<Node>();
                    for(Node n: trueDAG.getNodes())
                    {
                        if(truePAG.getNode(n.getName())==null)
                            latents.add(n);
                    }

                    Graph kciPAG = GraphUtils.loadGraphTxt(new File("kci_" + i + "_" + j + "_" + fe + ".txt"));
                    Graph multPAG = GraphUtils.loadGraphTxt(new File("mult_" + i + "_" + j + "_" + fe + ".txt"));

                    DataSet data = MixedUtils.loadDataSet2("SF_" + fe + "_data_NonMono.txt");
                    System.out.println(truePAG + "\n" + kciPAG + "\n" + multPAG + "\n" + trueDAG + "\n" + latents + "\n" + data);
                    double [][] kci_stats = MixedUtils.allEdgeStatsLatentNew(truePAG, kciPAG, trueDAG, latents,data);
                    double [][] mult_stats = MixedUtils.allEdgeStatsLatentNew(truePAG,multPAG,trueDAG,latents,data);
                    double lab = iarr[i];
                    for(int type = 0; type < 3; type++)
                    {
                        if(type==0)
                        {
                            kci_results.print(alp[j] + "\t" + lab + "\t" + fe +"\tCC\t");
                            mult_results.print(alp[j] + "\t" + lab + "\t" + fe +"\tCC\t");
                        }
                        else if (type==1)
                        {
                            kci_results.print(alp[j] + "\t" + lab + "\t" + fe +"\tCD\t");
                            mult_results.print(alp[j] + "\t" + lab + "\t" + fe +"\tCD\t");
                        }
                        else
                        {
                            kci_results.print(alp[j] + "\t" + lab + "\t" + fe +"\tDD\t");
                            mult_results.print(alp[j] + "\t" + lab + "\t" + fe +"\tDD\t");
                        }
                        kci_results.println(kci_stats[type][TPU] + "\t" + kci_stats[type][FPU]+"\t" + kci_stats[type][FNU]+"\t"+kci_stats[type][ETP]+"\t" + kci_stats[type][ETN] + "\t" + kci_stats[type][EFP] + "\t" + kci_stats[type][EFN] + "\t" + kci_stats[type][LTP]+"\t" + kci_stats[type][LFP] + "\t" + kci_stats[type][LFN]);
                        mult_results.println(mult_stats[type][TPU] + "\t" + mult_stats[type][FPU]+"\t" + mult_stats[type][FNU]+"\t"+mult_stats[type][ETP]+"\t" + mult_stats[type][ETN] + "\t" + mult_stats[type][EFP] + "\t" + mult_stats[type][EFN] + "\t" + mult_stats[type][LTP]+"\t" + mult_stats[type][LFP] + "\t" + mult_stats[type][LFN]);

                        kci_results.flush();
                        mult_results.flush();

                    }

                }
            }
        }

        kci_results.close();

        mult_results.close();
    }
}
