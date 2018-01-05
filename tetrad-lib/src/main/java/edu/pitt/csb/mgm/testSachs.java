package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.Knowledge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.performance.PerformanceTests;
import edu.cmu.tetrad.search.*;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

/**
 * Created by vinee_000 on 5/17/2016.
 */
public class testSachs {
    public static void main(String [] args) throws Exception {
        int[] vars = {2, 4, 7, 8};
        DataSet curr = MixedUtils.loadDataSet2("data_full_final_normalized.txt");

        //  DataSet da =curr.copy();
        //int inde = 2;
        //praf, pmek,plcg,pip2, pip3, p44/42, pakts473, pka, pkc, p38, pjnk, Dummy
        //String va = da.getVariable(vars[inde]).getName();
        //da.removeColumn(vars[inde]);
        //  IndependenceTest indte = new IndTestMultinomialAJ(curr,.01);
        //  CpcStable ppp = new CpcStable(indte);
        //  PrintStream pppout = new PrintStream("cpcstable_PKA_fullgraph.txt");
        //  pppout.println(ppp.search());
        PrintStream fciResults = new PrintStream("results/fci_results.txt");
        PrintStream mgmResults = new PrintStream("results/mgm_results.txt");
        PrintStream maxResults = new PrintStream("results/max_results.txt");
        PrintStream mfmResults = new PrintStream("results/mfm_results.txt");
        fciResults.println("Graph_ID\tUndirected Precision\tUndirected Recall\tHidden Variable Prec\tHidden Variable Rec\tOrientation Precision\tOrientation Recall\tExperimental Precision\tExperimental Recall");
        mgmResults.println("Graph_ID\tUndirected Precision\tUndirected Recall\tHidden Variable Prec\tHidden Variable Rec\tOrientation Precision\tOrientation Recall\tExperimental Precision\tExperimental Recall");
        maxResults.println("Graph_ID\tUndirected Precision\tUndirected Recall\tHidden Variable Prec\tHidden Variable Rec\tOrientation Precision\tOrientation Recall\tExperimental Precision\tExperimental Recall");
        mfmResults.println("Graph_ID\tUndirected Precision\tUndirected Recall\tHidden Variable Prec\tHidden Variable Rec\tOrientation Precision\tOrientation Recall\tExperimental Precision\tExperimental Recall");
        for (int index = 0; index < vars.length; index++) {
            double[] lambda = {.1, .1, .1};
            double alpha = .001;
            DataSet data = curr.copy();
            //praf, pmek,plcg,pip2, pip3, p44/42, pakts473, pka, pkc, p38, pjnk, Dummy
            String var = data.getVariable(vars[index]).getName();
            data.removeColumn(vars[index]);
           /* Knowledge k = new Knowledge();
            for (int i = 1; i < 9; i++) {
                k.addToTier(0, "Int" + i);
            }
            for (Node x : data.getVariables()) {
                if (!x.getName().startsWith("Int"))
                    k.addToTier(1, x.getName());

            }*/
            //TODO
            //Plan of attack
            //Get all children of removed node
            //if there is a bidirected edge between two of them: TP
            //if there is no edge between them: FN
            //if there is a bidirected edge elsewhere: FP
            Graph truth = GraphUtils.loadGraphTxt(new File("true_graph.txt"));
            //truth.getNode("plcg").setNodeType(NodeType.LATENT);
            truth.getNode(var).setNodeType(NodeType.LATENT);
            final DagToPag dagToPag = new DagToPag(truth);
            dagToPag.setCompleteRuleSetUsed(false);
            truth = dagToPag.convert();
            System.out.println("Finished Conversion");
            Graph realGraph = GraphUtils.loadGraphTxt(new File("true_graph.txt"));
            ArrayList<Node> latents = new ArrayList<Node>();
            latents.add(realGraph.getNode(var));
            realGraph.getNode(var).setNodeType(NodeType.LATENT);
            //System.out.println(data);
            //System.out.println(data.isMixed());
            /*int[] subSample = new int[1600];
            for (int i = 0; i < 1600; i++) {
                subSample[i] = (int) (Math.random() * 1600);
            }
            DataSet use = data.subsetRows(subSample);*/
            MGM m = new MGM(data, lambda);
            System.out.println("Finished initializing MGM");
            m.learnEdges(1000);
            System.out.println("Finished learning MGM");
            Graph init = m.graphFromMGM();
            IndependenceTest i = new IndTestMultinomialAJ(data, alpha);
            Fci f2 = new Fci(i);
            PrintStream out2 = new PrintStream("results/fci_hard_" + var + ".txt");
            //f2.setKnowledge(k);
            Graph fci_g = f2.search();

            out2.println(fci_g);
            out2.flush();
            PrintStream out3 = new PrintStream("results/mgm_fci_hard_" + var + ".txt");
            PrintStream out4 = new PrintStream("results/max_hard_" + var + ".txt");
            f2.setInitialGraph(init);
            System.out.println("Beginning fci search");
            Graph mgm_fci_g = f2.search();
            System.out.println("Ended FCI Search");
            out3.println(mgm_fci_g);
            out3.flush();

            FciMaxP f3 = new FciMaxP(i);
            //f3.setKnowledge(k);
            Graph fci_max_g = f3.search();
            out4.println(fci_max_g);
            out4.flush();

            FciMaxP f = new FciMaxP(i);
            f.setInitialGraph(init);
            //f.setKnowledge(k);

            Graph mfm_g = f.search();
            PrintStream output = new PrintStream("results/mfm_hard_" + var + ".txt");
            output.println(mfm_g);
            PrintStream out = new PrintStream("results/truePAG_hard_" + var + ".txt");
            out.println(truth);
            out.flush();

            output.flush();
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
            truth.addNode(fci_g.getNode("Dummy"));
            double[][] fci_counts = MixedUtils.allEdgeStatsLatentNew(truth, fci_g, realGraph, latents, data);
            double[][] mgm_fci_counts = MixedUtils.allEdgeStatsLatentNew(truth, mgm_fci_g, realGraph, latents, data);
            double[][] fci_max_counts = MixedUtils.allEdgeStatsLatentNew(truth, fci_max_g, realGraph, latents, data);
            double[][] mfm_counts = MixedUtils.allEdgeStatsLatentNew(truth, mfm_g, realGraph, latents, data);
            out3.close();
            out4.close();
            out.close();
            output.close();
            out2.flush();
            out2.close();//TP, FP, FN, Undetermined
            System.out.println(truth);
            System.out.println(fci_g);
            System.out.println(realGraph);
            System.out.println("Are the graphs");
            double[][] result_fci = MixedUtils.allEdgeStatsLatentPairwise(truth, fci_g, realGraph, latents, data);
            double[][] result_mgm = MixedUtils.allEdgeStatsLatentPairwise(truth, mgm_fci_g, realGraph, latents, data);
            double[][] result_max = MixedUtils.allEdgeStatsLatentPairwise(truth, fci_max_g, realGraph, latents, data);
            double[][] result_mfm = MixedUtils.allEdgeStatsLatentPairwise(truth, mfm_g, realGraph, latents, data);

            fciResults.println(var + "\t" + newPrec(fci_counts) + "\t" + newRec(fci_counts) + "\t" + newBidirPrec(fci_counts) + "\t" + newBidirRec(fci_counts) + "\t" + newOrPrec(result_fci) + "\t" + newOrRec(result_fci) + "\t" + newOrPrec(fci_counts) + "\t" + newOrRec(fci_counts));
            mgmResults.println(var + "\t" + newPrec(mgm_fci_counts) + "\t" + newRec(mgm_fci_counts) + "\t" + newBidirPrec(mgm_fci_counts) + "\t" + newBidirRec(mgm_fci_counts) + "\t" + newOrPrec(result_mgm) + "\t" + newOrRec(result_mgm) + "\t" + newOrPrec(mgm_fci_counts) + "\t" + newOrRec(mgm_fci_counts));
            maxResults.println(var + "\t" + newPrec(fci_max_counts) + "\t" + newRec(fci_max_counts) + "\t" + newBidirPrec(fci_max_counts) + "\t" + newBidirRec(fci_max_counts) + "\t" + newOrPrec(result_max) + "\t" + newOrRec(result_max) + "\t" + newOrPrec(fci_max_counts) + "\t" + newOrRec(fci_max_counts));
            mfmResults.println(var + "\t" + newPrec(mfm_counts) + "\t" + newRec(mfm_counts) + "\t" + newBidirPrec(mfm_counts) + "\t" + newBidirRec(mfm_counts) + "\t" + newOrPrec(result_mfm) + "\t" + newOrRec(result_mfm) +"\t" +  newOrPrec(mfm_counts) + "\t" + newOrRec(mfm_counts));

        }
    }
    public static double newOrPrec(double[][]arr)
    {
        return arr[0][3]/(arr[0][4]+arr[0][3]);
    }
    public static double newOrRec(double[][]arr)
    {
        return arr[0][3]/(arr[0][5]+arr[0][3]);
    }
    public static double newBidirPrec(double[][]arr)
    {
        return (arr[0][11]/(arr[0][11] + arr[0][12]));
    }
    public static double newBidirRec(double[][]arr)
    {
        return (arr[0][11])/(arr[0][13] + arr[0][11]);
    }
    public static double newPrec(double[][]arr)
    {
        return (arr[0][0]/(arr[0][0] + arr[0][1]));
    }
    public static double newRec(double[][]arr)
    {
        return (arr[0][0]/(arr[0][0] + arr[0][2]));
    }
    public static int getbidirFP(int[][]arr)
    {
        return arr[0][4] + arr[4][4]+arr[5][4]+arr[7][4];
    }
    public static int getbidirFN(int[][]arr)
    {
        return arr[6][0]+arr[6][3] +arr[6][5];
    }
    public static int getbidirAgg(int[][]arr)//Truth is undetermined, we give definitely bidirected
    {
        return arr[1][4] + arr[2][4] + arr[3][4];
    }
    public static int getbidirUnd(int[][]arr) //Truth is bidirected, we give possibly bidirected
    {
        return arr[6][1] + arr[6][2];
    }
    public static double getbidirRec(int[][] arr) {
        int row = 6;
        double sum = 0;
        //col 4 = correct
        for (int i = 0; i < 6; i++) {
            sum += arr[row][i];

        }
        return (double)arr[row][4] / (double)sum;
    }

    public static double getbidirPrec(int[][] arr) {
        int col = 4;
        double sum = 0;
        for (int i = 0; i < 8; i++) {
            sum += arr[i][col];
        }
        return (double)arr[6][col] / (double)sum;
    }
    public static int getbidirPrec2(int[][]arr) //returns false positives
    {
        int sum = 0;
        for(int i = 0; i < 8; i ++) {
            if (i != 6)
                sum += arr[i][4];
        }
        return sum;
    }
    public static int getbidirCorrect(int[][]arr)
    {
        return arr[6][4];
    }
    public static double getRec(int[][] arr) {
        double wrong = 0;
        double corr = 0;
        int second = 5;
        for (int i = 0; i < 7; i++) {
            wrong = wrong + arr[i][second];

        }
        for (int i = 0; i < 7; i++) {
            for (int j = 0; j < 5; j++) {
                corr += arr[i][j];

            }
        }
        return (double)corr / (double)(corr + wrong);
    }

    public static double getPrec(int[][] arr) {
        double corr = 0;
        double wrong = 0;
        for (int j = 0; j < 5; j++) {
            wrong += arr[7][j];
        }
        for (int i = 0; i < 7; i++) {
            for (int j = 0; j < 5; j++) {
                corr += arr[i][j];
            }
        }
        return (double)corr / (double)(corr + wrong);
    }

    public static double getAncestral(double[][] arr) {
        double wrong = 0;
        double sum = 0;
        for (int i = 0; i <= 7; i++) {
            if (i == 2) //A o-> B (B is not an ancestor of A)
            {
                for (int j = 2; j < 5; j++) {
                    sum += arr[i][j];
                }
                wrong += arr[i][5];
            } else if (i == 3) {
                wrong += arr[i][3];
                wrong += arr[i][5];
                sum += arr[i][4];
            } else if (i == 4) {
                sum += arr[i][2];
                sum += arr[i][3];
                wrong += arr[i][4];
                wrong += arr[i][5];

            } else if (i == 5) {
                wrong += arr[i][2];
                wrong += arr[i][3] * 2;
                sum += arr[i][4];
                wrong += arr[i][4];
                wrong += arr[i][5];

            } else if (i == 6) {
                sum += arr[i][2];
                sum += arr[i][3];
                wrong += arr[i][3];
                sum += arr[i][4] * 2;
                wrong += arr[i][4] * 2;
            } else if (i == 7) {
                wrong += arr[i][2];
                wrong += arr[i][3];
                wrong += arr[i][4] * 2;
            }
        }
        return sum / (sum + wrong);

    }

    public static double getAcc(int[][] arr) {
        //1,1 / 2,2 / 4,3 / 6,4
        double total = 0;
        double sum = 0;
        for (int i = 0; i < 7; i++) {
            for (int j = 0; j < 5; j++) {
                if ((i == 1 && j == 1) || (i == 2 && j == 2) || (i == 4 && j == 3) || (i == 6 && j == 4))
                    sum = sum + arr[i][j]*2;
                else if((i==1&&j==2) || (i==2 && j==1) || (i==2 && j==3) || (i==2 && j==4) || (i==3 && j ==1) || (i==3 && j ==4) || (i==4 && j ==2) || (i==4 && j ==4) || (i==6 && j ==2) || (i==6 && j ==3))
                    sum = sum + arr[i][j];

                total = total + arr[i][j]*2;

            }

        }
        return (double)sum / (double)total;

    }
}
