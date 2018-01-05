package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.simulation.ContinuousLinearGaussianSemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Endpoint;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.*;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by vinee_000 on 3/16/2017.
 */
public class crossValidation {
    public static void main(String [] args) throws Exception
    {
        PrintStream out = new PrintStream("pc_cv_1000.txt");
        out.println("Run\tTrue_Edge_Type\tProb_Adjacency(NonZero)\tProb_I_to_J\tProb_J_to_I\tProb_Undirected\tProb_Bidirected");

        PrintStream out2 = new PrintStream("ges_cv_1000.txt");
        out2.println("Run\tTrue_Edge_Type\tProb_Adjacency(NonZero)\tProb_I_to_J\tProb_J_to_I\tProb_Undirected\tProb_Bidirected");


        PrintStream out3 = new PrintStream("cpc_cv_1000.txt");
        out3.println("Run\tTrue_Edge_Type\tProb_Adjacency(NonZero)\tProb_I_to_J\tProb_J_to_I\tProb_Undirected\tProb_Bidirected");
        int numRuns = 25;
        for(int run =0; run < numRuns; run++) {
            Parameters p = new Parameters();
            ContinuousLinearGaussianSemSimulation c = new ContinuousLinearGaussianSemSimulation();
            c.simulate(p);
            DataSet d = c.getDataSet(0);
            double[] alpha = {.001, .005,.01,.03, .05,.08, .1, .2,.3,.4, .5};
            double[] penalty = {.01,.1,.5,1, 2, 4,7, 10, 20, 40,100};
            ArrayList<int[][]> pc = new ArrayList<int[][]>(); //stores the learned graph as an adjacency matrix
            ArrayList<int[][]> ges = new ArrayList<int[][]>(); //stores the learned graph as adjacency matrices
            ArrayList<int[][]> cpc = new ArrayList<int[][]>();
            for (int index = 0; index < alpha.length; index++) {
                int [][] pcAdj = new int[d.getNumColumns()][d.getNumColumns()];
                int [][] gesAdj = new int[d.getNumColumns()][d.getNumColumns()];
                int [][] cpcAdj = new int[d.getNumColumns()][d.getNumColumns()];
                IndependenceTest i = new IndTestFisherZ(d,alpha[index]);
                PcStable pcs = new PcStable(i);
                pcAdj = convert(pcs.search());
                Fgs2 f = new Fgs2(new SemBicScore(new CovarianceMatrix(d),penalty[index]));
                gesAdj = convert(f.search());
                CpcStable cpcS = new CpcStable(new IndTestFisherZ(d,alpha[index]));
                cpcAdj = convert(cpcS.search());
                pc.add(pcAdj);
                ges.add(gesAdj);
                cpc.add(cpcAdj);
            }
            Graph trueGraph = SearchGraphUtils.patternForDag(c.getTrueGraph());
                ArrayList<String> pcAdjDist = getAdjacencyDistribution(pc,trueGraph);

                ArrayList<String>gesAdjDist = getAdjacencyDistribution(ges,trueGraph);
for(String x: pcAdjDist)
{
    out.println(run + "\t" + x);
}
            for(String x: gesAdjDist)
            {
                out2.println(run + "\t" + x);
            }
            ArrayList<String>cpcAdjDist = getAdjacencyDistribution(cpc,trueGraph);
            for(String x: cpcAdjDist)
            {
                out3.println(run + "\t" + x);
            }

        }
        out.flush();
        out.close();
        out2.flush();
        out2.close();
        out3.flush();
        out3.close();
    }
    //0 not connected
    //1 i causes j
    //2 j causes i
    //3 undirected
    //4 bidirected edge between i and j

//This function produces a table of True_Adjacency(0 or 1) \t Probability of adjacency
//True_Edge_Type	Prob_Adjacency(NonZero)	Prob_I_to_J	Prob_J_to_I	Prob_Undirected	Prob_Bidirected"
    public static ArrayList<String> getAdjacencyDistribution(ArrayList<int[][]> adjacencies,Graph truth )
    {
        List<Node> nodes = truth.getNodes();
        int [][] orMat = convert(truth);
        ArrayList<String> output = new ArrayList<String>();
        for(int i = 0; i < nodes.size();i++) {

            for(int j = i + 1; j < nodes.size();j++)
            {
                Node one = nodes.get(i);
                Node two = nodes.get(j);
                int row = nodes.indexOf(one);
                int col = nodes.indexOf(two);
                int currOrientation = orMat[row][col];
                int [] guesses = new int[5];
                double count = 0;
                for(int [][] estOrientation: adjacencies)
                {
                    int orientation = estOrientation[row][col];
                    guesses[orientation]++;

                    count++;
                }
                String x = currOrientation + "\t" + (count-guesses[0])/count + "\t" + guesses[1]/count + "\t" + guesses[2]/count + "\t" + guesses[3]/count + "\t" + guesses[4]/count;
output.add(x);
            }
        }
        return output;
    }

    //this converts a graph into an orientation matrix, under the assumption of no latent variables
    public static int[][] convert(Graph g)
    {
        int [][] graphMat = new int[g.getNumNodes()][g.getNumNodes()];
        List<Node> allNodes = g.getNodes();
        for(int i = 0; i < g.getNodes().size();i++)
        {
            for(int j = i+1; j < g.getNodes().size();j++)
            {
                Edge e = g.getEdge(allNodes.get(i),allNodes.get(j));
                if(e==null)
                {
                    graphMat[i][j] = 0;
                    graphMat[j][i] = 0;
                }
                else {
                    graphMat[i][j] = edge2Num(e, allNodes.get(i), allNodes.get(j));
                    if(graphMat[i][j]==3 || graphMat[i][j] == 4)
                        graphMat[j][i] = graphMat[i][j];
                    else
                        graphMat[j][i] = -1*graphMat[i][j] + 3; // flip 2 and 1
                }

            }
        }


        return graphMat;
    }
    public static int edge2Num(Edge e,Node one, Node two)
    {
        if(e.getProximalEndpoint(one)==Endpoint.ARROW)
        {
            if(e.getProximalEndpoint(two)==Endpoint.ARROW)
            {
                return 4; // i <-> j
            }
            else
            {
                return 2; // j -> i
            }
        }
        else if(e.getProximalEndpoint(two)==Endpoint.ARROW)
        {
            return 1; // i -> j
        }
        else
            return 3; // i -- j
    }
}
