package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;
import edu.cmu.tetrad.search.SepsetMap;
import edu.pitt.csb.stability.StabilityUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.PrintStream;
import java.util.List;

/**
 * Created by vinee_000 on 8/31/2017.
 */
public class exploreMGM {
    public static void main(String [] args)throws Exception
    {
        int numVariables = 50;
        int sampleSize = 100;
        int numCategories = 4;
        int ns = 20;
        MixedLeeHastieSimulation c = new MixedLeeHastieSimulation();
        Parameters p = new Parameters();
        p.setValue("numMeasures", numVariables);
        NormalDistribution n = new NormalDistribution(2*numVariables,numVariables*.75);
        p.setValue("numEdges", (int)n.sample());
        p.setValue("sampleSize", sampleSize);
        p.setValue("numCategories",numCategories);
        int numRuns = 10;
        PrintStream out = new PrintStream("MGM_Stable_Exploration_" + numVariables + "_" + sampleSize + ".txt");
        PrintStream out2 = new PrintStream("PCS_Deletion_Exploration_" + numVariables + "_" + sampleSize + ".txt");
        PrintStream out3 = new PrintStream("Graph_Statistics_" + numVariables + "_" + sampleSize + ".txt");
        out.println("Stability\tTrue_Orientation\tType\tGraph_ID");
        out2.println("Edge_Type\tStability\tPCS_Orientation\tTrue_Orientation\tP-Value\tGraph_ID");
        out3.println("Graph_ID\tNum_Edges\tMax_Node_Degree\tAvg_Node_Degree\tStandard Deviation Node Degree\tCC_Colliders\tCD_Colliders\tDD_Colliders");
                //Loop through MGM edges
                //If the edge isn't in the PCS Graph, record the following:
                //CC,Stability,PCS Orientation (Collider, Not),True Orientation (Connected, Disconnected, Collider)
        for(int ii = 0; ii < numRuns;ii++) {
            c.simulate(p);
            double low = .1;
            double high = .8;
            double[] initLambdas = new double[40];
            for (int i = 0; i < 40; i++) {
                initLambdas[i] = i * (high - low) / 40 + low;
            }
            STEPS s = new STEPS(c.getDataSet(0), initLambdas, 0.05, ns);
            Graph g2 = s.runSteps();
            double[] lambda = s.lastLambda;
            IndependenceTest ind = new IndTestMultinomialAJ(c.getDataSet(0),.05);
            PcStable pc = new PcStable(ind);
            pc.setInitialGraph(g2);
            Graph pcs = pc.search();
            SepsetMap map = pc.getSepsets();
            int b = (int) Math.floor(10 * Math.sqrt(c.getDataSet(0).getNumRows()));
            int iterLimit = 1000;
            double[][] stabilities = new double[c.getDataSet(0).getNumColumns()][c.getDataSet(0).getNumColumns()];
            for (int i = 0; i < ns; i++) {
                boolean again = true;
                while (again) {
                    try {

                        DataSet temp = c.getDataSet(0).copy();
                        temp.permuteRows();
                        int[] removal = new int[b];
                        for (int j = 0; j < removal.length; j++)
                            removal[j] = j;
                        temp = temp.subsetRows(removal);
                        MGM m = new MGM(temp, lambda);
                        m.learnEdges(iterLimit);
                      Graph  g = m.graphFromMGM();
                        again = false;
                       // System.out.println(c.getDataSet(0));
                       // System.out.println(g);
                        for (int j = 0; j < c.getDataSet(0).getNumColumns(); j++) {
                            for (int k = 0; k < c.getDataSet(0).getNumColumns(); k++) {
                                if (g.isAdjacentTo(g.getNode(c.getDataSet(0).getVariable(j).getName()), g.getNode(c.getDataSet(0).getVariable(k).getName()))) {
                                    stabilities[j][k] += 1;
                                }

                            }
                        }
                    } catch (Exception e) {
                       // e.printStackTrace();
                    }
                }

            }

            
            for (int i = 0; i < numVariables; i++) {
                for (int j = i + 1; j < numVariables; j++) {
                    Graph a = c.getTrueGraph();
                    String val = marry(a, i, j, c.getDataSet(0));
                    if (stabilities[i][j] > ns / 2)
                        out.println(stabilities[i][j] + "\t" + val + "\t" + type(a,i,j,c.getDataSet(0)) + "\t" + ii);
                }
            }


            for(Edge e: g2.getEdges())
            {
                if(pcs.getEdge(pcs.getNode(e.getNode1().getName()),pcs.getNode(e.getNode2().getName()))==null)
                {
                    //CC,Stability,PCS Orientation (Collider, Not),True Orientation (Connected, Disconnected, Collider)
                    int i = c.getDataSet(0).getColumn(c.getDataSet(0).getVariable(e.getNode1().getName()));
                    int j = c.getDataSet(0).getColumn(c.getDataSet(0).getVariable(e.getNode2().getName()));
                    List<Node> l = map.get(ind.getVariable(e.getNode1().getName()),ind.getVariable(e.getNode2().getName()));
                   ind.isIndependent(ind.getVariable(e.getNode1().getName()),ind.getVariable(e.getNode2().getName()),l);
                   double pval = ind.getPValue();
                    out2.println(type(c.getDataSet(0),e) + "\t" + stabilities[i][j] + "\t" + orientation(pcs,e) + "\t" + torient(pcs,e,c.getTrueGraph()) + "\t" + pval + "\t" + ii);
                }
            }
            out3.println(ii + "\t" + c.getTrueGraph().getNumEdges() + "\t" + maxNode(c.getTrueGraph()) + "\t" + avgNode(c.getTrueGraph()) + "\t" + stdNode(c.getTrueGraph()) + "\t" + numColliders(c.getDataSet(0),c.getTrueGraph(),"CC") + "\t" + numColliders(c.getDataSet(0),c.getTrueGraph(),"CD") + "\t" + numColliders(c.getDataSet(0),c.getTrueGraph(),"DD"));
            out3.flush();
            out.flush();
            out2.flush();
        }
        out.flush();
        out.close();
        out2.flush();
        out2.close();
        out3.flush();
        out3.close();
    }
    private static String numColliders(DataSet d, Graph truth, String type)
    {
        int sum = 0;
        if(type.equals("CC"))
        {
            for(int i = 0; i < truth.getNodes().size();i++)
            {
                for(int j = i+1; j < truth.getNodes().size();j++)
                {
                    if(d.getVariable(truth.getNodes().get(i).getName())instanceof ContinuousVariable && d.getVariable(truth.getNodes().get(j).getName()) instanceof ContinuousVariable)
                    {
                       for(Node child: truth.getChildren(truth.getNodes().get(i)))
                       {
                           if(truth.isParentOf(truth.getNodes().get(j),child))
                               sum++;
                       }
                    }
                }
            }
        }
        if(type.equals("CD"))
        {
            for(int i = 0; i < truth.getNodes().size();i++)
            {
                for(int j = i+1; j < truth.getNodes().size();j++)
                {
                    Node one = d.getVariable(truth.getNodes().get(i).getName());
                    Node two = d.getVariable(truth.getNodes().get(j).getName());
                    if((one instanceof ContinuousVariable && two instanceof DiscreteVariable) || (one instanceof DiscreteVariable && two instanceof ContinuousVariable))
                    {
                        for(Node child: truth.getChildren(truth.getNodes().get(i)))
                        {
                            if(truth.isParentOf(truth.getNodes().get(j),child))
                                sum++;
                        }
                    }
                }
            }
        }
        if(type.equals("DD"))
        {
            for(int i = 0; i < truth.getNodes().size();i++)
            {
                for(int j = i+1; j < truth.getNodes().size();j++)
                {
                    if(d.getVariable(truth.getNodes().get(i).getName())instanceof DiscreteVariable && d.getVariable(truth.getNodes().get(j).getName()) instanceof DiscreteVariable)
                    {
                        for(Node child: truth.getChildren(truth.getNodes().get(i)))
                        {
                            if(truth.isParentOf(truth.getNodes().get(j),child))
                                sum++;
                        }
                    }
                }
            }
        }
        return Integer.toString(sum);
    }
    private static String torient(Graph pcs, Edge e, Graph truth)
    {
        Node x = truth.getNode(e.getNode1().getName());
        Node y = truth.getNode(e.getNode2().getName());
        Edge e2 = truth.getEdge(x,y);
        if(orientation(truth,x,y).equals("Collider"))
            return "Collider";
        else if(truth.isAdjacentTo(x,y))
            return "Connected";
        else
            return "Disconnected";


    }
    private static String orientation(Graph pcs,Node x, Node y)
    {
        List<Node> one =  pcs.getChildren(x);
        List<Node> two = pcs.getChildren(y);
        for(Node a: one)
        {
            if(two.contains(a))
            {
                if(pcs.getDirectedEdge(x,a)!=null && pcs.getDirectedEdge(y,a)!=null)
                    return "Collider";
            }
        }
        return "Non-Collider";
    }
    private static String orientation(Graph pcs, Edge e)
    {
        Node x = e.getNode1();
        Node y = e.getNode2();

        List<Node> one =  pcs.getChildren(x);
        List<Node> two = pcs.getChildren(y);
        for(Node a: one)
        {
            if(two.contains(a))
            {
                if(pcs.getDirectedEdge(x,a)!=null && pcs.getDirectedEdge(y,a)!=null)
                    return "Collider";
            }
        }
        return "Non-Collider";
    }
    private static String type(DataSet d, Edge e)
    {
        Node x = d.getVariable(e.getNode1().getName());
        Node y = d.getVariable(e.getNode2().getName());
        if(x instanceof ContinuousVariable && y instanceof ContinuousVariable)
        {
            return "CC";
        }
        else if(x instanceof DiscreteVariable && y instanceof DiscreteVariable)
        {
            return "DD";
        }
        else
            return "CD";

    }
    private static String type(Graph g, int i, int j, DataSet d)
    {
        if(d.getVariable(i) instanceof ContinuousVariable && d.getVariable(j) instanceof ContinuousVariable)
        {
            return "CC";
        }
        else if(d.getVariable(i) instanceof DiscreteVariable && d.getVariable(i)instanceof DiscreteVariable)
        {
            return "CD";
        }
        else
            return "DD";
    }
    private static String marry(Graph g, int i, int j, DataSet d)
    {
        Node x = g.getNode(d.getVariable(i).getName());
        Node y = g.getNode(d.getVariable(j).getName());
        
        if(g.getEdge(x,y)!=null)
            return "Connected";
        for(Node n: g.getChildren(x))
        {
            if(g.getChildren(y).contains(n))
                return "Collider Connection";
        }
        return "Disconnected";
    }
    public static double maxNode(Graph g)
    {
        double max = 0;
        for(Node n: g.getNodes())
        {
            if(g.getAdjacentNodes(n).size()>max)
                max = g.getAdjacentNodes(n).size();
        }
        return max;
    }
    public static double avgNode(Graph g)
    {
        int size = g.getNodes().size();
        double avg = 0;
        for(Node n:g.getNodes())
        {
            avg+= g.getAdjacentNodes(n).size();
        }
        return avg/size;
    }
    public static double stdNode(Graph g)
    {
        double avg = avgNode(g);
        double sum = 0;
        for(Node n:g.getNodes())
        {
            sum+=Math.pow((g.getAdjacentNodes(n).size()-avg),2);
        }
        sum= sum/(g.getNodes().size()-1);
        return Math.sqrt(sum);
    }
}
