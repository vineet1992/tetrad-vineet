package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

/**
 * Created by vinee_000 on 2/6/2017.
 */
public class ivy_analysis {
    public static void main(String [] args) throws Exception
    {
        String target = "BEST_PRED_FEV1_D_(CV2-CV1)";
        DataSet d = MixedUtils.loadDataSet2("pre_processed.txt");
        PrintStream out3 = new PrintStream("SCCOR_nonparanormal.txt");
        d.removeColumn(d.getVariable("SID"));
       DataSet cont =  MixedUtils.getContinousData(d);
        cont = DataUtils.getNonparanormalTransformed(cont);
       List<Node> nodes =  d.getVariables();
        for(int i = 0; i < nodes.size();i++)
        {
            if(cont.getVariable(nodes.get(i).getName())!=null)
            {
                d.removeColumn(nodes.get(i));
            }
        }
        for(Node x:cont.getVariables())
        {
            d.addVariable(x);
            int col = d.getNumColumns()-1;
           int contCol =  cont.getColumn(x);
            for(int i = 0; i < d.getNumRows();i++)
            {
                d.setDouble(i,col,cont.getDouble(i,contCol));
            }
        }
        System.out.println(d);
        for(Node x: d.getVariables())
        {
            if(x instanceof ContinuousVariable)
            {
                System.out.println(x + ":Continuous");
            }
            else if(x instanceof DiscreteVariable)
            {
                System.out.println(x + ":Discrete");
            }
        }
        out3.println(d);
        out3.flush();
        out3.close();

        double [] lambda = {0.2,0.2,0.2};
        MGM m = new MGM(d,lambda);
        m.learnEdges(1000);
        Graph g = m.graphFromMGM();
        PrintStream out2 = new PrintStream("ivy_graph_mgm.txt");
        out2.println(g);
        out2.flush();
        out2.close();
        HashSet<Node> dups = new HashSet<Node>();
        List<Node> temp = new ArrayList<Node>();
        List<Node> temp2 = new ArrayList<Node>();
        for(Node x: g.getAdjacentNodes(g.getNode(target)))
        {
            temp2.add(x);
            temp.add(d.getVariable(x.getName()));
            dups.add(d.getVariable(x.getName()));
            for(Node y: g.getAdjacentNodes(x))
            {
                if(dups.contains(d.getVariable(y.getName())))
                    continue;
                else
                {
                    temp2.add(y);
                    temp.add(d.getVariable(y.getName()));
                    dups.add(d.getVariable(y.getName()));
                }
            }
        }
        System.out.println(temp.size());
        System.out.println(temp + "\nTEMP2:\n" + temp2);
        d = d.subsetColumns(temp);
        IndependenceTest i = new IndTestMultinomialAJ(d,.1);
        FciMaxP f = new FciMaxP(i);
        System.out.println(g.subgraph(temp2));
        //f.setInitialGraph(g.subgraph(temp2));
        //f.setInitialGraph(g);
        PrintStream out = new PrintStream("ivy_graph.txt");
        Graph h = f.search();
        out.println(h);
        out.flush();
        out.close();
    }
}
