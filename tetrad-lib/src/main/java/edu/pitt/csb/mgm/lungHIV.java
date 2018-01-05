package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by vinee_000 on 3/1/2017.
 */
public class lungHIV {
    public static void main(String [] args)throws Exception
    {
        DataSet d = MixedUtils.loadDataSet2("lung_cleaned_pos_hiv.txt");
        d.removeColumn(d.getVariable("HIV_ST"));
        d.removeColumn(d.getVariable("A1HIV_ST"));

        //ONLY REMOVE THESE IF HIV NEGATIVE//
        //d.removeColumn(d.getVariable("D15PH"));
       // d.removeColumn(d.getVariable("D17PCP"));
        //d.removeColumn(d.getVariable("A21CHEMOTHERAPY"));
        /////////////////////////////////////////////////


int dv = 0;
int cv = 0;
        for(Node n: d.getVariables())
        {
            if(n instanceof DiscreteVariable)
            {
                dv++;
            }
            else
                cv++;

        }
        System.out.println("DV: " + dv + "\tCV: " + cv);
        IndependenceTest i = new IndTestMultinomialAJ(d,.05);
        List<Node> xxx = new ArrayList<Node>();
        xxx.add(i.getVariable("WT"));
        System.out.println(i.isIndependent(i.getVariable("BMI"),i.getVariable("HT"),xxx));
        System.out.println(i.getPValue());

       // System.out.println(MixedUtils.getContinousData(d));
        double [] lambda = {0.25,0.25,0.25};

        System.exit(0);
        MGM m = new MGM(d,lambda);
        m.learnEdges(1000);
        Graph g = m.graphFromMGM();
        PrintStream out = new PrintStream("mgm_graph_pos_NOVL.txt");
        out.println(g);

        PcStable p = new PcStable(i);
        p.setInitialGraph(g);
        out.flush();
        out.close();
        out = new PrintStream("stablepc_graph_pos_NOVL.txt");
        out.println(p.search());
        out.flush();
        out.close();
        out = new PrintStream("fcimax_graph_pos_NOVL.txt");
        FciMaxP f = new FciMaxP(i);
        f.setInitialGraph(g);
        out.println(f.search());
        out.flush();
        out.close();

    }
}
