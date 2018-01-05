package edu.pitt.csb.Statins;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintStream;

/**
 * Created by vinee_000 on 8/29/2017.
 */
public class Find_Neighbor {
    public static void main(String [] args) throws Exception
    {
        DataSet d = MixedUtils.loadDataSet2("expression_final.txt");
        DiscreteVariable d2 = new DiscreteVariable("Dummy",2);
        d.addVariable(d2);
        for(int i = 0; i < d.getNumRows();i++)
        {
            if(i==d.getNumRows()-1)
            {
                d.setInt(i,d.getColumn(d2),1);
            }
            else
            {
                d.setInt(i,d.getColumn(d2),0);
            }
        }
        System.out.println(d);
        String [] statins = {"Simvastatin","Lovastatin","Mevastatin","Atorvastatin"};
        for(int i = 0; i < 4;i++)
        {
            DataSet temp = d.copy();
            for(int j = 0; j < 4; j++)
            {
                if(j==i)
                    continue;
                temp.removeColumn(temp.getVariable(statins[j]));
            }
            double [] lam = {.2,.2,.2};
            MGM m = new MGM(temp,lam);
            m.learnEdges(1000);
            Graph g = m.graphFromMGM();
            PrintStream out = new PrintStream(statins[i] + "_neighbors.txt");
            for(Node n : g.getAdjacentNodes(g.getNode(statins[i])))
            {
                out.println(n.getName());

            }
            out.flush();
            out.close();
        }

    }
}
