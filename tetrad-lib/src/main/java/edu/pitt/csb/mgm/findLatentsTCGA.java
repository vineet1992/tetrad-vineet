package edu.pitt.csb.mgm;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Endpoint;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by vinee_000 on 10/30/2016.
 */
public class findLatentsTCGA {
    public static void main(String[]args) throws Exception
    {
        int numGraphs = 1;
        ArrayList<Graph> graphs = new ArrayList<Graph>();
        HashMap<String,Integer> edges = new HashMap<String,Integer>();
        for(int i = 0; i < numGraphs; i ++) {
            Graph curr = GraphUtils.loadGraphTxt(new File("fci.txt"));
            for(Edge e:curr.getEdges())
            {
                if(e.getEndpoint1()== Endpoint.ARROW && e.getEndpoint2()==Endpoint.ARROW)
                {
                    String one = e.getNode1().getName();
                    String two = e.getNode2().getName();
                    if(one.compareTo(two)>0)
                    {
                        if(edges.get(one+","+two)==null)
                            edges.put(one+","+two,1);
                        else
                            edges.put(one+","+two,edges.get(one+","+two)+1);
                    }
                    else
                    {
                        if(edges.get(two+","+one)==null)
                            edges.put(two+","+one,1);
                        else
                            edges.put(two+","+one,edges.get(two+","+one)+1);
                    }
                }
            }
        }
        PrintStream results = new PrintStream("latent_map_fci.txt");
        results.println("Gene_1\tGene_2\tScore");
        PrintStream out = new PrintStream("latents_fci.txt");
        for(String s: edges.keySet())
        {

            if(edges.get(s)==1 && !s.contains("Major_Diagnosis") && !s.contains("FEV1") && !s.contains("Cig"))
            {
                out.println(s.split(",")[0] +"\t"+ s.split(",")[1]);
            }
            results.println(s.split(",")[0]+"\t"+s.split(",")[1]+"\t"+edges.get(s));
        }
        out.flush();
        out.close();
        results.flush();
        results.close();
    }
}
