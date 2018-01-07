package edu.pitt.csb.Olja_Cancer_Analysis;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.HashMap;

/**
 * Created by vinee_000 on 10/14/2017.
 */
public class getResponseStabilities {
    public static void main(String [] args)throws Exception
    {
        String graphPath = "Dp_Full_No_Gender_IgG.txt";
        String stabPath = "Stabilities_Dp_Full_No_Gender_IgG.txt";
        //String outPath = "Stability_Dp_Full_IgG_No_Gender_Neighbors.txt";
        //String outPath2 = "Stability_Dp_Full_IgG_No_Gender_Non_Neighbors.txt";

        String outPath = "";
        String outPath2 = "";
        String target = "IgG_Ratio";

        int index = 0;
        while(index < args.length)
        {
            if(args[index].equals("-graph"))
            {
                graphPath = args[index+1];
                index+=2;
            }
            else if(args[index].equals("-stab"))
            {
                stabPath = args[index+1];
                index+=2;
            }
            else if(args[index].equals("-out"))
            {
                outPath = args[index+1];
                index+=2;
            }
            else if(args[index].equals("-out2"))
            {
                outPath2 = args[index]+1;
                index+=2;
            }
            else if(args[index].equals("-t"))
            {
                target = args[index+1];
                index+=2;
            }
        }
        if(outPath.equals(""))
            outPath = "Neighbors_" + graphPath;
        if(outPath2.equals(""))
            outPath2 = "Non_Neighbors_" + graphPath;
        Graph graph = GraphUtils.loadGraphTxt(new File(graphPath));
        BufferedReader b = new BufferedReader(new FileReader(stabPath));
        HashMap<String,Double> map = new HashMap<String,Double>();
        String [] names = b.readLine().split("\t");
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            if(line[0].equals(target))
            {
                for(int i = 1;i<line.length;i++)
                {
                    map.put(names[i-1],Double.parseDouble(line[i]));
                }
            }
        }
        b.close();
        PrintStream out = new PrintStream(outPath);
        for(Node n: graph.getAdjacentNodes(graph.getNode(target)))
        {
            out.println(n.getName() + "\tundir\t" + target + "\t"  + map.get(n.getName()));
        }
        out.flush();
        out.close();
        PrintStream out2 = new PrintStream(outPath2);
        for(String x:map.keySet())
        {
            if(map.get(x)>0)
            {
               out2.println(x + "\tundir\t" + target + "\t" + map.get(x));
            }
        }
        out2.flush();
        out2.close();
    }
}
