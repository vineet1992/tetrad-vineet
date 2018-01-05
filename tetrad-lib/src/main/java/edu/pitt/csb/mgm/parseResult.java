package edu.pitt.csb.mgm;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;

/**
 * Created by vinee_000 on 9/3/2016.
 */
public class parseResult {
    public static void main(String [] args) throws Exception
    {
        BufferedReader b = new BufferedReader(new FileReader("Atorva_tree.csv"));
        BufferedReader b2 = new BufferedReader(new FileReader("ratio_neighbors_rosuva.txt"));
        BufferedReader b3 = new BufferedReader(new FileReader("expression_computed_ratios_500.txt"));

        String [] genes = b3.readLine().replace("\"","").split("\t"); //remember to add one to index to get matlab index

        PrintStream p  = new PrintStream("rosuva_neighbors_clustered.txt");
        ArrayList<ArrayList<String>> clusters = new ArrayList<ArrayList<String>>();
        for(int i = 0; i < 19; i++)
        {
            clusters.add(new ArrayList<String>());
        }
        int count = 0;
        while(b.ready())
        {
            String [] line = b.readLine().split(",");

            for(String x:line)
            {
                String g = genes[Integer.parseInt(x)-1];
                clusters.get(count).add(g);
            }
            count++;
            b.readLine();
        }
        b.close();
        System.out.println(clusters);
        while(b2.ready())
        {
            ArrayList<Integer> ints = new ArrayList<Integer>();
            String x = b2.readLine();
            for(int i = 0; i < 19;i++)
            {

                if(clusters.get(i).contains(x))
                {
                    ints.add(i);
                }

            }
            p.print(x);
            for(Integer l:ints)
            {
                p.print("\t" + l);
            }
            p.println();
        }
p.flush();
        p.close();
        b2.close();
        b3.close();
    }
}
