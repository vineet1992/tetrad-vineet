package edu.pitt.csb.Olja_Cancer_Analysis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.HashMap;

/**
 * Created by vinee_000 on 10/15/2017.
 */
public class convertClustersToGeneList {
    public static void main(String [] args) throws Exception
    {
            String trecsFile = "trecs_clusters.txt";
            String outFile = "cluster";
            String deFile = "../DE_Genes.txt";
           BufferedReader b = new BufferedReader(new FileReader(trecsFile));
           int count = 0;
           while(b.ready())
           {
               b.readLine();
               count++;
           }
           HashMap[] maps = new HashMap[count];
           b = new BufferedReader(new FileReader(trecsFile));
           for(int j = 0; j < count;j++)
           {
               String [] line = b.readLine().split("\t");
               maps[j] = new HashMap<String,Integer>();
               int inArr = -1;
               A:for(int i = 0; i < line.length;i++)
               {
                   int p = inArray(line[i],maps);
                   if(p!=-1) {
                       inArr = p;
                       break A;
                   }
               }
               for(int i = 0; i < line.length;i++)
               {
                   if(inArr==-1)
                       maps[j].put(line[i],0);
                   else
                       maps[inArr].put(line[i],0);
               }
           }
           b.close();


           HashMap<String,Double> fc = new HashMap<String,Double>();
           BufferedReader b2 = new BufferedReader(new FileReader(deFile));
           b2.readLine();
           while(b2.ready())
           {
               String [] line = b2.readLine().split("\t");
               fc.put(line[0],Double.parseDouble(line[1]));
           }
           b2.close();
           for(int i = 0; i < maps.length;i++) {
               PrintStream out =new PrintStream(outFile + "_" + i + ".txt");
               for (Object s : maps[i].keySet()) {
                   out.println((String)s + "\t" + fc.get(s));
               }
               out.flush();
               out.close();
           }


    }
    public static int inArray(String x, HashMap[] maps)
    {
        for(int i = 0; i < maps.length;i++)
        {
            HashMap temp = maps[i];
            if(temp==null)
                continue;
            if(temp.get(x)!=null)
                return i;
        }
      return -1;
    }
}
