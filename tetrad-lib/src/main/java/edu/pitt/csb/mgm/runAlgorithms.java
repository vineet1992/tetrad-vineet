///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015 by Peter Spirtes, Richard Scheines, Joseph   //
// Ramsey, and Clark Glymour.                                                //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA //
///////////////////////////////////////////////////////////////////////////////

package edu.pitt.csb.mgm;

import cern.colt.matrix.DoubleMatrix2D;
import com.google.gson.JsonArray;
import com.google.gson.stream.JsonReader;
import edu.cmu.tetrad.cmd.TetradCmd;
import edu.cmu.tetrad.data.DataReader;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.performance.PerformanceTests;
import edu.cmu.tetrad.search.CpcStable;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcMax;
import edu.cmu.tetrad.search.PcStable;

import java.io.*;
import java.util.Arrays;
import java.util.List;

/**
 * Created by Vineet Raghu on 9/11/2017
 * The purpose of this class to create an interface to smoothly interact with MGM via several online platforms
 */
public class runAlgorithms {
    private static Graph trueGraph = null;
    private static DataSet d = null;
    private static double [] lambda = {.2,.2,.2};
    private static int count = 0;
    private static String alg = "None";
    private static double alpha = .05;
    private static String outputPath = "./output.txt";
    private static boolean outputSif = false;
    private static String sifPath = "";
    private static boolean runSteps = false;
    private static int ns = 20;
    private static double g = 0.05;
    private static boolean useKnowledge = false;
    private static String kFile = "";
    public static void main(String[] args) throws Exception{
        try {

            while(count < args.length)
            {
                if(args[count].equals("-k"))
                {
                    useKnowledge = true;
                    kFile = args[count+1];
                    count+=2;
                }
                if(args[count].equals("-steps"))
                {
                    runSteps = true;
                    count++;
                }
                if(args[count].equals("-g"))
                {
                    g = Double.parseDouble(args[count+1]);
                    count+=2;
                }
                if(args[count].equals("-ns"))
                {
                    ns = Integer.parseInt(args[count+1]);
                    count+=2;
                }
                if(args[count].equals("-sif"))
                {
                    outputSif = true;
                    sifPath = args[count+1];
                    count+=2;
                }
                if(args[count].equals("-d")) {
                    d = MixedUtils.loadDataSet2(args[count + 1]);
                    count+=2;
                }
                else if(args[count].equals("-o"))
                {
                    outputPath = args[count+1];
                    count+=2;
                }
                else if(args[count].equals("-l"))
                {
                    lambda[0] = Double.parseDouble(args[count+1]);
                    lambda[1] = Double.parseDouble(args[count+2]);
                    lambda[2] = Double.parseDouble(args[count+3]);
                    count += 4;
                }
                else if(args[count].equals("-alg"))
                {
                    alg = args[count+1];
                    count+=2;
                }
                else if(args[count].equals("-a"))
                {
                    alpha = Double.parseDouble(args[count+1]);
                    count+=2;
                }
                else
                {
                    throw new Exception("Unsupported Command Line Switch: " + args[count]);
                }
            }
            String [] algos = {"PCS","CPC","MAX","None"};
            boolean foundAl = false;
            for(String x:algos)
            {
                if(x.equals(alg))
                    foundAl = true;
            }
            if(!foundAl)
                throw new Exception("Unknown Algorithm: " + alg + ", Please use either PCS, CPC, MAX, or None");
            if(runSteps)
            {
                System.out.print("Running StEPS...");
                double low = .05;
                double high = .9;
                double[] initLambdas = new double[40];
                for (int i = 0; i < 40; i++) {
                    initLambdas[i] = i * (high - low) / 40 + low;
                }
                STEPS s = new STEPS(d,initLambdas,g,ns);
                s.runStepsArrayPar();
                double [] lbm = s.lastLambda;
                lambda = s.lastLambda;
                PrintStream out = new PrintStream(outputPath);
                out.println(lbm[0] + "\t" + lbm[1] + "\t" + lbm[2]);
                out.flush();
                out.close();
                System.out.println("Done");
            }

            MGM m = new MGM(d,lambda); //Create MGM object
            m.learnEdges(1000);//Use maximum 1000 iterations to learn the edges for the undirected MGM graph, stop searching if the edges in the graph don't change after 3 iterations
            Graph g = m.graphFromMGM(); //store the mgm graph
            Graph finalOutput = null;
            IKnowledge k = null;

            if(useKnowledge) {
                if (kFile == null) {
                    throw new IllegalStateException("No data file was specified.");
                }

                try {
                    File knowledgeFile = new File(kFile);

                    CharArrayWriter writer = new CharArrayWriter();

                    FileReader fr = new FileReader(knowledgeFile);
                    int i;

                    while ((i = fr.read()) != -1) {
                        writer.append((char) i);
                    }

                    DataReader reader = new DataReader();
                    char[] chars = writer.toCharArray();

                    k = reader.parseKnowledge(chars);
                }
                catch(Exception e)
                {
                    System.out.println("Unable to read knowledge file");
                    return;
                }
            }
            if(!alg.equals("None")) {
                IndependenceTest i = new IndTestMultinomialAJ(d, alpha); //Create an independence test object suitable for mixed data
                if(alg.equals("PCS")) {
                    System.out.print("Running PCS...");
                    PcStable p = new PcStable(i); //Create a PC Stable search (you can replace this with any search algorithm that you would want to use
                    p.setInitialGraph(g); //Restrict the search to only the edges chosen by MGM
                    if(useKnowledge)
                        p.setKnowledge(k);
                    finalOutput = p.search(); //Get final output from the pc-stable search
                    PrintStream out = new PrintStream(outputPath);
                    out.flush();
                    out.close();
                    System.out.println("Done");
                }
                else if(alg.equals("CPC"))
                {
                    System.out.print("Running CPC...");
                    CpcStable p = new CpcStable(i);
                    p.setInitialGraph(g);
                    if(useKnowledge)
                        p.setKnowledge(k);
                    finalOutput = p.search();
                    PrintStream out = new PrintStream(outputPath);
                    out.flush();
                    out.close();
                    System.out.println("Done");
                }
                else if(alg.equals("MAX"))
                {
                    System.out.print("Running PC-Max...");
                    PcMax p = new PcMax(i);
                    p.setInitialGraph(g);
                    if(useKnowledge)
                        p.setKnowledge(k);
                    finalOutput = p.search();
                    PrintStream out = new PrintStream(outputPath);
                    out.flush();
                    out.close();
                    System.out.println("Done");
                }
            }
            else
            {
                PrintStream out = new PrintStream(outputPath);
                finalOutput = g;
                out.println(g);
                out.flush();
                out.close();
            }
            if(outputSif)
            {
                PrintStream out = new PrintStream(sifPath);
                printSif(finalOutput,out);
            }





        } catch (IOException e){
            e.printStackTrace();
        }
    }


    public static void printSif(Graph g, PrintStream out)
    {
        for(Edge e:g.getEdges())
        {
            if(e.isDirected())
            {
                if(g.isParentOf(g.getNode(e.getNode1().getName()),g.getNode(e.getNode2().getName())))
                {
                    out.println(e.getNode1().getName() + "\tdir\t" + e.getNode2().getName());
                }
                else
                {
                    out.println(e.getNode2().getName() + "\tdir\t" + e.getNode2().getName());
                }
            }
            else
            {
                out.println(e.getNode1() + "\tundir\t" + e.getNode2());
            }
        }
        out.flush();
        out.close();
    }
    //Variable Name: "Undirected":[LIST OF UNDIRECTED EDGES] "Directed":LIST OF CHILDREN

}

