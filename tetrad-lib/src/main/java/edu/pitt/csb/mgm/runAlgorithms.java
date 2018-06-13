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
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.performance.PerformanceTests;
import edu.cmu.tetrad.search.*;
import edu.pitt.csb.stability.SearchWrappers;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by Vineet Raghu on 9/11/2017
 * The purpose of this class to create an interface to smoothly interact with MGM via several online platforms
 * This can further be used in general
 */
public class runAlgorithms {
    private static Graph trueGraph = null;
    private static String dataFile = "";
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
    private static int b = -1;
    private static double g = 0.05;
    private static boolean useKnowledge = false;
    private static String kFile = "";
    private static double penalty = 2;
    private static int maxNumDiscrete = 5;
    private static boolean runStars = false;
    private static int numParams = 20;
    private static boolean outputStabs = false;
    private static String stabPath = "";
    private static boolean outputOrients = false;
    private static String orientPath = "";
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
                else if(args[count].equals("-orientStabs"))
                {
                    outputOrients = true;
                    orientPath = args[count+1];
                    count+=2;
                }
                else if(args[count].equals("-stars"))
                {
                    if(args.length==count+1 || args[count+1].startsWith("-")) {
                        runStars = true;
                        count++;
                    }
                    else
                    {
                        numParams = Integer.parseInt(args[count+1]);
                        count+=2;
                    }

                }
                else if(args[count].equals("-b"))
                {
                    b = Integer.parseInt(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-maxCat"))
                {
                    //Maximum number of categories for a discrete variable
                    maxNumDiscrete = Integer.parseInt(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-steps"))
                {
                    runSteps = true;
                    count++;
                }
                else if(args[count].equals("-g"))
                {
                    g = Double.parseDouble(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-ns"))
                {
                    ns = Integer.parseInt(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-sif"))
                {
                    outputSif = true;
                    sifPath = args[count+1];
                    count+=2;
                }
                else if(args[count].equals("-stabs"))
                {
                    outputStabs = true;
                    stabPath = args[count+1];
                    count+=2;
                }
               else if(args[count].equals("-d")) {
                    dataFile = args[count + 1];
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
                else if(args[count].equals("-penalty"))
                {
                    penalty = Double.parseDouble(args[count+1]);
                    count+=2;
                }
                else
                {
                    throw new Exception("Unsupported Command Line Switch: " + args[count]);
                }
            }
            if(dataFile.equals(""))
            {
                System.err.println("No data file specified");
                System.exit(-1);
            }
            d = MixedUtils.loadDataSet2(dataFile,maxNumDiscrete);
            String [] algos = {"PCS","CPC","MAX","FGES","None"};
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
                STEPS s;
                if(b==-1)
                     s = new STEPS(d,initLambdas,g,ns);
                else
                    s = new STEPS(d,initLambdas,g,ns,b);
                s.runStepsPar();
                double [][] stabs = s.stabilities;
                double [] lbm = s.lastLambda;
                if(outputStabs)
                {
                    PrintStream stabOut = new PrintStream(stabPath);
                    edu.pitt.csb.mgm.runSteps.printStability(stabOut,d,stabs);
                }
                lambda = s.lastLambda;
                PrintStream out = new PrintStream(outputPath);
                out.println(lbm[0] + "\t" + lbm[1] + "\t" + lbm[2]);
                out.flush();
                out.close();
                System.out.println("Done");
            }

            Graph g = null;
            if(!alg.equals("FGES")) {
                MGM m = new MGM(d, lambda); //Create MGM object
                m.learnEdges(1000);//Use maximum 1000 iterations to learn the edges for the undirected MGM graph, stop searching if the edges in the graph don't change after 3 iterations
                g = m.graphFromMGM(); //store the mgm graph

            }
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
                    out.println(finalOutput);
                    out.flush();
                    out.close();
                    System.out.println("Done");
                }
                else if(alg.equals("FGES"))
                {
                    if(d.isMixed())
                    {
                        System.err.println("Can't apply FGES to mixed data");
                        System.exit(-1);
                    }
                    try{
                        System.out.print("Running StARS...");
                        double pLow = 0.01;
                        double pHigh = 15;

                        double [] penalties = new double[numParams];
                        for(int j = 0; j < penalties.length;j++)
                        {
                            penalties[j] = pLow + j*(pHigh-pLow)/numParams;
                        }
                        STARS strs;
                        if(b==-1)
                            strs = new STARS(d,penalties,0.05,ns,Algorithm.FGS);
                        else
                            strs = new STARS(d,penalties,0.05,ns,Algorithm.FGS,b);
                        penalty = strs.getAlpha(true);
                        System.out.println("STARS Chosen Parameter: " + penalty);
                        double [][] stab = strs.stabilities; //This is edge stability
                        if(outputStabs)
                        {
                            PrintStream stabOut = new PrintStream(stabPath);
                            edu.pitt.csb.mgm.runSteps.printStability(stabOut,d,stab);
                        }


                        System.out.println("Done");
                        System.out.print("Running FGES...");
                        Score s = new SemBicScore(new CovarianceMatrixOnTheFly(d),penalty);
                        Fgs2 fg = new Fgs2(s);
                        if(useKnowledge)
                            fg.setKnowledge(k);
                        finalOutput = fg.search();

                        PrintStream out = new PrintStream(outputPath);
                        out.println(finalOutput);
                        out.flush();
                        out.close();

                        if(outputOrients)
                        {
                            double[][][] orients = StabilityUtils.OrientationSearchPar(d, new SearchWrappers.FgsWrapper(penalty),ns,b);
                            outputOrient(orients,outputPath,finalOutput,orientPath);
                        }
                        System.out.println("Done");
                    }
                    catch(Exception e)
                    {
                        System.err.println("Error Running FGES, check for collinearity");
                        e.printStackTrace();
                        System.exit(-1);
                    }
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
                    out.println(finalOutput);
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
                    out.println(finalOutput);
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
    //                            outputOrient(orients,outputPath,finalOutput,orientPath);

    public static void outputOrient(double [][][] orients,String outputPath,Graph finalOutput,String orientPath) throws Exception
    {
        List<String> output= new ArrayList<String>();
        BufferedReader b = new BufferedReader(new FileReader(outputPath));
        for(int x =0 ; x < 3; x++)
            output.add(b.readLine());
        output.add(b.readLine() + " --> <-- --- <->");
        for(Edge e: finalOutput.getEdges())
        {
            boolean switched = false;
            int a1 = d.getColumn(d.getVariable(e.getNode1().getName()));
            int a2 = d.getColumn(d.getVariable(e.getNode2().getName()));

            if(a1 > a2)
                switched = true;
            String curr = b.readLine();
            int [] switchInts = {1,0,2,3};
            for(int x = 0; x < 4; x++)
            {
                if(switched)
                    curr = curr + " " + orients[switchInts[x]][a2][a1];
                else
                    curr = curr + " " + orients[x][a1][a2];
            }
            output.add(curr);
        }
        b.close();
        PrintStream out = new PrintStream(orientPath);

        for(String x:output)
            out.println(x);
        out.flush();
        out.close();

    }
    //Variable Name: "Undirected":[LIST OF UNDIRECTED EDGES] "Directed":LIST OF CHILDREN

}

