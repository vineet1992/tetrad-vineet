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
import edu.pitt.csb.BioInf_Paper.CPSS;
import edu.pitt.csb.stability.Bootstrap;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.SearchWrappers;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by Vineet Raghu on 9/11/2017
 * The purpose of this class is to create a smooth interface to run several causal discovery algorithms suitable for mixed discrete and continuous data
 * This class has been compiled into a .jar file for efficient use
 * TODO allow input of true graph in order to estimate typical orientation and adjacency problems
 * TODO Edit sif to work with PAGs
 * TODO Allow for bootstrapping of the entire causal process
 * TODO Allow choice of independence test
 */
public class runAlgorithms {
    private static Graph trueGraph = null; //TODO
    private static String dataFile = ""; //Data file to be used by the algorithms
    private static DataSet d = null; //The loaded dataset itself, must be specified by the user
    private static double [] lambda = {.2,.2,.2}; //Array of Lambda values for MGM
    private static int count = 0; //Upkeep variable to load parameters
    private static String alg = "None"; //Name of the algorithm to run
    private static double alpha = .05; //Independence test threshold for constraint-based causal discovery algorithms
    private static String outputPath = "./output.txt"; //Output file for the learned causal graph
    private static boolean outputSif = false; //Should we output an .sif format file
    private static String sifPath = ""; //Path to the .sif format file output
    private static boolean runSteps = false; //Should we run steps?
    private static int ns = 20; //Number of subsamples for StARS/StEPS
    private static int b = -1; //Subsample size, will be computed automatically later as 10*sqrt(sample size)
    private static double threshold = 0.05; //Stability threshold for StARS/StEPS
    private static boolean useKnowledge = false; //Should we use a tetrad knowledge file?
    private static String kFile = "";//Path to the tetrad knowledge file
    private static double penalty = 1; //Penalty discount/Structure Prior for Score for FGES
    private static int maxNumDiscrete = 5; //Maximum number of categories a variable can have to be considered a discrete variable
    private static boolean runStars = false; //Should we run StARS?
    private static boolean runMGM = false;// Should we run MGM to get an undirected skeleton?
    private static int numParams = 20; //Number of parameters to test for StARS/StEPS
    private static boolean runCPSS = false; //Should we run CPSS?
    private static double cpssBound = 0.05; //Error rate threshold for CPSS
    private static boolean runBootstrap = false; //Should we run bootstrapping?
    private static double bootBound = 50; //Number of times an orientation must appear to be added to the output
    private static int numBoots = 100; //Number of bootstrapped samples to run
    private static boolean outputStabs = false; //Should we output edge stabilities (StARS/StEPS)?
    private static String stabPath = ""; //Path to stability file
    private static boolean outputOrients = false; //Should we output orientation frequencies for CD algorithms or CPSS or Bootstrapping?
    private static String orientPath = ""; //Path to the orientation frequency file
    private static double paramLow = -1; //Lov value of parameter range to test for StEPS
    private static double paramHigh = -1; //High value of parameter range to test for StEPS
    private static double paramLowStars = -1;//Lov value of parameter range to test for StARS
    private static double paramHighStars = -1;//High value of parameter range to test for StARS
    private static Graph initGraph; //Initial Undirected Skeleton to use before running causal discovery
    public static void main(String[] args) throws Exception{
        try {

            //Interpet command line arguments
            while(count < args.length)
            {
                if(args[count].equals("-k")) //Knowledge file to use in Tetrad knowledge format
                {
                    useKnowledge = true;
                    kFile = args[count+1];
                    count+=2;
                }
                else if(args[count].equals("-orientStabs")) //File to output orientation stability
                {
                    outputOrients = true;
                    orientPath = args[count+1];
                    count+=2;
                }
                else if(args[count].equals("-stars")) //Should we run StARS? Only for causal (directed) search algorithms
                {
                    runStars = true;
                    if(args.length==count+1 || args[count+1].startsWith("-")) {

                        count++;
                    }
                    else //Optional: how many parameters should we test?
                    {
                        numParams = Integer.parseInt(args[count+1]);
                        count+=2;
                    }

                }
                else if(args[count].equals("-cpss"))
                {
                    runCPSS = true;
                    if(args.length==count+1||args[count+1].startsWith("-"))
                        count++;
                    else {
                        cpssBound = Double.parseDouble(args[count + 1]);
                        count += 2;
                    }
                }
                else if(args[count].equals("-initGraph"))
                {
                    initGraph = GraphUtils.loadGraphTxt(new File(args[count+1]));
                    count+=2;
                }
                else if(args[count].equals("-bootstrap"))
                {
                    runBootstrap = true;
                    count++;
                }
                else if(args[count].equals("-numBoots"))
                {
                    numBoots = Integer.parseInt(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-bootBound"))
                {
                    bootBound = Double.parseDouble(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-low")) //Low value of parameter range to test (for StARS/StEPS)
                {
                    paramLow = Double.parseDouble(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-lowStars"))
                {
                    paramLowStars = Double.parseDouble(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-high")) //High value of parameter range to test
                {
                    paramHigh = Double.parseDouble(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-highStars"))
                {
                    paramHighStars = Double.parseDouble(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-b")) //Subsample size for StARS/StEPS
                {
                    b = Integer.parseInt(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-mgm"))
                {
                    runMGM = true;
                    count++;
                }
                else if(args[count].equals("-maxCat")) //Maximum number of categories any discrete variable has
                {
                    //Maximum number of categories for a discrete variable
                    maxNumDiscrete = Integer.parseInt(args[count+1]);
                    if(maxNumDiscrete <= 0)
                    {
                        System.err.println("Maximum number of categories for discrete variable cannot be negative");
                        System.exit(-1);
                    }
                    count+=2;
                }
                else if(args[count].equals("-steps")) //Should steps be run?
                {
                    runSteps = true;
                    if(args.length==count+1 || args[count+1].startsWith("-")) {
                        count++;
                    }
                    else //Optional: how many parameters should we test?
                    {

                        numParams = Integer.parseInt(args[count+1]);
                        count+=2;
                    }
                }
                else if(args[count].equals("-g")) //Stability threshold for StEPS/StARS
                {
                    threshold = Double.parseDouble(args[count+1]);
                    count+=2;
                }
                else if(args[count].equals("-ns")) //Number of subsamples for StEPS/StARS
                {
                    ns = Integer.parseInt(args[count+1]);
                    if(ns <=0)
                    {
                        System.err.println("Number of Subsamples must be greater than zero: " + ns);
                        System.exit(-1);
                    }
                    count+=2;
                }
                else if(args[count].equals("-sif"))//Should we output a file in a format suitable for cytoscape
                {
                    outputSif = true;
                    sifPath = args[count+1];
                    count+=2;
                }
                else if(args[count].equals("-stabs")) //Should we output adjacency stabilities?
                {
                    outputStabs = true;
                    if(args.length==count+1 || args[count+1].startsWith("-"))
                    {
                        System.err.println("Must specify a file name to output the adjacency stabilities");
                        System.exit(-1);
                    }
                    stabPath = args[count+1];
                    count+=2;
                }
               else if(args[count].equals("-d"))  //Required: dataset filepath
               {
                    dataFile = args[count + 1];
                    count+=2;
                }
                else if(args[count].equals("-o")) //Path for output file for general tetrad graph output
                {
                    outputPath = args[count+1];
                    count+=2;
                }
                else if(args[count].equals("-l"))  //Lambda parameters for MGM
                {
                    lambda[0] = Double.parseDouble(args[count+1]);
                    lambda[1] = Double.parseDouble(args[count+2]);
                    lambda[2] = Double.parseDouble(args[count+3]);
                    for(int i = 0; i < lambda.length;i++)
                    {
                        if(lambda[i]<0 || lambda[i]>1)
                        {
                            System.err.println("Invalid value for lambda[" + i + "], " + lambda[i]);
                            System.exit(-1);
                        }
                    }
                    count += 4;
                }
                else if(args[count].equals("-alg")) //Which algorithm should be used
                {
                    alg = args[count+1];
                    count+=2;
                }
                else if(args[count].equals("-a")) //Alpha value for constraint-based causal discovery algorithms
                {
                    alpha = Double.parseDouble(args[count+1]);
                    if(alpha < 0 || alpha > 1)
                    {
                        System.err.println("Invalid value for alpha = " + alpha);
                        System.exit(-1);
                    }
                    count+=2;
                }
                else if(args[count].equals("-penalty")) //Penalty for score-based causal discovery algorithms
                {
                    penalty = Double.parseDouble(args[count+1]);
                    if(penalty < 0)
                    {
                        System.err.println("Invalid value for penalty: " + penalty);
                        System.exit(-1);
                    }
                    count+=2;
                }
                else
                {
                    throw new Exception("Unsupported Command Line Switch: " + args[count]);
                }
            }
            //Exit program if no data file is specified
            if(dataFile.equals(""))
            {
                System.err.println("Usage: java -jar causalDiscovery.jar -d <Data File>\n Please specify a data file");
                System.exit(-1);
            }
            //Load specified dataset from the file
            d = MixedUtils.loadDataSet2(dataFile,maxNumDiscrete);
            System.out.println("Loaded dataset: " + d + "\n Dataset is mixed? " + d.isMixed());

            //Exit program if no algorithm is specfiied
            if(alg.equals("None")&&!runMGM &&!runSteps)
            {
                System.out.println("No algorithm specified, and you are not running MGM\nTherefore no method will be run to analyze the data.\n Please use -runMGM to run MGM or -a <Algorithm Name> to run a causal discovery method");
                System.exit(-1);
            }
            String [] algos = {"FCI","FCI-MAX","PCS","CPC","MAX","FGES","None"};
            boolean foundAl = false;
            for(String x:algos)
            {
                if(x.equals(alg))
                    foundAl = true;
            }

            //Ensure that algorithm name is valid
            if(!foundAl)
                throw new Exception("Unknown Algorithm: " + alg + ", Please use either \"FCI\",\"FCI-MAX\",\"PCS\", \"CPC\", \"MAX\", \"FGES\" or \"None\" (MGM Only)");


            //Workflow for running StEPS to compute optimal Lambda parameters
            if(runSteps) {

                if (!d.isMixed())
                {
                    System.err.println("Cannot run StEPS on a dataset that isn't mixed discrete and continuous...exiting");
                    System.exit(-1);
                }


                //Compute range of lambda parameters based on user specified input
                double low = .05;
                double high = .9;
                if(paramLow!=-1)
                    low = paramLow;
                if(paramHigh!=-1)
                    high = paramHigh;
                System.out.print("Running StEPS with lambda range (" + low + "," + high + "), testing " + numParams + " params...");

                double[] initLambdas = new double[numParams];
                for (int i = 0; i < numParams; i++) {
                    initLambdas[i] = i * (high - low) / numParams + low;
                }
                STEPS s;


                //If size of subsamples are not given, these will be computed automatically
                if(b==-1)
                     s = new STEPS(d,initLambdas,threshold,ns);
                else
                    s = new STEPS(d,initLambdas,threshold,ns,b);
                s.setComputeStabs(outputStabs);


                //Run StEPS
                s.runStepsPar();

                //Get stabilities and optimal lambda
                double [][] stabs = s.stabilities;
                double [] lbm = s.lastLambda;



                //Output edge appearence stability to file
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


            //Workflow to run MGM based on optimal lambda or user specified lambda
            if(runMGM) {
                if(!d.isMixed())
                {
                    System.out.println("Dataset is not mixed continuous and discrete... cannot run MGM");
                    System.exit(-1);
                }
                if(runCPSS) //MGM with CPSS
                {
                    Bootstrap bs = new Bootstrap(d,cpssBound,50);
                    bs.setCPSS();
                    g = bs.runBootstrap(convert("MGM"),lambda);
                    if(outputOrients)
                        outputBootOrients(bs.getTabularOutput());
                }
                else if(runBootstrap) //MGM with Bootstrapping
                {
                    Bootstrap bs = new Bootstrap(d,bootBound,numBoots);
                    g = bs.runBootstrap(convert("MGM"),lambda);
                    if(outputOrients)
                        outputBootOrients(bs.getTabularOutput());
                }
                else //Vanilla MGM
                {
                    System.out.print("Running MGM with lambda params: " + Arrays.toString(lambda) + "...");
                    MGM m = new MGM(d, lambda); //Create MGM object
                    m.learnEdges(1000);//Use maximum 1000 iterations to learn the edges for the undirected MGM graph, stop searching if the edges in the graph don't change after 3 iterations
                    g = m.graphFromMGM(); //store the mgm graph
                    System.out.println("Done");
                    if(outputOrients)
                        System.out.println("For subsampled edge frequency information for MGM, please run StEPS");
                }
            }
            Graph finalOutput = null;
            IKnowledge k = null;
            //Incorporate prior knowledge file into cd algorithms if specified
            if(useKnowledge) {
                if (kFile == null) {
                    throw new IllegalStateException("No knowledge data file was specified.");
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

            //Run StARS to compute optimal parameter for specifieid cd algorithm
            if(runStars)
            {
                if(alg.equals("None"))
                {
                    System.err.println("Cannot run StARS without a specified directed causal discovery algorithm");
                    System.exit(-1);
                }
                double pLow = 0.0001;
                double pHigh = 0.8;

                //Tailor parameter range to algorithm being used
                if(alg.equals("FGES") && (d.isMixed() || d.isDiscrete()))
                {
                    pLow = 1;
                    pHigh = 10;
                }
                else if(alg.equals("FGES") && d.isContinuous())
                {
                    pLow = 0.01;
                    pHigh = 20;
                }

                //Or use user specified paramter range
                if(paramLowStars!=-1)
                    pLow = paramLowStars;
                if(paramHighStars!=-1)
                    pHigh = paramHighStars;

                //Compute paramter range
                double [] penalties = new double[numParams];
                for(int j = 0; j < penalties.length;j++)
                {
                    penalties[j] = pLow + j*(pHigh-pLow)/numParams;
                }


                STARS strs;
                Algorithm a = convert(alg);

                System.out.print("Running StARS for algorithm, " + a + " with " + numParams + " parameters in the range (" + pLow + "," + pHigh + ")...");
                if(b==-1)
                    strs = new STARS(d,penalties,threshold,ns,a);
                else
                    strs = new STARS(d,penalties,threshold,ns,a,b);
                //Run StARS to get optimal parameter, need to flip parameter range if FGS is being run (low to high)
                penalty = strs.getAlpha(a==Algorithm.FGS);
                System.out.println("STARS Chosen Parameter: " + penalty);


                //Output edge appearence stabilities (this will overwrite StEPS stabilities)
                double [][] stab = strs.stabilities; //This is edge stability
                if(outputStabs)
                {
                    PrintStream stabOut = new PrintStream(stabPath);
                    edu.pitt.csb.mgm.runSteps.printStability(stabOut,d,stab);
                }
            }

            //Run Regular CD Algorithm
            if(!runCPSS && !runBootstrap) {

                //Run any algorithm that isn't FGES
                if(!alg.equals("None") && !alg.equals("FGES")) {
                    System.out.print("Running " + alg + "...");
                    DataGraphSearch gs = Algorithm.algToSearchWrapper(convert(alg), new double[]{alpha});
                    if(g!=null)
                        gs.setInitialGraph(g);
                    else if(initGraph!=null)
                        gs.setInitialGraph(g);
                    if(k!=null)
                        gs.setKnowledge(k);
                    finalOutput = gs.search(d);
                    System.out.println("Done");
                }
                //Run FGES
                else if(alg.equals("FGES"))
                {
                    try{
                        System.out.print("Running FGES...");
                        Score s;
                        if(d.isMixed()) {
                            ConditionalGaussianScore s2 = new ConditionalGaussianScore(d);
                            s2.setStructurePrior(penalty);
                            s = s2;
                        }
                        else if(d.isContinuous()) {
                            s = new SemBicScore(new CovarianceMatrixOnTheFly(d), penalty);
                        }
                        else {
                            BDeuScore s2 = new BDeuScore(d);
                            s2.setStructurePrior(penalty);
                            s = s2;
                        }
                            Fges fg = new Fges(s);
                            if (useKnowledge)
                                fg.setKnowledge(k);
                            if(g!=null)
                                fg.setInitialGraph(g);
                            else if(initGraph!=null)
                                fg.setInitialGraph(initGraph);
                            finalOutput = fg.search();

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

                //Otherwise just save MGM output as the final output
                else
                {
                    finalOutput = g;
                }


                //This will output orientation stabilities for all edges that showed up at least once with some orientation
                //Could consider changing this to edges that showed up in the output graph only, but this will require a bit of post-processing
                //TODO
                if(outputOrients)
                {
                    Bootstrap bs = new Bootstrap(d,1,ns);
                    if(b<=0)
                        bs.setSubsample(0);
                    else
                        bs.setSubsample(b);
                    if(alg.equals("None"))
                    {
                        System.out.println("Cannot output orientation stability without a specified causal discovery algorithm");
                    }
                    else {
                        double param = 0;
                        if (alg.equals("FGES"))
                            param = penalty;
                        else
                            param = alpha;
                        bs.runBootstrap(convert(alg),new double[]{param});
                        outputBootOrients(bs.getTabularOutput());
                    }
                }
            }
            //Complimentary Pairs Stability Selection for all algorithms
            else if(runCPSS)
            {
                Bootstrap bs = new Bootstrap(d, cpssBound, 50);
                bs.setCPSS();
                if(!alg.equals("FGES") && !alg.equals("None")) {

                    finalOutput = bs.runBootstrap(convert(alg), new double[]{alpha});
                }
                else if(!alg.equals("None")) //Use penalty instead of alpha for FGES
                {
                    finalOutput = bs.runBootstrap(convert(alg), new double[]{penalty});
                }
                if(outputOrients)
                    outputBootOrients(bs.getTabularOutput());

            }
            //Bootstrapping, same workflow as CPSS (but no bs.setCPSS() line)
            else if(runBootstrap)
            {
                Bootstrap bs = new Bootstrap(d,bootBound,numBoots);
                if(!alg.equals("FGES") && !alg.equals("None")) {

                    finalOutput = bs.runBootstrap(convert(alg),new double[]{alpha});
                }
                else if(!alg.equals("None"))
                {
                    finalOutput = bs.runBootstrap(convert(alg),new double[]{penalty});
                }
                if(outputOrients)
                    outputBootOrients(bs.getTabularOutput());
            }

            if(outputSif)
            {
                PrintStream out = new PrintStream(sifPath);
                printSif(finalOutput,out);
            }
            PrintStream out = new PrintStream(outputPath);
            out.println(finalOutput);
            out.flush();
            out.close();




        } catch (IOException e){
            e.printStackTrace();
        }
    }


    public static void outputBootOrients(String o) throws Exception
    {
        PrintStream out = new PrintStream(orientPath);
        out.println(o);
        out.flush();
        out.close();
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




    //Helper method to convert a string to an Algorithm representation
    private static Algorithm convert(String alg)
    {
        Algorithm a;
        if(alg.equals("CPC"))
        {
            a = Algorithm.CPC;
        }
        else if(alg.equals("MAX"))
        {
            a = Algorithm.PCMAX;
        }
        else if(alg.equals("FCI"))
        {
            a = Algorithm.FCI;
        }
        else if(alg.equals("FCI-MAX"))
        {
            a = Algorithm.FCIMAX;
        }
        else if(alg.equals("FGES"))
        {
            a = Algorithm.FGS;
        }
        else if(alg.equals("None"))
        {
            a = Algorithm.MGM;
        }
        else
        {
            a = Algorithm.PCS;
        }
        return a;
    }
}

