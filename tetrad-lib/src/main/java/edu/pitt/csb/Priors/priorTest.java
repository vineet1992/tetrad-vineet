package edu.pitt.csb.Priors;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;
import edu.pitt.csb.stability.StabilityUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by vinee_000 on 4/26/2017.
 */
public class priorTest {
    public static void main(String [] args) throws Exception
    {
        boolean numEdgesRandom = true;
        boolean excludeUnreliable = true; //Should piMGM exclude priors below p-value threshold
        double reliabilityThreshold = 0.05; //piMGM will exclude priors with adjusted p-value < 0.05
        double amountPrior = .1;
        boolean reliable = false; //Are all priors reliable?
        boolean diffNumPrior = false; //Does each prior provide with the same number of edges?
        boolean correctEdges = false; //Determines whether or not we will use correct edges only for unreliable priors as well or not
        boolean pureRandom = false; //Only for different number of edges given by each prior, this sets the priors to be purely random with reliability computed after the fact
        int reliableExperts = 3; //How many priors are reliable?
        int numExperts = 5;
        int numLambdas = 40;
        int numVariables = 100;
        int numEdges = 75;
        int sampleSize = 500;
        int numSubsamples = 10;
        int numRuns = 15;
        int index = 0;
        int numCategories = 4;
        double gamma = 0.05;
        boolean saveData = true;
        boolean reuseData = true;
        boolean rerunAlgorithms = false;
        String directory = ".";
   //     String [] algs = {"mgm_one_steps","oracle_one","STEPS","oracle","mgm_priors","mgm_priors_split"};
   //     String [] algs = {"mgm_one_steps","STEPS"};
         String[] algs = {"mgm_priors","oracle_one","STEPS","oracle","mgm_one_steps"};
        while(index < args.length)
        {

            if(args[index].equals("-ex")) {
                numExperts = Integer.parseInt(args[index + 1]);
                index += 2;
            }
            else if(args[index].equals("-rd"))
            {
                reuseData = true;
                index++;
            }
           else if(args[index].equals("-d")) {
                directory = args[index + 1];
                index += 2;
            }
            else if(args[index].equals("-nc"))
            {
                numCategories = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-sv"))
            {
                saveData = true;
                index++;
            }
            else if(args[index].equals("-p"))
            {
                amountPrior = Double.parseDouble(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-l"))
            {
                numLambdas = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-v"))
            {
                numVariables = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-e"))
            {
                numEdges = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-s"))
            {
                sampleSize = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-ns"))
            {
                numSubsamples = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-r"))
            {
                numRuns = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-g"))
            {
                gamma = Double.parseDouble(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-re"))
            {
                reliableExperts = Integer.parseInt(args[index+1]);
                index+=2;
            }
        }

        double [] initLambdas = new double[numLambdas];
        double low = .05;
        double high = .9;
        double inc = (high-low)/(numLambdas-1);
        for(int i = 0; i < numLambdas;i++)
        {
            initLambdas[i] = i*inc + low;
        }

        double [][][] results = new double[numRuns][4][algs.length];

        File rFile = new File("Results");
        if(!rFile.isDirectory())
            rFile.mkdir();
        File gFile = new File("Graphs");
        File dFile = new File("Data");
        File pFile = new File("Priors");
        File estFile = new File("Estimated");
        File subFile = new File("Subsamples");
        if(!reliable)
        {
            File wFile = new File("Weights");
            if(!wFile.isDirectory())
                wFile.mkdir();
        }
        if(saveData) {
            if (!gFile.isDirectory())
                gFile.mkdir();
            if (!dFile.isDirectory())
                dFile.mkdir();
            if(!estFile.isDirectory())
                estFile.mkdir();
            if(!pFile.isDirectory())
                pFile.mkdir();
            if(!subFile.isDirectory())
                subFile.mkdir();
        }

        /*PrintStream pri = new PrintStream(directory + "/mgm_priors_" + amountPrior + "_" + numExperts + "_" + numVariables +  "_" + sampleSize + "_" + numSubsamples + ".txt");
        PrintStream step = new PrintStream(directory + "/STEPS_" + amountPrior + "_" + numExperts + "_" + numVariables + "_"  + sampleSize + "_" + numSubsamples + ".txt");
        PrintStream orc = new PrintStream(directory + "/oracle_" + amountPrior + "_" + numExperts + "_" + numVariables + "_" + sampleSize + "_" + numSubsamples + ".txt");
        PrintStream one = new PrintStream(directory + "/mgm_one_steps_" + amountPrior + "_" + numExperts + "_" + numVariables + "_" + sampleSize + "_" + numSubsamples + ".txt");
        PrintStream orcOne = new PrintStream(directory + "/oracle_one_" + amountPrior + "_" + numExperts + "_" + numVariables + "_" + sampleSize + "_" + numSubsamples + ".txt");*/
    NormalDistribution n = new NormalDistribution(numVariables,numVariables/2);
        PrintStream [] pri = new PrintStream [algs.length];
        for(int i = 0; i < algs.length;i++) {
            if(!reliable)
                pri[i] = new PrintStream("Results/" + algs[i] + "_" + amountPrior + "_" + reliableExperts + "_" + numExperts + "_" + numVariables + "_" + sampleSize + "_" + numSubsamples + ".txt");
            else
                pri[i] = new PrintStream("Results/" + algs[i] + "_" + amountPrior + "_" + numExperts + "_" + numVariables + "_" + sampleSize + "_" + numSubsamples + ".txt");
        }
        A:for(int i = 0; i < numRuns; i++) {
            MixedLeeHastieSimulation c = new MixedLeeHastieSimulation();
            System.out.println(i);
            if(numEdgesRandom)
            {
                numEdges = (int)n.sample();
                while(numEdges <(numVariables/4))
                    numEdges = (int)n.sample();
            }
            Parameters p = new Parameters();
            p.setValue("numMeasures", numVariables);
            p.setValue("numEdges", numEdges);
            p.setValue("sampleSize", sampleSize);
            p.setValue("numCategories",numCategories);
            c.simulate(p);
            TetradMatrix[] priors = new TetradMatrix[numExperts];
            DataSet [] subsamples = new DataSet[numSubsamples];
            if(reuseData)
            {
                boolean foundFile = false;
                boolean foundData = false;
                File f = new File("Graphs/Graph_" + i + "_" + numVariables + ".txt");
                if(f.exists()) {
                    foundFile = true;
                    c.setTrueGraph(GraphUtils.loadGraphTxt(f));
                    removeMoral(c.getTrueGraph());
                }
                f = new File("Data/Data_" + i + "_" + numVariables + "_" + sampleSize + ".txt");
                if(f.exists() && foundFile) {
                    foundData = true;
                    c.setDataSet(MixedUtils.loadDataSet2("Data/Data_" + i + "_" + numVariables + "_" + sampleSize + ".txt"), 0);
                }
                else
                    c.simulate(p);
                if(foundFile) {
                    for (int j = 0; j < numExperts; j++) {
                        if(reliable)
                         f = new File("Priors/Priors_" + i + "_" + numVariables + "_" + amountPrior + "_" + j + ".txt");
                        else
                            f = new File("Priors/Priors_" + i + "_" + numVariables + "_" + amountPrior + "_" + reliableExperts + "_" + j + ".txt");
                        if (f.exists()) {
                            priors[j] = new TetradMatrix(loadPrior(f, numVariables));
                        }
                    }
                }
                if(foundData) {
                    for (int j = 0; j < numSubsamples; j++) {
                        f = new File("Subsamples/Subsample_" + i + "_" + numVariables + "_" + sampleSize + "_" + j + ".txt");
                        if(f.exists())
                            subsamples[j] = MixedUtils.loadDataSet2(f.getAbsolutePath());
                    }
                }

            }

            Graph [] estimatedGraphs = new Graph[algs.length];

            if(!rerunAlgorithms) {
                for (int j = 0; j < algs.length; j++) {
                File f = new File("Estimated/" + algs[j] + "_" + i + "_" + numVariables + "_" + sampleSize + "_" + amountPrior + "_" + numExperts + ".txt");

                if(!reliable && algs[j].contains("priors"))
                    f = new File("Estimated/" + algs[j] + "_" + i + "_" + numVariables + "_" + sampleSize + "_" + amountPrior + "_" + reliableExperts + "_" + numExperts + ".txt");

                if (f.exists()) {

                        Graph est = GraphUtils.loadGraphTxt(f);
                        estimatedGraphs[j] = est;

                    }

                }
            }
            int [] edges = new int[numExperts];
            double [] reli = new double[numExperts];
           // System.out.println(c.getDataSet(0));
            c.setTrueGraph(moralize(c.getTrueGraph()));
           //  priors = simulatePrior(c.getTrueGraph(), amountPrior, numExperts,priors,c.getDataSet(0));
            if(reliable)
                priors = simulatePriorDim(c.getTrueGraph(), amountPrior, numExperts,priors,c.getDataSet(0));
            else if(!reliable && diffNumPrior) {

                priors = simulatePriorUnreli(c.getTrueGraph(), numExperts, c.getDataSet(0),edges,reli,pureRandom);
            }
          //  else if(reliable&&diffNumPrior)
            //    priors = simulatePriorDim(c.getTrueGraph(),numExperts,priors,c.getDataSet(0)); TODO
            else
                priors = simulatePriorUnreli(c.getTrueGraph(),amountPrior,numExperts,c.getDataSet(0),reliableExperts,correctEdges);
            boolean done = false;
            if (!diffNumPrior && correctEdges)
                 checkPriors(c.getTrueGraph(),c.getDataSet(0),priors,i,true);
            else
                checkPriors(c.getTrueGraph(),c.getDataSet(0),priors,i,false);

            while(!done)
            {
                System.out.println(i);
                try {
                    Graph mgmprior = null;
                    mgmPriors m = null;
                    boolean foundIt = false;
                    boolean foundSteps = false;
                    //String [] algs = {"mgm_one_steps","oracle_one","STEPS","oracle","mgm_priors"};
                    for(int pp = 0; pp < algs.length;pp++)
                    {
                        if((algs[pp].contains("oracle")||algs[pp].contains("prior")) && estimatedGraphs[pp]==null)
                            foundIt = true;
                        if((algs[pp].contains("STEPS")||algs[pp].contains("steps")) && estimatedGraphs[pp]==null)
                            foundSteps = true;
                    }

                    boolean nullSub = false;
                    for(int j = 0; j < subsamples.length;j++)
                    {
                        if(subsamples[j]==null)
                        {
                            nullSub = true;
                        }
                    }
                    if(nullSub) {
                        int b = (int) Math.floor(10 * Math.sqrt(c.getDataSet(0).getNumRows()));
                        if (b > c.getDataSet(0).getNumRows())
                            b = c.getDataSet(0).getNumRows() / 2;
                        int[][] samps = StabilityUtils.subSampleNoReplacement(c.getDataSet(0).getNumRows(), b, numSubsamples);
                        for (int j = 0; j < subsamples.length; j++) {
                            subsamples[j] = c.getDataSet(0).subsetRows(samps[j]);
                        }
                    }
                    STEPS s = new STEPS(c.getDataSet(0), initLambdas, gamma, subsamples);
                    Graph steps = null;
                    if(foundSteps)
                        steps = s.runStepsPar();
                    if(foundIt) {
                        m = new mgmPriors(numSubsamples, initLambdas, c.getDataSet(0), priors, c.getTrueGraph(),false,subsamples);
                        if(excludeUnreliable)
                            m.excludeUnreliablePriors(reliabilityThreshold);
                        mgmprior = m.runPriors();
                    }


                    for(int j = 0; j < algs.length;j++)
                    {
                        if(algs[j].equals("mgm_one_steps"))
                        {
                            if(estimatedGraphs[j]==null) {
                                double[] mgmOne = {s.origLambda, s.origLambda, s.origLambda};
                                MGM oneLambda = new MGM(c.getDataSet(0), mgmOne);
                                oneLambda.learnEdges(1000);
                                Graph oneLamb = oneLambda.graphFromMGM();
                                estimatedGraphs[j] = oneLamb;
                            }


                        }
                        else if(algs[j].equals("oracle_one"))
                        {
                            if(estimatedGraphs[j]==null) {
                                double[] lbOneOracle = {m.oracleAll, m.oracleAll, m.oracleAll};
                                MGM oneOracle = new MGM(c.getDataSet(0), lbOneOracle);
                                oneOracle.learnEdges(1000);
                                Graph finalOne = oneOracle.graphFromMGM();
                                estimatedGraphs[j] = finalOne;
                            }
                        }
                        else if(algs[j].equals("STEPS"))
                        {
                            if(estimatedGraphs[j]==null)
                            estimatedGraphs[j] = steps;
                        }
                        else if(algs[j].equals("oracle"))
                        {
                            if(estimatedGraphs[j]==null) {
                                double[] lb = {m.oracleCC, m.oracleCD, m.oracleDD};
                                MGM m2 = new MGM(c.getDataSet(0), lb);
                                m2.learnEdges(1000);
                                Graph oracle = m2.graphFromMGM();
                                estimatedGraphs[j] = oracle;
                            }
                        }
                        else if(algs[j].equals("mgm_priors"))
                        {
                            if(estimatedGraphs[j]==null)
                            estimatedGraphs[j] = mgmprior;
                            if(!reliable) {
                                PrintStream p3 = new PrintStream("Weights/Weight_" + algs[j] + "_" + i + "_" + numVariables + "_" + sampleSize + "_" + amountPrior + "_" + reliableExperts + "_" + numExperts + ".txt");
                                p3.println(Arrays.toString(m.expertWeights));
                                if(diffNumPrior)
                                {
                                    p3.println(Arrays.toString(reli));
                                    p3.println(Arrays.toString(edges));
                                    p3.println(Arrays.toString(m.normalizedTao));
                                    p3.println(Arrays.toString(m.pValues));
                                    p3.println(Arrays.toString(m.normalizedExpertWeights));
                                }
                                if(excludeUnreliable)
                                    p3.println(m.lastHavePrior.length);
                                p3.flush();
                                p3.close();
                            }
                        }
                        else if(algs[j].equals("mgm_priors_split"))
                        {
                            mgmPriors mtemp = new mgmPriors(numSubsamples, initLambdas, c.getDataSet(0), priors, c.getTrueGraph(),true,subsamples);
                            Graph mgmprior2 = mtemp.runPriors();
                            estimatedGraphs[j] = mgmprior2;
                            if(!reliable) {
                                PrintStream p3 = new PrintStream("Weights/Weight_" + algs[j] + "_" + i + "_" + numVariables + "_" + sampleSize + "_" + amountPrior + "_" + reliableExperts + "_" + numExperts + ".txt");
                                p3.println(Arrays.toString(mtemp.expertWeights));
                                if(diffNumPrior)
                                {
                                    p3.println(Arrays.toString(reli));
                                    p3.println(Arrays.toString(edges));
                                    p3.println(Arrays.toString(mtemp.normalizedTao));
                                    p3.println(Arrays.toString(mtemp.pValues));
                                    p3.println(Arrays.toString(mtemp.normalizedExpertWeights));
                                }
                                if(excludeUnreliable)
                                    p3.println(m.lastHavePrior.length);
                                p3.flush();
                                p3.close();
                            }
                        }
                        System.out.println(algs[j] + "," + mgmPriors.getF1(estimatedGraphs[j],c.getTrueGraph(),c.getDataSet(0),"CC"));
                        System.out.println(algs[j] + "," + mgmPriors.getF1(estimatedGraphs[j],c.getTrueGraph(),c.getDataSet(0),"CD"));
                        System.out.println(algs[j] + "," + mgmPriors.getF1(estimatedGraphs[j],c.getTrueGraph(),c.getDataSet(0),"DD"));
                        System.out.println(algs[j] + "," + mgmPriors.getF1(estimatedGraphs[j],c.getTrueGraph(),c.getDataSet(0),"All"));


                        results[i][0][j] = mgmPriors.getF1(estimatedGraphs[j],c.getTrueGraph(),c.getDataSet(0),"CC");
                        results[i][1][j] = mgmPriors.getF1(estimatedGraphs[j],c.getTrueGraph(),c.getDataSet(0),"CD");
                        results[i][2][j] = mgmPriors.getF1(estimatedGraphs[j],c.getTrueGraph(),c.getDataSet(0),"DD");
                        results[i][3][j] = mgmPriors.getF1(estimatedGraphs[j],c.getTrueGraph(),c.getDataSet(0),"All");

                    }


                    done = true;
                } catch (Exception e) {
                    e.printStackTrace();
                    if(!reuseData) {
                        removeMoral(c.getTrueGraph());
                        c.simulate(p);
                        moralize(c.getTrueGraph());
                    }
                }
            }

            if(saveData) {
                PrintStream p2 = new PrintStream("Graphs/Graph_" + i + "_" + numVariables +  ".txt");
                p2.println(c.getTrueGraph());
                p2.flush();
                p2.close();
                p2 = new PrintStream("Data/Data_" + i + "_" + numVariables + "_" + sampleSize + ".txt");
                p2.println(c.getDataSet(0));
                p2.flush();
                p2.close();
                for(int x = 0; x < algs.length;x++)
                {
                    p2 = new PrintStream("Estimated/" + algs[x] + "_" + i + "_" + numVariables +  "_" + sampleSize + "_" + amountPrior + "_" + numExperts + ".txt");
                    if(!reliable && algs[x].contains("priors"))
                        p2 = new PrintStream("Estimated/" + algs[x] + "_" + i + "_" + numVariables + "_" + sampleSize + "_" + amountPrior + "_" + reliableExperts + "_" + numExperts + ".txt");

                    p2.println(estimatedGraphs[x]);
                    p2.flush();
                    p2.close();
                }
                for(int j = 0; j < numExperts;j++) {
                    p2 = new PrintStream("Priors/Priors_" + i + "_" + numVariables + "_" + amountPrior + "_" + j + ".txt");
                    if(!reliable)
                        p2 = new PrintStream("Priors/Priors_"+ i + "_" + numVariables + "_" + amountPrior + "_" + reliableExperts + "_" + j + ".txt");
                   for(int k = 0; k < numVariables;k++)
                   {
                       for(int m = 0; m < numVariables;m++)
                       {
                           if(m==numVariables-1)
                            p2.println(priors[j].get(k,m));
                           else
                               p2.print(priors[j].get(k,m)+"\t");
                       }
                   }
                    p2.flush();
                    p2.close();

                }
                for(int j = 0; j < subsamples.length;j++)
                {
                    p2 = new PrintStream("Subsamples/Subsample_" + i + "_" + numVariables + "_" + sampleSize + "_" + j + ".txt");
                    p2.println(subsamples[j]);
                    p2.flush();
                    p2.close();
                }


                }

   //         System.out.println(steps + "\n" + mgmprior + "\n" + oracle);
        }
        for(int i = 0; i < pri.length;i++)
        {
            printData(pri[i],results,i);
        }
    }

    public static void printData(PrintStream p, double [][][] result,int ind)
    {

        p.println("Run\tCC\tCD\tDD\tAll");
        for(int i = 0; i < result.length;i++)
        {
            p.println(i + "\t" + result[i][0][ind] + "\t" + result[i][1][ind] + "\t" + result[i][2][ind] + "\t" + result[i][3][ind]);
        }
        p.flush();
        p.close();
    }
    public static String compareGraphs(Graph truth, Graph est)
    {
   //     System.out.println("Truth: " + truth);
     //   System.out.println("Est: " + est);
        String x = "\t---\tNo_Edge";
        x+="\n---\t";
        int total = 0;
        int count = 0;
        for(Edge e: truth.getEdges())
        {
            if(est.getEdge(est.getNode(e.getNode1().getName()),est.getNode(e.getNode2().getName()))!=null)
                count++;
            else
            total++;

        }
        x+= count + "\t" + total + "\n";
        x+= "No_Edge\t";
        count = 0;
        for(Edge e: est.getEdges())
        {
            if(truth.getEdge(truth.getNode(e.getNode1().getName()),truth.getNode(e.getNode2().getName()))==null)
                count++;
        }
        x+= count;

        return x;
    }
    public static TetradMatrix [] simulatePrior(Graph g, double ap, int numExperts,TetradMatrix [] currPrior,DataSet d)
    {
        Random rand = new Random();
       TetradMatrix [] temp = new TetradMatrix[numExperts];
        List<Node> list = g.getNodes();

       A: for(int i = 0; i < numExperts; i ++)
        {
            if(currPrior[i]!=null)
            {
                temp[i] = currPrior[i];
                continue A;
            }
            TetradMatrix curr = new TetradMatrix(g.getNumNodes(),g.getNumNodes());
            for(Edge e:g.getEdges())
            {
                if(rand.nextDouble() < ap)
                {
                    int x = d.getColumn(d.getVariable(e.getNode1().getName()));
                    int y = d.getColumn(d.getVariable(e.getNode2().getName()));
                    double randDoub = rand.nextDouble()*.3+.6;
                    curr.set(x,y,randDoub);
                    curr.set(y,x,randDoub);
                }
            }
            temp[i] = curr;
        }
        return temp;
    }



    //Generates unreliable priors, with differing amount of prior information from each expert HARD PRIOR (0-1) only
    public static TetradMatrix [] simulatePriorUnreli(Graph g, int numExperts, DataSet d, int [] edges, double [] reliability,boolean pureRandom)
    {
        TetradMatrix [] temp = new TetradMatrix [numExperts];
        for(int i = 0; i < numExperts;i++)
        {
            TetradMatrix curr = new TetradMatrix(g.getNumNodes(),g.getNumNodes());
            temp[i] = curr;
        }
        Random rand = new Random();
        int totalPossible = g.getNumNodes()*(g.getNumNodes()-1)/2;
        //Priors are purely random, with no relation to the underlying ground truth graph
        if(pureRandom)
        {
            for(int i = 0; i < numExperts;i++)
            {
                edges[i] = rand.nextInt((g.getNumNodes()*(g.getNumNodes()-1))/2);
                int numEdges = edges[i];
                int trueEdges = 0;
                while (numEdges > 0) {
                    int x = rand.nextInt(g.getNumNodes());
                    int y = rand.nextInt(g.getNumNodes());
                    if (x == y) {
                        continue;
                    }
                    if (x > y) {
                        int t = x;
                        x = y;
                        y = t;
                    }
                    if(temp[i].get(x,y)==0)
                    {
                        temp[i].set(x,y,1);
                        numEdges--;
                        if(g.getEdge(g.getNode(d.getVariable(x).getName()),g.getNode(d.getVariable(y).getName()))!=null)
                            trueEdges++;
                    }
                }

                reliability[i] = trueEdges/(double)edges[i];
            }
        }
        else {
            for (int i = 0; i < edges.length; i++) {
                edges[i] = rand.nextInt(g.getNumEdges() * 2);
                while (edges[i] > totalPossible || edges[i] == 0 || edges[i] < 10)
                    edges[i] = rand.nextInt(g.getNumEdges() * 2);
            }

            int[] trueEdges = new int[numExperts];
            int[] remaining = new int[numExperts];
            int iters = 0;
            for (int i = 0; i < numExperts; i++) {
                trueEdges[i] = rand.nextInt(g.getNumEdges() + 1);
                while ((trueEdges[i] == 0 || (edges[i] - trueEdges[i]) < 0) || edges[i] - trueEdges[i] > (totalPossible - g.getNumEdges()) || ((edges[i] / (double) trueEdges[i]) < 0.2 && iters < 1000)) {
                    trueEdges[i] = rand.nextInt(g.getNumEdges() + 1);
                    iters++;
                }
                remaining[i] = edges[i] - trueEdges[i];
                reliability[i] = trueEdges[i] / (double) edges[i];
            }


            for (int i = 0; i < numExperts; i++) {
                while (trueEdges[i] > 0 || remaining[i] > 0) {
                    int x = rand.nextInt(g.getNumNodes());
                    int y = rand.nextInt(g.getNumNodes());
                    if (x == y) {
                        continue;
                    }
                    if (x > y) {
                        int t = x;
                        x = y;
                        y = t;
                    }
                    if (g.getEdge(g.getNode(d.getVariable(x).getName()), g.getNode(d.getVariable(y).getName())) != null && temp[i].get(x, y) == 0 && trueEdges[i] > 0) {
                        temp[i].set(x, y, 1);
                        trueEdges[i]--;
                    } else if (g.getEdge(g.getNode(d.getVariable(x).getName()), g.getNode(d.getVariable(y).getName())) == null && temp[i].get(x, y) == 0 && remaining[i] > 0) {
                        temp[i].set(x, y, 1);
                        remaining[i]--;
                    }
                }
            }
        }
        return temp;
    }
    //Generates unreliable priors, currently generates the entire set of priors from scratch, can't deal with potentially generated priors TODO
    public static TetradMatrix [] simulatePriorUnreli(Graph g, double ap, int numExperts, DataSet d, int numReliable, boolean correctEdges) {
        Random rand = new Random();
        TetradMatrix[] temp = new TetradMatrix[numExperts];
        for (int i = 0; i < numExperts; i++)
        {
            TetradMatrix curr = new TetradMatrix(g.getNumNodes(), g.getNumNodes());
            temp[i] = curr;
        }

        if(correctEdges) //The priors only generate information about correct edges, reliable ones just give better information
        {
            List<Node> list = g.getNodes();
            for (Edge e : g.getEdges()) {

                if (rand.nextDouble() < ap) {
                    for (int i = 0; i < numExperts; i++) {

                        int x = d.getColumn(d.getVariable(e.getNode1().getName()));
                        int y = d.getColumn(d.getVariable(e.getNode2().getName()));
                        double randDoub = 0;
                        if (i < numReliable) //This is a reliable expert
                            randDoub = rand.nextDouble() * .3 + .6;
                        else //Unreliable experts gives a uniformly distributed number from 0 to 1
                            randDoub = rand.nextDouble();
                        temp[i].set(x, y, randDoub);
                        temp[i].set(y, x, randDoub);
                    }
                }
            }
        }
        else
        {
            //Sets the priors for reliable experts (only true edges)
            int count = 0;
            for(Edge e: g.getEdges())
            {
                if(rand.nextDouble() < ap)
                {
                    for(int i = 0; i < numExperts;i++)
                    {
                        if( i < numReliable)
                        {
                            int x = d.getColumn(d.getVariable(e.getNode1().getName()));
                            int y = d.getColumn(d.getVariable(e.getNode2().getName()));
                            double randDoub = rand.nextDouble()*.3+0.6;
                            temp[i].set(x,y,randDoub);
                            temp[i].set(y,x,randDoub);

                        }
                    }
                    count++;
                }
            }

            //Sets the priors for unreliable experts (purely random edges)
            for(int j = 0; j < count;j++)
            {
                int x = rand.nextInt(d.getNumColumns());
                int y = rand.nextInt(d.getNumColumns());
                while(x==y)
                    y = rand.nextInt(d.getNumColumns());
                for(int k = numReliable; k < numExperts;k++)
                {
                    double randDoub = rand.nextDouble()*0.3+0.6;
                    temp[k].set(x,y,randDoub);
                    temp[k].set(y,x,randDoub);
                }
            }
            //For the reliable experts, we randomly pick ap percent of true edges and set them to a probability
            //For the unreliable experts, we randomly pick (ap*num true edges) and set them to a probability

        }
        return temp;
    }
    public static TetradMatrix [] simulatePriorDim(Graph g, double ap, int numExperts,TetradMatrix [] currPrior,DataSet d)
    {
        Random rand = new Random();
        TetradMatrix [] temp = new TetradMatrix[numExperts];
        List<Node> list = g.getNodes();
        boolean addNew = true;
        for(int i = 0; i < numExperts;i++)
        {
            if(currPrior[i]==null) {
                TetradMatrix curr = new TetradMatrix(g.getNumNodes(), g.getNumNodes());
                temp[i] = curr;
            }
            else {
                addNew = false;
                temp[i] = currPrior[i];
            }

        }
        for(Edge e:g.getEdges()) {
            boolean oneExpert = false;
            boolean allExperts = true;
            for(int i = 0; i < numExperts;i++)
            {
                if(temp[i].get(d.getColumn(d.getVariable(e.getNode1().getName())),d.getColumn(d.getVariable(e.getNode2().getName())))!=0)
                {
                    oneExpert = true;
                }
                else
                    allExperts = false;


            }
            if(allExperts)
                continue;
            if(oneExpert)
            {
                for (int i = 0; i < numExperts; i++) {
                    int x = d.getColumn(d.getVariable(e.getNode1().getName()));
                    int y = d.getColumn(d.getVariable(e.getNode2().getName()));
                    if (temp[i].get(x,y) == 0)
                    {
                        double randDoub = rand.nextDouble() * .3 + .6;
                        temp[i].set(x, y, randDoub);
                        temp[i].set(y, x, randDoub);
                    }
                }
            }
            else {
                if (rand.nextDouble() < ap && addNew) {
                    for (int i = 0; i < numExperts; i++) {
                        int x = d.getColumn(d.getVariable(e.getNode1().getName()));
                        int y = d.getColumn(d.getVariable(e.getNode2().getName()));
                        double randDoub = rand.nextDouble() * .3 + .6;
                        temp[i].set(x, y, randDoub);
                        temp[i].set(y, x, randDoub);
                    }
                }
            }
        }
        return temp;

    }
    public static double [][] loadPrior(File data, int numVariables) throws Exception
    {
        double [][] temp = new double[numVariables][numVariables];
        BufferedReader b = new BufferedReader(new FileReader(data));
        for(int i = 0; i < numVariables;i++)
        {
            String [] line = b.readLine().split("\t");
            for(int j = 0; j < numVariables;j++)
            {
                temp[i][j] = Double.parseDouble(line[j]);
            }
        }
        return temp;
    }
    public static Graph moralize(Graph g)
    {
        Graph temp = g.subgraph(g.getNodes());

        for(Node n: g.getNodes())
        {
            for(Node child: g.getChildren(n))
            {
                for(Node parent: g.getParents(child))
                {
                    if(n==parent)
                        continue;
                    if(temp.getEdge(temp.getNode(n.getName()),temp.getNode(parent.getName()))==null)
                        temp.addEdge(new Edge(temp.getNode(n.getName()),temp.getNode(parent.getName()),Endpoint.TAIL, Endpoint.TAIL));
                }
            }
        }
        return temp;
    }
    public static void removeMoral(Graph g)
    {
        Set<Edge> ed = g.getEdges();
        for(Edge e:ed)
        {
            if(!e.isDirected())
            {
                g.removeEdge(e);
            }
        }
    }
    public static void checkPriors(Graph g, DataSet d, TetradMatrix [] priors,int index,boolean exit)
    {
        for(int i = 0; i < priors.length;i++)
        {
            int count = 0;
            for(int j = 0; j < d.getNumColumns();j++)
            {
                for(int k = 0; k < d.getNumColumns();k++)
                {

                    if(priors[i].get(j,k)!=0)
                    {
                        if(g.getEdge(g.getNode(d.getVariable(j).getName()),g.getNode(d.getVariable(k).getName()))==null)
                        {
                            if(exit) {
                                System.out.println("Violation in prior " + i + ", for run " + index);
                                System.exit(-1);
                            }
                            count++;
                        }
                    }
                }
            }
        }
    }
}
