package edu.pitt.csb.Priors;

import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;

/**
 * This is a packaged .jar file to run MGM Priors. The jar is expecting an expression dataset
 * and either a directory containing all prior sources
 * The prior sources can either be an N x N matrix of values corresponding to the variables in the
 * expression dataset, or a .sif style file, with src dest probability
 *
 * If the .sif files have only two columns, then it is assumed that the probability is 1 for each edge
 */

//TODO Change the loadPriors function here and in realDataPriorTest to treat missing value as null instead of 0
public class runPriors {

    private static boolean verbose = false;
    public static void main(String [] args) {

        String priorDirectory = "pathway_lists";
        boolean priorMatrices = true;
        String dataFile = "";
        String runName = "";
        boolean excludeUnreliable = false;
        double unreliableThreshold = 0.05;
        int ns = 20;
        int numLambdas = 40;
        double low = 0.05;
        double high = 0.95;
        boolean loocv = false;
        boolean makeScores = false;
        int index = 0;
        List<String> toRemove = new ArrayList<String>();
        try {
            while (index < args.length) {
                if (args[index].equals("-ns")) {
                    ns = Integer.parseInt(args[index + 1]);
                    index += 2;
                } else if (args[index].equals("-nl")) {
                    numLambdas = Integer.parseInt(args[index + 1]);
                    index += 2;
                } else if(args[index].equals("-llow")){
                    low = Double.parseDouble(args[index+1]);
                    index+=2;
                }else if (args[index].equals("-lhigh")){
                    high = Double.parseDouble(args[index+1]);
                    index+=2;
                } else if(args[index].equals("-run")) {
                    runName = args[index+1];
                    index+=2;
                }
                else if (args[index].equals("-priors")) {
                    priorDirectory = args[index + 1];
                    index += 2;
                } else if (args[index].equals("-sif")) {
                    priorMatrices = false;
                    index++;
                } else if (args[index].equals("-data")) {
                    dataFile = args[index + 1];
                    index += 2;
                } else if (args[index].equals("-ex")) {
                    excludeUnreliable = true;
                    index++;
                }
                else if(args[index].equals("-rm"))
                {
                    int count = index + 1;
                     while(count < args.length && !args[count].startsWith("-"))
                     {
                         toRemove.add(args[count]);
                         count++;
                     }
                     index = count;
                } else if (args[index].equals("-t")) {
                    unreliableThreshold = Double.parseDouble(args[index + 1]);
                    index += 2;
                } else if (args[index].equals("-loocv")) {
                    loocv = true;
                    index++;
                }
                else if(args[index].equals("-makeScores"))
                {
                    makeScores = true;
                    index++;
                }
                else if(args[index].equals("-v"))
                {
                    verbose = true;
                    index++;
                }
                else
                    index++;
            }
        }
        catch(ArrayIndexOutOfBoundsException e)
        {
            System.err.println("Command line arguments not specified properly at argument: " + args[index]);
            e.printStackTrace();
            System.exit(-1);
        }
        catch(NumberFormatException e)
        {
            System.err.println("Expected a number for element: " + args[index]);
            e.printStackTrace();
            System.exit(-1);
        }
        catch(Exception e)
        {
            System.err.println("Double check command line arguments");
            e.printStackTrace();
            System.exit(-1);
        }
        double [] initLambdas = new double[numLambdas];
        for(int i = 0; i < numLambdas;i++)
        {
            initLambdas[i] = low + i*(high-low)/numLambdas;
        }
        DataSet d = null;
        try {
            d = MixedUtils.loadDataSet2(dataFile);
        }
        catch(Exception e)
        {
            System.err.println("Error loading in data file");
            e.printStackTrace();
            System.exit(-1);
        }

        if(verbose)
        {
            System.out.println("Removing Variables... " + toRemove);
        }
        //System.out.println(toRemove);
        for(String s:toRemove)
        {
            d.removeColumn(d.getVariable(s));
        }
        boolean addedDummy = false;
        if(d.isContinuous())
        {
            System.out.print("Data is only continuous, adding a discrete variable...");
            Random rand = new Random();
            DiscreteVariable temp= new DiscreteVariable("Dummy",2);
            d.addVariable(temp);
            int column = d.getColumn(d.getVariable(temp.getName()));
            for(int i = 0; i < d.getNumRows();i++)
            {
                d.setInt(i,column,rand.nextInt(temp.getNumCategories()));
            }
            System.out.println("Done");
            addedDummy = true;
        }
        if(verbose)
        {
            System.out.println("Full Dataset: " + d);
            System.out.println("Is DataSet Mixed? " + d.isMixed());
        }
        try{
            File x = new File(runName);
            if(x.exists())
            {
                if(x.isDirectory())
                {
                    System.err.println("Please specify a run name that does not have an existing directory");
                    System.exit(-1);
                }
                else
                {
                    System.err.println("Please specify a run name that isn't a file");
                    System.exit(-1);
                }

            }
            x.mkdir();

            File f = new File(priorDirectory);
            if(!f.isDirectory())
            {
                System.err.println("Prior sources directory does not exist");
                System.exit(-1);
            }
            HashMap<Integer,String> fileMap = new HashMap<Integer,String>();
            int numPriors = f.listFiles().length;
            SparseDoubleMatrix2D[] priors = new SparseDoubleMatrix2D[numPriors];
            for(int i = 0;i < f.listFiles().length;i++)
            {
                fileMap.put(i,f.listFiles()[i].getName());
                String currFile = f.listFiles()[i].getPath();
                if(!priorMatrices)
                {
                    PrintStream out = new PrintStream("temp.txt");
                    createPrior(f.listFiles()[i],out,d.getVariableNames());
                    currFile = "temp.txt";
                }
                if(addedDummy)
                {
                    addLines(new File(currFile));
                    priors[i] = new SparseDoubleMatrix2D(realDataPriorTest.loadPrior(new File("temp_2.txt"),d.getNumColumns()));
                }
                else
                {
                    priors[i] = new SparseDoubleMatrix2D(realDataPriorTest.loadPrior(new File(currFile),d.getNumColumns()));
                }
            }
           /* if(verbose)
            {
                for(int i = 0; i < numPriors;i++)
                {
                    System.out.println(f.listFiles()[i].getName() + " Prior");
                    for(int j = 0; j < d.getNumColumns();j++)
                    {
                        System.out.print(d.getVariableNames().get(j) + " ");
                    }

                    System.out.println(priors[i]);
                }
            }*/


            File t = new File("temp.txt");
            t.deleteOnExit();
            t = new File("temp_2.txt");
            if(t.exists())
                t.deleteOnExit();

            int b = (int) Math.floor(10 * Math.sqrt(d.getNumRows()));
            if (b >= d.getNumRows())
                b = d.getNumRows() / 2;
            int [][] samps = new int[ns][];
            boolean done = false;
            int attempts = 10000;
            DataSet[] subsamples = new DataSet[ns];
            System.out.print("Generating subsamples and ensuring variance...");
            while(!done && attempts > 0) {
                done = true;
                if (loocv)
                    samps = StabilityUtils.generateSubsamples(d.getNumRows());
                else
                    samps = StabilityUtils.subSampleNoReplacement(d.getNumRows(), b, ns);

                for (int j = 0; j < ns; j++) {
                    subsamples[j] = d.subsetRows(samps[j]);
                    int col = checkForVariance(subsamples[j],d);
                    if(col!=-1)
                    {
                        if(loocv)
                        {
                            System.out.println("Can't perform Leave-one-out Cross Validation...leaving out sample " + j + " makes " + d.getVariable(col) + " have no variance");
                            System.exit(-1);
                        }
                        else {
                            attempts--;
                            done = false;
                        }
                    }
                }
            }

            subsamples = null;
            System.out.println("Done");
            System.out.print("Generating Lambda Params...");
            mgmPriors m = new mgmPriors(ns,initLambdas,d,priors,samps,verbose);
            System.out.println("Done");
            if(makeScores) {
                m.makeEdgeScores();
            }
            System.out.print("Running piMGM...");
            Graph g = m.runPriors();
            System.out.println("Done");
            System.out.print("Printing Results...");
            printAllResults(g,m,runName,fileMap);
            if(makeScores)
            {
                double [][] scores = m.edgeScores;
                printScores(scores,d,runName);
            }
            System.out.println("Done");
        }
        catch(Exception e)
        {
            System.err.println("Unknown Error");
            e.printStackTrace();
            System.exit(-1);
        }

    }


    public static synchronized int checkForVariance(DataSet d, DataSet full)
    {
        TetradMatrix t = d.getDoubleData();
        for(int i = 0; i < d.getNumColumns();i++)
        {
            if(d.getVariable(i)instanceof ContinuousVariable)
            {
                double [] curr = t.getColumn(i).toArray();
                curr = StatUtils.standardizeData(curr);
                double var = StatUtils.variance(curr);
                if(var <= 0.000001)
                    return i;

            }
            else
            {
                HashSet<Integer> cats = new HashSet<Integer>();
                for(int j = 0; j < full.getNumRows();j++)
                {
                    cats.add(full.getInt(j,i));
                }
                for(int j = 0; j < d.getNumRows();j++)
                {
                    cats.remove(d.getInt(j,i));
                }
                if(!cats.isEmpty())
                    return i;
            }
        }
        return -1;
    }

    public static void printScores(double [][] scores, DataSet d, String runName) throws Exception
    {
        PrintStream out = new PrintStream(runName + "/Edge_Scores.txt");
        for(int i = 0; i < d.getNumColumns();i++)
        {
            if(i==d.getNumColumns()-1)
                out.println(d.getVariable(i));
            else
                out.print(d.getVariable(i) + "\t");
        }
        for(int i = 0; i < scores.length;i++)
        {
            out.print(d.getVariable(i) + "\t");
            for(int j = 0; j < scores[i].length;j++)
            {
                if(j==scores[i].length-1)
                    out.println(scores[i][j]);
                else
                    out.print(scores[i][j] + "\t");
            }
        }
        out.flush();
        out.close();
    }

    public static void printAllResults(Graph g, mgmPriors m, String runName, HashMap<Integer,String> map) throws Exception
    {
        double[] weights = m.normalizedExpertWeights;
        double[] pValues = m.pValues;
        double[] normalizedTao = m.normalizedTao;
        double [] uncorrectedPVals = m.uncorrectedPValues;

        PrintStream out = new PrintStream(runName + "/Graph.txt");
        out.println(g);
        out.flush();
        out.close();
        out = new PrintStream(runName + "/Prior_Scores.txt");
        out.println("Name\tPrior_Weight\tCorrected p-Value\tUncorrected p-Value\tNormalized Deviance Score");
        for (int i = 0; i < weights.length; i++) {
            out.println(map.get(i) + "\t" + weights[i] + "\t" + pValues[i] + "\t" + uncorrectedPVals[i] + "\t" + normalizedTao[i]);
        }

        out.flush();
        out.close();

    }
    //Adds an extra column and row of zeroes to the file if there was a dummy variable added to the data
    public static void addLines(File f) throws Exception
    {
        PrintStream out = new PrintStream("temp_2.txt");
        BufferedReader b = new BufferedReader(new FileReader(f));

        int length = -1;
        while(b.ready())
        {
            String line = b.readLine();
            if(line.endsWith("\t")) {
                out.println(line + "0");
                length = line.split("\t").length;
            }
            else {
                length = line.split("\t").length+1;
                out.println(line + "\t0");
            }
        }
        for(int i = 0; i < length;i++)
        {
            if(i==length-1)
                out.println("0");
            else
                out.print("0\t");
        }
        out.flush();
        out.close();
        b.close();
    }

    public static void createPrior(File pathway, PrintStream out, List<String> vars) throws Exception
    {
        //loop through pathway file, and add elements to a double [][], then print it all out to the file
        double [][] prior = fileToPrior(pathway, vars);
        out.println(pathway.getName());
        for(int i = 0; i < prior.length;i++)
        {
            for(int j = 0; j < prior[i].length;j++)
            {
                if(j==prior[i].length-1)
                    out.println(prior[i][j]);
                else
                    out.print(prior[i][j] + "\t");
            }
        }
        out.flush();
        out.close();

    }
    public static double [][] fileToPrior(File pathway, List<String>vars) throws Exception
    {
        double [][] prior = new double[vars.size()][vars.size()];
        BufferedReader b = new BufferedReader(new FileReader(pathway.getAbsolutePath()));
        if(verbose)
        {
            System.out.println("Parsing prior for " + pathway.getName());
            System.out.println("Variables: " + vars);
        }
        b.readLine();//eat the title
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            int i = vars.indexOf(line[0]);
            int j = vars.indexOf(line[1]);
            double score = 1;
            if(line.length==3)
                score = Double.parseDouble(line[2]);
            if(i==-1||j==-1)
                continue;
            if(verbose)
            {
                System.out.println("Adding probability: " + score + " for edge " + line[0] + ":" + i + ", " + line[1] + ":" + j);
            }
            prior[i][j] = score;
            prior[j][i] = score;
        }
        b.close();
        return prior;
    }
}
