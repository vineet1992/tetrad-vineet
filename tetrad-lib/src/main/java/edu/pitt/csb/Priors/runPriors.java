package edu.pitt.csb.Priors;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

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
        int index = 0;
        try {
            while (index < args.length) {
                if (args[index].equals("-ns")) {
                    ns = Integer.parseInt(args[index + 1]);
                    index += 2;
                } else if (args[index].equals("-nl")) {
                    numLambdas = Integer.parseInt(args[index + 1]);
                    index += 2;
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
                } else if (args[index].equals("-t")) {
                    unreliableThreshold = Double.parseDouble(args[index + 1]);
                    index += 2;
                } else if (args[index].equals("-loocv")) {
                    loocv = true;
                    index++;
                }
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
            initLambdas[i] = i*(high-low)/numLambdas;
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
        d.removeColumn(d.getVariable("Response"));
        boolean addedDummy = false;
        if(d.isContinuous())
        {
            System.out.print("Data is only continuous, adding a discrete variable...");
            Random rand = new Random();
            DiscreteVariable temp= new DiscreteVariable("Dummy",3);
            d.addVariable(temp);
            int column = d.getColumn(d.getVariable(temp.getName()));
            for(int i = 0; i < d.getNumRows();i++)
            {
                d.setInt(i,column,rand.nextInt(temp.getNumCategories()));
            }
            System.out.println("Done");
            addedDummy = true;
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
            TetradMatrix [] priors = new TetradMatrix[numPriors];
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
                    priors[i] = new TetradMatrix(realDataPriorTest.loadPrior(new File("temp_2.txt"),d.getNumColumns()));
                }
                else
                {
                    priors[i] = new TetradMatrix(realDataPriorTest.loadPrior(new File(currFile),d.getNumColumns()));
                }
            }
            File t = new File("temp.txt");
            t.deleteOnExit();
            t = new File("temp_2.txt");
            if(t.exists())
                t.deleteOnExit();

            mgmPriors m = new mgmPriors(ns,initLambdas,d,priors);
            Graph g = m.runPriors();
            printAllResults(g,m,runName,fileMap);
        }
        catch(Exception e)
        {
            System.err.println("");
            System.exit(-1);
        }

    }


    public static void printAllResults(Graph g, mgmPriors m, String runName, HashMap<Integer,String> map) throws Exception
    {
        double[] weights = m.normalizedExpertWeights;
        double[] pValues = m.pValues;
        double[] normalizedTao = m.normalizedTao;

        PrintStream out = new PrintStream(runName + "/Graph.txt");
        out.println(g);
        out.flush();
        out.close();
        out = new PrintStream(runName + "/Prior_Scores.txt");
        out.println("Name\tPrior_Weight\tCorrected p-Value\tNormalized Deviance Score");
        for (int i = 0; i < weights.length; i++) {
            out.println(map.get(i) + "\t" + weights[i] + "\t" + pValues[i] + "\t" + normalizedTao[i]);
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
            prior[i][j] = score;
            prior[j][i] = score;
        }
        b.close();
        return prior;
    }
}
