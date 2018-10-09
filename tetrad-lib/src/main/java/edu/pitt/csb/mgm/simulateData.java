package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.graph.RandomForward;
import edu.cmu.tetrad.algcomparison.graph.RandomGraph;
import edu.cmu.tetrad.algcomparison.graph.RandomGraphUniform;
import edu.cmu.tetrad.algcomparison.graph.ScaleFree;
import edu.cmu.tetrad.algcomparison.simulation.*;
import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.File;
import java.io.PrintStream;
import java.util.List;



/**
 * Created by vinee_000 on 9/12/2017.
 */
public class simulateData {
    public static void main(String [] args)throws Exception {


       int numVariables = 50; //Number of variables in the dataset/graph
        int sampleSize = 1000; //Sample size of the dataset
        int numRuns = 100; //Number of dataset/graph pairs to generate
        int minCategories = 2; //Mininum number of categories for the discrete variables
        int maxCategories = 5; //Maximum number of categories for the discrete variables
        double percentDiscrete = 50; //Percent of variables that are discrete
        String directory = "."; //Directory to write files to
        double edgeCoefLow = 0.05; //Low range of edge coefficient, by default the values of these coefficients will be uniformly distributed from [edgeCoefLow to edgeCoefHigh]
        double edgeCoefHigh = 3; //High range of edge coefficients
        double varLow = 0.5;//Low range of variance, by default the variance will be uniformly distributed from [varLow to varHigh]
        double varHigh = 3; //High range of variance
        int numLatents = 0; //Number of latent variables in the graph
        int maxInDegree = 10; //Maximum number of parents for any child
        int maxOutDegree = 10;// Maximum number of children for any parent
        boolean connected = false;//Graphs must be connected (every node is reachable from every other node)
        Graph initGraph = null; //Generate num runs datasets using this as the graph
        boolean differentGraphs = false; //Should each run use a different graph?
        double mean = -1; //mean of the distribution to determine number of edges
        double sd = -1; //Standard deviation of that distribution
        String simType = "LH"; //Simulation Type, choices are: LH, CG, SEM, BAYES
        String graphType = "U"; //Graph Generation Type: U (Uniform), S (Scale-Free), F (Random-Forward)
        double alpha = 0.33; // probability of adding a new vertex and edge from v to existing vertex
        double beta = 0.33; //probability of adding an edge between two existing vertices
        double deltaIn = 0.1;  //larger value indicates that nodes with high indegree will be more likely to get more edges into
        double deltaOut = 0.1; //larger value indicates that nodes with high outdegree will be more likely to get more outgoing edges
        int index = 0;
        while (index < args.length) {
           if (args[index].equals("-d")) {
                directory = args[index + 1];
                index += 2;
            } else if(args[index].equals("-graph"))
           {
              initGraph = GraphUtils.loadGraphTxt(new File(args[index+1]));
              differentGraphs = false;
              index+=2;
           }  else if(args[index].equals("-diff")){
               differentGraphs = true;
               index++;
           } else if(args[index].equals("-connected")){
               connected = true;
               index++;
           } else if(args[index].equals("-nl")) {
               numLatents = Integer.parseInt(args[index+1]);
               index+=2;
           } else if(args[index].equals("-simulation")) {
                simType = args[index+1];
               index+=2;
           } else if(args[index].equals("-graphSim")) {
               graphType = args[index+1];
               index+=2;
           }  else if(args[index].equals("-maxIn")) {
               maxInDegree = Integer.parseInt(args[index+1]);
               index+=2;
           } else if(args[index].equals("-maxOut")) {
               maxOutDegree = Integer.parseInt(args[index+1]);
               index+=2;
           } else if (args[index].equals("-alpha")) {
               alpha = Double.parseDouble(args[index + 1]);
               index += 2;
           }else if (args[index].equals("-beta")) {
               beta = Double.parseDouble(args[index + 1]);
               index += 2;
           }else if (args[index].equals("-deltaIn")) {
               deltaIn = Double.parseDouble(args[index + 1]);
               index += 2;
           }else if (args[index].equals("-deltaOut")) {
               deltaOut = Double.parseDouble(args[index + 1]);
               index += 2;
           }else if (args[index].equals("-coefLow")) {
               edgeCoefLow = Double.parseDouble(args[index + 1]);
               index += 2;
           } else if (args[index].equals("-coefHigh")) {
               edgeCoefHigh = Double.parseDouble(args[index + 1]);
               index += 2;
           }else if (args[index].equals("-varLow")) {
               varLow = Double.parseDouble(args[index + 1]);
               index += 2;
           }else if (args[index].equals("-varHigh")) {
               varHigh = Double.parseDouble(args[index + 1]);
               index += 2;
           }else if (args[index].equals("-percentDisc")) {
               percentDiscrete = Double.parseDouble(args[index + 1]);
               index += 2;
           }else if (args[index].equals("-minC")) {
                minCategories = Integer.parseInt(args[index + 1]);
                index += 2;
            } else if (args[index].equals("-maxC")) {
               maxCategories = Integer.parseInt(args[index + 1]);
               index += 2;
           } else if (args[index].equals("-v")) {
                numVariables = Integer.parseInt(args[index + 1]);
                index += 2;
            } else if (args[index].equals("-er")) {
             mean = Double.parseDouble(args[index+1]);
             sd = Double.parseDouble(args[index+2]);
               index +=3 ;
           }else if (args[index].equals("-s")) {
                sampleSize = Integer.parseInt(args[index + 1]);
                index += 2;
            }  else if (args[index].equals("-r")) {
                numRuns = Integer.parseInt(args[index + 1]);
                index += 2;
            }
            else
           {

               throw new IllegalArgumentException("Commnd line switch " + args[index] + " not supported");

           }
        }


        if(mean <= 0 || sd <= 0)
        {
            mean = 2*numVariables; //mean of the distribution to determine number of edges
            sd = numVariables/2; //Standard deviation of that distribution
        }
       File direc = new File(directory);
        if(direc.isFile())
        {
            System.err.println("Cannot create a directory named: " + directory + ", file already exists");
            System.exit(-1);
        }
        if(graphType.equals("S") && (alpha+beta)>1)
        {
            System.err.println("Alpha + beta cannot be greater than one for a scale-free graph");
            System.exit(-1);
        }
        if(!direc.exists())
            direc.mkdir();
        File gFile = new File(directory + "/Graphs");
        File dFile = new File(directory + "/Data");
        File pFile = new File(directory + "/Parametric_Models/");
        if (!gFile.isDirectory())
            gFile.mkdir();
        if (!dFile.isDirectory())
            dFile.mkdir();
        if(!pFile.isDirectory())
            pFile.mkdir();


                Simulation c = stringToSimulation(simType, graphType);
                if(initGraph!=null)
                    c.setInitialGraph(initGraph);
                double meanDegree = mean*2/(double)numVariables;
                double devDegree = sd*2/(double)numVariables;



                Parameters p = new Parameters();
                p.set("numMeasures",numVariables);
                p.set("sampleSize",sampleSize);
                p.set("percentDiscrete",percentDiscrete);
                p.set("differentGraphs",differentGraphs);
                p.set("minCategories",minCategories);
                p.set("maxCategories",maxCategories);
                p.set("coefLow",edgeCoefLow);
                p.set("coefHigh",edgeCoefHigh);
                p.set("varLow",varLow);
                p.set("varHigh",varHigh);
                p.set("numLatents",numLatents);
                p.set("maxInDegree",maxInDegree);
                p.set("maxOutDegree",maxOutDegree);
                p.set("connected",connected);
                p.set("meanDegree",meanDegree);
                p.set("devDegree",devDegree);
                p.set("numRuns",numRuns);
                p.set("scaleFreeAlpha",alpha);
                p.set("scaleFreeBeta",beta);
                p.set("scaleFreeDeltaIn",deltaIn);
                p.set("scaleFreeDeltaOut",deltaOut);


                c.createData(p);
                for(int i = 0; i < numRuns;i++) {
                    PrintStream out = new PrintStream(directory + "/Graphs/Graph_" + i + ".txt");
                    Graph curr = c.getTrueGraph(i);
                    DataSet data = (DataSet)c.getDataModel(i);


                 out.println(curr);
                out.flush();
                out.close();
                out = new PrintStream(directory + "/Graphs/Adj_Mat_" + i + ".txt");
                printAdjMat(curr,out);
                out = new PrintStream(directory + "/Data/Data_" + i + ".txt");
                out.println(data);
                out.flush();
                out.close();


                if(percentDiscrete==0 || simType.equals("SEM")) {
                    ICovarianceMatrix cov = DataUtils.getCovMatrix(data);
                    out = new PrintStream(directory + "/Data/Covariance_" + i + ".txt");
                    printCovMat(cov, out);
                }
                out = new PrintStream(directory + "/Parametric_Models/Model_" + i + ".txt");
                out.println(c.getInstantiatedModel(i));
                out.flush();
                out.close();
            }
        }

    public static Simulation stringToSimulation(String sim, String graphType)
    {
        RandomGraph rg;
        if(graphType.equals("U"))
        {
            rg = new RandomGraphUniform();
        }else if(graphType.equals("S"))
        {
            rg = new ScaleFree();
        }else
        {
            rg = new RandomForward();
        }
        switch(sim){
            case("LH"):
                return new LeeHastieSimulation(rg);
            case("CG"):
                return new ConditionalGaussianSimulation(rg);
            case("SEM"):
                return new SemSimulation(rg);
            case("BAYES"):
                return new BayesNetSimulation(rg);
        }
        return null;
    }
    public static void printAdjMat(Graph g, PrintStream out)
    {
        List<Node> nodes = g.getNodes();
        for(int i = 0; i < nodes.size();i++)
        {
            for(int j = 0; j < nodes.size();j++)
            {
                if(g.getEdge(nodes.get(i),nodes.get(j))!=null) {
                    if(j==nodes.size()-1)
                        out.println(1);
                    else
                        out.print(1 + "\t");
                }
                else {
                    if(j==nodes.size()-1)
                        out.println(0);
                    else
                         out.print(0 + "\t");
                }
            }
        }
        out.flush();
        out.close();
    }
    public static void printCovMat(ICovarianceMatrix cov, PrintStream out)
    {
        for(int i = 0; i < cov.getDimension();i++)
        {
            for(int j = 0; j < cov.getDimension();j++)
            {
                    if(j==cov.getDimension()-1)
                        out.println(cov.getValue(i,j));
                    else
                        out.print(cov.getValue(i,j) + "\t");

            }
        }
        out.flush();
        out.close();
    }
}
