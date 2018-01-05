package edu.pitt.csb.mgm;

import edu.cmu.tetrad.algcomparison.simulation.ContinuousLinearGaussianSemSimulation;
import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
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


        boolean numEdgesRandom = false;
        boolean mixed = false;
        int numVariables = 25;
        int numEdges = 2*numVariables;
        int faithful = 0;
        int sampleSize = 1000;
        int numRuns = 100;
        int numCategories = 4;
        double percentDiscrete = 50;
        String directory = ".";
        int index = 0;
        while (index < args.length) {
           if (args[index].equals("-d")) {
                directory = args[index + 1];
                index += 2;
            } else if(args[index].equals("-f")) {
               faithful = 1;
               index++;
           }  else if (args[index].equals("-nc")) {
                numCategories = Integer.parseInt(args[index + 1]);
                index += 2;
            } else if (args[index].equals("-v")) {
                numVariables = Integer.parseInt(args[index + 1]);
                index += 2;
            } else if (args[index].equals("-e")) {
                numEdges = Integer.parseInt(args[index + 1]);
                index += 2;
            } else if (args[index].equals("-er")) {
             numEdgesRandom = true;
               index ++;
           }else if (args[index].equals("-s")) {
                sampleSize = Integer.parseInt(args[index + 1]);
                index += 2;
            }  else if (args[index].equals("-r")) {
                numRuns = Integer.parseInt(args[index + 1]);
                index += 2;
            } else if (args[index].equals("-m")) {
                mixed = true;
                index += 1;
            }
        }
        File gFile = new File("Graphs");
        File dFile = new File("Data");
        if (!gFile.isDirectory())
            gFile.mkdir();
        if (!dFile.isDirectory())
            dFile.mkdir();

       A: for(int i = 0; i < numRuns;i++) {

            System.out.println(i);


            if(mixed)
            {
                PrintStream out = new PrintStream(directory + "/Graphs/Graph_" + i + ".txt");
                MixedLeeHastieSimulation c = new MixedLeeHastieSimulation();
                Parameters p = new Parameters();
                p.setValue("numMeasures",numVariables);
                p.setValue("sampleSize",sampleSize);
                p.setValue("percentDiscreteForMixedSimulation",percentDiscrete);
                p.setValue("numCategories",numCategories);
                NormalDistribution n = new NormalDistribution(numVariables*2,numVariables/2);
                if(numEdgesRandom)
                     p.setValue("numEdges",(int)n.sample());
                else
                    p.setValue("numEdges",numEdges);

                c.simulate(p);
                out.println(c.getTrueGraph());
                out.flush();
                out.close();
                out = new PrintStream(directory + "/Graphs/Adj_Mat_" + i + ".txt");
                printAdjMat(c.getTrueGraph(),out);
                out = new PrintStream(directory + "/Data/Data_" + i + ".txt");
                out.println(c.getDataSet(0));
                out.flush();
                out.close();

            }
            else
            {
                PrintStream out = new PrintStream(directory + "/Graphs/Graph_" + i + ".txt");
                ContinuousLinearGaussianSemSimulation c = new ContinuousLinearGaussianSemSimulation();
                Parameters p = new Parameters();
                p.setValue("numMeasures",numVariables);
                p.setValue("sampleSize",sampleSize);
                p.setValue("percentDiscreteForMixedSimulation",percentDiscrete);
                p.setValue("numCategories",numCategories);
                NormalDistribution n = new NormalDistribution(numVariables,numVariables/4);
                if(numEdgesRandom)
                     p.setValue("numEdges",(int)n.sample());
                else
                    p.setValue("numEdges",numEdges);
                p.setValue("faithful",faithful);
                c.simulate(p);
                if(c.getDataSet(0)==null) {
                    i--;
                    continue A;
                }
                out.println(c.getTrueGraph());
                out.flush();
                out.close();
                out = new PrintStream(directory + "/Graphs/Adj_Mat_" + i + ".txt");
                printAdjMat(c.getTrueGraph(),out);
                out = new PrintStream(directory + "/Data/Covariance_" + i + ".txt");
                printCovMat(c.cov,out);
                out = new PrintStream(directory + "/Data/Data_" + i + ".txt");
                out.println(c.getDataSet(0));
                out.flush();
                out.close();

            }
        }
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
    public static void printCovMat(TetradMatrix cov, PrintStream out)
    {
        for(int i = 0; i < cov.rows();i++)
        {
            for(int j = 0; j < cov.columns();j++)
            {
                    if(j==cov.columns()-1)
                        out.println(cov.get(i,j));
                    else
                        out.print(cov.get(i,j) + "\t");

            }
        }
        out.flush();
        out.close();
    }
}
