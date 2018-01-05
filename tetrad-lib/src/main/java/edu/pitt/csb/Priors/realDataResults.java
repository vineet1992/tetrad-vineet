package edu.pitt.csb.Priors;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by vinee_000 on 12/11/2017.
 */
public class realDataResults {
    public static void main(String [] args)throws Exception
    {
        //Two metrics for each method: Output Two column data with PAM50 Genes and Score
        //Second data should be single file with precision and recall for edge recovery in each of the nine pathways
        //String [] methods = {"No_Prior","Irrelevant_Prior","Relevant_Prior"};
        String [] methods = {"No_Prior"};
        int numPathways = 9;
        PrintStream out2 = new PrintStream("Results/Pathway_Results.txt");
        PrintStream [] ps = new PrintStream[methods.length];
        out2.print("Pathway_Name\t");
        for(int i = 0; i < methods.length;i++)
        {
            ps[i] = new PrintStream("AUC_" + methods[i] + ".txt");
            ps[i].println("Edge\tScore\tTrue_Edge\tPathway");
            BufferedReader b = new BufferedReader(new FileReader("Results/Stabilities_" + methods[i] + ".txt"));
            HashMap<String,Double>tempStabs =  getStabs(b);
            PrintStream out = new PrintStream("Results/PAM50_Predictions_" + methods[i] + ".txt");
            out.println("Neighbor\tScore");
            Graph g = GraphUtils.loadGraphTxt(new File("Results/Graph_" + methods[i] + ".txt"));
            for(Node n: g.getAdjacentNodes(g.getNode("Subtype")))
            {
                out.println(n + "\t" + tempStabs.get(n.getName()));
            }
            out.flush();
            out.close();
            if(i!=methods.length-1)
                out2.print(methods[i] + "_PREC\t" + methods[i] + "_REC\t");
            else
                out2.println(methods[i] + "_PREC\t" + methods[i] + "_REC");

        }
        DataSet data = MixedUtils.loadDataSet2("genes_with_clinical.txt");
        for(int j = 0; j < numPathways;j++)
        {
            out2.print(j + "\t");
             double [][] pathway = loadPathway(data, new BufferedReader(new FileReader("prior_sources/Prior_" + j + ".txt")));
            Set<String> genes = getGenes(data,pathway); //Contains a list of the genes in this pathway
            List<String> trueEdges = getEdges(data,pathway);
            System.out.println(trueEdges);
            for(int i = 0; i < methods.length;i++) {

                Graph g = GraphUtils.loadGraphTxt(new File("Results/Graph_" + methods[i] + ".txt"));
                List<String> predictEdges = getEdges(genes, g);
                if(i==methods.length-1)
                    out2.println(getPrecision(predictEdges,trueEdges) + "\t" + getRecall(predictEdges,trueEdges));
                else
                    out2.print(getPrecision(predictEdges,trueEdges)+"\t" + getRecall(predictEdges,trueEdges) + "\t");


                BufferedReader b = new BufferedReader(new FileReader("Results/Stabilities_" + methods[i] + ".txt"));
                b.readLine();
                double [][] stabs = getAllStabs(b,data);
                printAUC(ps[i],predictEdges,trueEdges,j,stabs,data);
            }
        }
        for(int i = 0; i < ps.length;i++)
        {
            ps[i].flush();
            ps[i].close();
        }
        out2.flush();
        out2.close();

    }
    private static double [][] getAllStabs(BufferedReader b, DataSet data) throws Exception
    {
        double [][] stabs = new double[data.getNumColumns()][data.getNumColumns()];
        int x = 0;
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            for(int i = 1; i < line.length;i++)
            {
                stabs[x][i-1] = Double.parseDouble(line[i]);

            }
            x++;

        }
        return stabs;
    }
    private static void printAUC(PrintStream out, List<String> predicted, List<String> truth, int j, double [][] stabs, DataSet data)
    {
        for(int i = 0; i < predicted.size();i++)
        {
            if(truth.contains(predicted.get(i))) {
                String [] temp = predicted.get(i).split(",");

                out.println(predicted.get(i) + "\t" +stabs[data.getColumn(data.getVariable(temp[0]))][data.getColumn(data.getVariable(temp[1]))] + "\t" + 1 + "\t" + j);
            }
            else
            {
                String [] temp = predicted.get(i).split(",");
                if(truth.contains(temp[1] +"," + temp[0]))
                {
                    out.println(predicted.get(i) + "\t" + stabs[data.getColumn(data.getVariable(temp[1]))][data.getColumn(data.getVariable(temp[0]))] + "\t" + 1 + "\t" + j);

                }
                else
                {
                    out.println(predicted.get(i) + "\t" + stabs[data.getColumn(data.getVariable(temp[1]))][data.getColumn(data.getVariable(temp[0]))] + "\t" + 0 + "\t" + j);

                }
            }
        }
        for(int i = 0; i < truth.size();i++) //Now, whichever edges are not in the predicted set get a score of 0.0
        {
            String [] temp = truth.get(i).split(",");
            if(!predicted.contains(truth.get(i)) && !predicted.contains(temp[1] + "," + temp[0]))
            {
                out.println(truth.get(i) + "\t" + 0.0 + "\t" + 1 + "\t" + j);
            }

        }
        out.flush();
    }
    private static double getRecall(List<String>predictEdges,List<String>trueEdges)
    {
        double tp = 0;
        double fn = 0;
        for(int i = 0; i < trueEdges.size();i++)
        {
            if(predictEdges.contains(trueEdges.get(i)))
                tp++;
            else
            {
                String [] temp = trueEdges.get(i).split(",");
                if(predictEdges.contains(temp[1] + "," + temp[0]))
                    tp++;
                else
                    fn++;
            }
        }
        return tp/(tp+fn);
    }
    private static double getPrecision(List<String> predictEdges, List<String> trueEdges)
    {
        double tp = 0;
        double fp = 0;
        for(int i = 0; i < predictEdges.size();i++)
        {
            if(trueEdges.contains(predictEdges.get(i)))
                tp++;
            else
            {
                String [] temp = predictEdges.get(i).split(",");
                if(trueEdges.contains(temp[1]+","+temp[0]))
                    tp++;
                else
                    fp++;
            }
        }
        return tp/(tp+fp);
    }
    private static List<String> getEdges(Set<String> genes,Graph g)
    {
        List<String> temp = new ArrayList<String>();
        for(String x:genes)
        {
            for(Node n: g.getAdjacentNodes(g.getNode(x)))
            {
                if(genes.contains(n.getName()))
                {
                    if(!temp.contains(x + "," + n.getName()) && !temp.contains(n.getName() + "," + x))
                    temp.add(x + "," + n.getName());
                }
            }

        }
        return temp;
    }
    private static List<String> getEdges(DataSet data, double [][] pathway)
    {
        List<String> temp = new ArrayList<String>();
        for(int i = 0; i < pathway.length;i++)
        {
            for(int j = i+1; j < pathway[i].length;j++)
            {
                if(pathway[i][j]!=0)
                {
                    if(!temp.contains(data.getVariable(i).getName() + "," + data.getVariable(j).getName()) && !temp.contains(data.getVariable(i).getName() + "," + data.getVariable(j).getName()))
                    temp.add(data.getVariable(i).getName() + "," + data.getVariable(j).getName());
                }
            }
        }
        return temp;
    }
    private static Set<String> getGenes(DataSet data, double [][] pathway)
    {
        Set<String> temp = new HashSet<String>();
        for(int i = 0; i < pathway.length;i++)
        {
            for(int j = i+1; j < pathway[i].length;j++)
            {
                if(pathway[i][j]!=0)
                {
                    temp.add(data.getVariable(i).getName());
                    temp.add(data.getVariable(j).getName());
                }
            }
        }
        return temp;
    }
    private static double [][] loadPathway(DataSet data, BufferedReader b)throws Exception
    {
        b.readLine();
        double [][] prior = new double[data.getNumColumns()][data.getNumColumns()];
        int x = 0;
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            for(int j = 0; j < line.length;j++)
            {
                prior[x][j] = Double.parseDouble(line[j]);
            }
            x++;
        }
        b.close();
        return prior;


    }
    private static HashMap<String,Double> getStabs(BufferedReader b) throws Exception
    {
        HashMap<String,Double> temp = new HashMap<String,Double>();
        String [] line = b.readLine().split("\t");
        int index = -1;
        for(int i = 0; i < line.length;i++)
        {
            if(line[i].equals("Subtype"))
                index = i+1;
        }
        while(b.ready())
        {
           line = b.readLine().split("\t");
           temp.put(line[0], Double.parseDouble(line[index]));
        }
        return temp;
    }
}
