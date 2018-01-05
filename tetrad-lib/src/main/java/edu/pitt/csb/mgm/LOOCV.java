package edu.pitt.csb.mgm;


import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.Regression;
import edu.cmu.tetrad.regression.RegressionUtils;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;
import edu.cmu.tetrad.util.TetradMatrix;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
/**
 * Created by vinee_000 on 6/12/2016.
 */
public class LOOCV {
    public DataSet data;
    public String algo;
    public HashMap<String,Integer> geneRanks;
    public double stability;
    public double [] loss;
    private double [] lambda;
    private double alpha;
    private int respIndex;
    public LOOCV(DataSet data, String algorithm,double alpha,double [] lambda,int responseIndex)//data with rows as samples, columns as variables, name of algorithm, name of Y variable
    {
        this.lambda = lambda;
        this.data = data;
        this.alpha = alpha;
        algo = algorithm;
        respIndex = responseIndex;
    }
    //perform leave one out

    //for each i, set the row with i as the test set, the rest is the set we'll use
    //then we run MGM-PCSTABLE on the training dataset to produce a graph
    //take the first neighbors of response variable as the selected features
    //train an svm on the training dataset to estimate the loss on the testing dataset
    //save loss(i), selectedGenes(i), geneRanks,
    //at the end calculate stability on the selected genes
    public void run()
    {
        geneRanks = new HashMap<String,Integer>();
        int numRows = data.getNumRows();
        String [] selectedGenes = new String[numRows];
        loss = new double[numRows];
        try {
            PrintStream out = new PrintStream("result_" + alpha + "_b.txt");
            for (int i = 0; i < numRows; i++) {
                int[] rows = new int[1];
                rows[0] = i;
                DataSet train = data.copy();
                train.removeRows(rows);
                Graph g;
                if (algo.equals("MGM")) {
                    System.out.println("Initializing MGM");
                    MGM m = new MGM(train, lambda);
                    System.out.println("Running MGM");
                    m.learnEdges(150);
                    g = m.graphFromMGM();
                } else if (algo.equals("MGMPCS")) {
                    MGM m = new MGM(train, lambda);
                    m.learnEdges(90);
                    IndependenceTest ind = new IndTestMultinomialAJ(train, alpha);
                    PcStable p = new PcStable(ind);
                    p.setInitialGraph(m.graphFromMGM());
                    g = p.search();
                } else if (algo.equals("PCS")) {
                    IndependenceTest ind = new IndTestFisherZ(train, alpha);
                    PcStable p = new PcStable(ind);
                    g = p.search();
                } else if (algo.equals("PCS-Mixed")) {
                    IndependenceTest ind = new IndTestMultinomialAJ(train, alpha);
                    PcStable p = new PcStable(ind);
                    g = p.search();
                } else {
                    return;
                }
                List<Node> neighbors = g.getAdjacentNodes(train.getVariable(respIndex));
                if (neighbors.size() == 0) {
                    System.out.println("No neighbors of response");
                }
                for (Node nod : neighbors) {
                    if (geneRanks.get(nod.getName()) != null)
                        geneRanks.put(nod.getName(), geneRanks.get(nod.getName()) + 1);
                    else
                        geneRanks.put(nod.getName(), 1);
                    if (selectedGenes[i] == null)
                        selectedGenes[i] = Integer.toString(data.getColumn(nod));
                    else
                        selectedGenes[i] = selectedGenes[i] + "\t" + data.getColumn(nod);
                }
                System.out.println("Printing to Stream");
                out.println(i + "\t" + selectedGenes[i]);
                out.flush();
                System.out.println("Flushed stream");
                DataSet forSVM = train.subsetColumns(neighbors);

                //train on forSVM, test on just row i of original dataset with subsetted columns

            }
out.close();
            stability = calcStable(selectedGenes);
        }
        catch(Exception e) {
            e.printStackTrace();
            return;
        }
    }
    private double calcStable(String [] genes)//average pairwise intersection over union
    {
        int count = 0;
        double score = 0;
        for(int i = 0; i < genes.length;i++)
        {
            if(genes[i]==null)
                continue;
            String [] curr = genes[i].split("\t");
            for(int j = i+1; j < genes.length;j++)
            {
                if(genes[j]==null)
                    continue;
                String [] curr2 = genes[j].split("\t");

                count = count + 1;
                score = score + (intersect(curr,curr2)/(double)union(curr,curr2));
            }
        }
        return score/count;
    }
    private int intersect(String[] one, String [] two)
    {
        int found = 0;
        if(one.length ==0 || two.length==0)
        {
            return 0;
        }

        for (int i = 0; i < one.length; i++) {
            A:for(int j = 0;j < two.length;j++)
            {
                if(two[j].equals(one[i])) {

                    found++;
                    break A;
                }
            }
        }
        return found;
    }
    private int union(String [] one, String[] two)
    {
        ArrayList<String> temp = new ArrayList<String>();
        if(one.length==0)
            return two.length;
        if(two.length==0)
            return one.length;
        for(int i = 0; i < one.length;i++)
          {
             boolean found = false;
                for(String x : temp)
                {
                    if(x.equals(one[i]))
                    {
                        found = true;
                    }
                }
                if(!found)
                    temp.add(one[i]);
          }

        for(int i = 0; i < two.length;i++)
        {
            boolean found = false;
            for(String x : temp)
            {
                if(x.equals(two[i]))
                {
                    found = true;
                }
            }
            if(!found)
                temp.add(two[i]);
        }
        return temp.size();

    }
}
