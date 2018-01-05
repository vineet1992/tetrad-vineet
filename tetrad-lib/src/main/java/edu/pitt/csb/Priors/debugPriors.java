package edu.pitt.csb.Priors;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;

/**
 * Created by vinee_000 on 9/16/2017.
 */
public class debugPriors {
    public static double [][] numPriors;
    public static double nr;

  /*  public static void main(String [] args) throws Exception
    {
        int numVariables = 50;
        int numRuns = 20;
        Parameters p = new Parameters();
        int sampleSize = 500;
        TetradMatrix [] priors = new TetradMatrix[1];
        MixedLeeHastieSimulation c = new MixedLeeHastieSimulation();
        c.simulate(p);
        double ap = 0.5;
        priors = simulatePrior(c.getTrueGraph(),ap,1,priors);
        boolean [][] havePriors = new boolean[priors[0].rows()][priors[0].columns()];
        for(int i = 0; i < priors[0].rows();i++)
        {
            for(int j = 0; j < priors[0].columns();j++)
            {
                if(priors[0].get(i,j)!=0)
                    havePriors[i][j] = true;
            }
        }
        double [] lambda = {0.2,0.2,0.2};
        double [] lambdaWP = {0.9,0.2,0.9};
        double [] lambdaWP2 = {0.1,0.5,0.1};

        MGM_Priors m = new MGM_Priors(c.getDataSet(0),lambda,lambdaWP,havePriors);
        m.learnEdges(1000);
        Graph dense = m.graphFromMGM();
        m = new MGM_Priors(c.getDataSet(0),lambda,lambdaWP2,havePriors);
        m.learnEdges(1000);
        Graph sparse = m.graphFromMGM();
        double [] denseEd =  getPriorEdges(dense,c.getDataSet(0),havePriors);
        double [] sparseEd = getPriorEdges(sparse,c.getDataSet(0),havePriors);
        System.out.println(Arrays.toString(denseEd));
        System.out.println(Arrays.toString(sparseEd));
    }*/
    public static double[]getPriorEdges(Graph g, DataSet d, boolean [][] priors)
    {
        double [] temp = new double[4];
        int [] div = new int[4];
        for(int i = 0; i < d.getNumColumns();i++)
        {
            for(int j = i+1;j < d.getNumColumns();j++)
            {
                if(priors[i][j])
                {
                    temp[0]++;
                    Node datOne = d.getVariable(i);
                    Node datTwo = d.getVariable(j);
                    if(datTwo instanceof ContinuousVariable && datOne instanceof ContinuousVariable)
                    {
                        if(g.getEdge(g.getNode(datOne.getName()),g.getNode(datTwo.getName()))!=null)
                        {
                            temp[1]++;
                        }
                        div[1]++;

                    }
                    else if (datTwo instanceof DiscreteVariable && datOne instanceof DiscreteVariable)
                    {
                        if(g.getEdge(g.getNode(datOne.getName()),g.getNode(datTwo.getName()))!=null)
                        {
                            temp[3]++;
                        }
                        div[3]++;
                    }
                    else
                    {
                        if(g.getEdge(g.getNode(datOne.getName()),g.getNode(datTwo.getName()))!=null)
                        {
                            temp[2]++;
                        }
                        div[2]++;
                    }
                }
            }
        }
        for(int i = 1; i < 4;i++)
            temp[i] = temp[i]/div[i];
        return temp;
    }
    public static void main(String [] args)throws Exception
    {
        int numVariables = 50;
        int numRuns = 20;
        Parameters p = new Parameters();
        int numExperts = 5;
        int sampleSize = 500;
        double amountPrior = 0.6;
        PrintStream out = new PrintStream("full_data_" + amountPrior + ".txt");
        out.println("CC_F1\tCD_F1\tDD_F1\tALL_F1\tCC_F1_NP\tCD_F1_NP\tDD_F1_NP\tALL_F1_NP\tCC_F1_WP\tCD_F1_WP\tDD_F1_WP\tALL_F1_WP\tPercent_Prior_CC\tPercent_Prior_CD\tPercent_Prior_DD\tPrior_Overlap_CC\tPrior_Overlap_CD\tPrior_Overlap_DD\tNum_Edges_CC_NP\tNum_Edges_CD_NP\tNum_Edges_DD_NP\tEstimated_Edges_CC_NP\tEstimated_Edges_CD_NP\tEstimated_Edges_DD_NP\tNum_Edges_CC_WP\tNum_Edges_CD_WP\tNum_Edges_DD_WP\tEstimated_Edges_CC_WP\tEstimated_Edges_CD_WP\tEstimated_Edges_DD_WP");
        double [][] results = new double[numRuns][4];
        TetradMatrix [] priors = new TetradMatrix[numExperts];
        MixedLeeHastieSimulation c = new MixedLeeHastieSimulation();
        c.simulate(p);
        BufferedReader b = new BufferedReader(new FileReader("Results/mgm_priors" + "_" + amountPrior + "_" + numExperts + "_" + numVariables + "_" + sampleSize + "_10.txt"));
        b.readLine();
        for(int i = 0; i < numRuns ;i++) {
            System.out.println(i);
            nr = i;
            boolean foundFile = false;
            boolean foundData = false;
            //TODO load the graph and data, and save results into the results array
            File f = new File("Graphs/Graph_" + i + "_" + numVariables + ".txt");
            if (f.exists()) {
                foundFile = true;
                c.setTrueGraph(GraphUtils.loadGraphTxt(f));
            }
            f = new File("Data/Data_" + i + "_" + numVariables + "_" + sampleSize + ".txt");
            if (f.exists() && foundFile) {
                foundData = true;
                c.setDataSet(MixedUtils.loadDataSet2("Data/Data_" + i + "_" + numVariables + "_" + sampleSize + ".txt"), 0);
            } else
                c.simulate(p);
            if (foundFile && foundData) {
                for (int j = 0; j < numExperts; j++) {
                    f = new File("Priors/Priors_" + i + "_" + numVariables + "_" + amountPrior + "_" + j + ".txt");
                    if (f.exists()) {
                        priors[j] = new TetradMatrix(priorTest.loadPrior(f, numVariables));
                    }
                }
            }
                String [] line = b.readLine().split("\t");
                for(int j = 0; j < line.length-1;j++)
                {
                    out.print(Double.parseDouble(line[j+1]) + "\t");

                }
            Graph est = GraphUtils.loadGraphTxt(new File("Estimated/" + "mgm_priors" + "_" + i + "_" + numVariables + "_" + sampleSize + "_" + amountPrior + "_" + numExperts +  ".txt"));

            boolean[][] havePrior = findPrior(priors,c.getDataSet(0));
            for(int j = 0; j < 8;j++)
            {
                String type = "";
                if(j%4==0)
                    type = "CC";
                else if(j%4==1)
                    type="CD";
                else if(j%4==2)
                    type="DD";
                else
                    type="ALL";
                boolean WP = true;
                if(j/4==0)
                    WP = false;


                out.print(getF1(est,c.getTrueGraph(),c.getDataSet(0),type,WP,havePrior) + "\t");
            }
                double [][] result = getOverlap(c.getDataSet(0),c.getTrueGraph());
                //i = type of result
            //j = cc,cd,dd
            for(int k = 1; k >= 0;k--) {
                for (int m = 0; m < 3; m++)
                    out.print(result[k][m] + "\t");
            }

            double [][][] ne =  getNumEdges(est,c.getTrueGraph(),c.getDataSet(0),havePrior);
            for(int x = 0; x < 2;x++) {
                for (int k = 1; k >= 0; k--) {
                    for (int m = 0; m < 3; m++)
                        out.print(ne[x][k][m] + "\t");
                }
            }

            out.println();

        }
        b.close();
        out.flush();
        out.close();

        /*System.out.println(priors[0].toString());
        System.out.println(m.getTrueGraph());
        mgmPriors m2 = new mgmPriors(ns,initLambdas,m.getDataSet(0),priors);
        Graph g = m2.runPriors();
       // System.out.println("Oracle:" + m2.oracleAll);
      //  System.out.println("Oracle three lambdas: " + m2.oracleCC + "," + m2.oracleCD + "," + m2.oracleCD);
        System.out.println("DD:" + mgmPriors.getF1(g,m.getTrueGraph(),m.getDataSet(0),"DD"));
        System.out.println("CD:" + mgmPriors.getF1(g,m.getTrueGraph(),m.getDataSet(0),"CD"));
        System.out.println("Estimated: "  + g);
        System.out.println("Truth: " + m.getTrueGraph());
        System.out.println(m.getDataSet(0));*/
    }

    public static double [][][] getNumEdges(Graph est, Graph truth, DataSet data,boolean[][]havePrior)
    {
        double [][][] result = new double[2][2][3];
        for(int i = 0; i < data.getNumColumns();i++)
        {
            for(int j = i +1; j < data.getNumColumns();j++)
            {
                if(est.getEdge(est.getNode(data.getVariable(i).getName()),est.getNode(data.getVariable(j).getName()))!=null)
                {
                    if(data.getVariable(i)instanceof ContinuousVariable && data.getVariable(j) instanceof ContinuousVariable)
                    {
                        if(havePrior[i][j])
                        result[1][0][0]++;
                        else
                            result[0][0][0]++;
                    }
                    else if(data.getVariable(i) instanceof DiscreteVariable && data.getVariable(j) instanceof DiscreteVariable)
                    {
                        if(havePrior[i][j])
                            result[1][0][2]++;
                        else
                            result[0][0][2]++;
                    }
                    else
                    {
                        if(havePrior[i][j])
                            result[1][0][1]++;
                        else
                            result[0][0][1]++;
                    }
                }
                if(truth.getEdge(truth.getNode(data.getVariable(i).getName()),truth.getNode(data.getVariable(j).getName()))!=null)
                {
                    if(data.getVariable(i)instanceof ContinuousVariable && data.getVariable(j) instanceof ContinuousVariable)
                    {
                        if(havePrior[i][j])
                            result[1][1][0]++;
                        else
                            result[0][1][0]++;
                    }
                    else if(data.getVariable(i) instanceof DiscreteVariable && data.getVariable(j) instanceof DiscreteVariable)
                    {
                        if(havePrior[i][j])
                            result[1][1][2]++;
                        else
                            result[0][1][2]++;
                    }
                    else
                    {
                        if(havePrior[i][j])
                            result[1][1][1]++;
                        else
                            result[0][1][1]++;
                    }
                }

            }
        }
        return result;
    }
    public static double [][] getOverlap(DataSet d,Graph truth)
    {
        System.out.println(truth);
        int [][] num = new int[2][3];
        double [][] result = new double[2][3];
        for(int i = 0; i < numPriors.length;i++)
        {
            for(int j = i+1; j < numPriors[i].length;j++)
            {
                if(numPriors[i][j]!=0 && truth.getEdge(truth.getNode(d.getVariable(i).getName()),truth.getNode(d.getVariable(j).getName()))==null)
                {
                    System.out.println("Prior Information for incorrect Edge: " + d.getVariable(i) + "," + d.getVariable(j) + "," + numPriors[i][j]);
                }
                if(d.getVariable(i)instanceof ContinuousVariable &&d.getVariable(j) instanceof ContinuousVariable)
                {
                    result[0][0]+=numPriors[i][j];
                    if(numPriors[i][j]!=0) {
                        result[1][0]++;
                        num[0][0]++;
                    }
                    if(truth.getEdge(truth.getNode(d.getVariable(i).getName()),truth.getNode(d.getVariable(j).getName()))!=null)
                    num[1][0]++;
                }
                else if(d.getVariable(i)instanceof DiscreteVariable && d.getVariable(j)instanceof DiscreteVariable) {
                    result[0][2] += numPriors[i][j];
                    if (numPriors[i][j] != 0){
                        result[1][2]++;
                        num[0][2]++;
                    }
                    if(truth.getEdge(truth.getNode(d.getVariable(i).getName()),truth.getNode(d.getVariable(j).getName()))!=null)

                        num[1][2]++;
                }
                else
                {
                    result[0][1]+=numPriors[i][j];
                    if(numPriors[i][j]!=0) {
                        result[1][1]++;
                        num[0][1]++;
                    }
                    if(truth.getEdge(truth.getNode(d.getVariable(i).getName()),truth.getNode(d.getVariable(j).getName()))!=null)

                        num[1][1]++;
                }
            }
        }

        for(int i = 0; i < 2;i++) {
            for (int j = 0; j < 3; j++) {
                result[i][j] = result[i][j] / num[i][j];
            }
        }

        return result;
    }
    public static boolean[][] findPrior(TetradMatrix [] pInf,DataSet data) throws Exception
    {

        PrintStream temp2 = new PrintStream("temp.txt");
        boolean [][] temp = new boolean[pInf[0].rows()][pInf[0].columns()];
        numPriors = new double[pInf[0].rows()][pInf[0].columns()];
        for(int i = 0; i < pInf.length;i++)
        {

            TetradMatrix curr = pInf[i];

            for(int j = 0; j < data.getNumColumns();j++)
            {
                for(int k = j +1; k < data.getNumColumns();k++)
                {
                    if(curr.get(j,k)!=0) {
                        temp2.println(data.getVariable(j) + "\t" + data.getVariable(k) + "\t" + i + "\t" + nr);
                        numPriors[j][k]++;
                        temp[j][k] = true;
                        temp2.flush();
                    }
                }
            }
        }
        temp2.close();
        return temp;
    }
    public static double getF1(Graph est, Graph truth, DataSet d,String type,boolean WP,boolean[][]havePriors)
    {
        double tp = 0;
        double fp = 0;
        double fn = 0;
        if(type.equals("CC"))
        {
            for(Edge e:est.getEdges())
            {
                int i = d.getColumn(d.getVariable(e.getNode1().getName()));
                int j = d.getColumn(d.getVariable(e.getNode2().getName()));
                if (i > j)
                {
                    int x = i;
                    i = j;
                    j = x;
                }
                if(havePriors[i][j]!=WP)
                    continue;
                if(d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable)
                {
                    Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()),truth.getNode(e.getNode2().getName()));
                    if(temp!=null)
                    {
                        tp++;
                    }
                    else
                        fp++;
                }
            }
            for(Edge e:truth.getEdges())
            {
                int i = d.getColumn(d.getVariable(e.getNode1().getName()));
                int j = d.getColumn(d.getVariable(e.getNode2().getName()));
                if (i > j)
                {
                    int x = i;
                    i = j;
                    j = x;
                }
                if(havePriors[i][j]!=WP)
                    continue;
                if(d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable)
                {
                    Edge temp = est.getEdge(est.getNode(e.getNode1().getName()),est.getNode(e.getNode2().getName()));
                    if(temp==null)
                        fn++;
                }
            }
        }
        else if(type.equals("CD"))
        {
            for(Edge e:est.getEdges())
            {
                int i = d.getColumn(d.getVariable(e.getNode1().getName()));
                int j = d.getColumn(d.getVariable(e.getNode2().getName()));
                if (i > j)
                {
                    int x = i;
                    i = j;
                    j = x;
                }
                if(havePriors[i][j]!=WP)
                    continue;
                if((d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable) ||  (d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable))
                {
                    Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()),truth.getNode(e.getNode2().getName()));
                    if(temp!=null)
                    {
                        tp++;
                    }
                    else
                        fp++;
                }
            }
            for(Edge e:truth.getEdges())
            {
                int i = d.getColumn(d.getVariable(e.getNode1().getName()));
                int j = d.getColumn(d.getVariable(e.getNode2().getName()));
                if (i > j)
                {
                    int x = i;
                    i = j;
                    j = x;
                }
                if(havePriors[i][j]!=WP)
                    continue;
                if((d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable) ||  (d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable))
                {
                    Edge temp = est.getEdge(est.getNode(e.getNode1().getName()),est.getNode(e.getNode2().getName()));
                    if(temp==null)
                        fn++;
                }
            }
        }
        else if(type.equals("DD"))
        {
            for(Edge e:est.getEdges())
            {
                int i = d.getColumn(d.getVariable(e.getNode1().getName()));
                int j = d.getColumn(d.getVariable(e.getNode2().getName()));
                if (i > j)
                {
                    int x = i;
                    i = j;
                    j = x;
                }
                if(havePriors[i][j]!=WP)
                    continue;
                if(d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable)
                {
                    Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()),truth.getNode(e.getNode2().getName()));

                    if(temp!=null)
                    {
                        tp++;
                    }
                    else
                        fp++;
                }
            }
            for(Edge e:truth.getEdges())
            {
                int i = d.getColumn(d.getVariable(e.getNode1().getName()));
                int j = d.getColumn(d.getVariable(e.getNode2().getName()));
                if (i > j)
                {
                    int x = i;
                    i = j;
                    j = x;
                }
                if(havePriors[i][j]!=WP)
                    continue;
                if(d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable)
                {
                    Edge temp = est.getEdge(est.getNode(e.getNode1().getName()),est.getNode(e.getNode2().getName()));
                    if(temp==null)
                        fn++;
                }
            }
        }
        else
        {
            for(Edge e:est.getEdges())
            {
                int i = d.getColumn(d.getVariable(e.getNode1().getName()));
                int j = d.getColumn(d.getVariable(e.getNode2().getName()));
                if (i > j)
                {
                    int x = i;
                    i = j;
                    j = x;
                }
                if(havePriors[i][j]!=WP)
                    continue;
                Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()),truth.getNode(e.getNode2().getName()));
                if(temp!=null)
                {
                    tp++;
                }
                else
                    fp++;
            }
            for(Edge e:truth.getEdges())
            {
                int i = d.getColumn(d.getVariable(e.getNode1().getName()));
                int j = d.getColumn(d.getVariable(e.getNode2().getName()));
                if (i > j)
                {
                    int x = i;
                    i = j;
                    j = x;
                }
                if(havePriors[i][j]!=WP)
                    continue;
                Edge temp = est.getEdge(est.getNode(e.getNode1().getName()),est.getNode(e.getNode2().getName()));
                if(temp==null)
                    fn++;
            }
        }

        double prec = tp/(tp+fp);
        double rec = tp/(tp+fn);
        return (2*prec*rec)/(prec+rec);
    }
}
