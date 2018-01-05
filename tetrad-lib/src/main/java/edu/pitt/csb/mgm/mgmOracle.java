package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;

import java.io.File;

/**
 * Created by vinee_000 on 7/14/2017.
 */
public class mgmOracle {
    public mgmOracle()
    {
    }
    public static void main(String [] args) throws Exception
    {
        DataSet d = MixedUtils.loadDataSet2("C:/Users/vinee_000/Desktop/DataSet_50_100_101_D.txt");
        Graph truth = GraphUtils.loadGraphTxt(new File("C:/Users/vinee_000/Desktop/Graph_50_101_D.txt"));
        double [] lambdaValues = new double[40];
        for(int i = 0; i < lambdaValues.length;i++)
        {
            lambdaValues[i] = .1 + .7*i/40;
        }
      double [] params =   mgmOracle.getOracleParams(lambdaValues,d,truth);
        MGM m = new MGM(d,params);
        m.learnEdges(1000);
        Graph g = m.graphFromMGM();
        System.out.println(mgmOracle.getF1(g,truth,d,"CC"));
        System.out.println(mgmOracle.getF1(g,truth,d,"CD"));
        System.out.println(mgmOracle.getF1(g,truth,d,"DD"));


    }
    private static double getF1(Graph est, Graph truth, DataSet d,String type)
    {
        double tp = 0;
        double fp = 0;
        double fn = 0;
        if(type.equals("CC"))
        {
            for(Edge e:est.getEdges())
            {
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
            return -1;
        }

        double prec = tp/(tp+fp);
        double rec = tp/(tp+fn);
        return (2*prec*rec)/(prec+rec);
    }
    public static double[] getOracleParams(double[]init,DataSet data,Graph trueGraph) {

        double oracleCC = 0;
        double oracleCD = 0;
        double oracleDD = 0;
        double[][] numEdges = new double[init.length][3];
        int iterLimit = 1000;
        double[] bestF1 = new double[3];
        for (int i = 0; i < init.length; i++) //learn an MGM for each initial lambda and save the adjacency matrices
        {
            double[] lambda = {init[i], init[i], init[i]};
            MGM m = new MGM(data, lambda);
            m.learnEdges(iterLimit);
            Graph curr = m.graphFromMGM();
            if (trueGraph != null) {
                double F1CC = getF1(curr, trueGraph, data, "CC");
                if (F1CC > bestF1[0]) {
                    bestF1[0] = F1CC;
                    oracleCC = init[i];
                }
                double F1CD = getF1(curr, trueGraph, data, "CD");
                if (F1CD > bestF1[1]) {
                    bestF1[1] = F1CD;
                    oracleCD = init[i];
                }
                double F1DD = getF1(curr, trueGraph, data, "DD");
                if (F1DD > bestF1[2]) {
                    bestF1[2] = F1DD;
                    oracleDD = init[i];
                }
            }
        }

        double [] temp = {oracleCC,oracleCD,oracleDD};
        return temp;
    }
}
