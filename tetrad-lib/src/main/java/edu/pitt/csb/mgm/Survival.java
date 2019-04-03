package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.distribution.WeibullDistribution;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

/**
 * Created by vinee_000 on 2/1/2017.
 */
public class Survival {
    public static DataSet getSurvival(DataSet d, Graph g) {
        Node p = null;
        System.out.println(g);
        /*A:for(Node n:g.getNodes())
        {
            if((g.getChildren(n)==null || g.getChildren(n).size()<1) && g.getParents(n).size() > 1 && n instanceof DiscreteVariable)
            {
                List<Node> n2 = g.getParents(n);
                boolean dv = false;
                boolean cv = false;
                for(Node nn: n2)
                {
                    System.out.println(nn);
                    System.out.println(nn instanceof ContinuousVariable);
                    System.out.println(nn instanceof DiscreteVariable);
                    if(nn instanceof DiscreteVariable)
                        dv = true;
                    else
                        cv = true;
                }
                System.out.println(n2 + "," + dv + "," + cv);
                if(dv && cv) {
                    p = n;
                    break A;
                }
            }
        }*/



        p = g.getNode("Survival");

        int col = d.getColumn(d.getVariable(p.getName()));
        List<Node> temp = g.getParents(p);

        List<Node> parents = new ArrayList<Node>();
        for(Node n: temp)
        {
            parents.add(d.getVariable(n.getName()));
        }
        g.getNode(p.getName()).setName("Time");
        DataSet covariates = d.subsetColumns(parents);
        HashMap<String,Double> params = new HashMap<String,Double>();

        Random rand = new Random();
        for(int  i = 0; i < covariates.getNumColumns();i++)
        {
            String n = covariates.getVariable(i).getName();
            boolean cont = false;
            try {
                covariates.getDouble(0, i);
                cont = true;
            }
            catch(Exception e) {

            }

            if(cont)
            {
                double pos = rand.nextDouble();
                if(pos > 0.5)
                    pos = 1;
                else
                    pos = -1;

                double val = rand.nextDouble()+0.5;
                val = val * pos;
                params.put(n,val);
            }
            else
            {
                for(int j = 0; j < 3;j++)
                {
                    double pos = rand.nextDouble();
                    if(pos > 0.5)
                        pos = 1;
                    else
                        pos = -1;

                    double val = rand.nextDouble()+0.5;
                    val = val * pos;
                    params.put(n + "_" + j,val);
                }
            }

        }
       double beta = 2;
        double beta_0 = 1;
        double [][] survivalVar = new double[d.getNumRows()][2]; //first column is survival, second is censor
        WeibullDistribution w2 = new WeibullDistribution(beta,Math.exp(-1*beta*beta_0));
        for(int i = 0; i < d.getNumRows();i++)
        {
            double linpred = beta_0;
            for(int j = 0; j < covariates.getNumColumns();j++)
            {
                String n = covariates.getVariable(j).getName();
                if(params.get(n)==null)
                {
                    linpred += params.get(n + "_" + covariates.getInt(i,j));
                }
                else
                {
                    linpred+= params.get(n)*covariates.getDouble(i,j);
                }
            }
            System.out.println(linpred);
            double lambda = Math.exp(-1*beta*linpred);
            double time = w2.sample();
            WeibullDistribution w = new WeibullDistribution(beta,lambda);
            //Shape,Scale

          double prob =  w.cumulativeProbability(time);
            int factor = 1000;
            if(rand.nextDouble() < prob) //dead
            {
                survivalVar[i][0] = factor*time;
                survivalVar[i][1] = 0;
            }
            else
            {
                survivalVar[i][0] = factor*time;
                survivalVar[i][1]=1;
            }


        }
        d.removeColumn(col);
        List<String> categories = new ArrayList<String>();
        categories.add("0.0000");
        categories.add("1.0000");
        categories.add("2.0000");
        Node censor = new DiscreteVariable("Censor_Indicator",categories);
        Node survivalTime = new ContinuousVariable("Time");
        TetradMatrix t = new TetradMatrix(survivalVar);
        List<Node> cenList = new ArrayList<Node>();
        for(Node x:d.getVariables())
            cenList.add(x);
        cenList.add(survivalTime);
        cenList.add(censor);
        TetradMatrix data1 = d.getDoubleData();
        TetradMatrix finalOutput = new TetradMatrix(d.getNumRows(),d.getNumColumns()+2);
        for(int i = 0; i < d.getNumRows();i++)
        {
            for(int j = 0; j < d.getNumColumns();j++)
            {
                finalOutput.set(i,j,data1.get(i,j));
            }
        }
        for(int i = 0; i < d.getNumRows();i++)
        {
            for(int j = d.getNumColumns();j < d.getNumColumns()+2;j++)
            {
                finalOutput.set(i,j,t.get(i,j-d.getNumColumns()));
            }
        }


        //DataSet d = MixedUtils.makeMixedData()
        //col is the survival var, covariates has the covariates
        DataSet d3 = ColtDataSet.makeData(cenList,finalOutput);
        //change p's name to survival and break
        return d3;
    }
}
