package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.Knowledge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fci;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradVector;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * Created by vinee_000 on 5/16/2016.
 */
public class testStatin {
    public static void main(String[]args) throws Exception {

        PrintStream out3 = new PrintStream("second_neighbors_subtype.txt");
        Set<Node> s3 = new TreeSet<Node>();
        Graph g3 = GraphUtils.loadGraphTxt(new File("mfm_steps_classifying_subtype.txt"));
        for(Node n: g3.getParents(g3.getNode("Subtype")))
        {
            s3.add(n);
            for(Node j: g3.getParents(n))
            {
                s3.add(j);
            }
            for(Node x: g3.getChildren(n))
            {
                s3.add(x);
            }
        }
        for(Node n: g3.getChildren(g3.getNode("Subtype")))
        {
            s3.add(n);
            for(Node j: g3.getParents(n))
                s3.add(j);
            for(Node j: g3.getChildren(n))
                s3.add(j);
        }
        for(Node n: s3)
            out3.println(n.getName());
        out3.flush();
        out3.close();
        System.exit(0);

        DataSet data = MixedUtils.loadDataSet2("mfm_expression_dataset.txt");

        data.removeColumn(data.getVariable("Statin_Sensitivity"));
        data.removeColumn(data.getVariable("HER2_Status"));
        data.removeColumn(data.getVariable("ID"));
        data.removeColumn(data.getVariable("ER_Status"));
        data.removeColumn(data.getVariable("PR_Status"));
        data.removeColumn(data.getVariable("Gender"));
        data.removeColumn(data.getVariable("Metastasis"));
        System.out.println(data);
        System.out.println(data.getVariableNames());
        double [] lambda = new double[7];
        for(int i = 0; i < lambda.length;i++)
        {
            lambda[i] = .1*i + .06;
        }
        STEPS s = new STEPS(data,lambda,.01,3);
        Graph g = s.runSteps();
        FciMaxP f = new FciMaxP(new IndTestMultinomialAJ(data,.05));
        PrintStream out = new PrintStream ("mfm_steps_classifying_subtype.txt");
        f.setInitialGraph(g);
        out.println(f.search());
        System.exit(0);
       /* System.out.println(data.subsetColumns(rows));
        double[] lambda = {.55, .55, .55};
        int numSubsamples = 5;
        int subSampleSize = 420;
        for (int i = 0; i < numSubsamples; i ++) {
            int [] rows = new int[subSampleSize];
            for(int j = 0; j< subSampleSize ; j++)
            {
                rows[j] = j;
            }
          //  data.permuteRows();
        double [] stepsLambda = {.64,.55,.45,.35};
        STEPS s = new STEPS(data,stepsLambda,.05,5);
        Graph g = s.runSteps();
        IndependenceTest ii = new IndTestMultinomialAJ(data, .01);
        FciMaxP f = new FciMaxP(ii);
        PrintStream out = new PrintStream("steps.txt");
        f.setInitialGraph(g);
        out.println(f.search());
        out.flush();
        out.close();
            DataSet curr = data.copy();
           // curr = curr.subsetRows(rows);
            boolean again = false;
            while(!again) {
                try {
                    MGM m = new MGM(curr, lambda);
                    System.out.println("Initialized MGM");
                    //m.learnEdges(1000);
                    PrintStream out = new PrintStream("fci.txt");
                    IndependenceTest ii = new IndTestMultinomialAJ(curr, .001);
                    Knowledge k = new Knowledge();

                    k.addToTier(2,"Major_Diagnosis");
                    for(String x: curr.getVariableNames())
                    {
                        if(!x.equals("Major_Diagnosis"))
                            k.addToTier(1,x);
                    }

                    Fci f = new Fci(ii);
                    //f.setInitialGraph(m.graphFromMGM());
                   // f.setKnowledge(k);
                    Graph se = f.search();
                    System.out.println(se);
                    out.println(se);
                    out.flush();
                    out.close();

                    again = true;
                } catch (Exception e) {
                    e.printStackTrace();
                    data.permuteRows();
                    curr = data.copy();
                    curr = curr.subsetRows(rows);
                }
            }*/
    }
}
