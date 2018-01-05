package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataWriter;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphNode;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.sem.GeneralizedSemIm;
import edu.cmu.tetrad.sem.GeneralizedSemPm;
import edu.cmu.tetrad.util.RandomUtil;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

/**
 * Created by vinee_000 on 5/16/2017.
 */

/* This class generates networks and data for the FCI-MAX experiments*/

public class GenNetworksAndData {
    public static void main(String [] args)
    {

        int numSamples = 1000;
        int numGraphs = 7;
        int numCats = 3;
        int numVars = 50;
        NormalDistribution n2 = new NormalDistribution(750,200);
        int maxNumEdges = 1100;
        int addedCycleEdges = 20;
        int minCycleLength = 3;
        int maxDegree = 20;
        PrintWriter pw = null;
        String dir = "Graphs_50_Nodes";
        String dir2 = "Data_50_Nodes";
//        String dir = System.getProperty("user.dir") + "/" + args[0];
        //      System.out.println("Saving to " + dir);
        //String dir = "/Users/ajsedgewick/tetrad_mgm_runs/run1/";
        //String dir = "/Users/ajsedge/Dropbox/Research_Main/tetrad_stuff/run3/";
       /* try {
            FileOutputStream fos = new FileOutputStream(dir + "netsummary.txt");
            pw = new PrintWriter(fos);
        } catch (Throwable t) {
            t.printStackTrace();
        }*/

        //WARNING: this isn't enough to make output repeatable, need to figure that out
        RandomUtil.getInstance().setSeed(12345);


        //GRAPH:
        //for(graphTypes gt: graphTypes.values()) {
        String gt = "SF";
        Random rand = new Random();
        REP:
        for (int rep = 6; rep < numGraphs; rep++) {
            Graph graph;
            // if(gt == graphTypes.DAG) {
            List<Node> vars  = new ArrayList<>();
            for(int i = 0; i < numVars; i++) {
                vars.add(new GraphNode("X" + (i + 1)));
            }
            graph = GraphUtils.loadGraphTxt(new File(dir + "/DAG_" + rep + ".txt"));
           // graph = GraphUtils.randomDag(vars, 0, (int)n2.sample(), maxDegree, maxDegree, maxDegree, false);
            //else {if(gt == graphTypes.SF) { //no interconnection
            //graph = ScaleFreeGraphGenerator.randomSF(numVars,maxNumEdges*10,.45,0.1,.45,1.0,1.0, true);
            //}
            //double beta = rand.nextDouble();
             //graph = GraphUtils.scaleFreeGraph(numVars,0,.3,beta*.6+.1,1,1);
            //else if(gt == graphTypes.CYC) { //no interconnection
            //Dag dag = GraphUtils.randomDag(numVars, 0, maxNumEdges, maxDegree, maxDegree, maxDegree, false);
            //graph = GraphUtils.addCycles(dag, addedCycleEdges, maxDegree);
            //  graph = GraphUtils.cyclicGraph4(numVars, maxNumEdges);
            //} else {//SF2 more interconnection
            //    graph = ScaleFreeGraphGenerator.randomSF(numVars,maxNumEdges*10,.25,.5,.25,1.0,1.0, false);
            //}

            /*String prefix = "SF" + "_" + rep;
            File[] outfiles = new File(dir + "networks").listFiles();
            for (File _file : outfiles) {
                if (_file.getName().startsWith(prefix)) {
                    continue REP;
                }
            }*/


            //scale free nets have hubs in earlier variables, so randomize...

            //int numVars = variableNodes.size();

            int varsInCycles = 0;
            int numDisc = 0;
            // Graph graph = GraphUtils.loadGraphTxt(new File("SF_" + rep + "_graph.txt"));

            List<Node> variableNodes = graph.getNodes();
            //Collections.shuffle(variableNodes);
            HashMap<String, Integer> nodeDists = new HashMap<>();
            for (int i = 0; i < numVars; i++) {
                Node node = variableNodes.get(i);
                //if(GraphUtils.directedPathFromTo(graph, node, node) != null) {
                //    varsInCycles++;
                //}
                String nName = node.getName();
                int ind = Integer.parseInt(nName.substring(1,nName.length()));
                if (ind <= numVars / 2.0) {
                    nodeDists.put(nName, 0);
                } else {
                    nodeDists.put(nName, numCats);
                    numDisc++;
                }
            }
            graph = MixedUtils.makeMixedGraph(graph, nodeDists);
           // System.out.println(graph.getNumEdges());
            //GraphUtils.saveGraph(graph,new File("DAG_" + rep + ".txt"),false);
            //String paramTemplate = "Split(-.9,-.5,.5,.9)";
            //if(gt==graphTypes.DAG || gt==graphTypes.SF) {
            String paramTemplate = "Split(-2,-0.2,0.2,2)";
            //}
            //make sem with our dists
            GeneralizedSemPm semPm = MixedUtils.GaussianCategoricalPm(graph, paramTemplate);
            GeneralizedSemIm semIm = MixedUtils.GaussianCategoricalIm(semPm);

            //recursive works well with no cycles
            //DataSet ds = semIm.simulateDataRecursive(numSamples, false);
            DataSet ds = semIm.simulateDataAvoidInfinity(numSamples, false);
            //DataSet ds = semIm.simulateDataAJ(numSamples, false);
            try {
                PrintWriter p = new PrintWriter(new FileWriter(dir2 + "/DAG_" + rep + "_data_hard.txt"));
                p.println("/variables");
                for(Node n: variableNodes){
                    String nName = n.getName();
                    if(nodeDists.get(nName)==numCats){
                        p.println(nName + ":0.0000 1.0000 2.0000");
                    } else {
                        p.println(nName + ":Continuous");
                    }
                }
                p.println("/data");
                p.println(ds);

                p.flush();
                p.close();
            }
            catch(Exception e) {
                e.printStackTrace();
            }
            List<Node> dsVars = ds.getVariables();
                    //pw.println(prefix + "\t" + varsInCycles);

                    //PrintWriter impw = new PrintWriter(new File(dir + "im" ,prefix + "_im.txt"));
                    //impw.println(semIm);
                    //impw.close();

        }

    }
}
