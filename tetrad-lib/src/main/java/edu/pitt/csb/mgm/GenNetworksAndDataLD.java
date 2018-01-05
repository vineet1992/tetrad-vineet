package edu.pitt.csb.mgm;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataWriter;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.sem.GeneralizedSemIm;
import edu.cmu.tetrad.sem.GeneralizedSemPm;
import edu.cmu.tetrad.util.RandomUtil;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.*;
import java.util.*;

/**
 * Created by ajsedgewick on 8/27/15.
 * java -cp tetrad.jar GenNetworksAndDataLD run_dir
 * where run_dir must have directories named 'data', 'im', and 'networks' in it before
 * you run this
 */
public class GenNetworksAndDataLD {
    //private enum graphTypes {DAG}
    public static void main(String[] args) throws Exception {
        PrintStream out = new PrintStream("target_info.txt");
        out.println("Graph_Number\tTarget_Variable\tNeighborhood");
        int numSamples = 10000;
        int numGraphs = 1;
        int numCats = 3;
        int numVars = 100;
        int maxNumEdges = 200;
        int addedCycleEdges = 20;
        int minCycleLength = 3;
        int maxDegree = 10;
        PrintWriter pw = null;
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
        REP:
        for (int rep = 0; rep < numGraphs; rep++) {
            Graph graph;
            // if(gt == graphTypes.DAG) {
            List<Node> vars  = new ArrayList<>();
            for(int i = 0; i < numVars; i++) {
                vars.add(new GraphNode("X" + (i + 1)));
            }

            graph = GraphUtils.randomDag(vars, 0, maxNumEdges, maxDegree, maxDegree, maxDegree, false);
            //else {if(gt == graphTypes.SF) { //no interconnection
            //graph = ScaleFreeGraphGenerator.randomSF(numVars,maxNumEdges*10,.45,0.1,.45,1.0,1.0, true);
            //}
           // graph = GraphUtils.scaleFreeGraph(numVars,0,.45,0.1,1,1);
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

            A:while(true) {
                for (int j = 0; j < numVars; j++) {
                    Node curr = graph.getNode("X" + (j + 1));
                    if (graph.getAdjacentNodes(curr).size() > 5 && graph.getChildren(curr).size() > 1 && graph.getChildren(curr).size() < 5) {
                        String target = "X" + (j + 1);
                        out.println(rep + "\t" + target + "\t" + graph.getAdjacentNodes(curr));
                        break A;
                    }
                }
                graph = GraphUtils.scaleFreeGraph(numVars,0,.45,0.1,1,1);
            }
            /*******************************************************/
            //Add node with 5 discrete and 5 continuous parents here
            /*****************************************************/

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

            graph.addNode(new GraphNode("Survival"));
            Random rand = new Random();
            int nodesAdd = 0;
            while(nodesAdd < 5) {
                int rd = rand.nextInt(graph.getNumNodes());
                if (nodeDists.get("X" + rd) == 0 && !graph.getAdjacentNodes(graph.getNode("Survival")).contains(graph.getNode("X" + rd))) {
                    graph.addEdge(new Edge(graph.getNode("X" + rd), graph.getNode("Survival"), Endpoint.TAIL, Endpoint.ARROW));
                    nodesAdd++;
                }
            }

            nodesAdd = 0;
            while(nodesAdd < 5) {
                int rd = rand.nextInt(graph.getNumNodes());
                if (nodeDists.get("X" + rd) != 0 && !graph.getAdjacentNodes(graph.getNode("Survival")).contains(graph.getNode("X" + rd))) {
                    graph.addEdge(new Edge(graph.getNode("X" + rd), graph.getNode("Survival"), Endpoint.TAIL, Endpoint.ARROW));
                    nodesAdd++;
                }
            }

            //String paramTemplate = "Split(-.9,-.5,.5,.9)";
            //if(gt==graphTypes.DAG || gt==graphTypes.SF) {
            String paramTemplate = "Split(-1.5,-0.5,0.5,1.5)";
            //}
            //make sem with our dists
            GeneralizedSemPm semPm = MixedUtils.GaussianCategoricalPm(graph, paramTemplate);
            GeneralizedSemIm semIm = MixedUtils.GaussianCategoricalIm(semPm);

            //recursive works well with no cycles
            //DataSet ds = semIm.simulateDataRecursive(numSamples, false);
            DataSet ds = semIm.simulateDataAvoidInfinity(numSamples, false);
            DataSet ds2 = Survival.getSurvival(ds,graph);
            PrintStream out2 = new PrintStream("Survival_" + rep + "_graph.txt");
            out2.println(graph);
            //DataSet ds = semIm.simulateDataAJ(numSamples, false);
            try {
                PrintWriter p = new PrintWriter(new FileWriter("Survival_" + rep + "_data.txt"));
               /* p.println("/variables");
               /* for(Node n: variableNodes){
                    String nName = n.getName();
                    if(nodeDists.get(nName)==numCats){
                        p.println(nName + ":0.0000 1.0000 2.0000");
                    } else {
                        p.println(nName + ":Continuous");
                    }
                }
               for(Node n:ds2.getVariables())
               {
                   String nName = n.getName();
                   if(nodeDists.get(nName)!=null)
                   {
                       if(nodeDists.get(nName)==numCats)
                           p.println(nName + ":0.0000 1.0000 2.0000");
                       else
                           p.println(nName + ":Continuous");
                   }
                   else
                   {
                       if(nName.startsWith("Time"))
                           p.println(nName + ":Continuous");
                       else
                           p.println(nName + ":0.0000 1.0000 2.0000");
                   }
               }
                p.println("/data");
                DataWriter d= new DataWriter();*/
                p.println(ds2);

                //p.println(ds);
                p.flush();
                p.close();
            }
            catch(Exception e) {
                e.printStackTrace();
            }
            List<Node> dsVars = ds.getVariables();
            //int[] discInds = new int[numDisc];
            //int[] contInds = new int[numVars-numDisc];
            //numDisc = 0;
                /*for(int i=0; i < numVars;i++){
                    String nName = dsVars.get(i).getName();
                    if(nodeDists.get(nName).equals("Disc")){
                        discInds[numDisc] = i;
                        numDisc++;
                    } else{
                        contInds[i-numDisc] = i;
                    }
                }
                System.out.println("Mixed?: " + ds.isMixed());
                //DataSet dsDisc = ds.subsetColumns(discInds);
                //DataSet dsCont = ds.subsetColumns(contInds);
*/
             /*   try {
                    //File contFile = new File(dir + "data" ,prefix + "_data.C.txt");
                   // PrintWriter fulldpw = new PrintWriter(new File(dir + "data" ,prefix + "_data.txt"));

                    //DataWriter.writeRectangularData(ds, fulldpw, '\t');
                    //fulldpw.close();
                    //DataWriter.writeRectangularData(dsCont, new PrintWriter(new File(dir + "data" ,prefix + "_data.C.txt")), '\t');
                    //DataWriter.writeRectangularData(dsDisc, new PrintWriter(new File(dir + "data" ,prefix + "_data.D.txt")), '\t');
                   // GraphUtils.saveGraph(graph, new File(dir + "networks" , prefix + "_graph.txt"), false);
                    //pw.println(prefix + "\t" + varsInCycles);

                    //PrintWriter impw = new PrintWriter(new File(dir + "im" ,prefix + "_im.txt"));
                    //impw.println(semIm);
                    //impw.close();

                } catch (Throwable t)
                {
                    t.printStackTrace();
                }
            }
	    }

        pw.close();*/

        }

    }
}