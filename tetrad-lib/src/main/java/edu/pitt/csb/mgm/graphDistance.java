package edu.pitt.csb.mgm;

import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by vinee_000 on 1/23/2017.
 */
public class graphDistance
{
    //


    public static void intensity(String filename,Graph g1,String target) throws Exception
    {
        int maxDistance = 6;
        BufferedReader b = new BufferedReader(new FileReader("intensity_values.txt"));
        HashMap<String,Double> intensity = new HashMap<String,Double>();
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            intensity.put(line[0],Double.parseDouble(line[1]));
        }
        b.close();
        PrintStream out = new PrintStream(filename);
        out.println("Distance\tIntensity");
            for (Node n : g1.getNodes()) {
                if(intensity.get(n.getName())==null)
                    continue;
                if(n.getName().equals(target))
                    continue;
                List<List<Node>>  temp = GraphUtils.treks(g1,n,g1.getNode(target),maxDistance);
                if(temp==null)
                    out.println(maxDistance + "\t" + intensity.get(n.getName()));
                else
                {
                    double min = maxDistance;
                    for(List<Node> list: temp)
                    {
                        if(min > list.size()-1)
                            min = list.size()-1;
                    }
                    out.println(min + "\t" + intensity.get(n.getName()));
                }
            }
        out.flush();
        out.close();
    }
    public static double[] blanketStable(Graph g1, Graph g2, Graph g3, String target) throws Exception
    {

        Set<Node> graph1 = new HashSet<Node>(g1.getNodes());
        Set<Node> graph2 = new HashSet<Node>(g2.getNodes());
        Set<Node> graph3 = new HashSet<Node>(g3.getNodes());


        List<Node> n = g1.getAdjacentNodes(g1.getNode(target));
        Set<Node> allNodes = new HashSet<Node>(n);
        List<Node> temp = g1.getChildren(g1.getNode(target));
        for(Node k: temp)
        {
            List<Node> temp2 = g1.getParents(k);
            allNodes.addAll(temp2);
        }

        n = g2.getAdjacentNodes(g2.getNode(target));
        Set<Node> allNodes2 = new HashSet<Node>(n);
        temp = g2.getChildren(g2.getNode(target));
        for(Node k: temp)
        {
            List<Node> temp2 = g2.getParents(k);
            allNodes2.addAll(temp2);
        }

        n = g3.getAdjacentNodes(g3.getNode(target));
        Set<Node> allNodes3 = new HashSet<Node>(n);
        temp = g3.getChildren(g3.getNode(target));
        for(Node k: temp)
        {
            List<Node> temp2 = g3.getParents(k);
            allNodes3.addAll(temp2);
        }

        int count = 0;
        System.out.println("K=100 Graph:");
        for(Node nn:allNodes)
        {
            if(!nn.getName().toUpperCase().equals(nn.getName()))
                continue;
            for(Node k:allNodes2)
            {
                if(k.getName().equals(nn.getName()))
                {
                    System.out.println(nn.getName());
                    count++;
                }
            }
        }
        double [] results = new double[4];
        results[0] = count;
        System.out.println("K = 500 Graph:");
        count = 0;
        for(Node nn:allNodes)
        {
            if(!nn.getName().toUpperCase().equals(nn.getName()))
                continue;
            for(Node k:allNodes3)
            {
                if(k.getName().equals(nn.getName()))
                {
                    System.out.println(nn.getName());
                    count++;
                }
            }
        }

results[1] = count;

        count = 0;

        for(Node nn: graph1)
        {
            if(!nn.getName().toUpperCase().equals(nn.getName()))
                continue;
            A:for(Node nn2:graph2)
            {

                if(nn.getName().equals(nn2.getName())) {

                    count++;
                    break A;
                }
            }
        }
        results[2] = count;
        count = 0;

        for(Node nn: graph1)
        {
            if(!nn.getName().toUpperCase().equals(nn.getName()))
                continue;
           A: for(Node nn2:graph3)
            {
                if(nn.getName().equals(nn2.getName())) {
                    count++;
                    break A;
                }
            }
        }
        results[3] = count;

        return results;

    }
    public static void intensityDegree(String filename,Graph g1,String target) throws Exception
    {
        BufferedReader b = new BufferedReader(new FileReader("intensity_values.txt"));
        HashMap<String,Double> intensity = new HashMap<String,Double>();
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            intensity.put(line[0],Double.parseDouble(line[1]));
        }
        b.close();
        PrintStream out = new PrintStream(filename);
        out.println("Degree\tIntensity");
        for (Node n : g1.getNodes()) {
            if(intensity.get(n.getName())==null)
                continue;
           out.println(g1.getAdjacentNodes(n).size()+ "\t" + intensity.get(n.getName()));
        }
        out.flush();
        out.close();
    }
    public static double distance(Graph g1, Graph g2) throws Exception{
        PrintStream out = new PrintStream("temp1.txt");
        out.println(g1);
        out.flush();
        out.close();
        out = new PrintStream("temp2.txt");
        out.println(g2);
        out.flush();
        out.close();
            int distance = 0;
        Graph g3 = GraphUtils.loadGraphTxt(new File("temp1.txt"));
        Graph g4 = GraphUtils.loadGraphTxt(new File("temp2.txt"));
int nodeAdditions = 0;
        int nodeDeletions = 0;
        for(Node n:g3.getNodes())
        {
            if(g4.getNode(n.getName())==null)
            {
                nodeAdditions++;
                        distance++;
            }
        }
        List<Node> temp = new ArrayList<Node>();
        for(Node n: g4.getNodes())
        {
            if(g3.getNode(n.getName())==null)
            {
                temp.add(n);
                nodeDeletions++;
                distance++;
            }
        }

        HashMap<String,Integer> used = new HashMap<String,Integer>();
        int edgeChanges = 0;
        int edgeAdditions = 0;
        int edgeDeletions = 0;
        A:for(Edge e: g3.getEdges())
        {
            if(used.get(e.getNode1().getName()+","+e.getNode2().getName())==null)
            {
                used.put(e.getNode1().getName()+","+e.getNode2().getName(),2);
                used.put(e.getNode2().getName()+","+e.getNode1().getName(),2);
            }
            else
                continue A;

            if(g4.getEdge(g4.getNode(e.getNode1().getName()),g4.getNode(e.getNode2().getName()))==null)
            {
                edgeAdditions++;
                distance++;
            }
            else
            {
                Edge ed = g4.getEdge(g4.getNode(e.getNode1().getName()),g4.getNode(e.getNode2().getName()));
                if(ed.getNode1().getName().equals(e.getNode1().getName()))
                {
                    if(ed.getProximalEndpoint(ed.getNode1())!=e.getProximalEndpoint(e.getNode1()) || ed.getProximalEndpoint(ed.getNode2())!= e.getProximalEndpoint(e.getNode2()))
                    {
                        edgeChanges++;
                        distance++;
                    }
                }
                else
                {
                    if(ed.getProximalEndpoint(ed.getNode1())!=e.getProximalEndpoint(e.getNode2()) || ed.getDistalEndpoint(ed.getNode1())!=e.getDistalEndpoint(e.getNode2()))
                    {
                        edgeChanges++;
                        distance++;
                    }
                }
            }
        }
        used.clear();
        for(Edge e:g4.getEdges())
        {
            if(g3.getEdge(g3.getNode(e.getNode1().getName()),g3.getNode(e.getNode2().getName()))==null)
            {
                edgeDeletions++;
                distance++;
            }
        }
        System.out.println("Nodes_Added\tNodes_Deleted\tEdges_Added\tEdges_Deleted\tEdges_Changed");
        System.out.println(nodeAdditions + "\t" + nodeDeletions + "\t" + edgeAdditions + "\t" + edgeDeletions + "\t" + edgeChanges);
        File f = new File("temp1.txt");
        f.delete();
        f = new File("temp2.txt");
        f.delete();
        return distance;

    }

}
