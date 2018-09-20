package edu.pitt.csb.latents;

import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.DagToPag;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by vinee_000 on 10/25/2017.
 */
public class latentResults {
    public static void main(String [] args) throws Exception{
        int numRuns = 25;
        int numVars = 50;
        int sampleSize = 500;
        int numLatents = 10;
        int numSubsamples = 10;
        int numSubSets = 4;
        //String [] algs = {"FCI","MGM-FCI-MAX","Latent_FCI","Latent_MGM-FCI-MAX"};
        //String [] algs = {"Latent_FCI","Latent_MGM-FCI-MAX"};
        String [] algs = {"FCI","MGM-FCI-MAX","CPSS-FCI","CPSS-MGM-FCI-MAX","Bootstrap-FCI","Bootstrap-MGM-FCI-MAX"};
        PrintStream [] ps = new PrintStream[algs.length];
        PrintStream [] ps2 = new PrintStream[algs.length];
        /************************************************
         *
         * All True Latent Printing, can remove this later
         *************************************************/
            PrintStream out = new PrintStream("All_Latents.txt");


        for(int i = 0; i < algs.length;i++) {
            ps[i] = new PrintStream(algs[i] + "_Distribution.txt");
            ps[i].println("Run\tEdge Type\tStability\tTrue_Edge\tIdentifiable");
            if (!algs[i].contains("CPSS") && !algs[i].contains("Bootstrap")) {
                ps2[i] = new PrintStream(algs[i] + "_Orientations.txt");
                ps2[i].println("Run\tEdge Type\tStability\tTrue_Edge\tIdentifiable\tRule_1\tRule_2\tRule_1_Stability\tRule_2_Stability\tRule_1_Stability_Precise\tRule_2_Stability_Precise");
            }
        }
        for (int i = 0; i < numRuns; i++) {


            DataSet d = MixedUtils.loadDataSet2("Data/Data_" + i + "_" + numVars + "_" + sampleSize + "_" + numLatents + ".txt");
            Graph g = GraphUtils.loadGraphTxt(new File("Graphs/Graph_" + i + "_" + numVars + "_" + numLatents + ".txt"));
            for(int x = 1; x <= numVars;x++)
            {
                if(d.getVariable("X" + x)==null)
                    g.getNode("X" + x).setNodeType(NodeType.LATENT);
            }
            DagToPag pg = new DagToPag(g);
            pg.setCompleteRuleSetUsed(false);
            Graph truePag = pg.convert();


            /***********************************************************
             *
             * Remove this later, prints document with all latent variable edges and whether they are identifiable
             ************************************************************/
                ArrayList<LatentPrediction.Pair> list2 = LatentPrediction.getLatents(g, d, "CC");

                for (LatentPrediction.Pair p : list2) {
                    out.print(i + "\tCC\t");
                    Edge e = truePag.getEdge(truePag.getNode(p.one.getName()), truePag.getNode(p.two.getName()));
                    if (e != null && e.getEndpoint1() == Endpoint.ARROW && e.getEndpoint2() == Endpoint.ARROW)
                        out.println("T");
                    else
                        out.println("F");
                }

                list2 = LatentPrediction.getLatents(g, d, "CD");

                for (LatentPrediction.Pair p : list2) {
                    out.print(i + "\tCD\t");
                    Edge e = truePag.getEdge(truePag.getNode(p.one.getName()), truePag.getNode(p.two.getName()));
                    if (e != null && e.getEndpoint1() == Endpoint.ARROW && e.getEndpoint2() == Endpoint.ARROW)
                        out.println("T");
                    else
                        out.println("F");
                }

                list2 = LatentPrediction.getLatents(g, d, "DD");

                for (LatentPrediction.Pair p : list2) {
                    out.print(i + "\tDD\t");
                    Edge e = truePag.getEdge(truePag.getNode(p.one.getName()), truePag.getNode(p.two.getName()));
                    if (e != null && e.getEndpoint1() == Endpoint.ARROW && e.getEndpoint2() == Endpoint.ARROW)
                        out.println("T");
                    else
                        out.println("F");
                }

            out.flush();
            /*************************************************************
             *
             * End all latents printing
             ************************************************************/
            for (int j = 0; j < algs.length; j++) {
                //These maps contain the information for any edges that we predict to be latent variable edges
                HashMap<String,String>types = new HashMap<String,String>();
                HashMap<String,Double> stabilities = new HashMap<String,Double>();
                HashMap<String,String> truth = new HashMap<String,String>();
                HashMap<String,String>ident = new HashMap<String,String>();
                System.out.println(algs[j] + "\t" + i);


                BufferedReader b = new BufferedReader(new FileReader("Estimated/" + algs[j] + "_" + i + "_" + numVars + "_" + sampleSize + "_" + numLatents + "_" + numSubSets + ".txt"));
                while(b.ready())
                {
                    String [] line = b.readLine().split("\t");
                    if(line[0].equals("Orientations"))
                        break;
                    else if (line[0].equals("Splits")&&algs[j].contains("Latent"))
                        break;
                    ps[j].print(i + "\t");


                   //Get Edge type
                    if(d.getVariable(line[0])instanceof ContinuousVariable&& d.getVariable(line[1]) instanceof ContinuousVariable) {
                        types.put(line[0]+","+line[1],"CC");
                        ps[j].print("CC\t");
                    }
                    else if(d.getVariable(line[0])instanceof DiscreteVariable && d.getVariable(line[1]) instanceof DiscreteVariable) {
                        types.put(line[0]+","+line[1],"DD");
                        ps[j].print("DD\t");
                    }
                    else {
                        types.put(line[0]+","+line[1],"CD");
                        ps[j].print("CD\t");
                    }
                    stabilities.put(line[0]+","+line[1],Double.parseDouble(line[2]));
                    ps[j].print(line[2]+"\t");
                    ArrayList<LatentPrediction.Pair> list =  LatentPrediction.getLatents(g,d,"All");



                    //list has every single latent variable in it
                    //TODO Do i need to use the splits for anything here? Or do anything differently for Latents?
                    boolean found = false;
                        for (LatentPrediction.Pair p : list) {
                            if (p.one.getName().equals(line[0]) && p.two.getName().equals(line[1]))
                                found = true;
                            else if (p.two.getName().equals(line[0]) && p.one.getName().equals(line[1]))
                                found = true;
                        }
                        if(found) {
                            truth.put(line[0]+","+line[1],"T");
                            ps[j].print("T\t");
                        }
                        else {
                            truth.put(line[0]+","+line[1],"F");
                            ps[j].print("F\t");
                        }
                        Edge e = truePag.getEdge(truePag.getNode(line[0]),truePag.getNode(line[1]));
                        if(e!=null && e.getEndpoint1()== Endpoint.ARROW && e.getEndpoint2()==Endpoint.ARROW)
                        {
                           ident.put(line[0]+","+line[1],"T");
                            ps[j].println("T");
                        }
                        else {
                            ident.put(line[0]+","+line[1],"F");
                            ps[j].println("F");
                        }
                }
                if(algs[j].contains("Latent"))
                {
                    while(b.ready())
                    {
                        if(b.readLine().equals("Orientations"))
                            break;
                    }
                }
                HashMap<String,String> map = new HashMap<String,String>();

                System.out.println(types + "\n" + stabilities + "\n" + truth + "\n" + ident);
                while(b.ready())
                {
                    String [] line = b.readLine().split("\t");
                    String result = line[1];
                    for(int t = 2; t < line.length;t++)
                        result+=("\t"+line[t]);
                    map.put(line[0],result);
                }

                if(!algs[j].contains("CPSS") && !algs[j].contains("Bootstrap")) {
                    for (String x : map.keySet()) {
                        String[] nodes = x.split(",");


                        //If we get to this if statement, then we should have predicted this to be a latent variable edge at least in the final graph,
                        //but not necessarily above the threshold of tao for the full sample graph

                        if (map.get(nodes[1] + "," + nodes[0]) != null) {
                            String[] zero = map.get(x).split("\t");
                            String[] one = map.get(nodes[1] + "," + nodes[0]).split("\t");
                            HashMap<String, String> rules0 = new HashMap<String, String>();
                            HashMap<String, String> rules1 = new HashMap<String, String>();
                            String condSet0 = "";
                            String condSet1 = "";
                            for (int jj = 0; jj < zero.length; jj++) {
                                String[] temp = zero[jj].split(",");
                                String z = temp[1];
                                for (int k = 2; k < temp.length; k++)
                                    z += ("," + temp[k]);

                                rules0.put(temp[0], z);
                            }

                            for (int jj = 0; jj < one.length; jj++) {
                                String[] temp = one[jj].split(",");
                                String z = temp[1];
                                for (int k = 2; k < temp.length; k++)
                                    z += ("," + temp[k]);

                                rules1.put(temp[0], z);
                            }
                            double stab0 = 0;
                            double stab1 = 0;

                            if (rules0.get("FINAL") == null || rules1.get("FINAL") == null)
                                continue;

                            int rule0 = Integer.parseInt(rules0.get("FINAL").split(",")[0]);
                            int rule1 = Integer.parseInt(rules1.get("FINAL").split(",")[0]);
                            for (String key : rules0.keySet()) {
                                if (!key.equals("FINAL")) {
                                    String res = rules0.get(key).split(",")[0];
                                    if (Integer.parseInt(res) == rule0)
                                        stab0 += 1 / (double) numSubsamples;
                                } else {
                                    condSet0 = rules0.get(key);//res = 0,X40,X35
                                }
                            }
                            for (String key : rules1.keySet()) {
                                if (!key.equals("FINAL")) {
                                    String res = rules1.get(key).split(",")[0];
                                    if (Integer.parseInt(res) == rule1)
                                        stab1 += 1 / (double) numSubsamples;
                                } else {
                                    condSet1 = rules1.get(key);
                                }
                            }
                            double precStab0 = 0; //precise stability
                            double precStab1 = 0;
                            System.out.println(rules0 + "\t" + condSet0);
                            for (String key : rules0.keySet()) {
                                if (key.equals("FINAL"))
                                    continue;
                                String[] real = condSet0.split(",");
                                String[] curr = rules0.get(key).split(",");
                                if (rule0 != Integer.parseInt(curr[0]))
                                    continue;
                                ArrayList<String> realList = new ArrayList<String>();
                                ArrayList<String> currList = new ArrayList<String>();
                                for (int a = 1; a < real.length; a++)
                                    realList.add(real[a]);
                                for (int a = 1; a < curr.length; a++)
                                    currList.add(curr[a]);

                                int intersection = 0;
                                int union = 0;
                                Set<String> tempSet = new HashSet<String>();
                                for (int a = 0; a < realList.size(); a++) {
                                    if (currList.contains(realList.get(a)))
                                        intersection++;
                                    tempSet.add(realList.get(a));

                                }
                                for (int a = 0; a < currList.size(); a++)
                                    tempSet.add(currList.get(a));

                                union = tempSet.size();
                                if (realList.isEmpty() && currList.isEmpty())
                                    precStab0 += 1;
                                else
                                    precStab0 += intersection / (double) union;

                            }
                            for (String key : rules1.keySet()) {
                                if (key.equals("FINAL"))
                                    continue;
                                String[] real = condSet1.split(",");
                                String[] curr = rules1.get(key).split(",");
                                if (rule1 != Integer.parseInt(curr[0]))
                                    continue;
                                ArrayList<String> realList = new ArrayList<String>();
                                ArrayList<String> currList = new ArrayList<String>();
                                for (int a = 1; a < real.length; a++)
                                    realList.add(real[a]);
                                for (int a = 1; a < curr.length; a++)
                                    currList.add(curr[a]);

                                int intersection = 0;
                                int union = 0;
                                Set<String> tempSet = new HashSet<String>();
                                for (int a = 0; a < realList.size(); a++) {
                                    if (currList.contains(realList.get(a)))
                                        intersection++;
                                    tempSet.add(realList.get(a));

                                }
                                for (int a = 0; a < currList.size(); a++)
                                    tempSet.add(currList.get(a));

                                union = tempSet.size();
                                if (realList.isEmpty() && currList.isEmpty())
                                    precStab1 = 1;
                                else
                                    precStab1 += intersection / (double) union;

                            }
                            precStab0 = precStab0 / numSubsamples;
                            precStab1 = precStab1 / numSubsamples;
                            //intersection over union for each subsample * 1/numSubsamples
                            //sum all of those to get overall stability (if different rule than the full graph rule, contributes 0 to stability)
                            //How to compare final set of variables to other sets of variables of potentially different sizes with partial overlap


                            String edgeType = "";
                            double stability = 0;
                            String trueEdge = "";
                            String identifiable = "";
                            String o = nodes[0] + "," + nodes[1];
                            String o2 = nodes[1] + "," + nodes[0];
                            if (types.get(o) != null)
                                edgeType = types.get(o);
                            else if (types.get(o2) == null)
                                continue;
                            else
                                edgeType = types.get(o2);
                            if (stabilities.get(o) != null)
                                stability = stabilities.get(o);
                            else
                                stability = stabilities.get(o2);
                            if (truth.get(o) != null)
                                trueEdge = truth.get(o);
                            else
                                trueEdge = truth.get(o2);
                            if (ident.get(o) != null)
                                identifiable = truth.get(o);
                            else
                                identifiable = truth.get(o2);

                            ps2[j].println(i + "\t" + edgeType + "\t" + stability + "\t" + trueEdge + "\t" + identifiable + "\t" + rule0 + "\t" + rule1 + "\t" + stab0 + "\t" + stab1 + "\t" + precStab0 + "\t" + precStab1);


                        }
                    }
                    List<LatentPrediction.Pair> l = LatentPrediction.getLatents(g,d,"All");
                    for(int jj = 0; jj < l.size();jj++)
                    {
                        if(map.get(l.get(jj).one.getName() + "," + l.get(jj).two.getName())==null)
                        {

                            if(l.get(jj).one instanceof ContinuousVariable && l.get(jj).two instanceof ContinuousVariable)
                                ps2[j].println(i + "\t" + "CC" + "\t" + 0 + "\tT\t" + LatentPrediction.isIdentifiable(truePag,l.get(jj).one,l.get(jj).two) + "\t0\t0\t0\t0\t0\t0");
                            else if(l.get(jj).one instanceof DiscreteVariable && l.get(jj).two instanceof DiscreteVariable)
                                ps2[j].println(i + "\t" + "DD" + "\t" + 0 + "\tT\t" + LatentPrediction.isIdentifiable(truePag,l.get(jj).one,l.get(jj).two) + "\t0\t0\t0\t0\t0\t0");
                            else
                                ps2[j].println(i + "\t" + "CD" + "\t" + 0 + "\tT\t" + LatentPrediction.isIdentifiable(truePag,l.get(jj).one,l.get(jj).two) + "\t0\t0\t0\t0\t0\t0");

                        }

                    }


                }
                b.close();
                ps[j].flush();
            }
        }
        for(int i = 0; i < algs.length;i++)
            ps[i].close();
    }
}
