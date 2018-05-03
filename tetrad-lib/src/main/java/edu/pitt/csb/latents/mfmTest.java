package edu.pitt.csb.latents;

import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.algcomparison.simulation.Parameters;
import edu.cmu.tetrad.algcomparison.simulation.Simulation;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;
import edu.pitt.csb.stability.StabilityUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by vinee_000 on 4/26/2017.
 */
public class mfmTest {

    public static String [] scoreTypes = {"ADJ_PREC","ADJ_REC","OR_PREC","OR_REC","SHD"};
    public static String [] types = {"CC","CD","DD","All"};
    public static void main(String [] args) throws Exception
    {
        boolean numEdgesRandom = true;
        int numLambdas = 40;
        int numVariables = 50;
        int numLatents = 5;
        int numEdges = 75;
        int sampleSize = 200;
        int numSubsamples = 10;
        int numRuns = 11;
        int index = 0;
        int numCategories = 4;
        double gamma = 0.05;
        boolean saveData = true;
        boolean reuseData = true;
        boolean rerunAlgorithms = false;
        boolean includeBCCD = true;
        String directory = ".";
        double [] alphas = {0.001,0.01,0.03,0.05};
        String [] p_loc = {".25",".5",".75"};
        String [] p_edge = {".25",".5",".75"};
        String [] algs = {"MGM-FCI-MAX","FCI","CFCI","MGM-CFCI","MGM-FCI","FCI-MAX"};

      //   String [] algs = {"MGM-FCI-MAX","FCI","MGM-FCI"};

        while(index < args.length)
        {
            if(args[index].equals("-rd"))
            {
                reuseData = true;
                index++;
            }
            else if(args[index].equals("-alphas"))
            {
                index++;
                ArrayList<Double>alps = new ArrayList<Double>();
                while(index < args.length && !args[index].contains("-")) {
                    alps.add(Double.parseDouble(args[index]));
                    index++;
                }
            }
            else if(args[index].equals("-d")) {
                directory = args[index + 1];
                index += 2;
            }
            else if(args[index].equals("-nc"))
            {
                numCategories = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-sv"))
            {
                saveData = true;
                index++;
            }
            else if(args[index].equals("-l"))
            {
                numLambdas = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-v"))
            {
                numVariables = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-e"))
            {
                numEdges = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-s"))
            {
                sampleSize = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-ns"))
            {
                numSubsamples = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-r"))
            {
                numRuns = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-g"))
            {
                gamma = Double.parseDouble(args[index+1]);
                index+=2;
            }
            else if(args[index].equals("-nl"))
            {
                numLatents = Integer.parseInt(args[index+1]);
                index+=2;
            }
            else
            {
                System.out.println("Couldn't parse input argument " + args[index] + "...using default");
                index++;
            }
        }
        double [] initLambdas = new double[numLambdas];
        double low = .05;
        double high = .9;
        double inc = (high-low)/(numLambdas-1);
        for(int i = 0; i < numLambdas;i++)
        {
            initLambdas[i] = i*inc + low;
        }

        //Type: ["CC","CD","DD","All"] * Adjacency Precision, Adj Recall, Orientation Prec, Orientation Rec, SHD
        double [][][][] results = new double[numRuns][algs.length][alphas.length][20];
        double [][][][] bccdResults = new double[numRuns][p_loc.length][p_edge.length][20];

        File rFile = new File("Results");
        if(!rFile.isDirectory())
            rFile.mkdir();
        File gFile = new File("Graphs");
        File dFile = new File("Data");
        File estFile = new File("Estimated");
        File subFile = new File("Subsamples");
        if(saveData) {
            if (!gFile.isDirectory())
                gFile.mkdir();
            if (!dFile.isDirectory())
                dFile.mkdir();
            if(!estFile.isDirectory())
                estFile.mkdir();
            if(!subFile.isDirectory())
                subFile.mkdir();
        }

        /*PrintStream pri = new PrintStream(directory + "/mgm_priors_" + amountPrior + "_" + numExperts + "_" + numVariables +  "_" + sampleSize + "_" + numSubsamples + ".txt");
        PrintStream step = new PrintStream(directory + "/STEPS_" + amountPrior + "_" + numExperts + "_" + numVariables + "_"  + sampleSize + "_" + numSubsamples + ".txt");
        PrintStream orc = new PrintStream(directory + "/oracle_" + amountPrior + "_" + numExperts + "_" + numVariables + "_" + sampleSize + "_" + numSubsamples + ".txt");
        PrintStream one = new PrintStream(directory + "/mgm_one_steps_" + amountPrior + "_" + numExperts + "_" + numVariables + "_" + sampleSize + "_" + numSubsamples + ".txt");
        PrintStream orcOne = new PrintStream(directory + "/oracle_one_" + amountPrior + "_" + numExperts + "_" + numVariables + "_" + sampleSize + "_" + numSubsamples + ".txt");*/
        NormalDistribution n = new NormalDistribution(numVariables*2,numVariables/2);
        PrintStream [][] pri = new PrintStream [algs.length][alphas.length];
        for(int i = 0; i < algs.length;i++) {
            for(int j = 0; j < alphas.length;j++) {
                pri[i][j] = new PrintStream("Results/" + algs[i] + "_" + alphas[j] + "_" + numVariables + "_" + sampleSize + "_" + numSubsamples + "_" + numLatents + ".txt");
            }
        }
        A:for(int i = 0; i < numRuns; i++) {
            MixedLeeHastieSimulation c = new MixedLeeHastieSimulation();
            System.out.println(i);
            if(numEdgesRandom)
            {
                numEdges = (int)n.sample();
                while(numEdges <(numVariables/4))
                    numEdges = (int)n.sample();
            }
            Parameters p = new Parameters();
            p.setValue("numMeasures", numVariables);
            p.setValue("numEdges", numEdges);
            p.setValue("sampleSize", sampleSize);
            p.setValue("numCategories",numCategories);
            c.simulate(p);


            ///Allen this is the section where we remove latent variables
            Graph currPAG = addLatents(c,numLatents);


            DataSet [] subsamples = new DataSet[numSubsamples];
            if(reuseData)
            {
                boolean foundFile = false;
                boolean foundData = false;
                File f = new File("Graphs/Graph_" + i + "_" + numVariables + ".txt");
                if(f.exists()) {
                    foundFile = true;
                    c.setTrueGraph(GraphUtils.loadGraphTxt(f));
                }
                f = new File("Data/Data_" + i + "_" + numVariables + "_" + sampleSize + "_" + numLatents +  ".txt");
                if(f.exists() && foundFile) {
                    foundData = true;
                    c.setDataSet(MixedUtils.loadDataSet2("Data/Data_" + i + "_" + numVariables + "_" + sampleSize + "_" + numLatents + ".txt"), 0);
                    DataSet temp = c.getDataSet(0);
                    Graph gTemp = c.getTrueGraph();
                    //Specifying which latent variables are there, to get the correct current PAG
                    for(int x = 0;x < gTemp.getNodes().size();x++)
                    {
                        if(temp.getVariable(gTemp.getNodes().get(x).getName())==null)
                            gTemp.getNodes().get(x).setNodeType(NodeType.LATENT);
                    }
                    c.setTrueGraph(gTemp);
                    DagToPag p2 = new DagToPag(gTemp);
                    p2.setCompleteRuleSetUsed(false);
                    currPAG = p2.convert();
                }
                else if(foundFile) {
                    c.simulate(p);
                    currPAG = addLatents(c,numLatents);
                }
                if(foundData) {
                    for (int j = 0; j < numSubsamples; j++) {
                        f = new File("Subsamples/Subsample_" + i + "_" + numVariables + "_" + sampleSize + "_" + numLatents + "_" +  j + ".txt");
                        if(f.exists())
                            subsamples[j] = MixedUtils.loadDataSet2(f.getAbsolutePath());
                    }
                }

            }




            Graph [][] estimatedGraphs = new Graph[algs.length][alphas.length];

            if(!rerunAlgorithms) {
                for (int j = 0; j < algs.length; j++) {
                    for(int k = 0; k < alphas.length;k++) {
                        File f = new File("Estimated/" + algs[j] + "_" + alphas[k] + "_" + i + "_" + numVariables + "_" + sampleSize + "_" + numLatents + ".txt");

                        if (f.exists()) {
                            Graph est = GraphUtils.loadGraphTxt(f);
                            estimatedGraphs[j][k] = est;
                        }
                    }
                }
            }
                System.out.println(i);
                try {
                    Graph mgmSteps = null;
                    boolean foundSteps = false;
                    //String [] algs = {"mgm_one_steps","oracle_one","STEPS","oracle","mgm_priors"};
                    for(int pp = 0; pp < algs.length;pp++)
                    {
                        for(int k = 0; k < alphas.length;k++) {
                            if ((algs[pp].contains("MGM") && estimatedGraphs[pp][k] == null))
                                foundSteps = true;
                        }
                    }

                    boolean nullSub = false;
                    for(int j = 0; j < subsamples.length;j++)
                    {
                        if(subsamples[j]==null)
                        {
                            nullSub = true;
                        }
                    }
                    if(nullSub) {
                        int b = (int) Math.floor(10 * Math.sqrt(c.getDataSet(0).getNumRows()));
                        if (b > c.getDataSet(0).getNumRows())
                            b = c.getDataSet(0).getNumRows() / 2;
                        int[][] samps = StabilityUtils.subSampleNoReplacement(c.getDataSet(0).getNumRows(), b, numSubsamples);
                        for (int j = 0; j < subsamples.length; j++) {
                            subsamples[j] = c.getDataSet(0).subsetRows(samps[j]);
                        }
                    }
                    STEPS s = new STEPS(c.getDataSet(0), initLambdas, gamma, subsamples);
                    Graph steps = null;
                    if(foundSteps)
                        steps = s.runStepsPar();

                    for(int j = 0; j < algs.length;j++) {
                        for (int k = 0; k < alphas.length; k++) {
                            if (algs[j].equals("MGM-FCI-MAX")) {
                                if (estimatedGraphs[j][k] == null) {

                                    IndependenceTest ii = new IndTestMultinomialAJ(c.getDataSet(0), alphas[k]);
                                    FciMaxP f = new FciMaxP(ii);
                                    f.setInitialGraph(steps);
                                    estimatedGraphs[j][k] = f.search();
                                }


                            } else if (algs[j].equals("MGM-FCI")) {
                                if (estimatedGraphs[j][k] == null) {
                                    IndependenceTest ii = new IndTestMultinomialAJ(c.getDataSet(0), alphas[k]);
                                    Fci f = new Fci(ii);
                                    f.setInitialGraph(steps);
                                    estimatedGraphs[j][k] = f.search();
                                }
                            } else if (algs[j].equals("FCI")) {
                                if (estimatedGraphs[j][k] == null) {
                                    IndependenceTest ii = new IndTestMultinomialAJ(c.getDataSet(0), alphas[k]);
                                    Fci f = new Fci(ii);
                                    estimatedGraphs[j][k] = f.search();                                }
                            } else if (algs[j].equals("FCI-MAX")) {
                                if (estimatedGraphs[j][k] == null) {
                                    IndependenceTest ii = new IndTestMultinomialAJ(c.getDataSet(0), alphas[k]);
                                    FciMaxP f = new FciMaxP(ii);
                                    estimatedGraphs[j][k] = f.search();
                                }
                            } else if (algs[j].equals("CFCI")) {
                                if (estimatedGraphs[j][k] == null)
                                {
                                    IndependenceTest ii = new IndTestMultinomialAJ(c.getDataSet(0), alphas[k]);
                                    Cfci f = new Cfci(ii);
                                    estimatedGraphs[j][k] = f.search();

                                }
                            } else if (algs[j].equals("MGM-CFCI")) {
                                if (estimatedGraphs[j][k] == null)
                                {
                                    IndependenceTest ii = new IndTestMultinomialAJ(c.getDataSet(0), alphas[k]);
                                    Cfci f = new Cfci(ii);
                                    f.setInitialGraph(steps);
                                    estimatedGraphs[j][k] = f.search();

                                }
                            }
                           /* System.out.println(algs[j] + "," + mgmPriors.getF1(estimatedGraphs[j], c.getTrueGraph(), c.getDataSet(0), "CC"));
                            System.out.println(algs[j] + "," + mgmPriors.getF1(estimatedGraphs[j], c.getTrueGraph(), c.getDataSet(0), "CD"));
                            System.out.println(algs[j] + "," + mgmPriors.getF1(estimatedGraphs[j], c.getTrueGraph(), c.getDataSet(0), "DD"));
                            System.out.println(algs[j] + "," + mgmPriors.getF1(estimatedGraphs[j], c.getTrueGraph(), c.getDataSet(0), "All"));
*/

                            results[i][j][k] = getAllResults(estimatedGraphs[j][k], currPAG,c.getTrueGraph(), c.getDataSet(0));
                   //         done = true;
                        }
                    }
                    if(includeBCCD)
                    {
                        for(int j = 0; j < p_loc.length;j++)
                        {
                            for(int k = 0;k < p_edge.length;k++)
                            {

                                File tFile = new File("Estimated/BCCD_" + i + "_" + p_loc[j] + "_" + p_edge[k] + "_" + numVariables + "_" + sampleSize + "_" + numLatents + ".txt");
                                tFile.renameTo(new File("Estimated/BCCD_" +  "_" + p_loc[j] + "_" + p_edge[k] + "_" + i + "_" + numVariables + "_" + sampleSize + "_" + numLatents + ".txt"));
                                tFile = new File("Estimated/BCCD_" + "_" + p_loc[j] + "_" + p_edge[k] + "_" + i + "_" + numVariables + "_" + sampleSize + "_" + numLatents + ".txt");

                                Graph est = bccdToPAG(tFile,c.getDataSet(0));
                                bccdResults[i][j][k] = getAllResults(est,currPAG,c.getTrueGraph(),c.getDataSet(0));
                            }
                        }
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(-1);
                    if(!reuseData) {
                        c.simulate(p);
                    }
                }

            if(saveData) {
                PrintStream p2 = new PrintStream("Graphs/Graph_" + i + "_" + numVariables +  ".txt");
                p2.println(c.getTrueGraph());
                p2.flush();
                p2.close();
                p2 = new PrintStream("Data/Data_" + i + "_" + numVariables + "_" + sampleSize + "_" + numLatents + ".txt");
                p2.println(c.getDataSet(0));
                p2.flush();
                p2.close();
                for(int x = 0; x < algs.length;x++)
                {
                    for(int y = 0; y < alphas.length;y++) {
                        p2 = new PrintStream("Estimated/" + algs[x] + "_" + alphas[y] + "_" + i + "_" + numVariables + "_" + sampleSize + "_" + numLatents +  ".txt");
                        p2.println(estimatedGraphs[x][y]);
                        p2.flush();
                        p2.close();
                    }
                }

                for(int j = 0; j < subsamples.length;j++)
                {
                    p2 = new PrintStream("Subsamples/Subsample_" + i + "_" + numVariables + "_" + sampleSize + "_" + numLatents + "_" + j + ".txt");
                    p2.println(subsamples[j]);
                    p2.flush();
                    p2.close();
                }
            }
        }
        for(int i = 0; i < pri.length;i++)
        {
            for(int j = 0; j < pri[i].length;j++) {
                printData(pri[i][j], results, i,j);
            }
        }
        if(includeBCCD)
        {
            for(int i = 0; i < p_loc.length;i++)
            {
                for(int j = 0; j < p_edge.length;j++)
                {
                    PrintStream bccdPrints = new PrintStream("Results/BCCD_" + p_loc[i] + "_" + p_edge[j] + "_" + numVariables + "_" + sampleSize + "_" + numSubsamples + "_" + numLatents + ".txt");
                    printData(bccdPrints,bccdResults,i,j);
                }
            }
        }
    }
    //Need a printstream for each alpha as well
    public static void printData(PrintStream p, double [][][][] result,int algInd,int alpInd)
    {
        p.print("Run");
        for(int i = 0; i < scoreTypes.length;i++)
        {
            for(int j = 0; j < types.length;j++)
            {
                p.print("\t" + scoreTypes[i] + "_" + types[j]);
            }
        }
        p.println();

        //Result is run, algorithm, alpha, score Type
        for(int i = 0; i < result.length;i++)
        {
            p.print(i);
            for(int j = 0; j < scoreTypes.length*types.length;j++)
            {
               p.print("\t" + result[i][algInd][alpInd][j]);
            }
            p.println();
        }
        p.flush();
        p.close();
    }
    public static Graph bccdToPAG(File f, DataSet d)
    {
        try{
            BufferedReader b = new BufferedReader(new FileReader(f));
            int [][] adj = new int[d.getNumColumns()][d.getNumColumns()];
            int row = 0;
            while(b.ready())
            {
                String [] line = b.readLine().split("\t");
                for(int col = 0; col < line.length;col++)
                {
                    adj[row][col] = Integer.parseInt(line[col]);
                }
                row++;
            }
            return adjToPAG(adj,d);
        }
        catch(Exception e)
        {
            System.err.println("Couldn't load graph for bccd");
            e.printStackTrace();
            System.exit(-1);
        }
        return null;
    }
    public static Graph adjToPAG(int [][] adj, DataSet d)
    {
        Graph g = new EdgeListGraphSingleConnections(d.getVariables());
        for(int i = 0; i < adj.length;i++)
        {
            for(int j = i+1; j < adj[i].length;j++)
            {
                if(adj[i][j]!=0)
                {
                    //TODO Make sure this isn't reversed
                    g.addEdge(new Edge(g.getNode(d.getVariable(i).getName()),g.getNode(d.getVariable(j).getName()),toEndpoint(adj[i][j]),toEndpoint(adj[j][i])));
                }
            }
        }
        return g;
    }
    private static Endpoint toEndpoint(int i)
    {
        if(i==1)
            return Endpoint.TAIL;
        else if(i==2)
            return Endpoint.ARROW;
        else if(i==3)
            return Endpoint.CIRCLE;
        else
            return null;
    }
    public static String compareGraphs(Graph truth, Graph est)
    {
        //     System.out.println("Truth: " + truth);
        //   System.out.println("Est: " + est);
        String x = "\t---\tNo_Edge";
        x+="\n---\t";
        int total = 0;
        int count = 0;
        for(Edge e: truth.getEdges())
        {
            if(est.getEdge(est.getNode(e.getNode1().getName()),est.getNode(e.getNode2().getName()))!=null)
                count++;
            else
                total++;

        }
        x+= count + "\t" + total + "\n";
        x+= "No_Edge\t";
        count = 0;
        for(Edge e: est.getEdges())
        {
            if(truth.getEdge(truth.getNode(e.getNode1().getName()),truth.getNode(e.getNode2().getName()))==null)
                count++;
        }
        x+= count;

        return x;
    }





    public static double [] getAllResults(Graph est, Graph truth, Graph DAG, DataSet dat)
    {
        int TPU = 0;
        int FPU = 1;
        int FNU = 2;
        int ETP = 3;
        int EFP = 4;
        int EFN = 5;
        int ETN = 6;
        int SHD = 7;
        double [] res = new double[20];
        DagToPag p = new DagToPag(truth);
        p.setCompleteRuleSetUsed(false);
        Graph comp = p.convert();
        int count = 0;
        double [][] allResults = MixedUtils.newLatentScores(truth,est,DAG,dat,false);
        for(int j = 0; j < 5;j++)//Loop over types of scores, adj prec, adj rec, orientation prec, orientation rec, SHD
        {

            for (int i = 0; i < types.length; i++) {

                if(j==0)
                {
                    res[count] = allResults[i][TPU]/(allResults[i][TPU]+allResults[i][FPU]);
                    count++;
                }
                else if(j==1)
                {
                    res[count] = allResults[i][TPU]/(allResults[i][TPU]+allResults[i][FNU]);
                    count++;
                }
                else if(j==2)
                {
                    res[count] = allResults[i][ETP]/(allResults[i][ETP]+allResults[i][EFP]);
                    count++;
                }
                else if(j==3)
                {
                    res[count] = allResults[i][ETP]/(allResults[i][ETP]+allResults[i][EFN]);
                    count++;
                }
                else
                {
                    res[count] = allResults[i][SHD];
                    count++;
                }
                //Adj Prec for each edge type

            }
        }
        return res;
    }
    //Add latents to the dataset in c based on the graph in c and the specified number of latents
    //Each latent has at least two children
    public static Graph addLatents(MixedLeeHastieSimulation c, int numLatents)
    {
        Graph g = c.getTrueGraph(); //Get true graph
        boolean done = false;

        Graph p;
        int numNodesToRemove = numLatents; //The number of nodes we need to remove

        int totalNodes = g.getNumNodes();
        int numNodesRemoved = 0;
        Random rand = new Random();
        int count = 0;
        ArrayList<Integer> usedInts = new ArrayList<Integer>();
        ArrayList<Integer> removeInts = new ArrayList<Integer>();
        boolean reset = false;

        //Loop until enough latents are identified
        A:
        while (numNodesToRemove > 0) {
            count++;
            if (count == 1000) {
                count = 0;
                System.out.println("Stuck on run");
                reset = true;
                break A;
            }

            int nod = rand.nextInt(totalNodes);
            if (usedInts.contains(nod))
                continue;
            else
                usedInts.add(nod);
            if (g.getOutdegree(g.getNode("X" + (nod + 1))) > 1) {

                numNodesToRemove--;
                numNodesRemoved++;

                removeInts.add(nod);
                totalNodes--; //not sure if i need this

            }

        }
        if (reset)
            System.exit(0);
        //  continue;
        Collections.sort(removeInts);
        Collections.reverse(removeInts);

        ArrayList<Node> latents = new ArrayList<Node>();
        for (int curr = 0; curr < removeInts.size(); curr++) {
            latents.add(g.getNode("X" + (removeInts.get(curr) + 1)));
            g.getNode("X" + (removeInts.get(curr) + 1)).setNodeType(NodeType.LATENT);
        }
        usedInts.clear();

        c.setTrueGraph(g);
        DataSet temp = c.getDataSet(0);
        int [] datRemove = new int[numLatents];
        count = 0;
        for(int i = 0; i < temp.getNumColumns();i++)
        {
            if(g.getNode(temp.getVariable(i).getName()).getNodeType()==NodeType.LATENT) {
                datRemove[count] = i;
                count++;
            }
        }
        temp.removeCols(datRemove);
        c.setDataSet(temp,0);



        final DagToPag dagToPag = new DagToPag(g);
        dagToPag.setCompleteRuleSetUsed(false);
        p = dagToPag.convert(); //This is the graph we will compare against
        return p;
    }

}
