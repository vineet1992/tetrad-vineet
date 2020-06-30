package edu.pitt.csb.stability;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.jet.math.Functions;
import edu.cmu.tetrad.data.ColtDataSet;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.Fci;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.StatUtils;
import edu.pitt.csb.Priors.runPriors;
import edu.pitt.csb.latents.LatentPrediction;
import edu.pitt.csb.mgm.Algorithm;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.stability.DataGraphSearch;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

/**
 * Created by vinee_000 on 7/12/2018.
 * Class to run Complimentary Pairs Stability Selection with causal modeling algorithms
 *  Currently how this works is for bootstrap we only add an edge to the final graph if it showed up in a particular orientation greater than bound times
 *  For CPSS we first compute theta using average number of edges in the graph (edges are selections)
 *  But then we use our bound the same way as Bootstrapping
 *  TODO Make the bound match up with edge selection, and then orient by different strategies?
 *  Use bound as an initial graph and then run algorithm as normal on full data?
 *  Use bound to select edges and then orient according to most likely orientation?
 */
public class Bootstrap {
    private int B; //Number of bootstrap samples
    private int p;//total number of variables
    private double bound; //percentage of bootstrap samples cutoff
    private DataSet data; //Dataset
    private double [] lambda; //Array of lambda values for MGM
    private int [][] subs; //subsampled indicies
    private String tabularOutput; //Output of edge orientation frequencies in a tabular form
    private boolean cpss; //Should we use Complimentary pairs to create datasets?
    private boolean subsample; //Should we use subsampling to create datasets?
    private int b = 0; //Subsample size
    private boolean fullGraph;

    public Bootstrap (DataSet data, double bound, int B)
    {
        this.data = data;
        this.bound = bound;
        this.B = B;
        this.p = data.getNumColumns()*(data.getNumColumns()-1)/2;
    }
    public Bootstrap(DataSet data, double [] lambda, double bound, int B)
    {
        this.data = data;
        this.lambda = lambda;
        this.bound = bound;
        this.B = B;
        this.p = data.getNumColumns()*(data.getNumColumns()-1)/2;
    }

    public Bootstrap(DataSet data, double [] lambda, double bound, int [][]subs)
    {
        this.data = data;
        this.lambda = lambda;
        this.bound = bound;
        this.p = data.getNumColumns()*(data.getNumColumns()-1)/2;
        this.subs = subs;
        this.B = subs.length;
    }

    public String getTabularOutput() {return tabularOutput;}
    public void setFullGraph(){fullGraph=true;}

    public void setCPSS() {
        cpss = true;
        B = B*2;
    }
    public void setSubsample(int b)
    {
        subsample = true;
        this.b = b;
    }
    public static int[][] createSubs(DataSet data, int B)
    {
        return DataUtils.getBootstrapIndices(data,data.getNumRows(),B);
    }
    private int[][] createSubs(DataSet data)
    {
        subs = createSubs(data,B);
        return subs;
    }
    private int [][] createSubsSubsample(DataSet data)
    {
        if(b==0)
            this.b = StabilityUtils.getSubSize(data.getNumRows());
        return StabilityUtils.subSampleNoReplacement(data.getNumRows(),b,B);
    }
    private int[][] createSubsCPSS(DataSet data) {
        ArrayList<Integer>tempInds  = new ArrayList<Integer>();
        for(int i = 0; i < data.getNumRows();i++)
            tempInds.add(i);
        int [][] subs = new int[B][];
        for(int i = 0; i < B/2;i++)
        {
            Collections.shuffle(tempInds);
            int[] d1 = new int[tempInds.size() / 2];
            for (int j = 0; j < d1.length; j++) {
                d1[j] = tempInds.get(j);
            }
            int size2 = tempInds.size() / 2;
            if (tempInds.size() % 2 == 1)
                size2 = tempInds.size() / 2 + 1;
            int[] d2 = new int[size2];
            for (int j = d1.length; j < tempInds.size(); j++) {
                d2[j - d1.length] = tempInds.get(j);
            }

            subs[2*i] = d1;
            subs[2*i+1] = d2;
            if(runPriors.checkForVariance(data.subsetRows(d1),data)!=-1 || runPriors.checkForVariance(data.subsetRows(d2),data)!=-1)
                --i;
        }
        return subs;
    }


    //Right now Bootstrap is done in a directed sense, whereas CPSS is done fully undirected
    public Graph runBootstrap(final Algorithm a, final double [] params)
    {
        if(subs==null && !cpss && !subsample)
            subs = createSubs(data);
        else if(subs==null && cpss)
            subs = createSubsCPSS(data);
        else if(subs==null && subsample)
            subs = createSubsSubsample(data);

        final ArrayList<Graph> graphs = new ArrayList<Graph>();
        final double [] totalEdges = new double[B];
        final double [][][] theta = new double[8][data.getNumColumns()][data.getNumColumns()];
        if(cpss)
            System.out.println("Computing " + B + " Compliemntary Pairs Graphs in parallel... using " + a);
        else if(subsample)
            System.out.println("Computing " + B + " Subsampled Graphs in parallel... using " + a);
        else
            System.out.println("Computing " + B + " Bootstrap Graphs in parallel... using " + a);


        final ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

        class StabilityAction extends RecursiveAction {
            private int chunk;
            private int from;
            private int to;

            public StabilityAction(int chunk, int from, int to) {
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }

            private synchronized void addToGraphs(Graph g){graphs.add(g);}

            private synchronized DataSet subset(int [] x){return data.subsetRows(x);}
            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    A:for (int s = from; s < to; s++) {
                        DataSet data1 = subset(subs[s]);
                        System.out.println("Running sample " + s);
                        DataGraphSearch gs = Algorithm.algToSearchWrapper(a,params);
                        try {
                            Graph gt = gs.search(data1);
                            totalEdges[s] = gt.getNumEdges();
                            addToGraphs(gt);
                        }catch(Exception e)
                        {
                            continue A;
                        }
                    }

                    return;
                } else {
                    List<StabilityAction> tasks = new ArrayList<>();

                    final int mid = (to + from) / 2;

                    tasks.add(new StabilityAction(chunk, from, mid));
                    tasks.add(new StabilityAction(chunk, mid, to));

                    invokeAll(tasks);

                    return;
                }
            }

        }

        final int chunk = 5;

        pool.invoke(new StabilityAction(chunk, 0,B));



        for(int i = 0 ; i< graphs.size();i++)
        {
            Graph curr = graphs.get(i);
            for(int j = 0; j < data.getNumColumns();j++)
            {
                String node1 = data.getVariable(j).getName();
                for(int k = j+1; k < data.getNumColumns();k++)
                {
                    String node2 = data.getVariable(k).getName();
                    Edge currEdge = curr.getEdge(curr.getNode(node1),curr.getNode(node2));
                    theta[edgeToInt(currEdge,node1,node2)][j][k]++;
                    //TODO Debug this
                }
            }
            //Use stability utils orientation search
            //Include a function to convert double[][][] to a an output graph based on a bound
        }

        System.out.println("Done");
        //Partition the data B times and run MGM on each partition to get a count of each edge's appearence
        //Based on the average number of selected variables per run, compute theta
        //use theta, q, and the bound to get a value for tao
        //only keep the edges with probability greater than tao
        //Return the graph with these edges
        if(cpss) {
            double sum = StatUtils.sum(totalEdges);
            this.bound = CPSS.computeTao(p/sum,this.bound);
        }


        //Run this algorithm on the full graph and subset edges based on bootstrap probabilities
        Graph output;
        if(fullGraph)
        {
            DataGraphSearch gs = Algorithm.algToSearchWrapper(a,params);
            output = gs.search(data);
        }
        else
            output = createGraph(theta,data);
        tabularOutput = createStringOutput(output,theta,data);
        return output;
    }
    private String createStringOutput(Graph g, double [][][]theta, DataSet data)
    {
        String out = "Var 1\tInteraction\tVar 2\t---\t-->\t<--\t<->\to->\t<-o\to-o\tNo Edge\n";
        for(Edge e:g.getEdges())
        {
            int x = data.getColumn(data.getVariable(e.getNode1().getName()));
            int y = data.getColumn(data.getVariable(e.getNode2().getName()));
            String name1 = e.getNode1().getName();
            String name2= e.getNode2().getName();

            int ind1 = x;
            int ind2 = y;
            if(x > y)
            {
                name1 = e.getNode2().getName();
                name2 = e.getNode1().getName();

                int temp = ind2;
                ind2 = ind1;
                ind1 = temp;
            }
            out+= name1 +"\t" + edgeString(e) + "\t" + name2 +"\t";

            for(int i = 0; i < theta.length;i++)
            {
                if(i==theta.length-1)
                    out+=theta[i][ind1][ind2] +"\n";
                else
                    out+=theta[i][ind1][ind2] + "\t";
            }
        }
        return out;
    }
    private String edgeString(Edge e)
    {
        StringBuilder buf = new StringBuilder();

        Endpoint endptTypeA = e.getEndpoint1();
        Endpoint endptTypeB = e.getEndpoint2();


        if (endptTypeA == Endpoint.TAIL) {
            buf.append("-");
        } else if (endptTypeA == Endpoint.ARROW) {
            buf.append("<");
        } else if (endptTypeA == Endpoint.CIRCLE) {
            buf.append("o");
        }

        buf.append("-");

        if (endptTypeB == Endpoint.TAIL) {
            buf.append("-");
        } else if (endptTypeB == Endpoint.ARROW) {
            buf.append(">");
        } else if (endptTypeB == Endpoint.CIRCLE) {
            buf.append("o");
        }

        return buf.toString();
    }
    private Graph createGraph(double[][][]theta,DataSet data)
    {
        Graph temp = new EdgeListGraphSingleConnections(data.getVariables());
        for(int j = 0; j < data.getNumColumns();j++)
        {
            for(int k = j+1; k < data.getNumColumns();k++)
            {
                double max = 0;
                int index = -1;
                for(int i = 0; i < theta.length-1;i++)
                {
                    if(theta[i][j][k] > max) {
                        max = theta[i][j][k];
                        index = i;
                    }
                }
                if(max > bound)
                {
                    temp.addEdge(createEdge(data.getVariable(j),data.getVariable(k),index));
                }
            }
        }
        return temp;
    }
    private Edge createEdge(Node var1, Node var2, int type)
    {
        Endpoint e1 = null;
        Endpoint e2 = null;
        if(type==0 || type==1)
            e1 = Endpoint.TAIL;
        else if(type==2||type==3||type==5)
            e1 = Endpoint.ARROW;
        else
            e1 = Endpoint.CIRCLE;
        if(type==0||type==2)
            e2=Endpoint.TAIL;
        else if(type==1||type==3||type==4)
            e2=Endpoint.ARROW;
        else
            e2 = Endpoint.CIRCLE;
         return new Edge(var1,var2,e1,e2);
    }
    public ArrayList<LatentPrediction.Pair> runBootstrap(final String algName, double alpha)
    {
        final double tempAlpha = alpha;
        final double [] tempLambda = lambda;
        if(subs==null)
            subs = createSubs(data);
        final ArrayList<Graph> graphs = new ArrayList<Graph>();
        final DoubleMatrix2D edgeCounts = new SparseDoubleMatrix2D(data.getNumColumns(),data.getNumColumns());
        final int [] totalVars = new int[B];
        System.out.print("Computing " + B + " Graphs in parallel... using " + algName);



        final ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

        class StabilityAction extends RecursiveAction {
            private int chunk;
            private int from;
            private int to;

            public StabilityAction(int chunk, int from, int to) {
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }

            //could avoid using synchronized if we keep track of array of mats and add at end, but that needs lots of
            //memory

            private synchronized ArrayList<Integer> createTemp(ArrayList<Integer> inds) {
                ArrayList<Integer> temp = new ArrayList<Integer>();
                temp.addAll(inds);
                return temp;
            }

            private synchronized void addToGraphs(Graph g){graphs.add(g);}
            private synchronized void addToMat(Graph g, DataSet data, DoubleMatrix2D edgeCounts)
            {
                for(Edge e:g.getEdges())
                {
                    int x = data.getColumn(data.getVariable(e.getNode1().getName()));
                    int y = data.getColumn(data.getVariable(e.getNode2().getName()));
                    try {
                        edgeCounts.set(x, y, edgeCounts.get(x, y) + 1);
                        edgeCounts.set(y, x, edgeCounts.get(y, x) + 1);
                    }
                    catch(Exception e2)
                    {
                        System.out.println(e + "\t" + data + "\t" + g);
                        System.out.println(x + "\t" + y);
                        System.out.println(data.getVariable(e.getNode1().getName()));
                        System.out.println(data.getVariable(e.getNode2().getName()));
                        e2.printStackTrace();
                        System.exit(-1);
                    }
                }
            }

            private synchronized int getLatentCount(Graph g1)
            {
                int result = 0;
                for(Edge e:g1.getEdges())
                {
                    if(e.getEndpoint1()== Endpoint.ARROW && e.getEndpoint2()==Endpoint.ARROW)
                    {
                        result++;
                    }
                }
                return result;
            }
            private synchronized DataSet createData(){return new ColtDataSet(data.getNumColumns(),data.getVariables());}

            private synchronized DataSet subset(int [] x){return data.subsetRows(x);}
            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        DataSet data1 = subset(subs[s]);
                        Graph g1 = new EdgeListGraph(data.getVariables());
                        if(algName.equals("FCI") || algName.equals("MGM-FCI"))
                        {
                            IndependenceTest ind = new IndTestMultinomialAJ(data1,tempAlpha,true);
                            Fci f = new Fci(ind);
                            if(algName.equals("MGM-FCI"))
                            {
                                MGM m = new MGM(data1,tempLambda);
                                m.learnEdges(1000);
                                f.setInitialGraph(m.graphFromMGM());
                            }
                            g1 = f.search();
                        }
                        else
                        {
                            IndependenceTest ind = new IndTestMultinomialAJ(data1,tempAlpha,true);
                            FciMaxP f = new FciMaxP(ind);
                            MGM m = new MGM(data1,tempLambda);
                            m.learnEdges(1000);
                            f.setInitialGraph(m.graphFromMGM());
                            g1 = f.search();
                        }
                        //addToMat(g,data,edgeCounts);
                        addToGraphs(g1);
                    }

                    return;
                } else {
                    List<StabilityAction> tasks = new ArrayList<>();

                    final int mid = (to + from) / 2;

                    tasks.add(new StabilityAction(chunk, from, mid));
                    tasks.add(new StabilityAction(chunk, mid, to));

                    invokeAll(tasks);

                    return;
                }
            }

        }

        final int chunk = 2;

        pool.invoke(new StabilityAction(chunk, 0,B));
        for(int i = 0 ; i< graphs.size();i++)
        {
            Graph curr = graphs.get(i);
            for(Edge e: curr.getEdges())
            {
                if(e.getEndpoint1()==Endpoint.ARROW && e.getEndpoint2()==Endpoint.ARROW) {
                    int x = data.getColumn(data.getVariable(e.getNode1().getName()));
                    int y = data.getColumn(data.getVariable(e.getNode2().getName()));
                    edgeCounts.set(x, y, edgeCounts.get(x, y) + 1);
                    edgeCounts.set(y, x, edgeCounts.get(y, x) + 1);
                }
            }
        }
        System.out.println("Done");
        ArrayList<LatentPrediction.Pair> finalResult = new ArrayList<LatentPrediction.Pair>();
        for(int i = 0; i < edgeCounts.rows();i++)
        {
            for(int j = i+1; j < edgeCounts.columns();j++)
            {
                if(edgeCounts.get(i,j)/(B)>bound)
                {
                    finalResult.add(new LatentPrediction.Pair(data.getVariable(i),data.getVariable(j),edgeCounts.get(i,j)/(B)));
                }
            }
        }
        //Partition the data B times and run MGM on each partition to get a count of each edge's appearence
        //Based on the average number of selected variables per run, compute theta
        //use theta, q, and the bound to get a value for tao
        //only keep the edges with probability greater than tao
        //Return the graph with these edges
        return finalResult;
    }

    private int edgeToInt(Edge e, String name1, String name2)
    {

        if(e==null)
            return 7;//No edge

        boolean flip = false;

        if(e.getNode1().getName().equals(name2))
        {
            if(!e.getNode2().getName().equals(name1))
            {
                System.err.println("Inconsistent names in edgeToInt...");
                System.exit(-1);
            }
            flip = true;
        }

        if(e.getEndpoint1()==Endpoint.ARROW)
        {
            if(e.getEndpoint2()==Endpoint.TAIL)
            {
                return 2;// <-- this will never happen
            }
            else if(e.getEndpoint2()==Endpoint.ARROW)
            {
                return 3; //<->
            }
            else
            {
                return 5; //<-o this will never happen
            }
        }
        else if(e.getEndpoint1()==Endpoint.TAIL)
        {
            if(e.getEndpoint2()==Endpoint.ARROW) {
                if(flip)
                    return 2;
                else
                    return 1;//-->
            }
            else if(e.getEndpoint2()==Endpoint.TAIL)
                return 0; //---
        }
        else
        {
            if(e.getEndpoint2()==Endpoint.CIRCLE)
                return 6; //o-o
            else if(e.getEndpoint2()==Endpoint.ARROW) {
                if(flip)
                    return 5;
                else
                    return 4;//o->
            }
        }
        return -1;
    }

}