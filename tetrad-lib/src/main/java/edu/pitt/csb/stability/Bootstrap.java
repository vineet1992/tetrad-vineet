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
import edu.pitt.csb.Priors.runPriors;
import edu.pitt.csb.latents.LatentPrediction;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.stability.DataGraphSearch;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

/**
 * Created by vinee_000 on 7/12/2018.
 * Class to run Complimentary Pairs Stability Selection with causal modeling algorithms
 *
 */
public class Bootstrap {
    private int B; //Number of bootstrap samples
    private int p;//total number of variables
    private double bound; //percentage of bootstrap samples cutoff
    private DataSet data;
    private double [] lambda; //Array of lambda values for MGM
    private int [][] subs;


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

    public static int[][] createSubs(DataSet data, int B)
    {
        return DataUtils.getBootstrapIndices(data,data.getNumRows(),B);
    }
    private int[][] createSubs(DataSet data)
    {
        subs = createSubs(data,B);
        return subs;
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
    public double [][] addToMat(Graph g, DataSet data, double [][] edgeCounts)
    {
        for(Edge e:g.getEdges())
        {
            int x = data.getColumn(data.getVariable(e.getNode1().getName()));
            int y = data.getColumn(data.getVariable(e.getNode2().getName()));
            edgeCounts[x][y]++;
            edgeCounts[y][x]++;
        }
        return edgeCounts;

    }
}