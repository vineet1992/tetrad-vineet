package edu.pitt.csb.Pref_Div;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.*;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.StabilityUtils;

import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

import static org.apache.commons.collections4.CollectionUtils.intersection;
import static org.apache.commons.collections4.CollectionUtils.union;

/**
  Runs genomic Pref-Div algorithm with subsampling and clustering
    TODO Handle Discrete features in theory relevance and dissimilarity
 */
public class RunPrefDiv {

    private DataSet data;
    private ArrayList<Gene> genes;
    private float[] dissimilarity;
    private int numSubs = 20;
    private int topK = 50;
    private double accuracy = 0; //TODO what's the right value for this to get good clusters? 0?
    private double radius = 0.5; //TODO decide if we should use the radius or the max cluster size
    private int numAlphas = 20;
    private double alphaLimit = 1;
    private int [] diseaseIDs;
    private int maxIntDiscrete = 5;
    private double g = 0.05; //Stability threshold
    private ArrayList<Gene> lastGeneSet;
    private HashMap<Gene,List<Gene>> lastClusters;
    private String target;
    private boolean approxCorrelations = false;
    private boolean useClusterStability = false; //cluster stability, or gene wise stability

    private boolean loocv = false;
    private int [][] subs;



    public RunPrefDiv(float [] dissimilarity, ArrayList<Gene> genes, DataSet data,String target,boolean leaveOneOut)
    {
        this.data = data;
        this.genes = genes;
        this.dissimilarity = dissimilarity;
        this.target = target;
        this.loocv = leaveOneOut;
    }

    public HashMap<Gene,List<Gene>> getClusters()
    {
        return lastClusters;
    }
    public void setThreshold(double t)
    {
        g = t;
    }
    public void setApproxCorrelations(boolean a)
    {
        approxCorrelations = a;
    }
    public void setClusterStability(boolean x)
    {
        useClusterStability = x;
    }
    public void setMaxInt(int m)
    {
        this.maxIntDiscrete = m;
    }
    public void setAlphaLimit(double lim)
    {
        alphaLimit = lim;
    }
    public void setNS(int ns)
    {
        this.numSubs = ns;
    }
    public void setTopK(int k)
    {
        this.topK = k;
    }
    public void setAccuracy(double A)
    {
        this.accuracy = A;
    }
    public void setRadius(double r)
    {
        this.radius = r;
    }
    public void setNumAlphas(int na)
    {
        this.numAlphas = na;
    }
    public void setDiseases(int [] IDs)
    {
        this.diseaseIDs = IDs;
    }

    public ArrayList<Gene> runPD()
    {
        //Loop over alpha values
        //Call helper method to compute stability and result set for given alpha
        //If result is below threshold, then keep it
        //Otherwise increment by alphaStep and continue searching

        System.out.print("Generating Subsamples...");
        if(loocv)
            subs = StabilityUtils.generateSubsamples(data.getNumRows());
        else
            subs = StabilityUtils.generateSubsamples(numSubs,data.getNumRows());
        System.out.println("Done");
        double [] alphas = new double[numAlphas+1];
        for(int i = 0; i < alphas.length;i++)
        {
            alphas[i] = i*alphaLimit/(double)numAlphas;
        }

        for(int i = alphas.length-1; i >=0;i--)
        {
            System.out.print("Running Pref-Div for Alpha value: " + alphas[i] + " ...");
            double stab = stabilityPD(alphas[i]);

            System.out.println("Stability = " + (stab) + ", Done");
            if(1-stab < g) {
                System.out.println("Underneath threshold of " + g + ", Returning this set");
                return lastGeneSet;
            }
        }


        return lastGeneSet;

    }

    private double stabilityPD(final double alp)
    {

        final int numVars = data.getNumColumns();
        final List<List<Gene>> allGenes = new ArrayList<List<Gene>>();

        final List<HashMap<Gene,List<Gene>>> allClusters = new ArrayList<HashMap<Gene,List<Gene>>>();
        final ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

        class StabilityAction extends RecursiveAction {
            private int chunk;
            private int from;
            private int to;

            public StabilityAction(int chunk, int from, int to){
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }

            private synchronized void addToFullSet(ArrayList<Gene> x,int s){

                if(s < allGenes.size())
                    allGenes.add(s,x);
                else
                    allGenes.add(x);
            }

            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {

                        DataSet dataSubSamp = data.subsetRows(subs[s]);

                        ArrayList<Gene> curr = Functions.computeAllIntensities(genes,alp,dataSubSamp,target);
                        long time = System.nanoTime();
                        //PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alp,dataSubSamp,true);
                        //p.setCluster(true);
                        //ArrayList<Gene> result = p.diverset();
                        time = System.nanoTime()-time;
                       // System.out.println(s + " Computing all correlations " + time/Math.pow(10,9));

                        time = System.nanoTime();
                        PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alp,dataSubSamp,approxCorrelations);
                        p.setCluster(true);
                        ArrayList<Gene> result = p.diverset();
                        allClusters.add(p.clusters);
                        time = System.nanoTime()-time;
                       // System.out.println(s + " Computing on the fly " + time/Math.pow(10,9));


                       // ArrayList<Gene> result = p.diverset();
                        addToFullSet(result,s);
                        time = System.nanoTime()-time;
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

        final int chunk = 1000;

        pool.invoke(new StabilityAction(chunk, 0, subs.length));
        ArrayList<Gene> curr = Functions.computeAllIntensities(genes,alp,data,target);
        PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alp,data,approxCorrelations);
        p.setCluster(true);
        lastGeneSet = p.diverset();
        lastClusters = p.clusters;
        //do this elsewhere
        //thetaMat.assign(thetaMat.copy().assign(Functions.minus(1.0)), Functions.mult).assign(Functions.mult(-2.0));
        if(useClusterStability)
            return 1-clusterSim(allClusters);
        else
            return tanimotoSim(allGenes);
    }


    private void addToGeneList(HashMap<Gene,List<Gene>> map, List<List<Gene>> one)
    {
        for(Gene g: map.keySet())
        {
            ArrayList<Gene> temp = new ArrayList<Gene>();
            temp.add(g);
            if(map.get(g)!=null) {
                for (Gene x : map.get(g)) {
                    temp.add(x);
                }
            }
            one.add(temp);
        }
    }

    private double[][] computeCostMatrix(List<List<Gene>> one, List<List<Gene>> two)
    {
        double [][] costs = new double[one.size()][two.size()];
        for(int i = 0; i < one.size();i++)
        {
            for(int j = 0; j < two.size();j++)
            {
                costs[i][j] = 1-(intersection(one.get(i),two.get(j)).size()/(double)union(one.get(i),two.get(j)).size());
            }
        }
        return costs;
    }
    private double [][] normalizedCosts(double[][]costs)
    {
        double sum = 0;
        for(int i = 0; i < costs.length;i++)
        {
            for(int j = 0; j < costs[i].length;j++)
            {
                sum+= costs[i][j];
            }
        }

        for(int i = 0; i < costs.length;i++)
        {
            for(int j = 0; j < costs[i].length;j++)
            {
                costs[i][j]/= sum;
            }
        }

        return costs;

    }
    private double clusterSim(List<HashMap<Gene,List<Gene>>> maps)
    {
        double sim = 0;
        int count = 0;
        for(int i = 0; i < maps.size();i++)
        {
            for(int j = i + 1; j < maps.size();j++)
            {
                List<List<Gene>> one = new ArrayList<List<Gene>>();
                HashMap<Gene,List<Gene>> map = maps.get(i);
                addToGeneList(map,one);

                List<List<Gene>> two = new ArrayList<List<Gene>>();
                map = maps.get(j);
                addToGeneList(map,two);

                double [][] costs = computeCostMatrix(one,two);
                HungarianAlgorithm ha = new HungarianAlgorithm(costs);

                int [] result = ha.execute();

                costs = normalizedCosts(costs);
                for(int k = 0; k < result.length;k++)
                {
                    sim+=costs[k][result[k]];
                }
                count++;
            }
        }

        return sim/count;
    }

    private double tanimotoSim(List<List<Gene>> allGenes)
    {
        double sim = 0;
        int count = 0;
        for(int i = 0; i < allGenes.size();i++)
        {
            for(int j = i+1; j < allGenes.size();j++)
            {
                sim+=(intersection(allGenes.get(i),allGenes.get(j)).size()/(double)union(allGenes.get(i),allGenes.get(j)).size());
                count++;
            }
        }
        return sim/count;
    }

}
