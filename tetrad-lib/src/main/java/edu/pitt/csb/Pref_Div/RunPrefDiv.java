package edu.pitt.csb.Pref_Div;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.*;
import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.PrintStream;
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
    private int numAlphas = 20; //How many alpha values should we test? (Data-Theory tradeoff alpha)
    private double alphaLimit = 1; //The largest value of alpha we test
    private double lambdaLow = 0.05; //Lambda low limit
    private double lambdaHigh = 0.95; //Lambda high limit for MGM StEPS
    private double g = 0.05; //Stability threshold
    private ArrayList<Gene> lastGeneSet;
    private HashMap<Gene,List<Gene>> lastClusters;
    private String target; //Target variable
    private boolean approxCorrelations = false; //Compute some correlations and do the rest on the fly?
    private boolean useClusterStability = false; //cluster stability, or gene wise stability
    private boolean partialCorr = false; //Use partial correlations instead of correlation for continuous variables
    private double [][] lastStepsStabilities; //Last set of gene stabilities given by StEPS

    private ArrayList<String> stabPS;//PrintStream to output file where we have run\talpha\tstability
    private boolean includeAccuracy = false; //When writing the stability output to file, should we include an accuracy score for each stability data point
    private List<Node> targets = null; //The set of target genes (for computing accuracy score when ground truth is known)
    private int run; //purely for printing to the stability output file
    private boolean loocv = false;
    private int [][] subs;
    private boolean getCausalAccuracy = false;
    private Graph trueGraph;

    private boolean useStabilitySelection = false;

    public RunPrefDiv(float [] dissimilarity, ArrayList<Gene> genes, DataSet data,String target,boolean leaveOneOut)
    {
        this.data = data;
        this.genes = genes;
        this.dissimilarity = dissimilarity;
        this.target = target;
        this.loocv = leaveOneOut;
        this.stabPS = null;
    }

    public int [][] getSubs()
    {
        return subs;
    }

    public double[][]getLastStabilities(){return lastStepsStabilities;}
    public void setSubs(int [][] subs){this.subs = subs;}
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
    public void setStabilityOutput(boolean includeAccuracy){this.stabPS = new ArrayList<String>(); this.includeAccuracy=includeAccuracy;}
    public ArrayList<String> getStabilityOutput(){return stabPS;}
    public void setTargets(List<Node> targets){this.targets = targets;}
    public void setTrueGraph(Graph truth){this.trueGraph = truth;}
    public void setRun(int r){run = r;}
    public ArrayList<Gene> getLastGeneSet(){return lastGeneSet;}
    public void usePartialCorrelation(boolean pc){partialCorr = pc;}
    public void useStabilitySelection(){useStabilitySelection = true;}

    public ArrayList<Gene> runPD()
    {
        //Loop over alpha values
        //Call helper method to compute stability and result set for given alpha
        //If result is below threshold, then keep it
        //Otherwise increment by alphaStep and continue searching


        System.out.print("Generating Subsamples...");
        if(subs==null) {
            if (loocv)
                subs = StabilityUtils.generateSubsamples(data.getNumRows());
            else
                subs = StabilityUtils.generateSubsamples(numSubs, data.getNumRows());
        }
        System.out.println("Done");
        double [] alphas = new double[numAlphas+1];
        for(int i = 0; i < alphas.length;i++)
        {
            alphas[i] = i*alphaLimit/(double)numAlphas;
        }



        ArrayList<Gene> finalSet = null;


        for(int i = alphas.length-1; i >=0;i--)
        {
            if(useStabilitySelection)
            {

            }

            System.out.print("Running Pref-Div for Alpha value: " + alphas[i] + " ...");
            double stab = stabilityPD(alphas[i]);
            if(stabPS!=null) {
                if(includeAccuracy) {
                    if(getCausalAccuracy) {
                        Graph graph = learnGraph();
                        stabPS.add(run + "\t" + alphas[i] + "\t" + stab + "\t" + getAccuracy(targets, lastGeneSet)[1] + "\t" + getAccuracy(targets, lastGeneSet)[2] + "\t" + lastGeneSet + "\t" + targets + "\t" + lastClusters + "\t" + graphAccuracy(graph)[0] + "\t" + graphAccuracy(graph)[1] + "\t" + simulatePrefDivMGM.neighborAccuracy(graph,trueGraph,target)[0] + "\t" + simulatePrefDivMGM.neighborAccuracy(graph,trueGraph,target)[1]);
                    }else
                        stabPS.add(run + "\t" + alphas[i] + "\t" + stab + "\t" + getAccuracy(targets, lastGeneSet)[1] + "\t" + getAccuracy(targets, lastGeneSet)[2] + "\t" + lastGeneSet + "\t" + targets + "\t" + lastClusters);

                }else
                    stabPS.add(run + "\t" + alphas[i] + "\t" + stab + "\t" + lastClusters);
            }
            System.out.println("Stability = " + (stab) + ", Done");
            if(1-stab < g && stabPS==null) { //Continue to search alpha values if you want to output stabilities for them
                System.out.println("Underneath threshold of " + g + ", Returning this set");
                return lastGeneSet;
            }
            else if(1-stab < g && finalSet==null)
            {
                finalSet = lastGeneSet;
            }
        }


        if(finalSet!= null) {
            lastGeneSet = finalSet;
            return finalSet;
        }
        return lastGeneSet;

    }


    public double [] graphAccuracy(Graph gOut)
    {

        double tp = 0;
        double fp = 0;
        double fn = 0;
        for(Edge e:gOut.getEdges())
        {
            if(e.getNode1().getName().equals("Dummy")||e.getNode2().getName().equals("Dummy"))
                continue;
            if(trueGraph.getEdge(trueGraph.getNode(e.getNode1().getName()),trueGraph.getNode(e.getNode2().getName()))!=null)
                tp++;
            else if(trueGraph.getEdge(trueGraph.getNode(e.getNode2().getName()),trueGraph.getNode(e.getNode1().getName()))!=null)
                tp++;
            else
                fn++;

        }
        for(Edge e: trueGraph.getEdges())
        {
            if(gOut.getNode(e.getNode1().getName())==null || gOut.getNode(e.getNode2().getName())==null)
                continue;
            if(gOut.getEdge(gOut.getNode(e.getNode1().getName()),gOut.getNode(e.getNode2().getName()))==null && gOut.getEdge(gOut.getNode(e.getNode2().getName()),gOut.getNode(e.getNode1().getName()))==null)
                fp++;
        }
        return new double[]{tp/(tp+fp),tp/(tp+fn)};
    }
    public Graph learnGraph()
    {
        ArrayList<Node> names = new ArrayList<Node>();
        for(Gene x:lastGeneSet)
        {
            names.add(data.getVariable(x.symbol));
        }
        if(target!=null && !target.equals(""))
            names.add(data.getVariable(target));

        DataSet temp = data.subsetColumns(names);
        if(!temp.isMixed())
        {
            temp.addVariable(new DiscreteVariable("Dummy"));
            int col = temp.getColumn(temp.getVariable("Dummy"));
            Random rand = new Random();
            for(int i = 0; i < temp.getNumRows();i++)
            {
                temp.setInt(i,col,rand.nextInt(2));
            }
        }
        double [] lambdas = new double[numAlphas];
        for(int i = 0; i < numAlphas;i++)
        {
            lambdas[i] = ((i+1)*(lambdaHigh-lambdaLow)/numAlphas);
        }
        DataSet [] subsamples = new DataSet[subs.length];
        for(int i = 0; i < subs.length;i++)
        {
            subsamples[i] = temp.subsetRows(subs[i]);
        }
        STEPS s = new STEPS(temp,lambdas,g,subsamples);
        Graph output = s.runStepsPar();
        return output;
    }
    public Graph getCausalGraph(String target)
    {
        getCausalAccuracy = true;
        ArrayList<Gene> g2 = runPD();
        ArrayList<Node> names = new ArrayList<Node>();
        for(Gene x:g2)
        {
            names.add(data.getVariable(x.symbol));
        }
        if(target!=null && !target.equals(""))
             names.add(data.getVariable(target));

        DataSet temp = data.subsetColumns(names);
        if(!temp.isMixed())
        {
            temp.addVariable(new DiscreteVariable("Dummy"));
            int col = temp.getColumn(temp.getVariable("Dummy"));
            Random rand = new Random();
            for(int i = 0; i < temp.getNumRows();i++)
            {
                temp.setInt(i,col,rand.nextInt(2));
            }
        }
        double [] lambdas = new double[numAlphas];
        for(int i = 0; i < numAlphas;i++)
        {
            lambdas[i] = ((i+1)*(lambdaHigh-lambdaLow)/numAlphas);
        }
        DataSet [] subsamples = new DataSet[subs.length];
        for(int i = 0; i < subs.length;i++)
        {
            subsamples[i] = temp.subsetRows(subs[i]);
        }
        STEPS s = new STEPS(temp,lambdas,g,subsamples);
        Graph output = s.runStepsPar();
        lastGeneSet = g2;
        lastStepsStabilities = s.stabilities;
        return output;
    }
    private double [] getAccuracy(List<Node> targets, ArrayList<Gene> result)
    {
        double [] output = new double[3];
        double tp = 0;
        double fp = 0;
        double fn = 0;
        for(int i = 0; i < targets.size();i++)
        {
            boolean found = false;
            for(int j = 0; j < result.size();j++)
            {
                if(targets.get(i).getName().equals(result.get(j).symbol))
                    found = true;
            }
            if(found)
                tp++;
            else
                fn++;
        }
        for(int i = 0; i < result.size();i++)
        {
            boolean found = false;
            for(int j = 0; j < targets.size();j++)
            {
                if(targets.get(j).getName().equals(result.get(i).symbol))
                    found = true;
            }
            if(!found)
                fp++;
        }
        output[1] = tp/(tp+fp);
        output[2] = tp/(tp+fn);
        output[0] = 2*output[1]*output[2]/(output[1]+output[2]);
        return output;
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

                        ArrayList<Gene> curr = Functions.computeAllIntensities(genes,alp,dataSubSamp,target,partialCorr);
                        long time = System.nanoTime();
                        //PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alp,dataSubSamp,true);
                        //p.setCluster(true);
                        //ArrayList<Gene> result = p.diverset();
                        time = System.nanoTime()-time;
                       // System.out.println(s + " Computing all correlations " + time/Math.pow(10,9));

                        time = System.nanoTime();
                        Collections.sort(curr,Gene.IntensityComparator);
                        PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alp,dataSubSamp,approxCorrelations,partialCorr);
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
        ArrayList<Gene> curr = Functions.computeAllIntensities(genes,alp,data,target,partialCorr);
        Collections.sort(curr,Gene.IntensityComparator);
        PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alp,data,approxCorrelations,partialCorr);
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
