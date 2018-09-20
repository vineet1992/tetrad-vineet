package edu.pitt.csb.latents;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.Fci;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.pitt.csb.Pref_Div.Functions;
import edu.pitt.csb.Pref_Div.Gene;
import edu.pitt.csb.Priors.runPriors;
import edu.pitt.csb.latents.mirnaPrediction;
import edu.pitt.csb.mgm.*;
import edu.pitt.csb.stability.Bootstrap;
import edu.pitt.csb.stability.CPSS;
import edu.pitt.csb.stability.SearchWrappers;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.ThreadLocalRandom;

/**
 * Created by vinee_000 on 10/9/2017.
 * Runs modifications of FCI-flavor approaches to determine the variables likely to be confounded by latents
 */

//TODO Add in complementary pairs stability selection and compare the two approaches before finalizing a method for simulated data experiments
public class LatentPrediction {


    public Map<String,String> orientations;
    private DataSet data;
    private int nss;
    private double tao;
    private int ns;
    private int [][] subs;
    private int [][] lastSplit;
    private static double [] searchParams;
    private static double [] initAlphas;
    private static double [] initLambdas;
    private int numLambdas = 30; //Number of lambda values to test
    private int numAlphas = 30; //Number of alpha values to test
    private double lambdaLow = 0.05; //Range of lambda values to test
    private double lambdaHigh = 0.9; //Range of lambda values to test
    private double alphaLow = 0.001;//Range of alpha values to test
    private double alphaHigh = 0.5; // Range of alpha values to test
    private double stepsG = 0.05; //Stability threshold for StEPS
    private double starsG = 0.05; //Stability threshold for StARS
    private double useAlpha = -1;
    public Graph trueGraph;
    private int runNumber = 0;
    private boolean runCPSS; //To run or not to run CPSS?
    private boolean runBoot; //To run or not to run Bootstrapping?
    private double cpssAlpha = 0.05; //Alpha value for CPSS
    private double [] cpssLambda = {0.2,0.2,0.2}; //Lambda value for CPSS
    private int B = 50; //Number of Complimentary Pairs to generate
    private double bound; //CPSS bound
    private int [][] cpSubs; //Subsamples (Complimentary pairs), generated by CPSS
    public LatentPrediction(DataSet d,int numSubSets, double tao, int numSubsamples)
    {
        data = d;
        nss = numSubSets;
        this.tao = tao;
        this.ns = numSubsamples;
        searchParams = new double[4];
        initializeArrays();
        orientations = new HashMap<String,String>();
        this.subs = StabilityUtils.subSampleNoReplacement(data.getNumRows(),numSubsamples);
        while(!ensureVariance())
            subs = StabilityUtils.subSampleNoReplacement(data.getNumRows(),numSubsamples);

    }
    public LatentPrediction(DataSet d, int numSubSets, double tao, int [][] subs)
    {
        data = d;
        nss = numSubSets;
        this.tao = tao;
        this.subs = subs;
        this.ns = subs.length;
        searchParams = new double[4];
        orientations = new HashMap<String,String>();
        initializeArrays();

    }


    public void setCPSS(int B, double alpha, double [] lambda, double bound, int [][] cpSubs){this.cpSubs = cpSubs; runCPSS = true; this.bound = bound; this.B = B; this.cpssAlpha = alpha; this.cpssLambda = lambda; }
    public void setBootstrap(int B, double alpha, double [] lambda, double bound, int [][] bootSubs){this.cpSubs = bootSubs; runBoot = true; this.bound = bound; this.B = B; this.cpssAlpha = alpha; this.cpssLambda = lambda; }


    public void setStarsGamma(double starsGamma)
    {
        starsG = starsGamma;
    }
    public void setStepsGamma(double stepsGamma)
    {
        stepsG = stepsGamma;
    }

    public void setGraph(Graph g)
    {
        trueGraph = g;
    }
    public void setRunNumber(int i)
    {
        runNumber = i;
    }
    public int [][] getLastSplit()
    {
        return lastSplit;
    }
    private void initializeArrays()
    {
        initLambdas = new double[numLambdas];
        for (int i = 0; i < numLambdas; i++) {
            initLambdas[i] = i * (lambdaHigh - lambdaLow) / numLambdas + lambdaLow;
        }
        initAlphas = new double[numAlphas];
        for (int i = 0; i < numAlphas; i++) {
            initAlphas[i] = i * (alphaHigh - alphaLow) / numAlphas + alphaLow;
        }
    }
    //For FCI it's just [alpha]
    //For MGM flavors it's [labmda_1, lambda_2, lambda_3, alpha]
    public void setParams(double [] params)
    {
        searchParams = params;
    }


    public void setAlpha(double alp)
    {
        this.useAlpha = alp;
    }



    //Runs an algorithm through the full procedure -> gets the optimal parameter based on StaRS/StEPS
    //Runs the algorithm with this parameter across subsamples
    //Returns a hashmap of why edges were oriented (which FCI rule), and stability of each edge
    //Returns a list of edges that meet stability requirements
    //Algorithm options: "FCI","MGM-FCI","MGM-FCI-MAX"
    //Either generates subsamples or uses the subsamples present in the class already
    //Boolean says whether or not the edge must be in the graph generated by running on all of the samples in order to be in the output
    public ArrayList<Pair> runRegularAlgorithm(String algName, boolean requireFullEdge) throws Exception
    {
        return runRegularAlgorithm(algName,data,requireFullEdge);
    }

    //If MGM algorithm is used, then default lambdas are employed [0.2,0.2,0.2]
    private double getStarsAlpha(DataSet d, Algorithm a)
    {
        if(useAlpha > 0)
            return useAlpha;
        STARS s = new STARS(d,initAlphas,starsG,ns,a);
        s.setTrueGraph(trueGraph);
        s.setRun(runNumber);
        double [][] instab = s.runSTARS();
        for(int i = 0; i < instab.length;i++)
        {
            System.out.println(Arrays.toString(instab[i]));
        }
        return s.lastAlpha;
    }

    //For use with MGM algorithms, assuming that STEPS was already run to produce lambdas
    private double getStarsAlpha(DataSet d, Algorithm a, double [] lambda)
    {
        if(useAlpha > 0)
            return useAlpha;

        STARS s = new STARS(d,initAlphas,starsG,subs,a);
        s.setMGMLambda(lambda);
        s.setTrueGraph(trueGraph);
        s.setRun(runNumber);
        double [][] instab = s.runSTARS();
        for(int i = 0; i < instab.length;i++)
        {
            System.out.println(Arrays.toString(instab[i]));
        }
        return s.lastAlpha;
    }


    private double[] getStepsLambda(DataSet d)
    {
        STEPS s = new STEPS(d,initLambdas,stepsG,subs);
        s.runStepsPar();
        return s.lastLambda;
    }



    private synchronized void addToOrient(ConcurrentHashMap<String,String> tempOr)
    {
        for(String x:tempOr.keySet())
        {
            if(orientations.get(x)!=null)
            {
                orientations.put(x,orientations.get(x)+"\tFINAL," + tempOr.get(x));
            }
        }
    }

    //Private version of run regular algorithm, used for both the runAlgorithm subroutine, and the regular runRegularAlgorithm method
    private ArrayList<Pair> runRegularAlgorithm(String algName, DataSet d, boolean requireFullEdge) throws Exception
    {
        ArrayList<Pair> result = new ArrayList<Pair>();
        int b = (int)(10*Math.sqrt(d.getNumRows()));
        DoubleMatrix2D stabs;
        Graph gTotal = null;
        if(runCPSS)
        {
            CPSS cp = new CPSS(data,cpssLambda,bound,cpSubs);
            return cp.runCPSSLatent(algName,cpssAlpha);
        }
        if(runBoot)
        {
            Bootstrap bs = new Bootstrap(data,cpssLambda,bound,cpSubs);
            return bs.runBootstrap(algName,cpssAlpha);
        }
        else {
            if (algName.equals("FCI")) {
                searchParams[0] = getStarsAlpha(d, Algorithm.FCI);
                stabs = StabilityUtils.StabilitySearchParLatent(d, new SearchWrappers.FCIWrapper(searchParams), subs, orientations);
                IndependenceTest i = new IndTestMultinomialAJ(d, searchParams[0],true);
                Fci f = new Fci(i);
                gTotal = f.search();
                ConcurrentHashMap<String, String> tempOr = f.whyOrient;
                addToOrient(tempOr);

            } else if (algName.equals("MGM-FCI")) {
                double[] temp = getStepsLambda(d);
                searchParams[0] = temp[0];
                searchParams[1] = temp[1];
                searchParams[2] = temp[2];
                searchParams[3] = getStarsAlpha(d, Algorithm.MGMFCI, temp);
                stabs = StabilityUtils.StabilitySearchParLatent(d, new SearchWrappers.MGMFCIWrapper(searchParams), subs, orientations);
                MGM m = new MGM(data, temp);
                m.learnEdges(1000);
                IndependenceTest i = new IndTestMultinomialAJ(data, searchParams[3],true);
                Fci f = new Fci(i);
                f.setInitialGraph(m.graphFromMGM());
                gTotal = f.search();
                ConcurrentHashMap<String, String> tempOr = f.whyOrient;
                addToOrient(tempOr);
            } else if (algName.equals("MGM-FCI-MAX")) {
                System.out.println("Computing Optimal Lambda");
                double[] temp = getStepsLambda(d);
                System.out.println("Lambda: " + Arrays.toString(temp));

                searchParams[0] = temp[0];
                searchParams[1] = temp[1];
                searchParams[2] = temp[2];
                System.out.println("Computing Optimal Alpha...");
                searchParams[3] = getStarsAlpha(d, Algorithm.MGMFCIMAX, temp);
                System.out.println("Alpha: " + searchParams[3]);
                System.out.println("Computing Stable Edges...");
                stabs = StabilityUtils.StabilitySearchParLatent(d, new SearchWrappers.MFMWrapper(searchParams), subs, orientations);
                MGM m = new MGM(d, temp);
                m.learnEdges(1000);
                IndependenceTest i = new IndTestMultinomialAJ(d, searchParams[3],true);
                FciMaxP f = new FciMaxP(i);
                f.setInitialGraph(m.graphFromMGM());
                gTotal = f.search();
                ConcurrentHashMap<String, String> tempOr = (ConcurrentHashMap<String, String>) f.whyOrient;
                addToOrient(tempOr);

            } else {
                throw new Exception("Unrecognized Algorithm");
            }
            System.out.println("Graph for latent algorithm: " + gTotal);
            for (int i = 0; i < stabs.rows(); i++) {
                for (int j = i + 1; j < stabs.columns(); j++) {
                    if (d.getVariable(i) instanceof DiscreteVariable || d.getVariable(j) instanceof DiscreteVariable)
                        continue;
                    Edge e = gTotal.getEdge(gTotal.getNode(d.getVariable(i).getName()), gTotal.getNode(d.getVariable(j).getName()));
                    if (e != null && e.getEndpoint1() == Endpoint.ARROW && e.getEndpoint2() == Endpoint.ARROW) {
                        if (stabs.get(i, j) > tao) {
                            result.add(new Pair(d.getVariable(i), d.getVariable(j), stabs.get(i, j)));
                        }
                    } else if (!requireFullEdge && stabs.get(i, j) > tao) {
                        result.add(new Pair(d.getVariable(i), d.getVariable(j), stabs.get(i, j)));
                    }
                }
            }
        }
        System.out.println("List of latents for latent algorithm: " + result);
        return result;
    }


    private double getAlpha(final String algName)
    {
        if(algName.equals("FCI"))
            return searchParams[0];
        else if (algName.equals("MGM-FCI"))
            return searchParams[3];
        else if(algName.equals("MGM-FCI-MAX"))
            return searchParams[3];
        else
        {
            return -1;
        }
    }
    public static double getPrecision(ArrayList<Pair> latents, Graph trueGraph,DataSet dat,  String type)
    {
        ArrayList<Pair> trueLatents = getLatents(trueGraph,dat,type);
        System.out.println("True Latents: " + trueLatents);
        double tp = 0;
        double fp = 0;
        for(Pair p:latents)
        {
            String type2= "";
            if(dat.getVariable(p.one.getName()) instanceof ContinuousVariable && dat.getVariable(p.two.getName())instanceof ContinuousVariable)
                type2="CC";
            else if(dat.getVariable(p.one.getName())instanceof DiscreteVariable && dat.getVariable(p.two.getName())instanceof DiscreteVariable)
                type2="DD";
            else
                type2="CD";
            if(!type2.equals(type) && !type.equals("All"))
                continue;
            boolean found = false;
            for(Pair x: trueLatents)
            {
                if(x.one.getName().equals(p.one.getName()) && x.two.getName().equals(p.two.getName()))
                    found = true;
                if(x.two.getName().equals(p.one.getName()) && x.one.getName().equals(p.two.getName()))
                    found = true;
            }
            if(!found)
                fp++;
            else
                tp++;
        }
        return tp/(tp+fp);
    }
    public static ArrayList<Pair> getLatents(Graph trueGraph, DataSet dat, String type)
    {
        HashMap<String,String> already = new HashMap<String,String>();
        ArrayList<Pair> result = new ArrayList<Pair>();
        for(Node n: trueGraph.getNodes())
        {
            if(n.getNodeType()==NodeType.LATENT)
            {
                List<Node> children = trueGraph.getChildren(n);
                for (int i = 0; i < children.size(); i++) {
                    for (int j = i + 1; j < children.size(); j++) {
                        Node one = children.get(i);
                        Node two = children.get(j);
                        if(trueGraph.getNode(one.getName()).getNodeType()==NodeType.LATENT)
                            continue;
                        if(trueGraph.getNode(two.getName()).getNodeType()==NodeType.LATENT)
                            continue;
                        String type2 = "";
                        if (dat.getVariable(one.getName()) instanceof ContinuousVariable && dat.getVariable(two.getName()) instanceof ContinuousVariable)
                            type2 = "CC";
                        else if (dat.getVariable(one.getName()) instanceof DiscreteVariable && dat.getVariable(two.getName()) instanceof DiscreteVariable)
                            type2 = "DD";
                        else
                            type2 = "CD";
                        if (!type2.equals(type) && !type.equals("All"))
                            continue;
                        if(already.get(one.getName())==null || already.get(one.getName())!=two.getName())
                        {
                            if(already.get(two.getName())==null||already.get(two.getName())!=one.getName())
                            {
                                result.add(new Pair(one,two));
                                already.put(one.getName(),two.getName());
                            }
                        }
                    }
                }
            }
        }
        return result;
    }
    public static ArrayList<Pair> getLatents(Graph trueGraph, DataSet dat, String type, int [][] splits)
    {
        HashMap<String,String> already = new HashMap<String,String>();
        ArrayList<Pair> result = new ArrayList<Pair>();
        for(Node n: trueGraph.getNodes())
        {
            if(n.getNodeType()==NodeType.LATENT)
            {
                List<Node> children = trueGraph.getChildren(n);
                for (int i = 0; i < children.size(); i++) {
                    for (int j = i + 1; j < children.size(); j++) {
                        Node one = children.get(i);
                        Node two = children.get(j);
                        if(one.getNodeType()==NodeType.LATENT || two.getNodeType()==NodeType.LATENT)
                            continue;
                        String type2 = "";
                        if (dat.getVariable(one.getName()) instanceof ContinuousVariable && dat.getVariable(two.getName()) instanceof ContinuousVariable)
                            type2 = "CC";
                        else if (dat.getVariable(one.getName()) instanceof DiscreteVariable && dat.getVariable(two.getName()) instanceof DiscreteVariable)
                            type2 = "DD";
                        else
                            type2 = "CD";
                        if (!type2.equals(type) && !type.equals("All"))
                            continue;
                        if(already.get(one.getName())==null || already.get(one.getName())!=two.getName())
                        {
                            if(already.get(two.getName())==null||already.get(two.getName())!=one.getName())
                            {
                                int indexOne = dat.getColumn(dat.getVariable(one.getName()));
                                int indexTwo = dat.getColumn(dat.getVariable(two.getName()));
                                boolean found = false;
                                A:for(int k = 0; k < splits.length;k++)
                                {
                                    int score = 0;
                                    for(int m = 0; m < splits[k].length;m++)
                                    {
                                        if(splits[k][m] == indexOne)
                                            score+=1;
                                        if(splits[k][m]==indexTwo)
                                            score+=2;
                                    }
                                    if(score==3)
                                        found=true;
                                    else if(score!=0)
                                        break A;
                                }
                                if(found) {
                                    result.add(new Pair(one, two));
                                    already.put(one.getName(), two.getName());
                                }

                            }
                        }
                    }
                }
            }
        }
        return result;
    }

    public static double getRecall(ArrayList<Pair> latents, Graph trueGraph, DataSet dat, String type)
    {
        ArrayList<Pair> trueLatents = getLatents(trueGraph,dat,type);
       double tp = 0;
        double fn = 0;
        for(int i = 0; i < trueLatents.size();i++)
        {
            Pair currLatent = trueLatents.get(i);
            boolean found = false;
            for(int j = 0; j < latents.size();j++)
            {
                Pair currGuess = latents.get(j);
                if(currGuess.one.getName().equals(currLatent.one.getName()) && currGuess.two.getName().equals(currLatent.two.getName()))
                    found = true;
                if(currGuess.two.getName().equals(currLatent.one.getName())&& currGuess.one.getName().equals(currLatent.two.getName()))
                    found = true;
            }
            if(found)
                tp++;
            else
                fn++;
        }

        return tp/(tp+fn);
    }


    public static boolean isIdentifiable(Graph truePag, Node one, Node two)
    {
        Edge e = truePag.getEdge(truePag.getNode(one.getName()),truePag.getNode(two.getName()));
        if(e==null || e.getEndpoint1()!=Endpoint.ARROW || e.getEndpoint2()!=Endpoint.ARROW)
            return false;
        return true;
    }

    public static double getIdentifiableRecall(Graph PAG, ArrayList<Pair> latents, Graph trueDAG, DataSet d, String type, int [][] splits)
    {
        ArrayList<Pair> trueLatents = getLatents(trueDAG,d,type,splits);
        ArrayList<Pair> temp = new ArrayList<Pair>();
        for(Pair p:trueLatents)
        {
            Edge e = PAG.getEdge(PAG.getNode(p.one.getName()),PAG.getNode(p.two.getName()));
            if(e==null)
                continue;
            if(e.getEndpoint1()== Endpoint.ARROW && e.getEndpoint2()==Endpoint.ARROW)
                temp.add(p);
        }
        trueLatents = temp;
        temp = null;
        double tp = 0;
        double fn = 0;
        for(int i = 0; i < trueLatents.size();i++)
        {
            Pair currLatent = trueLatents.get(i);
            boolean found = false;
            for(int j = 0; j < latents.size();j++)
            {
                Pair currGuess = latents.get(j);
                if(currGuess.one.getName().equals(currLatent.one.getName()) && currGuess.two.getName().equals(currLatent.two.getName()))
                    found = true;
                if(currGuess.two.getName().equals(currLatent.one.getName())&& currGuess.one.getName().equals(currLatent.two.getName()))
                    found = true;
            }
            if(found)
                tp++;
            else
                fn++;
        }

        return tp/(tp+fn);
    }


    //Method to evaluate non-latent variable algorithms based upon how successfully they can recover latent connections that are possible to get
    //From the true PAG
    public static double getIdentifiableRecall(Graph PAG, ArrayList<Pair> latents, Graph trueDAG, DataSet d, String type)
    {
        ArrayList<Pair> trueLatents = getLatents(trueDAG,d,type);

        ArrayList<Pair> temp = new ArrayList<Pair>();
        for(Pair p:trueLatents)
        {
            Edge e = PAG.getEdge(PAG.getNode(p.one.getName()),PAG.getNode(p.two.getName()));
            if(e==null)
                continue;
            if(e.getEndpoint1()== Endpoint.ARROW && e.getEndpoint2()==Endpoint.ARROW)
                temp.add(p);
        }
        trueLatents = temp;
        temp = null;
        double tp = 0;
        double fn = 0;
        for(int i = 0; i < trueLatents.size();i++)
        {
            Pair currLatent = trueLatents.get(i);
            boolean found = false;
            for(int j = 0; j < latents.size();j++)
            {
                Pair currGuess = latents.get(j);
                if(currGuess.one.getName().equals(currLatent.one.getName()) && currGuess.two.getName().equals(currLatent.two.getName()))
                    found = true;
                if(currGuess.two.getName().equals(currLatent.one.getName())&& currGuess.one.getName().equals(currLatent.two.getName()))
                    found = true;
            }
            if(found)
                tp++;
            else
                fn++;
        }

        return tp/(tp+fn);

    }
    //Compute recall based on possibly identifiable latents
    public static double getRecall(ArrayList<Pair> latents, int [][] splits, Graph trueGraph,DataSet dat,  String type)
    {
        ArrayList<Pair> trueLatents = getLatents(trueGraph,dat,type,splits);
        double tp = 0;
        double fn = 0;
        for(int i = 0; i < trueLatents.size();i++)
        {
            Pair currLatent = trueLatents.get(i);
            boolean found = false;
            for(int j = 0; j < latents.size();j++)
            {
                Pair currGuess = latents.get(j);
                if(currGuess.one.getName().equals(currLatent.one.getName()) && currGuess.two.getName().equals(currLatent.two.getName()))
                    found = true;
                if(currGuess.two.getName().equals(currLatent.one.getName())&& currGuess.one.getName().equals(currLatent.two.getName()))
                    found = true;
            }
            if(found)
                tp++;
            else
                fn++;
        }

        return tp/(tp+fn);
    }
    //TODO not guaranteed to return all latent interactions (splitting by variables doesn't allow you to measure all n^2 interactions)
    public ArrayList<Pair> runAlgorithm(final String algName, final boolean prune, final boolean requireFullEdge) throws Exception
    {
        final int [][] subsets = splitDataByColumns(data,nss);
        lastSplit = subsets;
        final ArrayList<Pair> latents = new ArrayList<Pair>();

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

            private synchronized void add(ArrayList<Pair> temp, Pair p)
            {
                temp.add(p);
            }
            @Override
            protected void compute() {
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        DataSet temp = data.subsetColumns(subsets[s]);
                        ArrayList<Pair> origLatents = new ArrayList<Pair>();
                        try {
                            origLatents = runRegularAlgorithm(algName, temp,requireFullEdge);
                        }
                        catch(Exception e)
                        {
                            e.printStackTrace();
                            System.out.println("Unrecognized Algorithm Name");
                            System.exit(0);
                        }

                        IndependenceTest ind = new IndTestMultinomialAJ(data,getAlpha(algName),true);
                        System.out.println("Alpha:" + ind.getAlpha());
                        if(prune) {
                            B: for (Pair p : origLatents) {
                                for (int i = 0; i < data.getNumColumns(); i++) {
                                    if (temp.getVariable(data.getVariable(i).getName()) != null)
                                        continue;
                                    //for each latent, test if both variables are dependent on this third variable, and test
                                    // if variable one is independent of variable three given variable 2
                                    if (ind.isDependent(ind.getVariable(p.one.getName()), ind.getVariable(data.getVariable(i).getName())) && ind.isDependent(ind.getVariable(p.two.getName()), ind.getVariable(data.getVariable(i).getName()))) {
                                        if (ind.isIndependent(ind.getVariable(p.one.getName()), ind.getVariable(p.two.getName()), ind.getVariable(data.getVariable(i).getName())))
                                            continue B;
                                    }

                                }
                                add(latents, p);
                                continue B;
                            }
                        }
                        else
                        {
                            for(Pair p:origLatents)
                                add(latents,p);
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


        final int chunk = 1000;

        pool.invoke(new StabilityAction(chunk, 0, nss));

        /*******************TEMPORARY SERIAL VERSION FOR TESTING AND DEBUGGING***********************

        for(int s = 0; s < nss;s++) {
            DataSet temp = data.subsetColumns(subsets[s]);
            ArrayList<Pair> origLatents = new ArrayList<Pair>();
            try {
                origLatents = runRegularAlgorithm(algName, temp);
            } catch (Exception e) {
                e.printStackTrace();
                System.out.println("Unrecognized Algorithm Name");
                System.exit(0);
            }
            System.out.println(origLatents);
            IndependenceTest ind = new IndTestMultinomialAJ(data, getAlpha(algName));
            System.out.println("Alpha:" + ind.getAlpha());
            if (prune) {
                B:
                for (Pair p : origLatents) {
                    for (int i = 0; i < data.getNumColumns(); i++) {
                        if (temp.getVariable(data.getVariable(i).getName()) != null)
                            continue;
                        //for each latent, test if both variables are dependent on this third variable, and test
                        // if variable one is independent of variable three given variable 2
                        if (ind.isDependent(ind.getVariable(p.one.getName()), ind.getVariable(data.getVariable(i).getName())) && ind.isDependent(ind.getVariable(p.two.getName()), ind.getVariable(data.getVariable(i).getName()))) {
                            if (ind.isIndependent(ind.getVariable(p.one.getName()), ind.getVariable(p.two.getName()), ind.getVariable(data.getVariable(i).getName())))
                                continue B;
                        }

                    }
                    latents.add(p);
                    continue B;
                }
            } else {
                for (Pair p : origLatents)
                    latents.add(p);
            }
        }
        /*********************DELETE THIS SECTION WHEN DONE TESTING**********************************/
        return latents;
    }


    // Generates a random permutation of the column indices and converts this permutation to a 2d array as close to evenly sized subsets as posible
    private int [][] splitDataByColumns(DataSet d, int numSubSets)
    {
        //16 num vars
        //3 numSubSets
        //varsperset = 6
        //numWithExtra = 0
       // int numVars = d.getNumColumns();
        int [] discInds = MixedUtils.getDiscreteInds(d.getVariables());
        int [] contInds = MixedUtils.getContinuousInds(d.getVariables());
        int numVars = contInds.length;
        int varsPerSet = (int)Math.ceil(numVars/(double)numSubSets);
        int numWithExtra = numVars % numSubSets;
        int [][] result = new int[numSubSets][];

        int [] perm = new int[numVars];
        for(int i = 0; i < perm.length;i++)
            perm[i] = contInds[i];

        int count = 0;
        int resultCount = 0;
        perm = shuffleArray(perm);
        //7 vars
        //4 subsets
        //vars per set = 2
        if(numWithExtra==0)
        {
            for(int i = 0; i < numSubSets;i++)
            {
                result[i] = new int[varsPerSet + discInds.length];
                while (count < (i + 1) * (varsPerSet)) {
                    result[i][resultCount] = perm[count];
                    count++;
                    resultCount++;
                }
                for(int j = 0; j < discInds.length;j++)
                {
                    result[i][resultCount] = discInds[j];
                    resultCount++;
                }
                resultCount = 0;
            }

        }
        else {
            for (int i = 0; i < numWithExtra; i++) {
                result[i] = new int[varsPerSet + discInds.length];
                while (count < (i + 1) * (varsPerSet)) {
                    result[i][resultCount] = perm[count];
                    count++;
                    resultCount++;
                }
                for(int j = 0; j < discInds.length;j++)
                {
                    result[i][resultCount] = discInds[j];
                    resultCount++;
                }
                resultCount = 0;
            }
            for (int i = numWithExtra; i < numSubSets; i++) {
                int limit = count + varsPerSet - 1;
                result[i] = new int[varsPerSet - 1 + discInds.length];
                while (count < limit) {
                    result[i][resultCount] = perm[count];
                    count++;
                    resultCount++;

                }
                for(int j = 0; j < discInds.length;j++) {
                    result[i][resultCount] = discInds[j];
                    resultCount++;
                }
                resultCount = 0;
            }
        }

        for(int i = 0; i < result.length;i++)
            System.out.println(Arrays.toString(result[i]));

        return result;
    }
    private boolean ensureVariance()
    {
        int [][] sb = this.subs;
        for(int i = 0; i < subs.length;i++)
        {
            if(runPriors.checkForVariance(data.subsetRows(sb[i]),data)!=-1)
                return false;

        }
        return true;
    }
    private int [] shuffleArray(int[] x)
    {
            // If running on Java 6 or older, use `new Random()` on RHS here
            Random rnd = ThreadLocalRandom.current();
            for (int i = x.length - 1; i > 0; i--)
            {
                int index = rnd.nextInt(i + 1);
                // Simple swap
                int a = x[index];
                x[index] = x[i];
                x[i] = a;
            }
            return x;
    }


    public static class Pair
    {
        public Node one;
        public Node two;
        public int colIDOne;
        public int colIDTwo;
        public double stability;
        public Pair(Node x, Node y)
        {
            one = x;
            two = y;
        }
        public Pair(Node x, Node y, double stab)
        {
            one = x;
            two = y;
            stability = stab;
        }
        public String toString()
        {
            return ("[" + one.getName() + "," + two.getName() + "]" + Double.toString(stability));
        }
    }
}
