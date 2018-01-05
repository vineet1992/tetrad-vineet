package edu.pitt.csb.latents;

import cern.colt.matrix.DoubleMatrix2D;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.Fci;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.pitt.csb.Pref_Div.Functions;
import edu.pitt.csb.Pref_Div.Gene;
import edu.pitt.csb.latents.mirnaPrediction;
import edu.pitt.csb.mgm.*;
import edu.pitt.csb.stability.SearchWrappers;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.ThreadLocalRandom;

/**
 * Created by vinee_000 on 10/9/2017. Runs modifications of FCI-flavor approaches to determine the variables likely to be confounded by latents
 */
public class LatentPrediction {


    public Map<String,String> orientations;
    private DataSet data;
    private int nss;
    private double tao;
    private int ns;
    private DataSet[] subsamples;
    private int [][] lastSplit;
    private static double [] searchParams;
    private static double [] initAlphas;
    private static double [] initLambdas;
    private int numLambdas = 30;
    private int numAlphas = 30;
    private double lambdaLow = 0.05;
    private double lambdaHigh = 0.9;
    private double alphaLow = 0.001;
    private double alphaHigh = 0.5;
    private double stepsG = 0.05;
    private double starsG = 0.05;
    private double useAlpha = -1;

    public Graph trueGraph;
    private int runNumber = 0;
    public LatentPrediction(DataSet d,int numSubSets, double tao, int numSubsamples)
    {
        data = d;
        nss = numSubSets;
        this.tao = tao;
        this.ns = numSubsamples;
        searchParams = new double[4];
        initializeArrays();
        orientations = new HashMap<String,String>();
        //TODO Generate subsamples here


    }
    public LatentPrediction(DataSet d, int numSubSets, double tao, DataSet [] subs)
    {
        data = d;
        nss = numSubSets;
        this.tao = tao;
        this.subsamples = subs;
        this.ns = subs.length;
        searchParams = new double[4];
        orientations = new HashMap<String,String>();
        initializeArrays();

    }

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



    public ArrayList<Pair> runRegularAlgorithm(String algName, boolean requireFullEdge) throws Exception
    {
        Graph gTotal = null;
        ArrayList<Pair> result = new ArrayList<Pair>();
        int b = (int)(10*Math.sqrt(data.getNumRows()));
        DoubleMatrix2D stabs;
        if(subsamples==null) {
            if (algName.equals("FCI")) {
                searchParams[0] = getStarsAlpha(data, Algorithm.FCI);
                stabs = StabilityUtils.StabilitySearchParLatent(data, new SearchWrappers.FCIWrapper(searchParams), ns, b,orientations);
                IndependenceTest i = new IndTestMultinomialAJ(data,searchParams[0]);
                Fci f = new Fci(i);
                gTotal = f.search();
                ConcurrentHashMap<String,String>tempOr = f.whyOrient;
                for(String x:tempOr.keySet())
                {
                    if(orientations.get(x)!=null)
                    {
                        orientations.put(x,orientations.get(x)+"\tFINAL," + tempOr.get(x));
                    }
                }
            }
            else if (algName.equals("MGM-FCI")) {
                double [] temp = getStepsLambda(data);
                searchParams[0] = temp[0];
                searchParams[1] = temp[1];
                searchParams[2] = temp[2];
                searchParams[3] = getStarsAlpha(data,Algorithm.MGMFCI,temp);
               stabs = StabilityUtils.StabilitySearchParLatent(data, new SearchWrappers.MGMFCIWrapper(searchParams), ns, b,orientations);
                MGM m = new MGM(data,temp);
                m.learnEdges(1000);
                IndependenceTest i = new IndTestMultinomialAJ(data,searchParams[3]);
                Fci f = new Fci(i);
                f.setInitialGraph(m.graphFromMGM());
                gTotal = f.search();
               ConcurrentHashMap<String,String>tempOr = f.whyOrient;
                for(String x:tempOr.keySet())
                {
                    if(orientations.get(x)!=null)
                    {
                        orientations.put(x,orientations.get(x)+"\tFINAL," + tempOr.get(x));
                    }
                }
            }
            else if (algName.equals("MGM-FCI-MAX")) {
                double [] temp = getStepsLambda(data);

                searchParams[0] = temp[0];
                searchParams[1] = temp[1];
                searchParams[2] = temp[2];
                searchParams[3] = getStarsAlpha(data,Algorithm.MGMFCIMAX,temp);
                stabs = StabilityUtils.StabilitySearchParLatent(data,new SearchWrappers.MFMWrapper(searchParams),ns,b,orientations);
                MGM m = new MGM(data,temp);
                m.learnEdges(1000);
                IndependenceTest i = new IndTestMultinomialAJ(data,searchParams[3]);
                FciMaxP f = new FciMaxP(i);
                f.setInitialGraph(m.graphFromMGM());
                gTotal = f.search();
                ConcurrentHashMap<String,String>tempOr = f.whyOrient;
                for(String x:tempOr.keySet())
                {
                    if(orientations.get(x)!=null)
                    {
                        orientations.put(x,orientations.get(x)+"\tFINAL," + tempOr.get(x));
                    }
                }
            } else {
                throw new Exception("Unrecognized Algorithm");
            }
        }
        else
        {
            if (algName.equals("FCI")) {
                System.out.println("Computing Optimal Alpha...");
                searchParams[0] = getStarsAlpha(data,Algorithm.FCI);
                System.out.println("Alpha: " + searchParams[0]);
                System.out.println("Computing stable edges...");
                stabs = StabilityUtils.StabilitySearchParLatent(data, new SearchWrappers.FCIWrapper(searchParams), subsamples,orientations);
                IndependenceTest i = new IndTestMultinomialAJ(data,searchParams[0]);
                Fci f = new Fci(i);
                gTotal = f.search();
                ConcurrentHashMap<String,String>tempOr = f.whyOrient;
                for(String x:tempOr.keySet())
                {
                    if(orientations.get(x)!=null)
                    {
                        orientations.put(x,orientations.get(x)+"\tFINAL," + tempOr.get(x));
                    }
                }
            }
            else if (algName.equals("MGM-FCI")) {
                double [] temp = getStepsLambda(data);
                searchParams[0] = temp[0];
                searchParams[1] = temp[1];
                searchParams[2] = temp[2];
                searchParams[3] = getStarsAlpha(data,Algorithm.MGMFCI,temp);
                stabs = StabilityUtils.StabilitySearchParLatent(data, new SearchWrappers.MGMFCIWrapper(searchParams),subsamples,orientations);
                MGM m = new MGM(data,temp);
                m.learnEdges(1000);
                IndependenceTest i = new IndTestMultinomialAJ(data,searchParams[3]);
                Fci f = new Fci(i);
                f.setInitialGraph(m.graphFromMGM());
                gTotal = f.search();
                ConcurrentHashMap<String,String>tempOr = f.whyOrient;
                for(String x:tempOr.keySet())
                {
                    if(orientations.get(x)!=null)
                    {
                        orientations.put(x,orientations.get(x)+"\tFINAL," + tempOr.get(x));
                    }
                }
            }
            else if (algName.equals("MGM-FCI-MAX")) {
                System.out.println("Computing Optimal Lambda");
                double [] temp = getStepsLambda(data);
                System.out.println("Lambda: " + Arrays.toString(temp));

                searchParams[0] = temp[0];
                searchParams[1] = temp[1];
                searchParams[2] = temp[2];
                System.out.println("Computing Optimal Alpha");
                searchParams[3] = getStarsAlpha(data,Algorithm.MGMFCIMAX,temp);
                System.out.println("Alpha: " + searchParams[3]);
                System.out.println("Computing stable edges...");
                stabs = StabilityUtils.StabilitySearchParLatent(data,new SearchWrappers.MFMWrapper(searchParams),subsamples,orientations);
                MGM m = new MGM(data,temp);
                m.learnEdges(1000);
                IndependenceTest i = new IndTestMultinomialAJ(data,searchParams[3]);
                FciMaxP f = new FciMaxP(i);
                f.setInitialGraph(m.graphFromMGM());
                gTotal = f.search();
                Map<String,String>tempOr = f.whyOrient;
                for(String x:tempOr.keySet())
                {
                    if(orientations.get(x)!=null)
                    {
                        orientations.put(x,orientations.get(x)+"\tFINAL," + tempOr.get(x));
                    }
                }
            } else {
                throw new Exception("Unrecognized Algorithm");
            }
        }

        System.out.println("Graph for Regular Algorithm: " + gTotal);

        for(int i = 0; i < stabs.rows();i++)
        {
            for(int j = i+1; j < stabs.columns();j++)
            {
                if(data.getVariable(i)instanceof DiscreteVariable || data.getVariable(j) instanceof DiscreteVariable)
                    continue;
                Edge e = gTotal.getEdge(gTotal.getNode(data.getVariable(i).getName()),gTotal.getNode(data.getVariable(j).getName()));
                if(e!=null && e.getEndpoint1()==Endpoint.ARROW && e.getEndpoint2()==Endpoint.ARROW) {
                    if (stabs.get(i, j) > tao) {
                        result.add(new Pair(data.getVariable(i), data.getVariable(j), stabs.get(i, j)));
                    }
                }
                else if(!requireFullEdge && stabs.get(i,j) > tao)
                {
                    result.add(new Pair(data.getVariable(i),data.getVariable(j),stabs.get(i,j)));
                }
            }
        }
        System.out.println("Latent list for Regular Algorithm: " + result);
        return result;
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
    private double getStarsAlpha(DataSet d, Algorithm a,DataSet[]currSubSamples)
    {
        if(useAlpha > 0)
            return useAlpha;
        STARS s;
        //System.out.println(d);
        //System.out.println("Curr Sub Samples: " + Arrays.toString(currSubSamples));
        if(currSubSamples==null)
            s = new STARS(d,initAlphas,starsG,ns,a);
        else
            s = new STARS(d,initAlphas,starsG,currSubSamples,a);
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
        STARS s;
        if(subsamples==null)
            s = new STARS(d,initAlphas,starsG,ns,a);
        else
            s = new STARS(d,initAlphas,starsG,subsamples,a);
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

    private double getStarsAlpha(DataSet d, Algorithm a, double [] lambda,DataSet [] currSubSamples)
    {
        if(useAlpha > 0)
            return useAlpha;
        STARS s;
        if(currSubSamples==null)
            s = new STARS(d,initAlphas,starsG,ns,a);
        else
            s = new STARS(d,initAlphas,starsG,currSubSamples,a);
        s.setMGMLambda(lambda);
        s.setTrueGraph(trueGraph);
        s.setRun(runNumber);
        s.runSTARS();
        return s.lastAlpha;
    }
    private double[] getStepsLambda(DataSet d)
    {
        STEPS s;
        if(subsamples==null)
            s = new STEPS(d,initLambdas,stepsG,ns);
        else
            s = new STEPS(d,initLambdas,stepsG,subsamples);
        s.runStepsPar();
        return s.lastLambda;
    }
    private double[] getStepsLambda(DataSet d,DataSet[]currSubSamples)
    {
        STEPS s;
        if(currSubSamples==null)
            s = new STEPS(d,initLambdas,stepsG,ns);
        else
            s = new STEPS(d,initLambdas,stepsG,currSubSamples);
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

    //TODO Need to add gTotal and currSubSamples into all of the cases
    //Before using it
    private ArrayList<Pair> runRegularAlgorithm(String algName, DataSet d, boolean requireFullEdge) throws Exception
    {
        DataSet [] currSubSamples = new DataSet[subsamples.length];
        ArrayList<Pair> result = new ArrayList<Pair>();
        int b = (int)(10*Math.sqrt(d.getNumRows()));
        DoubleMatrix2D stabs;
        Graph gTotal = null;
        if(subsamples==null) {
            if (algName.equals("FCI")) {

                searchParams[0] = getStarsAlpha(d,Algorithm.FCI);
                stabs = StabilityUtils.StabilitySearchParLatent(d, new SearchWrappers.FCIWrapper(searchParams), ns, b,orientations);
            }
            else if (algName.equals("MGM-FCI")) {

                double [] temp = getStepsLambda(d);
                searchParams[0] = temp[0];
                searchParams[1] = temp[1];
                searchParams[2] = temp[2];
                searchParams[3] = getStarsAlpha(d,Algorithm.MGMFCI,temp);
                stabs = StabilityUtils.StabilitySearchParLatent(d, new SearchWrappers.MGMFCIWrapper(searchParams), ns, b,orientations);
            }
            else if (algName.equals("MGM-FCI-MAX")) {
                double [] temp = getStepsLambda(d);

                searchParams[0] = temp[0];
                searchParams[1] = temp[1];
                searchParams[2] = temp[2];
                searchParams[3] = getStarsAlpha(d,Algorithm.MGMFCIMAX,temp);
                stabs = StabilityUtils.StabilitySearchParLatent(d,new SearchWrappers.MFMWrapper(searchParams),ns,b,orientations);
            } else {
                throw new Exception("Unrecognized Algorithm");
            }
        }
        else
        {

            for(int i = 0; i < subsamples.length;i++)
            {
                currSubSamples[i] = subsamples[i].copy();
                int count = 0;
                int [] indsToKeep = new int[d.getNumColumns()];
                for(int j = 0; j < subsamples[i].getNumColumns();j++)
                {
                    if(d.getVariable(subsamples[i].getVariable(j).getName())!=null) {
                        indsToKeep[count] = j;
                        count++;
                    }
                }
                currSubSamples[i] = currSubSamples[i].subsetColumns(indsToKeep);
            }
            if (algName.equals("FCI")) {
                searchParams[0] = getStarsAlpha(d,Algorithm.FCI,currSubSamples);
                stabs = StabilityUtils.StabilitySearchParLatent(d, new SearchWrappers.FCIWrapper(searchParams), currSubSamples,orientations);
                IndependenceTest i = new IndTestMultinomialAJ(d,searchParams[0]);
                Fci f = new Fci(i);
                gTotal = f.search();
                ConcurrentHashMap<String,String>tempOr = f.whyOrient;
                addToOrient(tempOr);

            }
            else if (algName.equals("MGM-FCI")) {
                double [] temp = getStepsLambda(d);
                searchParams[0] = temp[0];
                searchParams[1] = temp[1];
                searchParams[2] = temp[2];
                searchParams[3] = getStarsAlpha(d,Algorithm.MGMFCI,temp);
                stabs = StabilityUtils.StabilitySearchParLatent(d, new SearchWrappers.MGMFCIWrapper(searchParams),currSubSamples,orientations);
            }
            else if (algName.equals("MGM-FCI-MAX")) {
                System.out.println("Computing Optimal Lambda");
                double [] temp = getStepsLambda(d,currSubSamples);
                System.out.println("Lambda: " + Arrays.toString(temp));

                searchParams[0] = temp[0];
                searchParams[1] = temp[1];
                searchParams[2] = temp[2];
                System.out.println("Computing Optimal Alpha...");
                searchParams[3] = getStarsAlpha(d,Algorithm.MGMFCIMAX,temp,currSubSamples);
                System.out.println("Alpha: " + searchParams[3]);
                System.out.println("Computing Stable Edges...");
                stabs = StabilityUtils.StabilitySearchParLatent(d,new SearchWrappers.MFMWrapper(searchParams),currSubSamples,orientations);
                MGM m = new MGM(d,temp);
                m.learnEdges(1000);
                IndependenceTest i = new IndTestMultinomialAJ(d,searchParams[3]);
                FciMaxP f = new FciMaxP(i);
                f.setInitialGraph(m.graphFromMGM());
                gTotal = f.search();
                ConcurrentHashMap<String,String>tempOr = (ConcurrentHashMap<String,String>)f.whyOrient;
                addToOrient(tempOr);

            } else {
                throw new Exception("Unrecognized Algorithm");
            }
        }

        System.out.println("Graph for latent algorithm: " + gTotal);
            for (int i = 0; i < stabs.rows(); i++) {
                for (int j = i + 1; j < stabs.columns(); j++) {
                    if(currSubSamples[0].getVariable(i) instanceof DiscreteVariable || currSubSamples[0].getVariable(j) instanceof DiscreteVariable)
                        continue;
                    Edge e = gTotal.getEdge(gTotal.getNode(currSubSamples[0].getVariable(i).getName()), gTotal.getNode(currSubSamples[0].getVariable(j).getName()));
                    if (e != null && e.getEndpoint1() == Endpoint.ARROW && e.getEndpoint2() == Endpoint.ARROW) {
                        if (stabs.get(i, j) > tao) {
                            result.add(new Pair(currSubSamples[0].getVariable(i), currSubSamples[0].getVariable(j), stabs.get(i, j)));
                        }
                    }
                    else if(!requireFullEdge && stabs.get(i,j)>tao)
                    {
                        result.add(new Pair(currSubSamples[0].getVariable(i),currSubSamples[0].getVariable(j),stabs.get(i,j)));
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
    //TODO Maintain all discrete variables in each split, just split the continuous variables into separate subsets (since this is what will be done in the real data)
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
                        System.out.println(data.getVariableNames());
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

                        IndependenceTest ind = new IndTestMultinomialAJ(data,getAlpha(algName));
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
