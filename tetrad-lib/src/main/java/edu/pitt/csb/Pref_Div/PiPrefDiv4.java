package edu.pitt.csb.Pref_Div;

import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.Priors.mgmPriors;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.stability.StabilityUtils;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.IntBuffer;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

/**
 * Created by vinee_000 on 10/19/2018.
 * This version uses a single file for each prior information source using both similarity and intensity
 * This version will also use only radius of similarity instead of both that and p-value (because some p-values are 0, making them impossible to discount by using prior information)
 */



public class PiPrefDiv4 {

    private DataSet data; //The current expression dataset to be analyzed, assumes that all variables are fair game for Pref-Div except the target
    private String target; //The target variable of interest
    private double [] initRadii; //Initial radius values to test
    private double lowRadii = 0.25; //Low range of radius values to test
    private double highRadii = 1; //High range of radius values to test
    private int numRadii = 40; //Number of radius values to test
    private int [][] subsamples; //subsamples for repeated Pref-Div
    private double normalEpsilon = 0.5;//Range around which to get probability for theta with prior information
    private int numGenes;//number of total variables in the data
    private double [] intTao;//Tao values for each intensity source
    private double [] simTao;//Tao values for each similarity source
    private boolean LOOCV = false; //Should we do leave-one-out cross validation for the data-theory tradeoff selection
    private Map<Gene,List<Gene>> lastCluster;//The last set of clusters from the run of selectGenes()
    private boolean verbose = false; //Should we print full debugging information?
    private int K = 5; //Number of genes to select
    private boolean clusterByCorrs = true;//Do we use all correlations to cluster the genes or standard Pref-Div?
    private double [] lastSimilarityWeights;
    private boolean outputScores; //Should we output parameter scores to file?
    private PrintStream scoreStream; //Printstream to output parameter scores to file
    private double lastRadius; //Last Radius chosen via best score
    private int numFolds; //Number of folds for cross-validation
    private boolean useStabilitySelection; //Should stability selection be used to select genes in PD?
    private DataSet summarizedData; //Summarized dataset using clusters as variables
    private boolean parallel; //Should we compute stabilities in parallel?
    private RunPrefDiv.ClusterType ctype; //Which clustering method should be used to produce summarized data?
    private boolean pdStability; //Should we run Pref-Div when computing stabilities or just use correlationa with p-value threshold?
    private boolean partialCorrs; //Should we use partial correlations?
    private double lastRadiusWP; //Last Radius WP selected
    private boolean [] iPriors; //Contains which features have prior intensity information
    private boolean [] dPriors; //Contains which features have prior dissimiliarty information
    private boolean saveMemory = false; //Should we conserve disc space and use the slower version of PiPrefDiv
    public double lastThreshold = 0;
    public double lastThresholdWP = 0;


    //Constructor that uses default parameters for both initial parameter ranges
    public PiPrefDiv4(DataSet data, String target,int K)
    {
        this.K = K;
        this.data = data;
        this.numGenes = data.getNumColumns()-1;
        this.target = target;
        this.initRadii = initRadiiDefault();
    }


    //Constructor that specifies the same number of parameters for radius and k
    public PiPrefDiv4(DataSet data, String target, int K,int numParams)
    {

        this.K = K;
        this.numGenes = data.getNumColumns()-1;
        this.data = data;
        this.target = target;
        this.numRadii = numParams;
        this.initRadii = initRadiiDefault();
        this.numFolds = 5;
    }


    //Constructor with initial ranges both specified
    public PiPrefDiv4(DataSet data, String target, int K, double [] radii)
    {
        this.data = data;
        this.numGenes = data.getNumColumns()-1;
        this.target = target;
        this.initRadii = radii;
        this.numRadii = radii.length;
        this.K = K;
        this.numFolds = 5;
    }


    public void setSubsamples(int[][] subs){this.subsamples = subs;}
    public Map<Gene,List<Gene>> getLastCluster()
    {
        return lastCluster;
    }
    public double[] getLastSimilarityWeights(){return lastSimilarityWeights;}
    public double[] getLastRadius(){return new double[]{lastRadius,lastRadiusWP};}
    public double getLastRadiusWP(){return lastRadiusWP;}
    public void setOutputScores(PrintStream out){outputScores=true; scoreStream=out;}
    public void setUseStabilitySelection(boolean b){useStabilitySelection = b;}
    public void setParallel(boolean b){parallel = b;}
    public void setClusterType(RunPrefDiv.ClusterType c){ctype = c;}
    public DataSet getSummarizedData(){return summarizedData;}
    public void setPdStability(boolean p){pdStability = p;}
    public void setPartialCorrs(boolean pc){partialCorrs = pc;}
    public void saveMemory(){saveMemory = true;}
    public boolean[] getIntensityPriors(){return iPriors;}
    public boolean[] getDissimilarityPriors(){return dPriors;}
    public double [] getLastThreshold(){return new double[]{lastThreshold,lastThresholdWP};}
    public double [] getTestedThresholds(){return new double[]{0,0.5,1};}
    public double [] getLastIntensityWeights(){return new double[10];}
    public double getLastThresholdWP(){return lastThresholdWP;}


    //Set numbe of folds
    public void setNumFolds(int nf){numFolds = nf;}

    //Setter for whether or not to do leave-one-out cross validation
    public void setLOOCV(boolean loocv)
    {
        LOOCV = loocv;
    }

    //Sets verbose to true for full debugging information
    public void setVerbose()
    {
        verbose = true;
    }

    public double[] getTestedRadii()
    {
        return initRadii;
    }
    //Create default radius of similarity range to examine
    private double[] initRadiiDefault()
    {
        double [] radii = new double[numRadii];

        for(int i = 0; i < numRadii;i++)
        {
            radii[i] = lowRadii + (highRadii-lowRadii)*i/(double)numRadii;
        }

        return radii;
    }

    /***Constricts the testing range of radius and p-value thresholds
     *
     * @param init initial range of the parameter
     * @param data Full Dataset to compute independence tests and p-values
     * @param threshold Whether or not this is p-value threshold or absolute correlation
     * @return double [] with a constricted range of the parameter
     */

    //TODO instead of running this on the full dataset, we want the average # across subsamples
    private double[] constrictRange(double[]init, DataSet data,boolean threshold)
    {
        int[] num = new int[init.length];
        for(int x = 0; x < subsamples.length; x++) {

            DataSet curr = data.subsetRows(subsamples[x]);
            ArrayList<Gene> temp = createGenes(data,target,false);
            Gene tgt = new Gene(temp.size());
            tgt.symbol="Target";
            temp.add(tgt);

            float[] corrs = Functions.computeAllCorrelations(temp, curr, partialCorrs, false, threshold, 1);

            for (int i = 0; i < init.length; i++) {
                if (!threshold)
                    num[i] += countLessThan(corrs, init[i]) / (double)init.length;
                else
                    num[i] += countPLessThan(corrs, init[i])/(double)init.length;
            }
        }
        TetradMatrix t = new TetradMatrix(init.length,2);
        for(int i = 0; i < t.rows();i++) {

            t.set(i,0,init[i]);
            t.set(i,1,num[i]);
        }
        double[] limits = mgmPriors.getLimit(t);
        double [] realLimits = new double[init.length];
        for(int i = 0; i < init.length;i++)
            realLimits[i] = mgmPriors.map(init[i],init,limits);

        return realLimits;
        /***Loop through radii, for each radii, how many edges exist***/

    }


    private int countPLessThan(float[] pvals, double t)
    {
        int num = 0;
        for(int i = 0; i < pvals.length;i++)
        {
            double curr = Math.log(pvals[i]);
            if(curr < t)
                num++;
        }
        return num;
    }
    private int countLessThan(float [] c, double t)
    {
        int result = 0;
        for(int i = 0; i < c.length;i++)
        {
            if(1-c[i] < t)
                result++;
        }
        return result;
    }


    public double [][] evaluatePriors(boolean boot, int numSamples,String iFile, String [] dFile,boolean useCausalGraph)
    {

        //Generate subsamples of the data
        System.out.print("Generating subsamples...");
        if(subsamples==null)
        {
            subsamples = genSubsamples(boot,numSamples,data,LOOCV);
        }
        System.out.println("Done");

        System.out.print("Constricting Parameter Range...");
        initRadii = constrictRange(initRadii,data,false);
        System.out.println("Done");

        //Compute stability of each gene selection and each gene-gene relationship
        System.out.print("Computing stability across radii and threshold values...");
        int [] sums = computeSums(); //i -> Radii, j -> TopK, k -> geneID
        System.out.println("Done");

        System.out.print("Computing intensity weights...");
        //Compute intensity weights
        double [] intWeights = getIntensityWeights(sums,iFile);
        System.out.println("Done");
        if(verbose)
        {
            //TODO get header from iFile to attach source to weights
            System.out.println("Intensity weights are: " + Arrays.toString(intWeights));
        }

        System.out.print("Computing similarity weights...");
        //Compute similarity weights for each prior knowledge source
        double [] simWeights = getSimilarityWeights(sums,dFile);

        double[][] allWeights = new double[2][];
        allWeights[0] = intTao;
        System.out.println("Intensity Tao: " + Arrays.toString(allWeights[0]));
        allWeights[1] = simTao;
        System.out.println("Similarity Tao: " + Arrays.toString(allWeights[1]));
        return allWeights;
    }


    /***
     *
     * @param boot Should we bootstrap instead of subsampling?
     * @param numSamples How many subsamples should we draw?
     * @param useCausalGraph Should we use a causal graph?
     * @return
     */
    public ArrayList<Gene> selectGenes(boolean boot, int numSamples,boolean useCausalGraph)
    {
        if(data.isContinuous()) {
            System.out.print("Standardizing Data...");
            data = DataUtils.standardizeData(data);
            System.out.println("Done");
        }

        //Generate subsamples of the data
        System.out.print("Generating subsamples...");
        if(subsamples==null)
        {
            subsamples = genSubsamples(boot,numSamples,data,LOOCV);
        }
        System.out.println("Done");

        //Constrict parameter ranges
        System.out.print("Constricting Parameter Range...");
        initRadii = constrictRange(initRadii,data,false);

        System.out.println("Done");

        if(verbose)
        {
            System.out.println("Selected radii range: (" + initRadii[0] + "," + initRadii[initRadii.length-1] + ")");

        }


        //Compute stability of each gene selection and each gene-gene relationship
        System.out.print("Computing stability across radii and threshold values...");


        double[] scores = computeScores();


        double maxRadii = -1;
        double bestScore = -1;
        for(int i = 0; i < scores.length;i++)
        {
                if(scores[i] > bestScore)
                {
                    bestScore = scores[i];
                    maxRadii = initRadii[i];

                }
        }

        lastRadius = maxRadii;

        if(verbose)
        {
            System.out.println("Best Radii NP: " + lastRadius);
        }


        System.out.println("Done");





        System.out.print("Running Pref-Div with optimal parameters...");
        //Run Pref-Div with optimal parameters


        ArrayList<Gene> temp = createGenes(data,target,false);
        ArrayList<Gene> meanGenes = Functions.computeAllIntensities(temp,1,data,target,partialCorrs,false,false,1);

        float [] meanDis = Functions.computeAllCorrelations(meanGenes,data,partialCorrs,false,false,1);

        Collections.shuffle(meanGenes);
        Collections.sort(meanGenes,Gene.IntensityComparator);



        RunPrefDiv rpd;
        rpd = new RunPrefDiv(meanDis,meanGenes,data,target,LOOCV);

        rpd.setTopK(K);
        rpd.setAccuracy(0);
        rpd.setNS(subsamples.length);
        rpd.setNumFolds(numFolds);
        rpd.setCausalGraph(useCausalGraph);
        rpd.useCrossValidation(false);
        rpd.useStabilitySelection(useStabilitySelection);
        rpd.setClusterType(ctype);
        if(clusterByCorrs)
            rpd.clusterByCorrs();
        rpd.setRadius(lastRadius);
        rpd.setSubs(subsamples);
        rpd.setPThreshold(1);

        ArrayList<Gene> top = rpd.runPD();
        HashMap<Gene,List<Gene>> map = rpd.getClusters();
        summarizedData = rpd.getSummarizedData();
        lastCluster = map;
        System.out.println("Done");

        if(verbose)
        {
            System.out.println("Selected top K features\n" + top);
            System.out.println("All clusters\n" + map);
        }


        return top;
    }

    /***Full procedure to select genes based on prior information Pref-Div
     //Input: boot (should we do bootstrapping (true) or subsampling (false))
     //Input: numSamples (number of samples to subsample or bootstrap)
     //Input: path to theory Intensity File -> iFile (should have a header)
     //Input: path to theory Dissimilarity File -> dFile (no header or rownames for these files)****/
    public ArrayList<Gene> selectGenes(boolean boot, int numSamples, String [] dFile,boolean useCausalGraph)
    {
        if(data.isContinuous()) {
            System.out.print("Standardizing Data...");
            data = DataUtils.standardizeData(data);
            System.out.println("Done");
        }

        dPriors = new boolean[numGenes * (numGenes + 1) / 2];

        //Generate subsamples of the data
        System.out.print("Generating subsamples...");
        if(subsamples==null)
        {
            subsamples = genSubsamples(boot,numSamples,data,LOOCV);
        }
        System.out.println("Done");

        //Constrict parameter ranges
        System.out.print("Constricting Parameter Range...");
        initRadii = constrictRange(initRadii,data,false);

        System.out.println("Done");

        if(verbose)
        {
            System.out.println("Selected radii range: (" + initRadii[0] + "," + initRadii[initRadii.length-1] + ")");
        }

        //Compute stability of each gene selection and each gene-gene relationship
        System.out.print("Computing stability across radii and threshold values...");


        double[] scores = computeScores(dFile);

        int [] inds = getMaxScoreSeparate(scores);


/***The method works from here below TODO debug everything above this (why can't we select great parameters?)***/


         double bestRadii = initRadii[inds[0]];
         double bestRadiiNP = initRadii[inds[1]];

        lastRadius = bestRadiiNP;
        lastRadiusWP = bestRadii;


        System.out.println("Done");

        if(verbose)
        {
            System.out.println("Best Radii WP/NP: " + bestRadii + "," + bestRadiiNP);
        }


        System.out.print("Running Pref-Div with optimal parameters...");
        //Run Pref-Div with optimal parameters

        ArrayList<Gene> temp = createGenes(data,target,false);

        iPriors = loadIPrior(dFile);
        boolean [] tempPrior = new boolean[numGenes*(numGenes-1)/2];
        for(int i = numGenes; i < dPriors.length;i++)
        {
            tempPrior[i-numGenes] = dPriors[i];
        }
        ArrayList<Gene> meanGenes = Functions.computeAllIntensities(temp,1,data,target,false,1,1,iPriors);
        //If you aren't within radius of the target then you are shrunk to zero
        for(int i = 0; i < meanGenes.size();i++)
        {
            if(iPriors[i] && 1-meanGenes.get(i).intensityValue>lastRadiusWP)
                meanGenes.get(i).intensityValue = 0;
            else if(!iPriors[i] && 1-meanGenes.get(i).intensityValue>lastRadius)
                meanGenes.get(i).intensityValue=0;
        }
        float [] meanDis = Functions.computeAllCorrelations(meanGenes,data,false,1,1,tempPrior);


        Collections.shuffle(meanGenes);
        Collections.sort(meanGenes,Gene.IntensityComparator);

        RunPrefDiv rpd;
        rpd = new RunPrefDiv(meanDis,meanGenes,data,target,LOOCV);


        rpd.setWithPrior(tempPrior);



        rpd.setAllParams(bestRadiiNP,bestRadii,1,1);

        rpd.setTopK(K);
        rpd.setAccuracy(0);
        rpd.setNS(subsamples.length);
        rpd.setNumFolds(numFolds);
        rpd.setCausalGraph(useCausalGraph);
        rpd.useStabilitySelection(useStabilitySelection);
        rpd.useCrossValidation(false);
        rpd.setClusterType(ctype);
        rpd.setSubs(subsamples);
        if(clusterByCorrs)
            rpd.clusterByCorrs();

        ArrayList<Gene> top = rpd.runPD();
        HashMap<Gene,List<Gene>> map = rpd.getClusters();
        summarizedData = rpd.getSummarizedData();
        lastCluster = map;
        System.out.println("Done");

        if(verbose)
        {
            System.out.println("Selected top K features\n" + top);
            System.out.println("All clusters\n" + map);
        }


        return top;
    }



    private boolean [] loadIPrior(String [] dFile)
    {
        boolean [] temp = new boolean[numGenes];
        try{
            for(int i = 0; i < dFile.length;i++)
            {
                BufferedReader b = new BufferedReader(new FileReader(dFile[i]));
                String [] line = b.readLine().split("\t");
                for(int j = 1; j < line.length;j++)
                {
                    if(Double.parseDouble(line[j])>=0)
                        temp[j-1] = true;
                }

            }
            return temp;
        }
        catch(Exception e)
        {
            System.err.println("Could not load intensity prior for file run of Pref-Div because of this error:");
            e.printStackTrace();
        }
        return null;

    }


    private int[] getMaxScoreSeparate(double [] avg)
    {
        int [] inds = new int[2];
        double maxScore = -1;
        double maxScoreNP = -1;
        for(int i = 0; i < avg.length/2;i++)
        {

                if(avg[i]>maxScore)
                {
                    inds[0] = i;
                    maxScore = avg[i];
                }
                if(avg[i+avg.length/2]>maxScoreNP)
                {
                    inds[1] = i;
                    maxScoreNP = avg[i+avg.length/2];
                }
        }
        return inds;
    }


    /****
     *
     * @param n Number of times the edge was selected across subsamples for this parameter
     * @return The instability value for this number
     */
    private double getG(int n)
    {
        double f = n/(double)subsamples.length;
        return 4*f*(1-f);
    }

    //Get theta value given that the gene showed up n times for these parameters, and sum times for all parameters
    private double getTheta(int sum, int n)
    {
        double P = sum/(double)(numRadii*subsamples.length);
        //Given binomial with probability P, what's the probability the gene showed up n times in q subsamples
        BinomialDistribution b = new BinomialDistribution(subsamples.length,P);
        return b.probability(n);
    }

    //Overloaded version of getTheta for those selections where we have information from at least one prior information source
    private double getTheta(float uPosterior, float varPosterior, int n)
    {
        NormalDistribution b;
        try {
            b = new NormalDistribution(uPosterior, varPosterior);
            return b.probability(n - normalEpsilon, n + normalEpsilon);
        } catch (Exception e) //If there is no variance, then only give probability if we match the posterior mean exactly
        {
            if (varPosterior == 0) {
                //           System.out.println("Temp: " + temp.get(ii,jj) + ", Posterior: " + u_posterior.get(ii,jj));
                if (n == (int)uPosterior)
                    return 1;
                else
                    return 0;
            } else {
                System.out.println(uPosterior + "\t" + varPosterior + "\t" + n);
                throw e;
            }
        }
    }




    /**
     * @param sums, Array of number of times each gene and gene-gene pair appeared
     * @param uPost, Posterior mean for each gene OR each gene-gene pair appearence
     * @param varPost, Posterior variance for each gene OR each gene-gene pair appearence
     * @return Scores for each pair of radii, threshold params. Either for intensity or dissimilarity depending upon input
     */
    private double[] getScoresSeparate(int [] sums, float [] uPost, float [] varPost)
    {
        double [] score = new double[numRadii*2];

        for (int i = 0; i < numRadii; i++) {
            int [] all = new int[sums.length];

                int [] curr = null;
                if(saveMemory) {

                    curr = getCounts(i, true);

                }
                else {
                    curr = getCountsFromFile(i,true);

                }

                for(int x = 0; x < curr.length;x++)
                    all[x] = curr[x];

                for(int g = 0; g < uPost.length;g++) //Loop through each gene
                {
                    double G = getG(curr[g]); //Clusts[i][j][g] is number of times gene was selected in the subsamples for this parameter setting
                    if (uPost[g] < 0)//No Prior information, contributes only stability to the score
                    {
                        double theta = getTheta(sums[g],curr[g]); //sums is the total number of times this gene was selected across all parameter settings
                        score[i + numRadii] += theta*(1-G);
                    }
                    else
                    {
                        if(varPost[g]<0)
                            System.out.println(g);
                        double theta = getTheta(uPost[g],varPost[g],curr[g]);
                        score[i]+=theta*(1-G);
                    }
                }

            }

        normalizeScore(score);
        return score;
    }

    /**
     *
     * @param i, index of initRadii to use for computations
     * @param needCorrs, Do we need to perform correlation computation?
     * @return An array of the number of times each gene was selected and each gene-gene pair was clustered
     *
     */
    private int[] getCounts(int i, boolean needCorrs)
    {
        int [] curr = new int[(numGenes*(numGenes+1))/2];
        for(int k = 0; k < subsamples.length;k++)
        {
            ArrayList<Gene> temp = createGenes(data,target,true);
            DataSet currData = data.subsetRows(subsamples[k]);
            float [] corrs = Functions.computeAllCorrelations(temp,currData,partialCorrs,false,false,1);
            addGeneConnections(corrs,curr,initRadii[i]);

        }
        return curr;
    }


    /**
     *
     * @param i, index of initRadii to use for computations
     * @param needCorrs, Do we need to perform correlation computation?
     * @return An array of the number of times each gene was selected and each gene-gene pair was clustered
     * When needCorrs is true the array has both number of times each gene was selected and each gene-gene pair
     *
     */
    private int[] getCountsFromFile(int i, boolean needCorrs)
    {
        int fullSize = numGenes*(numGenes+1)/2;
        int [] curr = new int[fullSize];

        try{
            File iFile =new File("Correlations_" + i +  ".txt");
            BufferedInputStream corrReader = new BufferedInputStream(new FileInputStream(iFile));
            int [] array = Functions.readAsIntArray(corrReader,fullSize);


            curr = array;
            array = null;
            corrReader.close();
            if(needCorrs)
                iFile.deleteOnExit();


        }
        catch(Exception e)
        {
            e.printStackTrace();
            System.err.println("Unable to load correlation values from file");
        }

        return curr;
    }


    /**
     *
     * @param i, index of initRadii to use for computations
     * @param corrs, correlation matrix between genes
     * @param read, boolean indicating whether or not there is already data in the output files to read
     * @return An array of the number of times each gene was selected and each gene-gene pair was clustered
     *
     */
    private int[] getCounts(int i, float [] corrs, boolean read)
    {

        int[] curr = new int[(numGenes*(numGenes+1))/2];


        for (int x = 0; x < corrs.length; x++) {
            if (1 - corrs[x] < initRadii[i])
                curr[x]++;
        }

        /***EXPERIMENTAL***/
        if(!saveMemory) {
            try {

                File f = new File("Correlations_" + i + ".txt");
                int[] array = new int[curr.length];
                if (read) {
                    BufferedInputStream corrReader = new BufferedInputStream(new FileInputStream("Correlations_" + i +  ".txt"));
                    array = Functions.readAsIntArray(corrReader,curr.length);
                }
                for (int k = 0; k < array.length; k++)
                    array[k] += curr[k];

                BufferedOutputStream b = new BufferedOutputStream(new FileOutputStream(f, false));
                Functions.writeAsByteArray(array,b);
                array = null;
            } catch (Exception e) {
                e.printStackTrace();
                System.err.println("Unable to print correlation/intensity to file");
                System.exit(-1);
            }
        }
        return curr;
    }



    private void normalizeScore(double [] score)
    {
        double max = 0;
        for (int i = 0; i < score.length/2; i++) {
                if (score[i] > max)
                    max = score[i];
        }
        for (int i = 0; i < score.length/2; i++) {
                score[i]/=max;

        }

        for (int i = score.length/2; i < score.length; i++) {
                if (score[i] > max)
                    max = score[i];
        }
        for (int i = score.length/2; i < score.length; i++) {
                score[i]/=max;

        }


    }


    /***
     *
     * @param sums How many times did each edge appear across the subsamples radii combinations?
     * @param var mixture distribution variance estimates for each edge
     * @return Posterior variance estimate for each edge
     */
    private float[] getVarPosterior(int [] sums, float [] var)
    {
        float [] varPost = new float[var.length];
        double denom = (double)(numRadii*subsamples.length);

        for(int i = 0; i < var.length;i++)
        {
            if(var[i]<0)
                varPost[i] = -1;
            else
            {
                double p = sums[i] / denom;
                double currVar = p * (1 - p) * subsamples.length;
                varPost[i] =(float)( currVar * var[i]/(currVar+var[i]));
            }
        }

        return varPost;
    }

    /***
     *
     * @param sums Number of times each "edge" appeared across all parameter subsample combinations
     * @param genes Proportion of times the prior knowledge sources expected the edge to appear (mean/num subsamples) of the mixture distribution
     * @param vars Variance of the mixture distribution for each edge
     * @return Posterior mean for each "edge" combining data and theory estimates
     */
    private float[] getMeanPosterior(int[]sums,float[]genes,float [] vars )
    {
        float [] uPost = new float[vars.length];
        for(int i = 0; i < vars.length;i++)
        {
            if (genes[i] < 0)
                uPost[i] = -1;
            else
            {
                //Probability estimate from the data (what percent of times do we expect edge to show up
                double p = sums[i] / (double) (numRadii*subsamples.length);
                //Variance estimate for the edge from the data
                double currVar = p * (1 - p) * subsamples.length;
                uPost[i] = (float) (((sums[i] / (double) (numRadii)) * vars[i] + currVar * genes[i]*subsamples.length) / (currVar + vars[i]));

            }
        }

        return uPost;
    }

    //Overloaded method to handle computation of variance more memory efficiently for similarities!

    /***
     *
     * @param weights Weight of each prior information source (w in the paper)
     * @param tao Tao (unreliability) of each prior information source
     * @param mean w*phi for each source (this is what percent of times did the priors expect this edge to show across subsamples)
     * @param dFile Array of filenames
     * @return An array of the variance value for the mixture distribution for each edge
     */
    private float[] getVarMixture(double [] weights, double [] tao, float []mean, String [] dFile)
    {
        float[] vars = new float[numGenes*(numGenes+1)/2];
        float [] totalWeight = new float[mean.length];
        for(int j = 0; j < dFile.length;j++)
        {
            float [] prior = Functions.loadTheoryMatrix(dFile[j],false,numGenes+1);
            for(int i = 0; i < mean.length;i++)
            {
                if(prior[i]>0) {
                    vars[i] += weights[j] * (tao[j] * tao[j] + Math.pow((mean[i]*subsamples.length - prior[i]*subsamples.length), 2));
                    totalWeight[i] += weights[j];
                }
            }
        }
        for(int i= 0; i < mean.length;i++) {
            if(totalWeight[i]==0)
                vars[i] = -1;
            else
                vars[i] /= totalWeight[i];
        }
        return vars;
    }




    //This is written as a separate function due to the offset in computing the sums and due to memory requirements preventing loading all sources at once
    private double [] getSimilarityWeights(int[] sums, String [] dFile)
    {
        //Same procedure as the intensity weights except don't load all sources at once for memory concerns
        double [] tao = new double[dFile.length];
        for(int i = 0; i < dFile.length;i++)
        {
            float [][] temp = new float[1][(numGenes*(numGenes+1)/2)];
            temp[0] = Functions.loadTheoryMatrix(dFile[i],false,numGenes+1);
            whichPriors(temp,dPriors);

            getPhi(temp[0],subsamples.length);
            tao[i] = getTao(temp,sums)[0];
        }
        simTao = tao;
        double [] alpha = mgmPriors.getAlpha(tao);
        return mgmPriors.getWeights(alpha);

    }

    private double [] getIntensityWeights(int[]sums, String iFile)
    {
        //Load each intensity source from the iFile
        float [][] sources = Functions.loadIntensityData(iFile,false);

        //Identify for which genes we have prior information
        whichPriors(sources,iPriors);
        //Compute phi matrix for each source (expected number of appearences of each gene)
        for(int i = 0; i< sources.length;i++)
        {
            getPhi(sources[i],subsamples.length);
        }

        //Sum up differences between phi matrix and actual appearances (absolute differences) to compute tao
        double [] tao = getTao(sources,sums);
        intTao = tao;

        //Compute alpha and then weight based on the piMGM paper
        double [] alpha = mgmPriors.getAlpha(tao);
        return mgmPriors.getWeights(alpha);

    }

    private void whichPriors(float[][]sources,boolean [] prior)
    {
        for(int j =0; j < sources[0].length;j++) //Loop through genes
        {
            for(int i = 0; i < sources.length;i++) //loop through sources
            {
                if(sources[i][j] >= 0)
                    prior[j] = true;
            }
        }
    }


    //TODO if no intensity priors given by a source, then NA is reported, need to deal with this
    private double [] getTao(float[][]phi, int[]counts)
    {
        double [] tao = new double[phi.length];
        //For each gene, compute deviation and sum these, then divide by number of things you summed
        for(int i = 0; i < phi.length;i++) //Loop over sources
        {
            int numPriors = 0;
            for(int j = 0; j<phi[i].length;j++)
            {
                if(phi[i][j]>=0) {
                    numPriors++;
                    tao[i] += Math.abs(phi[i][j] - counts[j] / (double) (numRadii));
                }
            }
            if(numPriors == 0)
                tao[i] = -1;
            else
                tao[i] /= numPriors;
        }
        return tao;
    }

    private void getPhi(float[] in, int numSubs)
    {
        for(int i = 0; i < in.length;i++)
        {
            if(in[i]>=0.0)
                in[i] = in[i]*numSubs;
        }
    }

    /***
     *
     * @return int [] consisting of the appearence frequency of each gene and gene-gene pair across all radii,threshold combinations
     */
    private int[]computeSums()
    {
        int [] sums = new int[numGenes + (numGenes)*(numGenes-1)/2];
        for(int i = 0; i < initRadii.length;i++)
        {
                int [] curr = getCounts(i,true);
                for(int k = 0; k < curr.length;k++)
                {
                    sums[k]+=curr[k];
                }
        }
        return sums;
    }


    /***
     * This function computes scores when there is no prior information available
     * @return a 3-d array of scores where dim1 is 2 (for intensity and similarity), dim2 is numThresholds, dim3 is numRadii
     */
    private double [] computeScores()
    {
        int [] sums = new int[(numGenes)*(numGenes+1)/2];
        boolean read = false;
        for(int k = 0; k < subsamples.length;k++)
        {
            ArrayList<Gene> temp = createGenes(data,target,true);
            DataSet currData = data.subsetRows(subsamples[k]);

            //Threshold intensities and correlations after the fact
            float [] corrs = Functions.computeAllCorrelations(temp,currData,partialCorrs,false,false,1);

            for(int i = 0; i < initRadii.length;i++)
            {
                    int [] curr = getCounts(i,corrs,read);

                    for(int x = 0; x < curr.length;x++)
                    {
                        sums[x]+=curr[x];
                    }
            }
            read = true;
        }

        if(verbose)
        {
            System.out.println("Computed counts for intensities and correlations: " + Arrays.toString(sums));
        }





        float [] uPost = new float[(numGenes*(numGenes+1))/2];
        for(int i = 0; i < uPost.length;i++)
            uPost[i] = -1;
        double [] score = getScoresSeparate(sums,uPost,null);


        double [] scoreResult = new double[numRadii];

        for(int i = 0; i < numRadii;i++)
        {
                scoreResult[i] = score[i+numRadii];
        }

        return scoreResult;
    }


    /****
     *
     * @param dFile Array of files, each one having a matrix of dissimilarity information
     * @return a 3-d array of scores where dim1 is 2, dim2 is numThresholds, and dim3 is numRadii
     */
    private double[] computeScores(String[] dFile)
    {
        int [] sums = new int[(numGenes)*(numGenes+1)/2];
        boolean read = false;
        for(int k = 0; k < subsamples.length;k++)
        {
            ArrayList<Gene> temp = createGenes(data,target,true);
            DataSet currData = data.subsetRows(subsamples[k]);

            //Threshold intensities and correlations after the fact
            float [] corrs = Functions.computeAllCorrelations(temp,currData,partialCorrs,false,false,1);

            for(int i = 0; i < initRadii.length;i++)
            {

                    int [] curr = getCounts(i,corrs,read);

                    for(int x = 0; x < curr.length;x++)
                    {
                        sums[x]+=curr[x];
                    }
            }
            read = true;
        }

        if(verbose)
        {
            System.out.println("Computed counts for intensities and correlations: " + Arrays.toString(sums));
        }

        System.out.print("Computing similarity weights...");
        //Compute similarity weights for each prior knowledge source
        double [] simWeights = getSimilarityWeights(sums,dFile);
        lastSimilarityWeights = simWeights;
        System.out.println("Done");
        if(verbose)
        {
            System.out.println("Similarity weights: " + Arrays.toString(simWeights));
        }


        System.out.print("Computing scores for similarities...");
        //Load each prior knowledge source, weighted by reliability
        float [] meanDis = loadTheoryFiles(dFile,simWeights,numGenes);


        //Get variance of mixture distribution for each edge
        float [] varDis = getVarMixture(simWeights,simTao,meanDis,dFile);



        //Get posterior distributions of predicted dissimilarity and intensity from prior sources
        float[] uPost = getMeanPosterior(sums,meanDis, varDis);
        float [] varPost = getVarPosterior(sums, varDis);

        //Compute posterior probabilities and stability of each gene, synthesize these results into a score for each parameter setting
        double[] scoresSim;
        scoresSim = getScoresSeparate(sums,uPost,varPost);
        System.out.println("Done");

        return scoresSim;
    }



    private synchronized void addGeneConnections(float [] corrs, int [] curr, double radii)
    {
        for(int i = 0; i < corrs.length;i++)
        {
            if(1-corrs[i]<radii) //TODO Is this correct?
                curr[i]++;

        }
    }


    public static synchronized ArrayList<Gene> createGenes(DataSet data,String target, boolean includeTarget)
    {
        ArrayList<Gene> gList = new ArrayList<Gene>();
        int ID = 0;
        for(int i = 0; i < data.getNumColumns();i++)
        {
            if(!includeTarget && data.getVariable(i).getName().equals(target))
                continue;
            Gene g = new Gene(ID);
            ID++;
            g.symbol = data.getVariable(i).getName();
            gList.add(g);
        }
        return gList;
    }

    //Give an array of filenames for all theory files, load each one into a row of the matrix to get all theory information in a single matrix
    //This computes the phi function from the paper
    public static float[] loadTheoryFiles(String [] dFile, double [] weights, int numGenes)
    {
        float [] result = new float[numGenes*(numGenes+1)/2];
        float [] totalWeight = new float[result.length];
        //Load in each theory source
        for(int i = 0; i < dFile.length;i++)
        {
            float [] temp = Functions.loadTheoryMatrix(dFile[i],false,numGenes+1);

            //For each "edge" or relationship
            for(int j = 0; j < temp.length;j++)
            {
                //If the source gives information about this relationship
                if(temp[j]>0)
                {
                    result[j] += temp[j]*weights[i];
                    totalWeight[j] += weights[i];
                }
            }

        }
        //Normalize for edges where not all sources give information
        for(int i = 0; i < result.length;i++) {
            if(totalWeight[i]==0)
                result[i] = -1;
            else
                result[i] /= totalWeight[i];
        }
        return result;

    }

    //Generate subsampled indices for the given dataset
    //Input: bootstrap ( should we bootstrap to generate samples or should we use subsampling)
    //Input: numSamples (number of samples to test)
    //Output: an array of integers specifying the samples that will be used in the repeated Pref-Div runs
    public static int [][] genSubsamples(boolean bootstrap, int numSamples,DataSet data,boolean LOOCV)
    {
        if(bootstrap) {
            return DataUtils.getBootstrapIndices(data,data.getNumRows(),numSamples);
        }
        else if(LOOCV)
        {
            return StabilityUtils.generateSubsamples(data.getNumRows());
        }
        else
        {
            return StabilityUtils.subSampleNoReplacement(data.getNumRows(),numSamples);
        }
    }



}
