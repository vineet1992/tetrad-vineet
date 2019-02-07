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
 * Automatically uses separate params, and uses memory efficient methodology
 * Need to parallelize correlation testing and correlation computing
 */

//TODO Scores seem to be completely invariant across threshold values, regardless of which radius is chosen


public class PiPrefDiv3 {

    private DataSet data; //The current expression dataset to be analyzed, assumes that all variables are fair game for Pref-Div except the target
    private String target; //The target variable of interest
    private double [] initRadii; //Initial radius values to test
    private double lowRadii = 0; //Low range of radius values to test
    private double highRadii = 1; //High range of radius values to test
    private int numRadii = 100; //Number of radius values to test

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
    private double [] lastIntensityWeights;
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


    //Constructor that uses default parameters for both initial parameter ranges
    public PiPrefDiv3(DataSet data, String target,int K)
    {
        this.K = K;
        this.data = data;
        this.numGenes = data.getNumColumns()-1;
        this.target = target;
        this.initRadii = initRadiiDefault();
    }


    //Constructor that specifies the same number of parameters for radius and k
    public PiPrefDiv3(DataSet data, String target, int K,int numParams)
    {
        this.K = K;
        this.numGenes = data.getNumColumns()-1;
        this.data = data;
        this.target = target;
        this.numRadii = numParams;
        this.initRadii = initRadiiDefault();
        this.numFolds = 5;
    }


    //Construct with number of radii and k parameters to test
    public PiPrefDiv3(DataSet data, String target,int K, int numRadii, int numThreshold)
    {
        this.K = K;
        this.data = data;
        this.numGenes = data.getNumColumns()-1;
        this.target = target;
        this.numRadii = numRadii;
        this.initRadii = initRadiiDefault();
    }


    //Constructor with initial ranges both specified
    public PiPrefDiv3(DataSet data, String target, int K, double [] radii, double [] threshold)
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
    public double[] getLastIntensityWeights(){return lastIntensityWeights;}
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
     * @return double [] with a constricted range of the parameter
     */

    //TODO instead of running this on the full dataset, we want the average # across subsamples
    private double[] constrictRange(double[]init, DataSet data)
    {
        int[] num = new int[init.length];
        for(int x = 0; x < subsamples.length; x++) {

            DataSet curr = data.subsetRows(subsamples[x]);
            ArrayList<Gene> temp = createGenes(data,target);
            Gene tgt = new Gene(temp.size());
            tgt.symbol="Target";
            temp.add(tgt);

            float[] corrs = Functions.computeCorrelationsOnly(temp, curr,false,false);

            for (int i = 0; i < init.length; i++) {
                    num[i] += countLessThan(corrs, init[i]) / (double)init.length;

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
        initRadii = constrictRange(initRadii,data);
        System.out.println("Done");

        //Compute stability of each gene selection and each gene-gene relationship
        System.out.print("Computing stability across radii values...");
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


    /***Full procedure to select genes based on prior information Pref-Div
     //Input: boot (should we do bootstrapping (true) or subsampling (false))
     //Input: numSamples (number of samples to subsample or bootstrap)
     //Input: path to theory Intensity File -> iFile (should have a header)
     //Input: path to theory Dissimilarity File -> dFile (no header or rownames for these files)****/
    public ArrayList<Gene> selectGenes(boolean boot, int numSamples, String iFile, String [] dFile,boolean useCausalGraph)
    {

        iPriors = new boolean[numGenes];
        dPriors = new boolean[numGenes * (numGenes - 1) / 2];

        //Generate subsamples of the data
        System.out.print("Generating subsamples...");
        if(subsamples==null)
        {
            subsamples = genSubsamples(boot,numSamples,data,LOOCV);
        }
        System.out.println("Done");

        //Constrict parameter ranges
        System.out.print("Constricting Parameter Range...");
        initRadii = constrictRange(initRadii,data);

        System.out.println("Done");

        //Compute stability of each gene selection and each gene-gene relationship
        System.out.print("Computing stability across radii and threshold values...");


        double[][] scores = computeScores(iFile,dFile);


        System.out.print("Merging scores to get top parameters...");
        //Find best radii and topK values and then Run Pref-Div with cross-validation to get final Output!!
        double [] avgScore = getAvgScore(scores[1],scores[0]);

        if(verbose)
        {
            try {

                System.out.println("Printing all scores...");
                for (int i = 0; i < initRadii.length*2; i++) {
                    System.out.print(initRadii[i % initRadii.length] + "\t");
                }

                for (int i = 0; i < initRadii.length*2; i++) {
                    System.out.print(avgScore[i] + "\t");
                }
                System.out.println();


                if(outputScores)
                {
                    for (int i = 0; i < initRadii.length*2; i++) {
                        scoreStream.print(initRadii[i % initRadii.length] + "\t");
                    }

                    scoreStream.println();
                    for (int i = 0; i < initRadii.length*2; i++) {
                        scoreStream.print(avgScore[i] + "\t");
                    }
                    scoreStream.flush();
                }
            }catch(Exception e)
            {
                System.err.println("Couldn't write scores to file");
            }
        }


        int [] inds = getMaxScoreSeparate(avgScore);



        double bestRadii = initRadii[inds[0]];
        double bestRadiiNP = initRadii[inds[1]];
        lastRadius = bestRadiiNP;
        lastRadiusWP = bestRadii;
        System.out.println("Done");

        if(verbose)
        {
            System.out.println("Best Radii: " + bestRadii);
        }




        avgScore = null;


        System.out.print("Running Pref-Div with optimal parameters...");
        //Run Pref-Div with optimal parameters


        ArrayList<Gene> temp = createGenes(data,target);
        //TODO Make this a funciton that doesn't use a p-value threshold at all
        ArrayList<Gene> meanGenes = Functions.computeAllIntensities(temp,1,data,target,false,1,1,iPriors);
        Collections.sort(meanGenes,Gene.IntensityComparator);

        //TODO this one too
        float [] meanDis = Functions.computeAllCorrelations(meanGenes,data,false,1,1,dPriors);

        RunPrefDiv rpd;
        rpd = new RunPrefDiv(meanDis,meanGenes,data,target,LOOCV);
        rpd.setWithPrior(dPriors);
        rpd.setAllParams(bestRadiiNP,bestRadii,1,1);

        rpd.setTopK(K);
        rpd.setAccuracy(0);
        rpd.setNS(subsamples.length);
        rpd.setNumFolds(numFolds);
        rpd.setCausalGraph(useCausalGraph);
        rpd.useStabilitySelection(useStabilitySelection);
        rpd.setClusterType(ctype);
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


    /***
     *
     * @param avg An array of the average (from intensity/similarity) for each parameter (first half is WP, second half is NP)
     * @return int [] where index 0 is the index of the maximum no prior score, and index 1 is the index of the maximum with prior score
     */
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


    /***
     *
     * @param s1 Array of intensity scores
     * @param s2 Array of similarity scores
     * @return An array of the average score for each parameter
     */
    private double[] getAvgScore(double [] s1, double [] s2)
    {
        double [] result = new double[s1.length];
        for(int i = 0; i < s1.length;i++)
        {
            result[i] = (s1[i] + s2[i])/2;
        }
        return result;
    }

    //Stability of showing up n times across subsamples.length subsamples
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

    private float[] convertListToMatrix(ArrayList<Gene> genes)
    {
        float [] result = new float[genes.size()];
        for(int i = 0; i < genes.size();i++) {
            result[i] = (float) genes.get(i).theoryIntensity;
        }
        return result;
    }


    /**
     * @param sums, Array of number of times each gene and gene-gene pair appeared
     * @param uPost, Posterior mean for each gene OR each gene-gene pair appearence
     * @param varPost, Posterior variance for each gene OR each gene-gene pair appearence
     * @return Scores for each radii param. Either for intensity or dissimilarity depending upon input (first half is WP, second half is NP)
     */
    private double[] getScoresSeparate(int [] sums, float [] uPost, float [] varPost)
    {
        boolean needCorrs = false;
        double [] score = new double[numRadii*2];
        int offset = 0;
        if(uPost.length>numGenes) {
            offset = numGenes;
            needCorrs = true;
        }

        for (int i = 0; i < numRadii; i++) {

                int [] curr = null;
                if(saveMemory) {

                    curr = getCounts(i, needCorrs);

                }
                else {
                    curr = getCountsFromFile(i, needCorrs);

                }


                for(int g = 0; g < uPost.length;g++) //Loop through each gene
                {
                    double G = getG(curr[g+offset]); //Clusts[i][j][g] is number of times gene was selected in the subsamples for this parameter setting
                    if (uPost[g] < 0)//No Prior information, contributes only stability to the score
                    {
                        double theta = getTheta(sums[g+offset],curr[g+offset]); //sums is the total number of times this gene was selected across all parameter settings
                        score[i + numRadii] += theta*(1-G);
                    }
                    else
                    {
                        if(varPost[g]<0)
                            System.out.println(g);
                        double theta = getTheta(uPost[g],varPost[g],curr[g+offset]);
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
        int [] curr = new int[numGenes + (numGenes*(numGenes-1))/2];

        /***Loop through subsamples and compute number of times each gene/gene-gene pair are selected***/
        for(int k = 0; k < subsamples.length;k++)
        {
            ArrayList<Gene> temp = createGenes(data,target);
            DataSet currData = data.subsetRows(subsamples[k]);


            //TODO Can replace both all intensities and all correlations methods to improve speed since threshold is unnecessary
            temp = Functions.computeAllIntensities(temp,1,currData,target,partialCorrs,false,false,1);

            if (needCorrs || pdStability) {
                float [] corrs = Functions.computeAllCorrelations(temp,currData,partialCorrs,false,false,1);


                if(pdStability) {
                    Collections.sort(temp,Gene.IntensityComparator);
                    PrefDiv pd = new PrefDiv(temp, K, 0, initRadii[i], corrs);
                    pd.setCluster(clusterByCorrs);
                    ArrayList<Gene> topGenes = pd.diverset();
                    addFoundGenes(topGenes, curr,initRadii[i]);
                }
                else
                {
                    addFoundGenes(temp,curr,initRadii[i]);
                }
                addGeneConnections(corrs,curr,initRadii[i]);
            }
            else
            {
                addFoundGenes(temp,curr,initRadii[i]);
            }
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
        int fullSize = numGenes + (numGenes*(numGenes-1)/2);
        int [] curr = new int[fullSize];

        try{
            File iFile =new File("Correlations_" + i + ".txt");
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
     * @param j, index of initThreshold to use for computation
     * @param temp, List of genes for computation
     * @param corrs, correlation matrix between genes
     * @return An array of the number of times each gene was selected and each gene-gene pair was clustered
     *
     */
    private int[] getCounts(int i, ArrayList<Gene> temp, float [] corrs, boolean read)
    {

        int[] curr;

        /***Create curr which holds the number of times this gene/edge was selected across subsamples***/
        if(corrs==null)
            curr = new int[numGenes];
        else
            curr = new int[numGenes + (numGenes*(numGenes-1))/2];


        /***Did the current gene meet the specified radii requirements***/
        for(int x = 0; x < temp.size();x++)
        {
            if(1-temp.get(x).intensityValue < initRadii[i])
                curr[x]++;
        }


        /***If we're using correlations, did the current correlation meet the specified radii requirement?***/
        if(corrs!=null) {
            for (int x = 0; x < corrs.length; x++) {
                if (1 - corrs[x] < initRadii[i])
                    curr[x + temp.size()]++;
            }
        }

        /***EXPERIMENTAL save memory by writing out curr to file and re-reading it when computing this subsample with another radii value***/
        if(!saveMemory) {
            try {

                File f = new File("Correlations_" + i + ".txt");
                int[] array = new int[curr.length];
                if (read) {
                    BufferedInputStream corrReader = new BufferedInputStream(new FileInputStream("Correlations_" + i + ".txt"));
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


    /**
     *
     * @param score Array of scores (first half are WP and second half are NP)
     */
    private void normalizeScore(double [] score)
    {

        /***Normalize WP scores (first half of array)***/


        /**Identify the max value**/
        double max = 0;
        for (int i = 0; i < score.length/2; i++) {

                if (score[i] > max)
                    max = score[i];
        }

        /***Divide all by the max value**/
        for (int i = 0; i < score.length/2; i++) {

                score[i]/=max;

        }

        /***Normalize NP scores***/

        /***Identify the maximum value***/
        for (int i = score.length/2; i < score.length; i++) {
                if (score[i] > max)
                    max = score[i];

        }

        /***Divide each score by the maximum value***/
        for (int i = score.length/2; i < score.length; i++) {
                score[i]/=max;

        }


    }

    //Get posterior variance based on computed counts and variance, modified to work for both intensity and similarity
    private float[] getVarPosterior(int [] sums, float [] var)
    {
        float [] varPost = new float[var.length];
        int offset = 0;
        double denom = (double)(numRadii*subsamples.length);
        if(var.length>numGenes)
            offset = numGenes;
        for(int i = 0; i < var.length;i++)
        {
            if(var[i]<0)
                varPost[i] = -1;
            else
            {
                double p = sums[i + offset] / denom;
                double currVar = p * (1 - p) * subsamples.length;
                varPost[i] =(float)( currVar * var[i]/(currVar+var[i]));
            }
        }

        return varPost;
    }

    //Get Posterior mean for intensity values from theory
    //Clusts is the number of times each gene and gene-gene connection appeared in radii i, threshold j, and subsample k
    //genes is the mixture mean predicted by the prior information sources
    //vars is the mixture variance from the prior information sources
    private float[] getMeanPosterior(int[]sums,float[]genes,float [] vars )
    {
        float [] uPost = new float[vars.length];
        int offset = 0; //For intensity posterior means
        if(uPost.length>numGenes) //Then we are doing similarity posterior means
            offset = numGenes;
        for(int i = 0; i < vars.length;i++)
        {
            if (genes[i] < 0)
                uPost[i] = -1;
            else
            {
                double p = sums[i+offset] / (double) (numRadii * subsamples.length);
                double currVar = p * (1 - p) * subsamples.length;
                uPost[i] = (float) (((sums[i+offset] / (double) (numRadii)) * vars[i] + currVar * genes[i]) / (currVar + vars[i]));

            }
        }

        return uPost;
    }

    //Overloaded method to handle computation of variance more memory efficiently for similarities!
    private float[] getVarMixture(double [] weights, double [] tao, float []mean, String [] dFile)
    {
        float[] vars = new float[numGenes*(numGenes-1)/2];
        float [] totalWeight = new float[mean.length];
        for(int j = 0; j < dFile.length;j++)
        {
            float [] prior = Functions.loadTheoryMatrix(dFile[j],false,numGenes);
            for(int i = 0; i < mean.length;i++)
            {
                if(prior[i]>0) {
                    vars[i] += weights[j] * (tao[j] * tao[j] + Math.pow((mean[i] - prior[i]) * subsamples.length, 2));
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
    //Compute the Variance mixture of the theory sources that provide intensity information for each gene
    private float [] getVarMixture(double [] weights, double [] tao, float[] genes, String iFile)
    {
        float[][] priors = Functions.loadIntensityData(iFile,false);
        float[] vars = new float[numGenes];

        for(int i = 0; i < genes.length;i++)
        {
            double sum = 0;
            boolean found = false;
            for(int j = 0; j < priors.length;j++)
            {
                if(priors[j][i]>=0)
                    found = true;
            }
            if(found)
            {
                double weight = 0;
                for(int j = 0; j < priors.length;j++)
                {
                    if(priors[j][i]>=0) {
                        sum += weights[j] * (tao[j] * tao[j] + Math.pow((genes[i] - priors[j][i]) * subsamples.length, 2));
                        weight += weights[j];
                    }
                }
                vars[i] = (float)(sum/weight);

            }
            else
                vars[i] = -1;
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
            float [][] temp = new float[1][(numGenes*(numGenes-1)/2)];
            temp[0] = Functions.loadTheoryMatrix(dFile[i],false,numGenes);
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

    //Load Gene Data from the specified intensity file
    public static ArrayList<Gene> loadGenes(String iFile,double[]iWeights)
    {
        return Functions.loadGeneData(iFile,false,iWeights);
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



    private double[][] computeScores(String iFile, String[] dFile)
    {
        int [] sums = new int[numGenes + (numGenes)*(numGenes-1)/2];
        boolean read = false;
        for(int k = 0; k < subsamples.length;k++)
        {
            ArrayList<Gene> temp = createGenes(data,target);
            DataSet currData = data.subsetRows(subsamples[k]);

            //TODO This method unnecessarily computes p-values too can make faster by changing
            temp = Functions.computeAllIntensitiesWithP(temp,1,currData,target,partialCorrs);


            //Threshold intensities and correlations after the fact
            float [] corrs = Functions.computeAllCorrelations(temp,currData,partialCorrs,false,false,1);

            for(int i = 0; i < initRadii.length;i++)
            {
                    int [] curr = getCounts(i,temp,corrs,read);

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



        System.out.print("Computing intensity weights...");
        //Compute intensity weights
        double [] intWeights = getIntensityWeights(sums,iFile);
        lastIntensityWeights = intWeights;
        System.out.println("Done");
        if(verbose)
        {
            //TODO get header from iFile to attach source to weights
            System.out.println("Intensity weights are: " + Arrays.toString(intWeights));
        }

        //Get Mixture Distributions of predicted dissimilarity and intensity from the prior sources

        System.out.print("Loading genes with weighted intensities...");
        //Load theory intensity file -> for mean intensity from theory
        ArrayList<Gene> meanGenes = loadGenes(iFile,intWeights);
        System.out.println("Done");
        if(verbose)
        {
            System.out.println("All Genes: " + meanGenes);
        }

        float [] means = convertListToMatrix(meanGenes);
        System.out.println(Arrays.toString(means));

        System.out.print("Computing scores for each radii, k combination from intensities...");
        //Compute variance of mixture distribution of prior knowledge sources
        float [] varIntense = getVarMixture(intWeights, intTao, means,iFile);

        //Get posterior distributions of predicted dissimilarity and intensity from prior sources
        float[] uPostInt = getMeanPosterior(sums, means, varIntense);
        float [] varPostInt = getVarPosterior(sums, varIntense);

        System.out.println(Arrays.toString(uPostInt));
        System.out.println(Arrays.toString(varPostInt));
        //Compute posterior probabilities and stability of each gene, synthesize these results into a score for each parameter setting


        //First half is no prior information, second half is with prior information
        double [] scoresInt = getScoresSeparate(sums,uPostInt,varPostInt);

        System.out.println("Done");



        if(verbose)
        {
            try {

                System.out.println("Printing all scores...");
                for (int i = 0; i < initRadii.length*2; i++) {
                    System.out.print(initRadii[i % initRadii.length] + "\t");
                }

                for (int i = 0; i < initRadii.length*2; i++) {
                    System.out.print(scoresInt[i] + "\t");
                }
                    System.out.println();


                if(outputScores)
                {
                    for (int i = 0; i < initRadii.length*2; i++) {
                        scoreStream.print(initRadii[i % initRadii.length] + "\t");
                    }

                    scoreStream.println();
                    for (int i = 0; i < initRadii.length*2; i++) {
                        scoreStream.print(scoresInt[i] + "\t");
                    }
                    scoreStream.flush();
                }
            }catch(Exception e)
            {
                System.err.println("Couldn't write scores to file");
            }
        }

        //Set not needed variables to null to save memory
        uPostInt = null;
        varPostInt = null;
        varIntense = null;
        intWeights = null;
        intTao = null;




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


        if(verbose)
        {
            try {

                System.out.println("Printing all scores...");
                for (int i = 0; i < initRadii.length*2; i++) {
                    System.out.print(initRadii[i % initRadii.length] + "\t");
                }

                for (int i = 0; i < initRadii.length*2; i++) {
                    System.out.print(scoresSim[i] + "\t");
                }
                System.out.println();


                if(outputScores)
                {
                    for (int i = 0; i < initRadii.length*2; i++) {
                        scoreStream.print(initRadii[i % initRadii.length] + "\t");
                    }

                    scoreStream.println();
                    for (int i = 0; i < initRadii.length*2; i++) {
                        scoreStream.print(scoresSim[i] + "\t");
                    }
                    scoreStream.flush();
                }
            }catch(Exception e)
            {
                System.err.println("Couldn't write scores to file");
            }
        }

        /***Create an array to store all scores (intensity and similarity)***/
        double[][]scores = new double[2][scoresInt.length];
        scores[0] = scoresInt;
        scores[1] = scoresSim;

        return scores;
    }

    private synchronized void addFoundGenes(ArrayList<Gene> top, int [] curr, double radii) {
        for (int i = 0; i < top.size(); i++)
        {
            if(1-top.get(i).intensityValue<radii)
                curr[top.get(i).ID]++;
        }

    }

    private synchronized void addGeneConnections(float [] corrs, int [] curr, double radii)
    {
        for(int i = 0; i < numGenes;i++)
        {
            for (int j = i+1; j < numGenes; j++) {
                int temp = Functions.getIndex(i,j,numGenes);
                if(1-corrs[temp]<radii) //TODO Is this correct?
                    curr[temp + numGenes]++;

            }
        }
    }


    public static synchronized ArrayList<Gene> createGenes(DataSet data, String target)
    {
        ArrayList<Gene> gList = new ArrayList<Gene>();
        int ID = 0;
        for(int i = 0; i < data.getNumColumns();i++)
        {
            if(data.getVariable(i).getName().equals(target))
                continue;
            Gene g = new Gene(ID);
            ID++;
            g.symbol = data.getVariable(i).getName();
            gList.add(g);
        }
        return gList;
    }

    //Give an array of filenames for all theory files, load each one into a row of the matrix to get all theory information in a single matrix
    public static float[] loadTheoryFiles(String [] dFile, double [] weights, int numGenes)
    {
        float [] result = new float[numGenes*(numGenes-1)/2];
        float [] totalWeight = new float[result.length];
        for(int i = 0; i < dFile.length;i++)
        {
            float [] temp = Functions.loadTheoryMatrix(dFile[i],false,numGenes);
            for(int j = 0; j < temp.length;j++)
            {
                if(temp[j]>0)
                {
                    result[j] += temp[j]*weights[i];
                    totalWeight[j] += weights[i];
                }
            }

        }
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