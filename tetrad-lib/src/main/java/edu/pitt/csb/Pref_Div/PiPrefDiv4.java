package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.Pref_Div.Comparisons.ComparablePD;
import edu.pitt.csb.Priors.mgmPriors;
import edu.pitt.csb.stability.StabilityUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.util.CombinatoricsUtils;

import java.io.*;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

/**
 * Created by vinee_000 on 10/19/2018.
 * This version uses a single file for each prior information source using both similarity and intensity
 * This version will also use only radius of similarity instead of both that and p-value (because some p-values are 0, making them impossible to discount by using prior information)
 */



public class PiPrefDiv4 implements ComparablePD {


    private int numSamples = 1000; //Number of null distribution samples to generate
    private DataSet data; //The current expression dataset to be analyzed, assumes that all variables are fair game for Pref-Div except the target
    private String target; //The target variable of interest
    private double [] initRadii; //Initial radius values to test
    private double lowRadii = 0.3; //Low range of radius values to test
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
    private DataSet summarizedData; //Summarized dataset using clusters as variables
    private boolean parallel; //Should we compute stabilities in parallel?
    private RunPrefDiv.ClusterType ctype; //Which clustering method should be used to produce summarized data?
    private boolean partialCorrs; //Should we use partial correlations?
    private double lastRadiusWP; //Last Radius WP selected
    private boolean [] iPriors; //Contains which features have prior intensity information
    private boolean [] dPriors; //Contains which features have prior dissimiliarty information
    private boolean [] dPriorsWT; // Contains which features have prior dissimilarity information (excludes the target variable)
    private boolean saveMemory = false; //Should we conserve disc space and use the slower version of PiPrefDiv
    private List<Gene> lastSelected; //Last set of selected genes
    private double wpCutoff = 0.5; //Cutoff for inclusion into the withprior set
    private Random rand = new Random(42);
    private boolean quiet = false; //Suppress all output
    private boolean computeP = false; //Should p-values be computed?



    private int targetIndex = 0; //Index of the target variable
    private double [] pvals; //P-value for each prior information source
    private double [] adjustP; //Adjusted p-value for each prior information source
    private double [] normTao; //Tao value normalized by the mean of the null distribution


    /***TODO Fix up and remove this when done***/

    public double [] constrictCorrs;


    private boolean profiling = false;

    private double lastNormTao = -1;


    //Constructor that uses default parameters for both initial parameter ranges
    public PiPrefDiv4(DataSet data, String target,int K)
    {
        this.K = K;
        this.data = data;
        this.numGenes = data.getNumColumns()-1;
        this.target = target;
        this.initRadii = initRadiiDefault();
        this.targetIndex = data.getColumn(data.getVariable(target));
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
        this.targetIndex = data.getColumn(data.getVariable(target));

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
        this.targetIndex = data.getColumn(data.getVariable(target));

    }


    public void setQuiet(){this.quiet = true;}
    public void setSubsamples(int[][] subs){this.subsamples = subs;}
    public Map<Gene,List<Gene>> getLastCluster()
    {
        return lastCluster;
    }
    public double[] getLastWeights(){return lastSimilarityWeights;}
    public double[] getLastTao(){return simTao;}
    public double[] getLastRadius(){return new double[]{lastRadius,lastRadiusWP};}
    public double getLastRadiusWP(){return lastRadiusWP;}
    public void setOutputScores(PrintStream out){outputScores=true; scoreStream=out;}
    public void setParallel(boolean b){parallel = b;}
    public void setClusterType(RunPrefDiv.ClusterType c){ctype = c;}
    public DataSet getSummarizedData(){return summarizedData;}
    public void setPartialCorrs(boolean pc){partialCorrs = pc;}
    public void saveMemory(){saveMemory = true;}
    public boolean[] getDissimilarityPriors(){return dPriorsWT;}
    public List<Gene> getLastSelected(){return lastSelected;}
    public void setWPCutoff(double wp){wpCutoff = wp;}
    public double [] getPValues(){return pvals;}
    public double [] getAdjustP(){return adjustP;}
    public void computePValue(boolean pvals){computeP = pvals;}
    public void profile(){profiling = true;}
    public double [] getNormTao()
    {
        if(normTao==null)
            normTao = new double[getPValues().length];
        return normTao;
    }


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

    private double[] constrictRange(double[]init, DataSet data,boolean threshold)
    {
        int[] num = new int[init.length];

        for(int x = 0; x < subsamples.length; x++) {

            DataSet curr = data.subsetRows(subsamples[x]);
            ArrayList<Gene> temp = createGenes(data,target,false);
            Gene tgt = new Gene(temp.size());
            tgt.symbol=target;
            temp.add(tgt);



            float [] corrs;
            if(parallel)
            {
                corrs = Functions.computeAllCorrelationsPar(temp,curr);
            }else
            {
                corrs = Functions.computeAllCorrelations(temp, curr, partialCorrs, false, threshold, 1);
            }


            for (int i = 0; i < init.length; i++) {
                if (!threshold)
                    num[i] += countLessThan(corrs, init[i]) / (double)subsamples.length;
                else
                    num[i] += countPLessThan(corrs, init[i])/(double)subsamples.length;
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
     * @param radiiNP Radius of similarity for relationships with no prior knowledge
     * @return List of selected genes based on this radius
     */
    public ArrayList<Gene> selectGenes(double radiiNP)
    {

        /***Create genes list based upon the order of the columns in the dataset***/
        ArrayList<Gene> temp = createGenes(data,target,false);


        long time = System.nanoTime();
        /***Add correlation to the target variable in the intensity value and fold change fields of each gene***/
        ArrayList<Gene> meanGenes = Functions.computeAllIntensities(temp,1,data,target,partialCorrs,false,false,1);


        /***Order of correlation is based upon the order of the genes in the array list (needs to be the same as the data order)***/
        float [] meanDis = Functions.computeAllCorrelations(meanGenes,data,partialCorrs,false,false,1);

        if(profiling)
        {
            time = System.nanoTime()-time;
            System.out.println("Computing Correlations/Intensities to select; " + time/Math.pow(10,9));
        }



        time = System.nanoTime();
        /***Shuffle and sort the list of genes***/
        Collections.shuffle(meanGenes,rand);
        Collections.sort(meanGenes,Gene.IntensityComparator);


        /***Create RunPrefDiv wrapper object and set all parameters***/
        RunPrefDiv rpd;
        rpd = new RunPrefDiv(meanDis,meanGenes,data,LOOCV);

        rpd.setTopK(K);
        rpd.setAccuracy(0);
        rpd.setNS(subsamples.length);
        rpd.setClusterType(ctype);
        if(clusterByCorrs)
            rpd.clusterByCorrs();
        rpd.setRadius(radiiNP);
        rpd.setSubs(subsamples);


        /***Run Pref-Div and get selected features and clusters***/
        ArrayList<Gene> top = rpd.runPD(target);
        Map<Gene,List<Gene>> map = rpd.getClusters();
        summarizedData = rpd.getSummarizedData();
        lastCluster = map;

        /*if(verbose)
        {
            System.out.println("Selected top K features\n" + top);
            System.out.println("All clusters\n" + map);
        }*/
        lastSelected = top;
        if(profiling)
        {
            time = System.nanoTime()-time;
            System.out.println("Shuffling, sorting, Pref-Div; " + time/Math.pow(10,9));
        }

        return top;
    }


    /***
     *
     * @param radiiNP Radius for relationships that do not have prior knowledge
     * @param radiiWP Radius for relationships with prior knowledge
     * @param dFile Array of filenames for similarity matrices
     * @return List of selected genes
     */
    public ArrayList<Gene> selectGenes(double radiiNP, double radiiWP, String [] dFile)
    {

        /****Create gene list based upon the order of columns in the data***/
        ArrayList<Gene> temp = createGenes(data,target,false);



        /***Load which genes have prior info about relationship to target***/
        if(iPriors==null)
            iPriors = loadIPrior(dFile);


        /***Load which gene-gene relationships have prior information (based upon the order in the dataset)***/
        if(dPriorsWT==null)
            dPriorsWT = getPriorNoTarget();


        /***Give intensity values to each gene in "temp"***/
        ArrayList<Gene> meanGenes = Functions.computeAllIntensities(temp,1,data,target,false,1,1,iPriors);


        /***Compute correlation for each pair of genes***/
        float [] meanDis = Functions.computeAllCorrelations(meanGenes,data,false,1,1,dPriorsWT);


        /***Shuffle and sort the list of genes based upon intensity value***/
        Collections.shuffle(meanGenes,rand);

        Collections.sort(meanGenes,Gene.IntensityComparator);

        /***If you aren't within radius of the target then you are shrunk to zero***/
        for(int i = 0; i < meanGenes.size();i++)
        {
            if(iPriors[meanGenes.get(i).ID] && 1-meanGenes.get(i).intensityValue>radiiWP)
                meanGenes.get(i).intensityValue = 0;
            else if(!iPriors[meanGenes.get(i).ID] && 1-meanGenes.get(i).intensityValue>radiiNP)
                meanGenes.get(i).intensityValue=0;
        }




        Collections.sort(meanGenes,Gene.IntensityComparator);



        /***TODO Temp prior is in the order of the dataset, but meanDis is in the order of the gene list***/


        RunPrefDiv rpd;
        rpd = new RunPrefDiv(meanDis,meanGenes,data,LOOCV);


        rpd.setWithPrior(dPriorsWT);


        rpd.setAllParams(radiiNP,radiiWP);

        rpd.setTopK(K);
        rpd.setAccuracy(0);
        rpd.setNS(subsamples.length);
        rpd.setClusterType(ctype);
        rpd.setSubs(subsamples);
        if(clusterByCorrs)
            rpd.clusterByCorrs();

        ArrayList<Gene> top = rpd.runPD(target);
        Map<Gene,List<Gene>> map = rpd.getClusters();
        summarizedData = rpd.getSummarizedData();
        lastCluster = map;

        /*if(verbose)
        {
            System.out.println("Selected top K features\n" + top);
            System.out.println("All clusters\n" + map);
        }*/
        lastSelected = top;
        return top;
    }

    /***
     *
     * @param boot Should we bootstrap instead of subsampling?
     * @param numSamples How many subsamples should we draw?
     * @return
     */
    public ArrayList<Gene> selectGenes(boolean boot, int numSamples)
    {
        long time = System.nanoTime();
        if(data.isContinuous()) {
            if(!quiet)
                System.out.print("Standardizing Data...");
            data = DataUtils.standardizeData(data);

            if(!quiet)
                System.out.println("Done");
        }
        if(profiling)
        {
            time = System.nanoTime()-time;
            System.out.println("Standardizing Data: " + time/Math.pow(10,9));
        }

        //Generate subsamples of the data
        if(!quiet)
            System.out.print("Generating subsamples...");
        if(subsamples==null)
        {
            subsamples = genSubsamples(boot,numSamples,data,LOOCV);
        }
        if(!quiet)
            System.out.println("Done");


        time = System.nanoTime();
        //Constrict parameter ranges
        if(!quiet)
         System.out.print("Constricting Parameter Range...");
        initRadii = constrictRange(initRadii,data,false);

        if(profiling)
        {
            time = System.nanoTime()-time;
            System.out.println("Constricting Range; " + time/Math.pow(10,9));
        }

        if(!quiet)
            System.out.println("Done");

        if(verbose)
        {
            System.out.println("Selected radii range: (" + initRadii[0] + "," + initRadii[initRadii.length-1] + ")");

        }


        if(!quiet)
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


        if(!quiet)
            System.out.println("Done");





        //Run Pref-Div with optimal parameters

        return selectGenes(lastRadius);

    }



    /***Full procedure to select genes based on prior information Pref-Div
     //Input: boot (should we do bootstrapping (true) or subsampling (false))
     //Input: numSamples (number of samples to subsample or bootstrap)
     //Input: path to theory Intensity File -> iFile (should have a header)
     //Input: path to theory Dissimilarity File -> dFile (no header or rownames for these files)****/
    public ArrayList<Gene> selectGenes(boolean boot, int numSamples, String [] dFile)
    {
        long time = System.nanoTime();
        if(data.isContinuous()) {
            if(!quiet)
                System.out.print("Standardizing Data...");
            data = DataUtils.standardizeData(data);

            if(!quiet)
                System.out.println("Done");
        }

        if(profiling)
        {
            time = System.nanoTime()-time;
            System.out.println("Standardizing Data: " + time/Math.pow(10,9));
        }

        dPriors = new boolean[numGenes * (numGenes + 1) / 2];

        //Generate subsamples of the data
        if(!quiet)
            System.out.print("Generating subsamples...");
        if(subsamples==null)
        {
            subsamples = genSubsamples(boot,numSamples,data,LOOCV);
        }

        if(!quiet)
            System.out.println("Done");

        time = System.nanoTime();

        //Constrict parameter ranges
        if(!quiet)
            System.out.print("Constricting Parameter Range...");
        initRadii = constrictRange(initRadii,data,false);

        if(profiling)
        {
            time = System.nanoTime()-time;
            System.out.println("Constricting Range; " + time/Math.pow(10,9));
        }


        if(!quiet)
            System.out.println("Done");

        if(verbose)
        {
            System.out.println("Selected radii range: (" + initRadii[0] + "," + initRadii[initRadii.length-1] + ")");
        }

        if(!quiet)
            //Compute stability of each gene selection and each gene-gene relationship
        System.out.print("Computing stability across radii and threshold values...");


        double[] scores = computeScores(dFile);

        int [] inds = getMaxScoreSeparate(scores);

         double bestRadii = initRadii[inds[0]];
         double bestRadiiNP = initRadii[inds[1]];

        lastRadius = bestRadiiNP;
        lastRadiusWP = bestRadii;

        if(!quiet)
            System.out.println("Done");

        if(verbose)
        {
            System.out.println("Best Radii WP/NP: " + bestRadii + "," + bestRadiiNP);
        }

        //Run Pref-Div with optimal parameters

        return selectGenes(bestRadiiNP,bestRadii,dFile);

    }


    /***
     *
     * @return A boolean [] specifying which gene-gene relationships had prior information, excluding the target varialbe
     * This is useful for the final run of Pref-Div after selecting parameters
     */
    private boolean[] getPriorNoTarget()
    {
        boolean [] tempPrior = new boolean[numGenes* (numGenes-1)/2];
        int allCount = 0; //Counter for the prior info including the target
        int newCount = 0; //Where we are in the newly created array
        //Loop over all genes and the target variable
        for(int i = 0; i < numGenes + 1;i++)
        {
            for(int j = i+1; j < numGenes + 1;j++)
            {
                if(i!=targetIndex && j!=targetIndex) {
                    tempPrior[newCount] = dPriors[allCount];
                    newCount++;
                }
                allCount++;
            }
        }
        return tempPrior;

    }
    /***
     *
     * @param dFile Array of filenames that contain prior information files
     * @return A boolean [] specifying whether or not each gene has prior information about its relationship to the target variable
     */
    private boolean [] loadIPrior(String [] dFile)
    {
        boolean [] temp = new boolean[numGenes];
        try{
            for(int i = 0; i < dFile.length;i++)
            {
                BufferedReader b = new BufferedReader(new FileReader(dFile[i]));

                //Read all lines up to the one about the target variable
                for(int j = 0; j < targetIndex;j++)
                    b.readLine();


                String [] line = b.readLine().split("\t");

                int count = 0;
                for(int j = 0; j < line.length;j++)
                {
                    if(j==targetIndex)
                        continue;

                    //If this gene has information then mark it
                    if(Double.parseDouble(line[j])> 0)
                        temp[count] = true;
                    count++;
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
        for(int i = avg.length/2-1; i >= 0;i--)
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

    /**
     *
     * @param sum Total number of times this feature appeared
     * @param n Number of times the feature appeared according to this parameter value
     * @param qj num subsamples * number of parameters tested
     * @return The concordance between this parameter's prediction and all parameter's predictions
     */
    private double getTheta(int sum, int n,double qj)
    {
        double choose = CombinatoricsUtils.binomialCoefficient(subsamples.length,n);
        return choose*Math.pow(sum/qj,n)*Math.pow(1-(sum/qj),subsamples.length-n);
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
    private double[] getScoresSeparatePar(final int [] sums, final float [] uPost, final float [] varPost)
    {


        final double [] score = new double[numRadii*2];

        long readingTime = 0;
        final double qj = numRadii*subsamples.length;


        long time = System.nanoTime();

        for (int i = 0; i < numRadii; i++) {

            int [] curr;
            if(saveMemory) {

                curr = getCounts(i, true);

            }
            else {
                curr = getCountsFromFile(i,true);


            }



            final int[] curr2 = curr;

            final ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

            final int radIndex = i;

            class StabilityAction extends RecursiveAction {
                private int chunk;
                private int from;
                private int to;

                public StabilityAction(int chunk, int from, int to){
                    this.chunk = chunk;
                    this.from = from;
                    this.to = to;
                }



                @Override
                protected void compute(){
                    if (to - from <= chunk) {
                        for (int s = from; s < to; s++) {
                            double G = getG(curr2[s]); //Clusts[i][j][g] is number of times gene was selected in the subsamples for this parameter setting
                            if (uPost[s] < 0)//No Prior information, contributes only stability to the score
                            {


                                double theta = getTheta(sums[s], curr2[s], qj); //sums is the total number of times this gene was selected across all parameter settings
                                int index = radIndex+numRadii;
                                addToScore(score,index,theta*(1-G));

                            } else {

                                if (varPost[s] < 0)
                                    System.out.println(s);
                                double theta = getTheta(uPost[s], varPost[s], curr2[s]);

                                addToScore(score,radIndex,theta*(1-G));

                            }
                        }
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

            final int chunk = uPost.length/Runtime.getRuntime().availableProcessors();
            StabilityAction sa = new StabilityAction(chunk,0, uPost.length);
            pool.invoke(sa);


        }
        normalizeScore(score);

        if(profiling)
        {
            readingTime += System.nanoTime()-time;
        }
        return score;
    }

    private synchronized void addToScore(double [] score,int index, double value)
    {
        score[index]+=value;
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

        /**Temporary for debugging**/
        double [] stabs = new double[numRadii];
        double [] match = new double[numRadii];

        long readingTime = 0;
        double qj = numRadii*subsamples.length;


        long time = System.nanoTime();

        for (int i = 0; i < numRadii; i++) {

                int [] curr = null;
                if(saveMemory) {

                    curr = getCounts(i, true);

                }
                else {
                    curr = getCountsFromFile(i,true);


                }



            for(int g = 0; g < uPost.length;g++) //Loop through each gene
                {

                    double G = getG(curr[g]); //Clusts[i][j][g] is number of times gene was selected in the subsamples for this parameter setting
                    if (uPost[g] < 0)//No Prior information, contributes only stability to the score
                    {


                        double theta = getTheta(sums[g],curr[g],qj); //sums is the total number of times this gene was selected across all parameter settings


                        //b.probability(curr[g]);


                        score[i + numRadii] +=  (theta * (1-G));


                    }
                    else
                    {

                        if(varPost[g]<0)
                            System.out.println(g);
                        double theta = getTheta(uPost[g],varPost[g],curr[g]);
                        score[i]+=(theta *(1-G));
                        stabs[i]+=(1-G);
                        match[i]+=theta;

                    }
                }


        }


        normalizeScore(score);

        if(profiling)
        {
            readingTime += System.nanoTime()-time;
        }
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


        /***Shrink anyone that didn't meet radius threshold**/
        for (int x = 0; x < corrs.length; x++) {
            if (1 - corrs[x] < initRadii[i])
                curr[x]++;
        }

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
                /***TODO is this if statement right?***/
                if(prior[i]>0 && mean[i] > 0) {
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
        pvals = new double[dFile.length];
        normTao = new double[dFile.length];

        long pval = 0;
        long load = 0;
        long phiTao = 0;
        for(int i = 0; i < dFile.length;i++)
        {
            float [][] temp = new float[1][(numGenes*(numGenes+1)/2)];
            long time = System.nanoTime();
            temp[0] = Functions.loadTheoryMatrix(dFile[i],false,numGenes+1);
            whichPriors(temp,dPriors);
            load+= (System.nanoTime()-time);


            time = System.nanoTime();
            getPhi(temp[0],subsamples.length);
            tao[i] = getTao(temp,sums)[0];
            phiTao += (System.nanoTime()-time);


            time = System.nanoTime();

            System.out.println(dFile[i]);
            if(computeP) {
                pvals[i] = getPValue(tao[i], temp, sums);
                normTao[i] = lastNormTao;
            }
            pval  += (System.nanoTime()-time);



        }

        if(profiling)
        {
            System.out.println(Arrays.toString(pvals));

            System.out.println("loading: " + load/Math.pow(10,9) + ", phi/tao:" + phiTao/Math.pow(10,9) + ", P-Val:" + pval/Math.pow(10,9));
        }
        simTao = tao;
        double [] alpha = mgmPriors.getAlpha(tao);
        if(computeP)
            adjustP = mgmPriors.adjustPValues(pvals);
        return mgmPriors.getWeights(alpha);

    }

    /***
     *
     * @param tao The current tao value to query the distribution against
     * @param temp The current prior information source's info
     * @param sums The data-driven frequency estimates
     * @return The p-value of this tao against a set of randomly shuffled temp's
     */
    private double getPValue(double tao, float [][] temp, int [] sums)
    {
        float [][] t2 = new float[1][temp[0].length];
        for(int i = 0; i < temp[0].length;i++)
            t2[0][i] = temp[0][i];
        double [] hist = new double[numSamples];

        long time = 0;
        long shuffle = 0;
        long hisTime = 0;
        Random rgen = new Random(42);

        for(int i = 0; i < numSamples;i++)
        {
            /***Compute tao***/

            if(i==0)
                System.out.println(getTao(t2,sums)[0]);
            time = System.nanoTime();
            hist[i] = getTao(t2,sums)[0];
            hisTime += (System.nanoTime()-time);

            time = System.nanoTime();
            /***Shuffle values in the array***/
           t2[0] = RandomizeArray(t2[0],rgen);

           shuffle+= (System.nanoTime()-time);



        }
        time = System.nanoTime();
        Arrays.sort(hist);
        System.out.println(tao);
        System.out.println(Arrays.toString(hist));

        lastNormTao = tao/StatUtils.mean(hist);

        if (profiling)
        {
            System.out.println("Sorting: " + (System.nanoTime()-time)/Math.pow(10,9) + ", Shuffling: " + shuffle/Math.pow(10,9) + ", Hist Tao: " + hisTime/Math.pow(10,9));
        }

        int index =0;
        while(index < hist.length && tao > hist[index])
        {
            index++;
        }
        return index/(double)numSamples;
    }



    public static float[] RandomizeArray(float[] array, Random rgen){

        for (int i = 0; i < array.length; i++)
        {
            int index = rgen.nextInt(i + 1);
            // Simple swap
            float a = array[index];
            array[index] = array[i];
            array[i] = a;
        }

        return array;
    }

    private static ArrayList<Integer> getIndices(int length)
    {
        ArrayList<Integer> list = new ArrayList<Integer>();
        for(int i = 0; i < length;i++)
            list.add(i);
        return list;
    }


    private static HashMap<Float,Integer> getArrayCounts(float [] array)
    {
        HashMap<Float,Integer> counts = new HashMap<Float,Integer>();
        for(int i = 0; i < array.length;i++)
        {
            if(array[i]<=0)
                continue;
            if(counts.get(array[i])==null)
                counts.put(array[i],1);
            else
                counts.put(array[i],counts.get(array[i])+1);

        }
        return counts;
    }


    private static float[] RandomizeArray2(HashMap<Float,Integer>map,ArrayList<Integer> indices, Random rgen){
        float [] result = new float[indices.size()];
        Collections.shuffle(indices,rgen);
        int count = 0;
        for(Float f:map.keySet())
        {
            int times = map.get(f);
            for(int i = count; i < count+ times;i++)
            {
                result[indices.get(i)] = f;
            }
            count+=times;
        }
        return result;

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
                if(sources[i][j] > 0)
                    prior[j] = true;
            }
        }
    }


    private double [] getTao(float[][]phi, int[]counts)
    {
        double [] tao = new double[phi.length];
        //For each gene, compute deviation and sum these, then divide by number of things you summed
        for(int i = 0; i < phi.length;i++) //Loop over sources
        {
            int numPriors = 0;
            for(int j = 0; j<phi[i].length;j++)
            {
                if(phi[i][j]> 0) {
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
            if(in[i]> 0.0)
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
        long corrTime = 0;
        long readTime = 0;
        for(int k = 0; k < subsamples.length;k++)
        {

            ArrayList<Gene> temp = createGenes(data,target,true);
            DataSet currData = data.subsetRows(subsamples[k]);

            //Threshold intensities and correlations after the fact
            float [] corrs;
            long time = System.nanoTime();

            if(parallel)
            {
                corrs = Functions.computeAllCorrelationsPar(temp,currData);
            }
            else
            {
                corrs = Functions.computeAllCorrelations(temp,currData,partialCorrs,false,false,1);
            }
            corrTime += System.nanoTime()-time;

            time = System.nanoTime();

            for(int i = 0; i < initRadii.length;i++)
            {
                    int [] curr = getCounts(i,corrs,read);

                    for(int x = 0; x < curr.length;x++)
                    {
                        sums[x]+=curr[x];
                    }
            }
            readTime += System.nanoTime()-time;
            read = true;
        }

        if(verbose)
        {
            System.out.println("Computed counts for intensities and correlations: ");
        }

        if(profiling)
        {
            System.out.println("Correlations across subsamples: " + corrTime/Math.pow(10,9));
            System.out.println("Reading correlations across subsamples " + readTime/Math.pow(10,9));
        }



        long time = System.nanoTime();
        float [] uPost = new float[(numGenes*(numGenes+1))/2];
        for(int i = 0; i < uPost.length;i++)
            uPost[i] = -1;

        double [] score;
        if(parallel)
        {
            score = getScoresSeparatePar(sums,uPost,null);
            constrictCorrs = score;
        }
        else
        {
            score = getScoresSeparate(sums,uPost,null);
            constrictCorrs = score;

        }


        double [] scoreResult = new double[numRadii];

        for(int i = 0; i < numRadii;i++)
        {
                scoreResult[i] = score[i+numRadii];
        }

        if(profiling)
        {
            time = System.nanoTime()-time;
            System.out.println("Posterior Distribution: " + time/Math.pow(10,9));
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


        long time = System.nanoTime();
        for(int k = 0; k < subsamples.length;k++)
        {
            ArrayList<Gene> temp = createGenes(data,target,true);
            DataSet currData = data.subsetRows(subsamples[k]);

            //Threshold intensities and correlations after the fact
            float [] corrs;
            if(parallel)
            {
                corrs = Functions.computeAllCorrelationsPar(temp,currData);

            }else
            {
                corrs = Functions.computeAllCorrelations(temp,currData,partialCorrs,false,false,1);

            }

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
            System.out.println("Computed counts for intensities and correlations: ");
        }

        if(profiling)
        {
            time = System.nanoTime()-time;
            System.out.println("Correlation computation with counts: " + time/Math.pow(10,9));
        }

        if(!quiet)
            System.out.print("Computing similarity weights...");
        //Compute similarity weights for each prior knowledge source
        time = System.nanoTime();
        double [] simWeights = getSimilarityWeights(sums,dFile);

        if(profiling)
        {
            time = System.nanoTime()-time;
            System.out.println("Similarity Weights Computation " + time/Math.pow(10,9));
        }

        lastSimilarityWeights = simWeights;
        if(!quiet)
            System.out.println("Done");
        if(verbose)
        {
            System.out.println("Similarity weights: " + Arrays.toString(simWeights));
        }

        if(!quiet)
            System.out.print("Computing scores for similarities...");


        time = System.nanoTime();
        //Load each prior knowledge source, weighted by reliability
        float [] meanDis = loadTheoryFiles(dFile,simWeights,numGenes);


        //Get variance of mixture distribution for each edge
        float [] varDis = getVarMixture(simWeights,simTao,meanDis,dFile);


        //Load Intensity prior using non -1 entries from meandis
        iPriors = loadIPrior(meanDis);

        //Load dissimilarity prior using non -1 entries from meandis
        dPriorsWT = loadPriorsNoTarget(meanDis);

        //Get posterior distributions of predicted dissimilarity and intensity from prior sources
        float[] uPost = getMeanPosterior(sums,meanDis, varDis);
        float [] varPost = getVarPosterior(sums, varDis);

        if(profiling)
        {
            time = System.nanoTime()-time;
            System.out.println("Load theory sources and get posterior distribution; " + time/Math.pow(10,9));
        }

        //Compute posterior probabilities and stability of each gene, synthesize these results into a score for each parameter setting
        double[] scoresSim;


        time = System.nanoTime();
        if(parallel) {
            scoresSim = getScoresSeparatePar(sums, uPost, varPost);
            constrictCorrs = scoresSim;
        }
        else {
            scoresSim = getScoresSeparate(sums, uPost, varPost);
            constrictCorrs = scoresSim;
        }
        if(!quiet)
            System.out.println("Done");

        if(profiling)
        {
            time = System.nanoTime()-time;
            System.out.println("Compute Scores (Normal Distribution) " + time/Math.pow(10,9));
        }

        return scoresSim;
    }


    /***
     *
     * @param means Mean similarity scores from the prior knowledge sources
     * @return boolean [] where true means that some source gave information about gene-gene relationship i
     */
    private boolean [] loadPriorsNoTarget(float [] means)
    {
        boolean [] prior = new boolean[numGenes*(numGenes-1)/2];
        for(int i = numGenes;i<means.length;i++)
        {
            if(means[i]>0)
                prior[i-numGenes] = true;
        }
        return prior;
    }

    /***
     *
     * @param meanDis Mean similarity scores from the prior knowledge sources
     * @return boolean [] where true means that some source gave information
     */
    private boolean[] loadIPrior(float [] meanDis)
    {
        boolean [] prior = new boolean[numGenes];
        for(int i = 0; i < numGenes;i++)
        {
            if(meanDis[i] >0)
                prior[i] = true;
        }
        return prior;
    }


    private synchronized void addGeneConnections(float [] corrs, int [] curr, double radii)
    {
        for(int i = 0; i < corrs.length;i++)
        {
            if(1-corrs[i]<radii)
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
    private float[] loadTheoryFiles(String [] dFile, double [] weights, int numGenes)
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

            /**TODO Make sure this is right, no prior information if it says the edge probability doesn't exist***/
            if(result[i] < wpCutoff)
                result[i] = -1;
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
