package edu.pitt.csb.Pref_Div;

import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.Priors.mgmPriors;
import edu.pitt.csb.stability.StabilityUtils;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by vinee_000 on 10/19/2018.
 */
public class PiPrefDiv {

    private DataSet data; //The current expression dataset to be analyzed, assumes that all variables are fair game for Pref-Div except the target
    private String target; //The target variable of interest
    private double [] initRadii; //Initial radius values to test
    private int [] initTopK; //Initial Top K values to test
    private double lowRadii = 0.01; //Low range of radius values to test
    private double highRadii = 0.99; //High range of radius values to test
    private int numRadii = 40; //Number of radius values to test
    private int numK = 40; //Number of top-K values to test
    private int lowK = 100; //Low range of top-K values to test
    private int highK = 2000; //High range of top-K values to test
    private int [][] subsamples; //subsamples for repeated Pref-Div
    private double normalEpsilon = 0.5;//Range around which to get probability for theta with prior information
    private int numGenes;//number of total variables in the data
    private double [] intTao;//Tao values for each intensity source
    private double [] simTao;//Tao values for each similarity source
    private boolean LOOCV = false; //Should we do leave-one-out cross validation for the data-theory tradeoff selection
    private Map<Gene,List<Gene>> lastCluster;//The last set of clusters from the run of selectGenes()
    private boolean verbose = false; //Should we print full debugging information?

    //Constructor that uses default parameters for both initial parameter ranges
    public PiPrefDiv(DataSet data, String target)
    {
        this.data = data;
        this.numGenes = data.getNumColumns()-1;
        this.target = target;
        this.initRadii = initRadiiDefault();
        this.initTopK = initTopKDefault();
    }


    //Constructor that specifies the same number of parameters for radius and k
    public PiPrefDiv(DataSet data, String target, int numParams)
    {
        this.numGenes = data.getNumColumns()-1;
        this.data = data;
        this.target = target;
        this.numRadii = numParams;
        this.numK = numParams;
        this.initRadii = initRadiiDefault();
        this.initTopK = initTopKDefault();
    }


    //Construct with number of radii and k parameters to test
    public PiPrefDiv(DataSet data, String target, int numRadii, int numK)
    {
        this.data = data;
        this.numGenes = data.getNumColumns()-1;
        this.target = target;
        this.numRadii = numRadii;
        this.numK = numK;
        this.initRadii = initRadiiDefault();
        this.initTopK = initTopKDefault();
    }


    //Constructor with initial ranges both specified
    public PiPrefDiv(DataSet data, String target, double [] radii, int [] k)
    {
        this.data = data;
        this.numGenes = data.getNumColumns()-1;
        this.target = target;
        this.initRadii = radii;
        this.numRadii = radii.length;
        this.initTopK = k;
        this.numK = k.length;
    }

    public Map<Gene,List<Gene>> getLastCluster()
    {
        return lastCluster;
    }

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


    //Create default range of number of variables to select
    private int [] initTopKDefault()
    {
        int [] topK = new int[numK];

        for(int i = 0; i < numK;i++)
        {
            topK[i] = lowK + (int)((highK-lowK)*i/(double)numK);
        }
        return topK;
    }

    //Full procedure to select genes based on prior information Pref-Div
    //Input: boot (should we do bootstrapping (true) or subsampling (false))
    //Input: numSamples (number of samples to subsample or bootstrap)
    //Input: path to theory Intensity File -> iFile (should have a header)
    //Input: path to theory Dissimilarity File -> dFile (no header or rownames for these files)
    public ArrayList<Gene> selectGenes(boolean boot, int numSamples, String iFile, String [] dFile,boolean useCausalGraph)
    {

        //Generate subsamples of the data
        System.out.print("Generating subsamples...");
        if(subsamples==null)
        {
            subsamples = genSubsamples(boot,numSamples);
        }
        System.out.println("Done");

        //Compute stability of each gene selection and each gene-gene relationship
        System.out.print("Computing stability across radii and k values...");
        int[][][] clusts = computeStabs(); //i -> Radii, j -> TopK, k -> geneID
        System.out.println("Done");

        int [] sums = getFullCounts(clusts);
        System.out.println(Arrays.toString(sums));

        System.out.print("Computing intensity weights...");
        //Compute intensity weights
        double [] intWeights = getIntensityWeights(sums,iFile);
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
        float[] uPostInt = getMeanPosterior(clusts, means, varIntense);
        float [] varPostInt = getVarPosterior(clusts, varIntense);

        System.out.println(Arrays.toString(uPostInt));
        System.out.println(Arrays.toString(varPostInt));
        //Compute posterior probabilities and stability of each gene, synthesize these results into a score for each parameter setting
        double[][] scoresInt = getScores(clusts,uPostInt,varPostInt); //i -> radii, j -> topK
        System.out.println("Done");

        if(verbose)
        {
            System.out.println("Printing all scores...");
            for(int i = 0; i < scoresInt.length;i++)
            {
                System.out.print(initRadii[i] + "\t");
                for(int j = 0; j < scoresInt[i].length;j++)
                {
                    System.out.print(initTopK[j] + ": " + scoresInt[i][j] + "\t");
                }
                System.out.println();
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
        System.out.println("Done");
        if(verbose)
        {
            System.out.println("Similarity weights: " + Arrays.toString(simWeights));
        }


        System.out.print("Computing scores for similarities...");
        //Load each prior knowledge source, weighted by reliability
        float [] meanDis = loadTheoryFiles(dFile,simWeights);


        //Get variance of mixture distribution for each edge
        float [] varDis = getVarMixture(simWeights,simTao,meanDis,dFile);



        //Get posterior distributions of predicted dissimilarity and intensity from prior sources
        float[] uPost = getMeanPosterior(clusts,meanDis, varDis);
        float [] varPost = getVarPosterior(clusts, varDis);

        //Compute posterior probabilities and stability of each gene, synthesize these results into a score for each parameter setting
        double[][] scoresSim = getScores(clusts,uPost,varPost); //i -> radii, j -> topK
        System.out.println("Done");


        if(verbose)
        {
            System.out.println("Printing all similarity scores...");
            for(int i = 0; i < scoresInt.length;i++)
            {
                System.out.print(initRadii[i] + "\t");
                for(int j = 0; j < scoresInt[i].length;j++)
                {
                    System.out.print(initTopK[j] + ": " + scoresInt[i][j] + "\t");
                }
                System.out.println();
            }
        }

        //TODO How to combine scoresInt with scoresSim to get top parameters??


        System.out.print("Merging scores to get top parameters...");
        //Find best radii and topK values and then Run Pref-Div with cross-validation to get final Output!!
        double [][] avgScore = getAvgScore(scoresSim,scoresInt);

        int[]inds = getMaxScore(avgScore);


        double bestRadii = initRadii[inds[0]];
        int bestK = initTopK[inds[1]];
        System.out.println("Done");

        if(verbose)
        {
            System.out.println("Best Radii: " + bestRadii + "\nBest K: " + bestK);
        }



        //Set not needed variables to null to save memory
        simWeights = null;
        varDis = null;
        uPost = null;
        varPost = null;
        scoresSim = null;
        scoresInt = null;
        avgScore = null;


        System.out.print("Running Pref-Div with optimal parameters...");
        //Run Pref-Div with optimal parameters

        RunPrefDiv rpd = new RunPrefDiv(means,meanGenes,data,target,LOOCV);
        rpd.setRadius(bestRadii);
        rpd.setTopK(bestK);
        rpd.setAccuracy(0);
        rpd.setNumAlphas(numK);
        rpd.setNS(subsamples.length);
        rpd.setCausalGraph(useCausalGraph);
        rpd.useStabilitySelection();

        ArrayList<Gene> top = rpd.runPD();
        HashMap<Gene,List<Gene>> map = rpd.getClusters();
        lastCluster = map;
        System.out.println("Done");

        if(verbose)
        {
            System.out.println("Selected top K features\n" + top);
            System.out.println("All clusters\n" + map);
        }
        return top;
    }

    private int[] getMaxScore(double [][] avg)
    {
        int x = 0;
        int y = 0;
        double maxScore = -1;
        for(int i = 0; i < avg.length;i++)
        {
            for(int j = 0; j < avg[i].length;j++)
            {
                if(avg[i][j]>maxScore)
                {
                    x = i;
                    y = j;
                    maxScore = avg[i][j];
                }
            }
        }
        return new int[]{x,y};
    }

    private double[][] getAvgScore(double [][] s1, double [][] s2)
    {
        double [][] result = new double[s1.length][s1[0].length];
        for(int i = 0; i < s1.length;i++)
        {
            for(int j = 0; j < s1[i].length;j++)
            {
                result[i][j] = (s1[i][j] + s2[i][j])/2;
            }
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
        double P = sum/(double)(numRadii*numK*subsamples.length);
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

    private double[][] getScores(int [][][] clusts, float [] uPost, float [] varPost)
    {
        int [] sums = getFullCounts(clusts);
        double [][] score = new double[numRadii][numK];
        int offset = 0;
        if(uPost.length>numGenes)
            offset = numGenes;
        for (int i = 0; i < numRadii; i++) {
            for (int j = 0; j < numK; j++) {
                for(int g = 0; g < uPost.length;g++) //Loop through each gene
                {
                    double G = getG(clusts[i][j][g+offset]); //Clusts[i][j][g] is number of times gene was selected in the subsamples for this parameter setting
                    if (uPost[g] < 0)//No Prior information, contributes only stability to the score
                    {
                        double theta = getTheta(sums[g+offset],clusts[i][j][g+offset]); //sums is the total number of times this gene was selected across all parameter settings
                        score[i][j] += theta*(1-G);
                    }
                    else
                    {
                        double theta = getTheta(uPost[g],varPost[g],clusts[i][j][g+offset]);
                        score[i][j]+=theta*(1-G);
                    }
                }
            }
        }
        return score;
    }


    //Get posterior variance based on computed counts and variance, modified to work for both intensity and similarity
    private float[] getVarPosterior(int [][][] clusts, float [] var)
    {
        float [] varPost = new float[var.length];
        int [] sums = getFullCounts(clusts);
        int offset = 0;
        if(var.length>numGenes)
            offset = numGenes;
        for(int i = 0; i < var.length;i++)
        {
            if(var[i]<0)
                varPost[i] = -1;
            else
            {
                double p = sums[i + offset] / (double) (numRadii * numK * subsamples.length);
                double currVar = p * (1 - p) * subsamples.length;
                varPost[i] =(float)( currVar * var[i]/(currVar+var[i]));
            }
        }

        return varPost;
    }

    //Get Posterior mean for intensity values from theory
    private float[] getMeanPosterior(int[][][]clusts,float[]genes,float [] vars )
    {
        int [] sums = getFullCounts(clusts);
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
                double p = sums[i+offset] / (double) (numRadii * numK * subsamples.length);
                double currVar = p * (1 - p) * subsamples.length;
                uPost[i] = (float) (((sums[i+offset] / (double) (numRadii * numK)) * vars[i] + currVar * genes[i]) / (currVar + vars[i]));

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
            float [] prior = Functions.loadTheoryMatrix(dFile[j],false);
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
                    sum+=weights[j]*(  tao[j]*tao[j]+Math.pow((genes[i]-priors[j][i])*subsamples.length,2) );
                    weight+=weights[j];
                }
                vars[i] = (float)(sum/weight);

            }
            else
                vars[i] = -1;
        }
        return vars;
    }
    private int[] getFullCounts(int [][][]clusts)
    {
        int [] res = new int[clusts[0][0].length];
        for(int i = 0; i < res.length;i++)
        {
            for(int j = 0; j < clusts.length;j++)
            {
                for(int k = 0; k < clusts[j].length;k++)
                {
                    res[i]+=clusts[j][k][i];
                }
            }
        }
        return res;
    }


    //This is written as a separate function due to the offset in computing the sums and due to memory requirements preventing loading all sources at once
    private double [] getSimilarityWeights(int[] sums, String [] dFile)
    {
        //Same procedure as the intensity weights except don't load all sources at once for memory concerns
        double [] tao = new double[dFile.length];
        for(int i = 0; i < dFile.length;i++)
        {
            float [][] temp = new float[1][(numGenes*(numGenes-1)/2)];
            temp[0] = Functions.loadTheoryMatrix(dFile[i],false);

            getPhi(temp[0],subsamples.length);
            tao[i] = getTao(temp,sums)[0];
        }
        double [] alpha = mgmPriors.getAlpha(tao);
        return mgmPriors.getWeights(alpha);

    }

    private double [] getIntensityWeights(int[]sums, String iFile)
    {
        //Load each intensity source from the iFile
        float [][] sources = Functions.loadIntensityData(iFile,false);


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
                    tao[i] += Math.abs(phi[i][j] - counts[j] / (double) (numRadii * numK));
                }
            }
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
    private ArrayList<Gene> loadGenes(String iFile,double[]iWeights)
    {
        return Functions.loadGeneData(iFile,false,iWeights);
    }


    //Compute the set of pref-div genes and their clusters for each setting of radii, k, and subsample
    //Output: Array of size radii tested x k values tested x (number of genes + (numGenes*(numGenes-1))/2)
    //For each 1-D array (the last dimension), the first numGenes boxes are how many times gene i showed up in these subsamples
    //the rest of the boxes denotes how many times genes i and j were in the same cluster
    //TODO Need to parallelize this stability search
    private int[][][] computeStabs()
    {
        ArrayList<Gene> temp = createGenes();
        int[][][] result = new int[initRadii.length][initTopK.length][numGenes + (numGenes)*(numGenes-1)/2];
        for(int i = 0; i < initRadii.length;i++)
        {
            for(int j = 0; j < initTopK.length;j++)
            {
                int [] curr = new int[numGenes+ (numGenes*(numGenes-1)/2)];
                for(int k = 0; k < subsamples.length;k++)
                {
                    DataSet currData = data.subsetRows(subsamples[k]);
                    float [] corrs = Functions.computeAllCorrelations(temp,currData,false,true);
                    temp = Functions.computeAllIntensities(temp,1,currData,target,false,false,true);
                    Collections.sort(temp,Gene.IntensityComparator);
                    PrefDiv pd = new PrefDiv(temp,initTopK[j],0,initRadii[i], corrs);
                    pd.setCluster(true);
                    ArrayList<Gene> topGenes = pd.diverset();
                    HashMap<Gene,List<Gene>> clusts = pd.clusters;
                    addFoundGenes(topGenes,curr);
                    addClusterGenes(clusts,curr);
                }
                result[i][j] = curr;
            }
        }
        return result;
    }

    private void addFoundGenes(ArrayList<Gene> top, int [] curr) {
        for (int i = 0; i < top.size(); i++)
        {
            curr[top.get(i).ID]++;
        }

    }

    private void addClusterGenes(HashMap<Gene,List<Gene>> map, int [] curr)
    {
        for(Gene g:map.keySet())
        {
            List<Gene> list = map.get(g);
            int x = g.ID;
            if(list!=null) {
                for (int i = 0; i < list.size(); i++) {
                    int y = list.get(i).ID;
                    curr[Functions.getIndex(x, y, numGenes) + numGenes]++;
                }
            }
        }
    }


    private ArrayList<Gene> createGenes()
    {
        ArrayList<Gene> gList = new ArrayList<Gene>();
        for(int i = 0; i < numGenes;i++)
        {
            if(data.getVariable(i).getName().equals(target))
                continue;
            Gene g = new Gene(i);
            g.symbol = data.getVariable(i).getName();
            gList.add(g);
        }
        return gList;
    }

    //Give an array of filenames for all theory files, load each one into a row of the matrix to get all theory information in a single matrix
    private float[] loadTheoryFiles(String [] dFile, double [] weights)
    {
        float [] result = new float[(data.getNumColumns()*data.getNumColumns()-1)/2];
        float [] totalWeight = new float[result.length];
        for(int i = 0; i < dFile.length;i++)
        {
            float [] temp = Functions.loadTheoryMatrix(dFile[i],false);
            for(int j = 0; j < temp.length;j++)
            {
                if(temp[j]>0)
                {
                    result[j] += temp[j]*weights[i];
                    totalWeight[j] += weights[i];
                }
            }

        }
        for(int i = 0; i < result.length;i++)
            result[i] /= totalWeight[i];
        return result;

    }

    //Generate subsampled indices for the given dataset
    //Input: bootstrap ( should we bootstrap to generate samples or should we use subsampling)
    //Input: numSamples (number of samples to test)
    //Output: an array of integers specifying the samples that will be used in the repeated Pref-Div runs
    private int [][] genSubsamples(boolean bootstrap, int numSamples)
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
