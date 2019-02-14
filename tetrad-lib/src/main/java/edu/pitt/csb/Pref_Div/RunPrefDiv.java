package edu.pitt.csb.Pref_Div;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.*;
import edu.cmu.tetrad.algcomparison.simulation.MixedLeeHastieSimulation;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.*;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.*;
import edu.pitt.csb.stability.CrossValidationSets;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.StabilityUtils;
import org.apache.commons.math3.analysis.function.Logistic;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;

import java.io.*;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

import static org.apache.commons.collections4.CollectionUtils.intersection;
import static org.apache.commons.collections4.CollectionUtils.union;

/**
  Runs genomic Pref-Div algorithm with subsampling and clustering
    TODO Handle Discrete features in theory relevance and dissimilarity
 TODO Switch to median accuracy instead of mean for cross validation?
    TODO Clean up this entire file, reduce some of the unnecessary bloat, and test the reading from file for accuracy
 */
public class RunPrefDiv {

    public enum ClusterType{
        PCA,MEAN,MEDIAN,NONE
    }

    private DataSet data;
    private Random rand;
    private int B = 50; //Number of Subsamples for CPSS
    private ArrayList<Gene> genes;
    private float[] dissimilarity;
    private int numSubs = 20; //Number of subsamples
    private int numFolds = 10; //Number of folds for internal cross validation
    private int topK = 50;
    private double accuracy = 0;
    private double radius = 0.5;
    private int numAlphas = 20; //How many alpha values should we test? (Data-Theory tradeoff alpha)
    private double alphaLimit = 1; //The largest value of alpha we test
    private double lambdaLow = 0.05; //Lambda low limit
    private double lambdaHigh = 0.95; //Lambda high limit for MGM StEPS
    private double g = 0.05; //Stability threshold
    private double pThresh = 0.05; //p-value threshold
    private double pThreshWP = 0.05; //p-value threshold for those which we have prior information (only if separate params is true)
    private double radiusWP = 0.5; //Radius for those which we have prior information (only if separate params is true)

    private boolean [] withPrior; //Says which edges we have prior information for

    private boolean separateParams = false; //Do we treat edges where we have prior information differently?
    private ArrayList<Gene> lastGeneSet;
    private HashMap<Gene,List<Gene>> lastClusters;
    private String target; //Target variable
    private boolean useClusterStability = false; //cluster stability, or gene wise stability
    private boolean partialCorr = false; //Use partial correlations instead of correlation for continuous variables
    private double [][] lastStepsStabilities; //Last set of gene stabilities given by StEPS

    private boolean loocv = false;
    private int [][] subs;
    private Graph trueGraph;
    private DataSet summarizedData; //The dataset you get using the prefidv method on the current dataset with the current clustering procedure


    private boolean useCrossValidation = true; //Should we use cross-validation to determine alpha?
    private boolean useStabilitySelection = false; //Should we use stability selection when running Pref-Div?
    private boolean useCausalGraph = false; //Should we use a causal graph when determining predictive neighbors for cross-validation
    private boolean usePThreshold = false;//Should we use a p-value threshold to shrink insignificant p-value correlations to 0
    private boolean clusterByCorrs = false; //Should we cluster purely by significant correlations instead of using Pref-Div
    private ClusterType clusterType; //Which clustering method should we use (Options are PCA, Median, Mean, None)




    private double newTime = 0;
    private double oldTime = 0;
    private int diffs = 0;


    public RunPrefDiv(float [] dissimilarity, ArrayList<Gene> genes, DataSet data,String target,boolean leaveOneOut)
    {
        this.rand = new Random();
        this.data = data;
        this.genes = genes;
        this.dissimilarity = dissimilarity;
        this.target = target;
        this.loocv = leaveOneOut;
        this.clusterType = ClusterType.NONE;
    }


    public void setWithPrior(boolean [] wp)
    {
        withPrior = wp;
        separateParams = true;
        useCrossValidation = false;
    }
    public void setAllParams(double radius, double radiusWP, double pThresh, double pThreshWP)
    {
        this.radius = radius;
        this.radiusWP = radiusWP;
        this.pThresh = pThresh;
        this.pThreshWP = pThreshWP;
        this.usePThreshold = true;
        useCrossValidation = false;
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
    public void setPThreshold(double p){pThresh =  p; usePThreshold=true;}
    public void setThreshold(double t)
    {
        g = t;
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
    public void setNumFolds(int nf){this.numFolds = nf;}
    public void setTrueGraph(Graph truth){this.trueGraph = truth;}
    public ArrayList<Gene> getLastGeneSet(){return lastGeneSet;}
    public void usePartialCorrelation(boolean pc){partialCorr = pc;}
    public void useCrossValidation(boolean b){useCrossValidation = b;}
    public void useStabilitySelection(boolean b){useStabilitySelection = b;}
    public void setCausalGraph(boolean cg){useCausalGraph = cg;}
    public void clusterByCorrs(){clusterByCorrs=true;}
    public void setB(int B){this.B = B;}
    public void setClusterType(ClusterType c){clusterType = c;}
    public DataSet getSummarizedData(){return summarizedData;}





    public ArrayList<Gene> runPD()
    {


        if(subs==null) {
            System.out.print("Generating Subsamples...");

            if (loocv)
                subs = StabilityUtils.generateSubsamples(data.getNumRows());
            else if(useCrossValidation)
                subs = StabilityUtils.generateSubsamples(numFolds, data.getNumRows());
            else
                subs = StabilityUtils.generateSubsamples(numSubs,data.getNumRows());
            System.out.println("Done");

        }
        double [] alphas = new double[numAlphas+1];
        for(int i = 0; i < alphas.length;i++)
        {
            alphas[i] = i*alphaLimit/(double)numAlphas;
        }


        /***run PD across B = 100 complementary pairs, and get those genes which appear(in the diverset) in those pairs most frequently
        //Then run StEPS with these selected genes and use the neighborhood of the target to get a predicted accuracy
        //Repeat this for all alphas and choose the one with the best predicted accuracy. Thus, stability comes from the selection of genes
        //Accuracy is used to select the right data-theory tradeoff parameter***/


        if(useCrossValidation)
        {
            System.out.println("Running cross-validation version of Pref-Div");

            List<List<Gene>> allGenes = new ArrayList<List<Gene>>();
            List<HashMap<Gene,List<Gene>>> allClusters = new ArrayList<>();
            double [][] accuracies = new double[alphas.length][numFolds];
            double [][] stabilities = new double[alphas.length][numFolds];
            int [][] trainSubs = new int[1][1];
            int maxAlpha = -1;
            double maxAcc = Double.MAX_VALUE;
            for(int i = 0; i < alphas.length;i++)
            {

                if(useStabilitySelection)
                System.out.println("Computing accuracy for " + alphas[i]);


                //Get the current training dataset and testing dataset
                A:for(int j = 0; j < subs.length;j++) {
                    Arrays.sort(subs[j]);
                    DataSet test = data.subsetRows(subs[j]);
                    DataSet train = getTrainingData(data,subs[j]);



                    /***Generate subsamples for the training dataset for internal CV by StEPS the first time through***/
                    if(j==0) {
                        int b = (int)(10*Math.sqrt(train.getNumRows()));
                        if(b>=train.getNumRows())
                            b = train.getNumRows()/2;
                        trainSubs = StabilityUtils.subSampleNoReplacement(train.getNumRows(),b,numSubs);
                    }

                    /*** Run Pref-Div with or without Stability Selection ***/
                    if(useStabilitySelection)
                        System.out.println("Running PD for alpha " + alphas[i] + ", and cv fold: " + j);


                    double stab;
                    if(useStabilitySelection)
                        stab = pdCPSS(alphas[i],train);
                    else
                        stab = runPD(alphas[i],train);
                    ArrayList<Gene> genes = lastGeneSet;
                    List<Node> dNeighbors = new ArrayList<Node>();


                    /***Temporarily add target variable to list of genes, so that they are included in the summarized dataset***/
                    ArrayList<Gene> sumTemp = (ArrayList<Gene>)genes.clone();
                    Gene xx = new Gene(-1);
                    xx.symbol = "Target";
                    sumTemp.add(xx);


                    /**Use summarization method specified along with selected genes to subset down to a dataset for causal analysis or straight prediciton**/
                    DataSet [] summarized = summarize(train,test,sumTemp,lastClusters,clusterType);


                    /***Construct Regression Dataset to test this value of alpha***/
                    train = summarized[0];
                    test = summarized[1];

                    /***Run Pref-Div and then use StEPS to select regression features***/
                    if(useCausalGraph)
                    {
                        //Create a dataset with only those variables selected by PD
                       DataSet temp = summarized[0];
                        if(!temp.isMixed())
                        {
                            temp.addVariable(new DiscreteVariable("Dummy"));
                            int col = temp.getColumn(temp.getVariable("Dummy"));
                            Random rand = new Random();
                            for(int x = 0; x < temp.getNumRows();x++)
                            {
                                temp.setInt(x,col,rand.nextInt(2));
                            }
                        }
                        DataSet[] subsampled = new DataSet[numSubs];
                        for (int k = 0; k < subsampled.length; k++) {
                            subsampled[k] = temp.subsetRows(trainSubs[k]);
                        }
                        double[] lambdas = new double[numAlphas];
                        for (int k = 0; k < lambdas.length; k++) {
                            lambdas[k] = lambdaLow + (lambdaHigh - lambdaLow) * k / numAlphas;
                        }
                        STEPS s = new STEPS(temp, lambdas, g, trainSubs);
                        System.out.print("Running StEPS...");
                        s.runStepsArrayPar();
                        System.out.println("Done");
                        List<Node> neighbors = s.lastGraph.getAdjacentNodes(s.lastGraph.getNode(target));
                        dNeighbors = new ArrayList<Node>();
                        for(Node n:neighbors)
                        {
                            if(n.getName().equals("Dummy"))
                                continue;
                            dNeighbors.add(temp.getVariable(n.getName()));
                        }
                    }
                    else
                    {
                        for(int k = 0; k< summarized[0].getNumColumns();k++)
                        {
                            if(summarized[0].getVariable(k).getName().equals("Dummy"))
                                continue;
                            if(summarized[0].getVariable(k).getName().equals(target))
                                continue;
                            dNeighbors.add(summarized[0].getVariable(k));
                        }
                    }



                    /***Randomly select a gene to be connected if none are connected in the causal graph?***/
                    if(dNeighbors.size()==0)
                    {
                        int r = rand.nextInt(train.getNumColumns());
                        while(train.getVariable(r).getName().equals("Dummy")||train.getVariable(r).getName().equals(target))
                            r = rand.nextInt(train.getNumColumns());
                        dNeighbors.add(train.getVariable(r));
                    }


                    /***Do a linear or logistic regression here depending upon the type of the target on the test set
                    //With the variables in the neighborhood of the target as the features***/

                    if(test.getVariable(target)instanceof ContinuousVariable)
                    {
                        if(useStabilitySelection)
                            System.out.print("Testing accuracy via regression...");
                        RegressionDataset rd = new RegressionDataset(train);
                        try {
                            RegressionResult res = rd.regress(train.getVariable(target), dNeighbors);
                            double RMSE = testRegression(res,test,dNeighbors);
                            accuracies[i][j] = RMSE;
                            if(useStabilitySelection)
                            System.out.println(RMSE + " Done");
                        }catch(Exception e)
                        {
                            continue A;
                        }


                    }
                    else
                    {
                       System.out.println("Currently discrete targets are not supported, please bug Vineet to fix this ASAP");
                       System.exit(-1);
                    }
                    stabilities[i][j]= stab;
                }


                /***Keep the alpha that gives the best accuracy***/
                if(useStabilitySelection)
                    System.out.println("Accuracy for alpha: " + alphas[i] + ", " + accuracies[i] + ", Stability: " + stabilities[i]);


                if(StatUtils.median(accuracies[i]) < maxAcc)
                {
                    maxAcc = StatUtils.median(accuracies[i]);
                    maxAlpha = i;
                }

                double stab;
                /***Use this alpha for the full dataset and see what genes you get***/
                if(useStabilitySelection)
                    stab = pdCPSS(alphas[i],data);
                else
                    stab = runPD(alphas[i],data);
                ArrayList<Gene> genes = lastGeneSet;
                if(useStabilitySelection)
                System.out.println("Gene Set for " + alphas[i] + ": " + genes);
                allGenes.add(genes);
                allClusters.add(lastClusters);

            }


            /***Construct summarized dataset based on selected genes and clusters***/
            lastGeneSet = (ArrayList<Gene>)allGenes.get(maxAlpha);
            lastClusters = allClusters.get(maxAlpha);


            ArrayList<Gene> sumTemp = (ArrayList<Gene>)lastGeneSet.clone();
            Gene xx = new Gene(-1);
            xx.symbol=target;
            sumTemp.add(xx);

            summarizedData = summarize(data,data,sumTemp,lastClusters,clusterType)[0];

        }

        /***If we aren't using cross validation then no alpha value is required, just run PrefDiv as is***/
        /***This will only be used if we are running PiPrefDiv***/
        else {
            PrefDiv p;
            if(withPrior==null)
            {
                p = new PrefDiv(genes,topK,accuracy,radius,dissimilarity);
            }
            else
            {
                 p = new PrefDiv(genes,topK,accuracy,radius,radiusWP,withPrior,dissimilarity);
            }
            p.setCluster(clusterByCorrs);
            lastGeneSet = p.diverset();
            lastClusters = p.clusters;

            summarizedData = summarize(data,data,(ArrayList<Gene>)lastGeneSet.clone(),lastClusters,clusterType)[0];

        }
        return lastGeneSet;

    }


    /*** Summarizes dataset dat, subsetted by the clusters given in clusters based on the cluster type "type"***/
    public static DataSet summarizeData(DataSet dat,ArrayList<Gene> genes, Map<Gene,List<Gene>> clusters, ClusterType type )
    {
        if(type==ClusterType.NONE)
        {
            Gene temp = new Gene(-1);
            temp.symbol="Target";
            genes.add(temp);
            return subset(dat,genes);
        }
        else
        {
            /***Produce one column at a time via subsetting by genes in cluster i and then running PCA on the resulting matrix and taking dimension 1***/

            double[][] td = new double[clusters.keySet().size() + 1][dat.getNumRows()];
            int col = 0;
            List<Node> nodes = new ArrayList<Node>();
            for(Gene g:clusters.keySet()) {
                List<Gene> currGenes = clusters.get(g);
                if(!currGenes.contains(g))
                    currGenes.add(g);

                DataSet temp = subset(dat,(ArrayList<Gene>)currGenes);

                double [] res;
                if(type==ClusterType.PCA) {
                    res = PCA(temp);
                }
                else if(type==ClusterType.MEAN)
                {
                    res = mean(temp);
                }
                else
                {
                    res = median(temp);
                }

                td[col] = res;

                col++;
                String name = "";
                for(Gene x:currGenes)
                    name+= "," + x.symbol;
                nodes.add(new ContinuousVariable(name));
            }
            /***Transpose summarized data to prepare for conversion to Dataset***/
            td = new TetradMatrix(td).transpose().toArray();

            for(int i = 0; i < td.length;i++)
            {
                td[i][td[i].length-1] = dat.getDouble(i,dat.getColumn(dat.getVariable("Target")));
            }
            nodes.add(new ContinuousVariable("Target"));
            DataSet finalData = new BoxDataSet(new DoubleDataBox(td),nodes);
            return finalData;
        }
    }

    private static DataSet subset(DataSet data, ArrayList<Gene> genes)
    {
        int [] cols = new int[genes.size()];
        for(int i = 0; i < genes.size();i++)
        {
            cols[i] = data.getColumn(data.getVariable(genes.get(i).symbol));
        }
        return data.subsetColumns(cols);
    }
    private DataSet[] summarize(DataSet train,DataSet test, ArrayList<Gene> genes, Map<Gene,List<Gene>> clusters, ClusterType type)
    {

        if(type==ClusterType.NONE)
        {
            return new DataSet[]{subset(train,genes),subset(test,genes)};
        }
        else
        {
            /***Produce one column at a time via subsetting by genes in cluster i and then running PCA on the resulting matrix and taking dimension 1***/

            double[][] trainData = new double[clusters.keySet().size() + 1][train.getNumRows()];
            double[][] testData = new double[clusters.keySet().size()+1][test.getNumRows()];
            int col = 0;
            List<Node> nodes = new ArrayList<Node>();
            for(Gene g:clusters.keySet()) {
                List<Gene> currGenes = clusters.get(g);
                currGenes.add(g);

                DataSet temp = subset(train,(ArrayList<Gene>)currGenes);
                DataSet temp2 = subset(test,(ArrayList<Gene>) currGenes);

                double [] res;
                double [] res2;
                if(clusterType==ClusterType.PCA) {
                    res = PCA(temp);
                    res2 = PCA(temp2);
                }
                else if(clusterType==ClusterType.MEAN)
                {
                    res = mean(temp);
                    res2 = mean(temp2);
                }
                else
                {
                    res = median(temp);
                    res2 = median(temp2);
                }

                trainData[col] = res;
                testData[col] = res2;

                col++;
                String name = g.symbol;
                for(Gene x:currGenes)
                    name+= "," + x.symbol;
                nodes.add(new ContinuousVariable(name));
            }
            /***Transpose summarized data to prepare for conversion to Dataset***/
            trainData = new TetradMatrix(trainData).transpose().toArray();
            testData = new TetradMatrix(testData).transpose().toArray();


            for(int i = 0; i < trainData.length;i++)
            {
                trainData[i][trainData[0].length-1] = train.getDouble(i,train.getColumn(train.getVariable("Target")));
            }

            for(int i = 0; i < testData.length;i++)
            {
                testData[i][testData[0].length-1] = test.getDouble(i,test.getColumn(test.getVariable("Target")));
            }
            nodes.add(new ContinuousVariable("Target"));
            DataSet trainFinal = new BoxDataSet(new DoubleDataBox(trainData),nodes);
            DataSet testFinal = new BoxDataSet(new DoubleDataBox(testData),nodes);
            return new DataSet[]{trainFinal,testFinal};
        }
    }

    private static double [] mean(DataSet temp)
    {
        double [] res = new double[temp.getNumRows()];
        for(int i = 0; i < res.length;i++)
        {
            for(int j = 0; j < temp.getNumColumns();j++)
            {
                res[i]+=temp.getDouble(i,j);
            }
            res[i] /= temp.getNumColumns();
        }
        return res;
    }
    private static double [] median(DataSet temp)
    {
        double [] res = new double[temp.getNumRows()];
        for(int i = 0; i < res.length;i++)
        {
            double [] curr = new double[temp.getNumColumns()];
            for(int j = 0; j < temp.getNumColumns();j++)
            {
                curr[j] = temp.getDouble(i,j);
            }
            res[i] = StatUtils.median(curr);
        }
        return res;
    }



    private static double [] PCA(DataSet temp)
    {
            RealMatrix realMatrix = MatrixUtils.createRealMatrix(temp.getDoubleData().toArray());

            //create covariance matrix of points, then find eigen vectors
            //see https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues

            Covariance covariance = new Covariance(realMatrix);
            RealMatrix covarianceMatrix = covariance.getCovarianceMatrix();
            EigenDecomposition ed = new EigenDecomposition(covarianceMatrix);
            double[] weights = ed.getEigenvector(0).toArray();
            double[] result = new double[temp.getNumRows()];
            for (int i = 0; i < temp.getNumRows(); i++) {
                double res = 0;
                for (int j = 0; j < temp.getNumColumns(); j++) {
                    res += weights[j] * temp.getDouble(i, j);
                }
                result[i] = res;
            }
            return result;
    }


    //TODO Use all genes significantly correlated to those selected as clusters
    //Compute PCA or Median or whatever of selected genes
    //Use the constructed genes for downstream prediction computation

    //Input: alpha (data-theory tradeoff),
    //Input: train (Dataset)
    //Output: Stability score, either cluster similarity between runs of PD or how similar selected genes are
    //Output: lastGeneSet, topK most stable genes (those that appear most often in repeateded runs of PD)
    //Output: lastClusters, clusters selected using PD with current alpha value (no subsampling)
    private double pdCPSS(final double alpha, final DataSet train)
    {

        System.out.println("Running CPSS");
        //Now allInds has the complementary pairs, and we can go ahead and perform stability selection

        int nr = (int)Math.floor(train.getNumRows()/2);
        final ArrayList<Integer> listInds = new ArrayList<Integer>();
        for(int i = 0; i < 2*nr;i++)
        {
            listInds.add(i);
        }
        final int [][][] allInds = new int[B][2][nr];
        for(int i = 0; i < B;i++)
        {
            Collections.shuffle(listInds);
            for(int k = 0; k < nr;k++)
            {
                allInds[i][0][k] = listInds.get(k);
            }
            for(int k = nr;k < 2*nr;k++)
            {
                allInds[i][1][k-nr] = listInds.get(k);
            }

        }

        final List<List<Gene>> allGenes = new ArrayList<List<Gene>>();
        final HashMap<Gene,Integer> geneCount = new HashMap<Gene,Integer>();

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


            private synchronized void addCounts(List<Gene> x, Map<Gene,Integer>gc)
            {
                for(int i = 0; i < x.size();i++)
                {
                    if(gc.get(x.get(i))!=null)
                        gc.put(x.get(i),gc.get(x.get(i))+1);
                    else
                        gc.put(x.get(i),1);
                }
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

                        DataSet dataSubSamp = train.subsetRows(allInds[s/2][s%2]);

                        ArrayList<Gene> curr;
                        if(usePThreshold)
                        {
                            if(separateParams)
                                curr = Functions.computeAllIntensities(genes,alpha,dataSubSamp,target,partialCorr,pThresh,pThreshWP,withPrior);
                            else
                                curr = Functions.computeAllIntensities(genes,alpha,dataSubSamp,target,partialCorr,false,false,pThresh);
                        }
                        else
                        {
                            curr = Functions.computeAllIntensities(genes,alpha,dataSubSamp,target,partialCorr,false,false,1);
                        }
                        long time = System.nanoTime();
                        //PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alp,dataSubSamp,true);
                        //p.setCluster(true);
                        //ArrayList<Gene> result = p.diverset();
                        time = System.nanoTime()-time;
                        // System.out.println(s + " Computing all correlations " + time/Math.pow(10,9));

                        time = System.nanoTime();
                        Collections.sort(curr,Gene.IntensityComparator);
                        PrefDiv p;
                        if(usePThreshold)
                        {
                            if(separateParams)
                            {
                                float[] corrs = Functions.computeAllCorrelations(curr, dataSubSamp, partialCorr, pThresh,pThreshWP,withPrior);
                                p = new PrefDiv(curr, topK, accuracy, radius, radiusWP,withPrior,corrs);
                            }
                            else {
                                float[] corrs = Functions.computeAllCorrelations(curr, dataSubSamp, partialCorr, false, false, pThresh);
                                p = new PrefDiv(curr, topK, accuracy, radius, corrs);
                            }

                        }
                        else
                        {
                           p = new PrefDiv(curr,topK,accuracy,radius,dissimilarity,alpha,dataSubSamp,false,partialCorr);
                        }
                        //PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alpha,dataSubSamp,approxCorrelations,partialCorr);
                        p.setCluster(clusterByCorrs);
                        ArrayList<Gene> result = p.diverset();
                        allClusters.add(p.clusters);

                        time = System.nanoTime()-time;
                        // System.out.println(s + " Computing on the fly " + time/Math.pow(10,9));


                        // ArrayList<Gene> result = p.diverset();
                        addToFullSet(result,s);
                        addCounts(result,geneCount);
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

        pool.invoke(new StabilityAction(chunk, 0, B*2));

        //Extract the stable genes using the counts hashmap (topK stable genes)
        System.out.println("Gene appearence map: " + geneCount);
        lastGeneSet = stableSet(geneCount);
        System.out.println("Top K Genes: " + lastGeneSet);


        if(usePThreshold)
        {
            ArrayList<Gene> curr;
            if(separateParams)
                curr = Functions.computeAllIntensities(genes,alpha,train,target,partialCorr,pThresh,pThreshWP,withPrior);
            else
                curr = Functions.computeAllIntensities(genes,alpha,train,target,partialCorr,false,false,pThresh);
            Collections.sort(curr,Gene.IntensityComparator);
            float [] corrs;
            if(separateParams)
                corrs = Functions.computeAllCorrelations(curr,train,false,pThresh,pThreshWP,withPrior);
            else
                corrs = Functions.computeAllCorrelations(curr,train,false,false,false,pThresh);

            PrefDiv p;
            if(separateParams)
                p = new PrefDiv(curr,topK,accuracy,radius,radiusWP,withPrior,corrs);
            else
                p = new PrefDiv(curr,topK,accuracy,radius,corrs);
            p.setCluster(clusterByCorrs);
            p.diverset();
            lastClusters = p.clusters;
        }
        else
        {
            ArrayList<Gene> curr = Functions.computeAllIntensities(genes,alpha,train,target,partialCorr,false,false,1);
            Collections.sort(curr,Gene.IntensityComparator);
            PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,dissimilarity,alpha,train,false,partialCorr);
            //PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alpha,data,approxCorrelations,partialCorr);
            p.setCluster(clusterByCorrs);
            p.diverset();
            lastClusters = p.clusters;

        }

        //thetaMat.assign(thetaMat.copy().assign(Functions.minus(1.0)), Functions.mult).assign(Functions.mult(-2.0));
        if(useClusterStability)
            return 1-clusterSim(allClusters);
        else
            return tanimotoSim(allGenes);
    }



    //Input: alpha (data-theory tradeoff),
    //Input: train (Dataset)
    //Input: k Subsample index, if less than 0 then just compute the correlations here
    //Output: always -1, just a dummy variable to maintain congruence with other version of PD
    //Output: lastGeneSet, topK genes selected by PD
    //Output: lastClusters, clusters selected using PD with current alpha value

    //TODO If this will be a useful approach, then change this to compute all correlations and output to a file beforehand
    //TODO In here you can just read as is
    private double runPD(final double alpha, final DataSet train)
    {
            ArrayList<Gene> curr = Functions.computeAllIntensities(genes,alpha,train,target,partialCorr,false,false,1);
            Collections.sort(curr,Gene.IntensityComparator);
            PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,dissimilarity,alpha,train,false,partialCorr);
            //PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alpha,data,approxCorrelations,partialCorr);
            p.setCluster(clusterByCorrs);
            lastGeneSet = p.diverset();
            lastClusters = p.clusters;
        return -1;
    }


    private ArrayList<Gene> stableSet(Map<Gene,Integer> map)
    {
        int [] counts = new int[map.keySet().size()];
        List<Gene> genes = new ArrayList<Gene>();
        int i = 0;
        for(Gene g:map.keySet())
        {
            genes.add(g);
            counts[i] = map.get(g);
            i++;
        }
        Comparator c = new GeneComparator(counts,clone(genes));
        Collections.sort(genes,c);

        ArrayList<Gene> result = new ArrayList<Gene>();
        for(int j = 0; j < topK;j++)
        {
            result.add(genes.get(j));
        }
        return result;

    }
    private ArrayList<Gene> clone(List<Gene> l)
    {
        ArrayList<Gene> res = new ArrayList<Gene>();
        for(int i =0; i < l.size();i++)
            res.add(l.get(i));
        return res;
    }
    private class GeneComparator implements Comparator<Gene>
    {
        private int [] counts;
        private List<Gene> gs;

        public GeneComparator(int [] c, List<Gene>g)
        {
            counts = c;
            gs = g;
        }


        @Override
        public int compare(Gene g1, Gene g2)
        {
            if(counts[gs.indexOf(g1)] > counts[gs.indexOf(g2)])
                return -1;
            else if(counts[gs.indexOf(g1)]==counts[gs.indexOf(g2)])
                return 0;
            else
                return 1;
        }
    }
    private double testRegression(RegressionResult res, DataSet test, List<Node> dNeighbors)
    {
        double [] predictions = new double[test.getNumRows()];
        double [] real = new double[test.getNumRows()];
        for(int k = 0; k < test.getNumRows();k++)
        {
            double [] testVec = new double[dNeighbors.size()];
            for(int x = 0; x < dNeighbors.size();x++)
            {
                testVec[x] = test.getDouble(k,test.getColumn(dNeighbors.get(x)));
            }

            real[k] = test.getDouble(k,test.getColumn(test.getVariable(target)));
            predictions[k] = res.getPredictedValue(testVec);
        }

        return RMSE(predictions,real);
    }

    private double RMSE(double [] preds, double [] truth)
    {
        int T = preds.length;
        double sum = 0;
        for(int i = 0; i < preds.length;i++)
        {
            sum+=Math.pow((preds[i]-truth[i]),2);
        }
        return Math.sqrt(sum/T);
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

        STEPS s = new STEPS(temp,lambdas,g,subs);
        Graph output = s.runStepsPar();
        return output;
    }
    public Graph getCausalGraph(String target)
    {
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

        STEPS s = new STEPS(temp,lambdas,g,subs);
        Graph output = s.runStepsPar();
        lastGeneSet = g2;
        lastStepsStabilities = s.stabilities;
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

                        ArrayList<Gene> curr = Functions.computeAllIntensities(genes,alp,dataSubSamp,target,partialCorr,false,false,1);
                        long time = System.nanoTime();
                        //PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alp,dataSubSamp,true);
                        //p.setCluster(true);
                        //ArrayList<Gene> result = p.diverset();
                        time = System.nanoTime()-time;
                       // System.out.println(s + " Computing all correlations " + time/Math.pow(10,9));

                        time = System.nanoTime();
                        Collections.sort(curr,Gene.IntensityComparator);
                        PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,dissimilarity,alp,dataSubSamp,false,partialCorr);
                        //PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alp,dataSubSamp,approxCorrelations,partialCorr);
                        p.setCluster(clusterByCorrs);
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
        ArrayList<Gene> curr = Functions.computeAllIntensities(genes,alp,data,target,partialCorr,false,false,1);
        Collections.sort(curr,Gene.IntensityComparator);
        PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,dissimilarity,alp,data,true,partialCorr);
        //PrefDiv p = new PrefDiv(curr,topK,accuracy,radius,PrefDiv.findTopKIntensity(curr,topK),dissimilarity,alp,data,approxCorrelations,partialCorr);
        p.setCluster(clusterByCorrs);
        lastGeneSet = p.diverset();
        lastClusters = p.clusters;

        //do this elsewhere
        //thetaMat.assign(thetaMat.copy().assign(Functions.minus(1.0)), Functions.mult).assign(Functions.mult(-2.0));
        if(useClusterStability)
            return 1-clusterSim(allClusters);
        else
            return tanimotoSim(allGenes);
    }


    /***
     *
     * @param subs the set of indices for the testing data
     * @return the Dataset consisting of all indices not in subs
     */
    private static DataSet getTrainingData(DataSet data, int[]subs)
    {
        int[] trainInds = new int[data.getNumRows() - subs.length];
        int tempCount = 0;
        int subCount = 0;
        for (int k = 0; k < data.getNumRows(); k++) {
            if (subs[subCount]!=k) {
                trainInds[tempCount] = k;
                tempCount++;
            }
            else
            {
                subCount++;
                if(subCount==subs.length)
                    subCount=0;
            }
        }
        return data.subsetRows(trainInds);
    }


    private static void addToGeneList(HashMap<Gene,List<Gene>> map, List<List<Gene>> one)
    {
        for(Gene g: map.keySet())
        {
            ArrayList<Gene> temp = new ArrayList<Gene>();
            temp.add(g);
            if(map.get(g)!=null) {
                for (Gene x : map.get(g)) {
                    if(!x.symbol.equals(g.symbol))
                     temp.add(x);
                }
            }
            one.add(temp);
        }
    }

    private static double[][] computeCostMatrix(List<List<Gene>> one, List<List<Gene>> two)
    {
        double [][] costs = new double[one.size()][two.size()];
        for(int i = 0; i < one.size();i++)
        {
            for(int j = 0; j < two.size();j++)
            {
                costs[i][j] = 1-(intersect(one.get(i),two.get(j)).size()/(double)union(one.get(i),two.get(j)).size());
            }
        }
        return costs;
    }

    private static List<Gene> intersect(List<Gene> one, List<Gene>two)
    {
        List<Gene> result = new ArrayList<Gene>();
        for(int i = 0; i < one.size();i++)
        {
            for(Gene g:two)
            {
                if(one.get(i).symbol.equals(g.symbol))
                    result.add(g);
            }
        }
        return result;
    }
    private static List<String> union(List<Gene>one, List<Gene>two)
    {
        List<String> result = new ArrayList<String>();
        for(int i = 0; i < one.size();i++)
        {
            if(!result.contains(one.get(i).symbol))
                result.add(one.get(i).symbol);
        }
        for(int i = 0; i < two.size();i++)
        {
            if(!result.contains(two.get(i).symbol))
                result.add(two.get(i).symbol);
        }
        return result;
    }
    private static double [][] normalizedCosts(double[][]costs)
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

    public static double clusterSim(Map<Gene,List<Gene>> m1, Map<Gene,List<Gene>> m2)
    {
        List<List<Gene>> one = new ArrayList<List<Gene>>();
        Map<Gene,List<Gene>> map = m1;
        addToGeneList((HashMap<Gene,List<Gene>>)map,one);

        List<List<Gene>> two = new ArrayList<List<Gene>>();
        map = m2;
        addToGeneList((HashMap<Gene,List<Gene>>)map,two);

        double [][] costs = computeCostMatrix(one,two);
        HungarianAlgorithm ha = new HungarianAlgorithm(costs);




        int [] result = ha.execute();
        //costs = normalizedCosts(costs);
        double sim = 0;
        for(int k = 0; k < result.length;k++)
        {
            sim+=costs[k][result[k]];
        }
        return 1-(sim/result.length);
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
