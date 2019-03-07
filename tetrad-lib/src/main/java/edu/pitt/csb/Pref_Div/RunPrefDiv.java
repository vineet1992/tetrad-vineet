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

import static edu.pitt.csb.Pref_Div.Comparisons.ClusterSim.tanimotoSim;
import static org.apache.commons.collections4.CollectionUtils.intersection;
import static org.apache.commons.collections4.CollectionUtils.union;

/**
  Runs genomic Pref-Div algorithm with subsampling and clustering
    TODO Handle Discrete features in theory relevance and dissimilarity

 */
public class RunPrefDiv {

    public enum ClusterType{
        PCA,MEAN,MEDIAN,NONE
    }

    private DataSet data;
   private ArrayList<Gene> genes;
    private float[] dissimilarity;
    private int numSubs = 20; //Number of subsamples
    private int topK = 50;
    private double accuracy = 0;
    private double radius = 0.5;
    private int numAlphas = 20; //How many alpha values should we test? (Data-Theory tradeoff alpha)
    private double alphaLimit = 1; //The largest value of alpha we test
   private double radiusWP = 0.5; //Radius for those which we have prior information (only if separate params is true)

    private boolean partialCorr = false;
    private boolean [] withPrior; //Says which edges we have prior information for

    private ArrayList<Gene> lastGeneSet;
    private Map<Gene,List<Gene>> lastClusters;

    private boolean loocv = false;
    private int [][] subs;
    private DataSet summarizedData; //The dataset you get using the prefidv method on the current dataset with the current clustering procedure



    private boolean clusterByCorrs = false; //Should we cluster purely by significant correlations instead of using Pref-Div
    private ClusterType clusterType; //Which clustering method should we use (Options are PCA, Median, Mean, None)

    public RunPrefDiv(float [] dissimilarity, ArrayList<Gene> genes, DataSet data,boolean leaveOneOut)
    {
        this.data = data;
        this.genes = genes;
        this.dissimilarity = dissimilarity;
        this.loocv = leaveOneOut;
        this.clusterType = ClusterType.NONE;
    }


    public void setWithPrior(boolean [] wp)
    {
        withPrior = wp;
    }
    public void setAllParams(double radius, double radiusWP)
    {
        this.radius = radius;
        this.radiusWP = radiusWP;

    }


    public int [][] getSubs()
    {
        return subs;
    }
    public void setSubs(int [][] subs){this.subs = subs;}
    public Map<Gene,List<Gene>> getClusters()
    {
        return lastClusters;
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
    public void usePartialCorrelation(boolean pc){partialCorr = pc;}
    public void clusterByCorrs(){clusterByCorrs=true;}
    public void setClusterType(ClusterType c){clusterType = c;}
    public DataSet getSummarizedData(){return summarizedData;}





    public ArrayList<Gene> runPD()
    {


        if(subs==null) {
            System.out.print("Generating Subsamples...");

            if (loocv)
                subs = StabilityUtils.generateSubsamples(data.getNumRows());
            else
                subs = StabilityUtils.generateSubsamples(numSubs,data.getNumRows());
            System.out.println("Done");

        }
        /***If we aren't using cross validation then no alpha value is required, just run PrefDiv as is***/
        /***This will only be used if we are running PiPrefDiv***/

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

        return lastGeneSet;

    }


    /*** Summarizes dataset dat, subsetted by the clusters given in clusters based on the cluster type "type"***/
    public static DataSet summarizeData(DataSet dat,ArrayList<Gene> genes, Map<Gene,List<Gene>> clusters, ClusterType type )
    {
        if(type==ClusterType.NONE)
        {
            ArrayList<Gene> tempList = new ArrayList<Gene>(genes);
            Gene temp = new Gene(-1);
            temp.symbol="Target";
            tempList.add(temp);
            return subset(dat,tempList);
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
                    name+= "|" + x.symbol;

                name = name.substring(1);
                nodes.add(new ContinuousVariable(name));
            }
            /***Transpose summarized data to prepare for conversion to Dataset***/
            td = new TetradMatrix(td).transpose().toArray();

            for(int i = 0; i < td.length;i++)
            {
                td[i][td[i].length-1] = dat.getDouble(i,dat.getColumn(dat.getVariable("Target")));
            }
            nodes.add(new ContinuousVariable("Target"));


            DataSet finalData = new ColtDataSet(td.length,nodes);
            for(int i = 0; i < td.length;i++)
            {
                for(int j = 0; j < td[i].length;j++)
                {
                    finalData.setDouble(i,j,td[i][j]);
                }
            }
            return finalData;
        }
    }

    private static List<Gene> removeDuplicates(List<Gene> genes)
    {
        List<Gene> result = new ArrayList<Gene>();
        Map<String,Integer> set = new HashMap<String,Integer>();
        for(Gene g:genes)
        {
            if(set.get(g.symbol)==null)
            {
                result.add(g);
                set.put(g.symbol,0);
            }
        }

        return result;
    }
    private static DataSet subset(DataSet data, ArrayList<Gene> genes)
    {
        //First get only unique elements of the list
        List<Gene> newList = removeDuplicates(genes);

        int [] cols = new int[newList.size()];
        for(int i = 0; i < newList.size();i++)
        {
            cols[i] = data.getColumn(data.getVariable(newList.get(i).symbol));
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








}
