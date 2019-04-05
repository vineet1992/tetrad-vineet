package edu.pitt.csb.Pref_Div.Comparisons;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.pitt.csb.Pref_Div.Gene;
import edu.pitt.csb.Pref_Div.PiPrefDiv4;
import edu.pitt.csb.Pref_Div.RunPrefDiv;

import java.util.*;

/**
 * Created by vinee_000 on 2/21/2019.
 * The purpose of this class is to produce output statistics for Pi-PrefDiv on simulated datasets
 */
public class PrefDivComparator {

    private ComparablePD est; //Estimated PiPrefDiv results
    private String target; //Target variable of interest
    private Graph trueGraph; //True Data Generating Graph
    private Map<String,List<Integer>> trueClusters; //True Clusters

    public PrefDivComparator(ComparablePD pd, String target, Graph trueGraph,Map<String,List<Integer>> trueClusters)
    {
        this.est = pd;
        this.target = target;
        this.trueGraph = trueGraph;
        this.trueClusters = trueClusters;
    }


    /***
     *
     * @param numClusters Total number of clusters
     * @return Friendly cluster accuracy looks first at the best matched cluster for each estimated cluster and then takes average I/U
     */
    public double getFriendlyCluster(int numClusters)
    {
        HashMap<Gene,List<Gene>> real = createRealCluster(numClusters);

        /***Find Accuracy mapping using cost matrix***/

        ClusterSim cs = new ClusterSim(est.getLastCluster(),real);
        return cs.friendlyClusterSim();
    }


    /***
     *
     * @return The number of clusters represented by the clusters selected by the PrefDiv object
     */
    public int getNumClusters()
    {
        /**Estimated clusters**/
        Map<Gene,List<Gene>> estClusters = this.est.getLastCluster();

        /***Unique clusters touched***/
        HashSet<Integer> set = new HashSet<Integer>();


        /***Load up all of the unique clusters touched***/
        for(Gene g:estClusters.keySet())
        {
            List<Integer> temp = trueClusters.get(g.symbol);
            set.addAll(temp);

            for(Gene gPrime: estClusters.get(g))
            {
                temp = trueClusters.get(gPrime.symbol);
                set.addAll(temp);
            }
        }

        /***Return this set***/
        return set.size();

    }

    private Set<Integer> getTrueClustIndices()
    {


        Set<Integer> result = new HashSet<Integer>();
        Set<String> parents = new HashSet<String>();

        for(Node n: trueGraph.getParents(trueGraph.getNode(target)))
        {
            parents.add(n.getName());
        }


        for(String s: trueClusters.keySet())
        {
            if(parents.contains(s))
            {
                result.addAll(trueClusters.get(s));
            }
        }
        return result;
    }
    /***
     *
     * @param numReal Number of causally connected clusters
     * @return The number of true clusters touched by
     */
    public int getTrueClusters(int numReal)
    {
        /**Estimated clusters**/
        Map<Gene,List<Gene>> estClusters = this.est.getLastCluster();

        /***Real clusters***/
        Map<Gene,List<Gene>> realClusters = createRealCluster(numReal);

        /***Convert representation type***/
        Map<String,List<Integer>> trueClusters = convertRepresentation(realClusters);

        /***Which cluster indices are causally related?***/
        Set<Integer> trueIndices = getTrueClustIndices();

        //X34 -> 1,2,3

        HashSet<Integer> covered = new HashSet<Integer>();


        /***Loop over estimated clusters***/
        for(Gene g:estClusters.keySet())
        {
            List<Gene> temp = estClusters.get(g);
            temp.add(g);

            List<Integer> indices = new ArrayList<Integer>();

            for(Gene gene: temp)
            {
                /***If this gene isn't in any of the true causal clusters then this will be null***/
                if(trueClusters.get(gene.symbol)!=null)
                    indices.addAll(trueClusters.get(gene.symbol));
            }

            /***Keep only the causal cluster indices***/
            indices.retainAll(trueIndices);


            covered.addAll(indices);
        }

        return covered.size();


    }

    /***
     *
     * @param numClusters Total number of clusters
     * @return cluster accuracy metric via Hungarian algorithm
     */
    public double getClusterAccuracy(int numClusters)
    {
        HashMap<Gene,List<Gene>> real = createRealCluster(numClusters);
        ClusterSim cs = new ClusterSim(est.getLastCluster(),real);
        return cs.clusterSim();
    }


    /***
     * Helper function to convert gene to list of gene to a string to list of integer representation of the clusters
     * @return
     */
    private Map<String,List<Integer>> convertRepresentation(Map<Gene,List<Gene>> input)
    {

        Map<String,List<Integer>> result = new HashMap<String,List<Integer>>();
        int x = 0;
        for(Gene g:input.keySet())
        {
            List<Integer> temp = new ArrayList<Integer>();
            if(result.get(g.symbol)!=null)
            {
                temp = result.get(g.symbol);
            }
            temp.add(x);
            result.put(g.symbol,temp);
            for(Gene gPrime: input.get(g))
            {
                temp = new ArrayList<Integer>();
                if(result.get(gPrime.symbol)!=null)
                {
                    temp = result.get(gPrime.symbol);
                }
                temp.add(x);
                result.put(gPrime.symbol,temp);

            }
            x++;
        }

        return result;

    }
    /**
     *
     * @param numClusters Total number of clusters
     * @return Converts true clusters to a format where it is Gene to List of Gene only for clusters which affect the target
     */
    private HashMap<Gene,List<Gene>> createRealCluster(int numClusters)
    {
        /***Create list of genes representing each actual pathway***/
        List<List<Gene>> temp = new ArrayList<List<Gene>>();
        for(int i = 0; i < numClusters;i++)
            temp.add(new ArrayList<Gene>());
        int x = 0;

        /***Loop through the true clusters and convert them to a single list of integers for each cluster***/
        for(String s: trueClusters.keySet())
        {
            Gene g = new Gene(x);
            g.symbol = s;
            List<Integer> currPaths = trueClusters.get(s);
            for(int i = 0; i < currPaths.size();i++)
            {
                temp.get(currPaths.get(i)).add(g);
            }
            x++;
        }

        boolean [] clustsToKeep = new boolean[numClusters];
        /***Identify clusters with at least one member causally connected to the target***/
        for(Node n: trueGraph.getAdjacentNodes(trueGraph.getNode(target)))
        {
            List<Integer> causalCluster = trueClusters.get(n.getName());
            for(Integer cause: causalCluster)
                clustsToKeep[cause] = true;

        }

        /****Remove Pathways from true list if they aren't causally connected***/

        for(int i = clustsToKeep.length-1;i >= 0; i--)
        {
            if(!clustsToKeep[i])
                temp.remove(i);
        }

        /***Convert actual clusters to Gene -> List<Gene> mapping***/
        HashMap<Gene,List<Gene>> real = new HashMap<Gene,List<Gene>>();
        for(int i = 0; i < temp.size();i++)
        {
            List<Gene> curr = temp.get(i);
            Gene key = curr.get(0);
            curr.remove(0);
            real.put(key,curr);
        }
        return real;
    }


    /***
     *
     * @return What percent of the true causal clusters are represented by the result
     * NOTE this does not handle genes being in multiple clusters (arbitrarily chooses their first one)
     */
    public double getRepresentAccuracy()
    {

        int trueLength = 0; /***How many clusts are represented in the truth?***/
        ArrayList<Integer> trueClusts = new ArrayList<Integer>();


        /***Get true cluster assignments***/
        Node tNode = trueGraph.getNode(target);
        for(Node n: trueGraph.getAdjacentNodes(tNode))
        {
            trueClusts.add(trueClusters.get(n.getName()).get(0));
        }

        trueLength = trueClusts.size();
        ArrayList<Integer> estClusts = new ArrayList<Integer>();

        /***Get estimated cluster assignments***/
        for(Gene n: est.getLastSelected())
        {
            estClusts.add(trueClusters.get(n.symbol).get(0));
        }

        /***Compute how many clusters in the truth are represented in the est***/

        int count = 0;
        for(int i = 0; i < estClusts.size();i++)
        {
            if(trueClusts.contains(estClusts.get(i))) {
                trueClusts.remove(estClusts.get(i));
                count++;
            }
        }
        return count/(double)trueLength;

    }


    /***
     *
     * @param type ACC or F1 (use ACC if the number of selected features matches the true number of causal features) else use F1
     * @return The accuracy or F1 score of the selected features against the truth
     */
    public double getFeatAccuracy(String type)
    {
        double acc = 0;
        double tp = 0;
        double fp = 0;
        double fn = 0;
        List<Gene> selected = est.getLastSelected();
        List<Node> correct = trueGraph.getAdjacentNodes(trueGraph.getNode(target));
        for(int i = 0; i < selected.size();i++)
        {
            boolean found = false;
            for(Node n: correct)
            {
                if(n.getName().equals(selected.get(i).symbol))
                {
                    acc++;
                    tp++;
                    found = true;
                }

            }
            if(!found)
            {
                fp++;
            }
        }
        if(type.equals("ACC"))
            return acc/correct.size();

        for(int i = 0; i < correct.size();i++)
        {
            Node curr = correct.get(i);
            boolean found = false;

            for(int j = 0; j < selected.size();j++)
            {
                if(selected.get(j).symbol.equals(curr.getName()))
                    found = true;
            }
            if(!found)
                fn++;
        }

        double prec = tp/(tp+fp);
        double rec = tp/(tp+fn);

        return (2*(prec*rec)/(prec+rec));

    }

    /***
     *
     * @param train Training set where feature selection was run
     * @param test Testing set to compute prediction accuracy
     * @return RMSE error of predicting the target variable on the test samples using the train and test sets specified
     */
    public double getPredictionAccuracy(DataSet train, DataSet test)
    {
        //Extract last set of selected genes
            List<Gene> selected = est.getLastSelected();




            //Create regression dataset and use selected features as regressors
            RegressionDataset rd = new RegressionDataset(train);
            List<Node> regressors = new ArrayList<Node>();
            if(selected!=null) {
                for (int i = 0; i < selected.size(); i++)
                    if(train.getVariable(selected.get(i).symbol)!=null)
                        regressors.add(train.getVariable(selected.get(i).symbol));
            }
            //If there were no selected features (e.g. if causal modeling was used intermediately, then use all of them)
            if(selected==null || regressors.size()==0)
            {
                for (int i = 0; i < train.getNumColumns(); i++) {
                    if(train.getVariable(i).getName().equals("Target"))
                        continue;
                    regressors.add(train.getVariable(i));
                }
            }

            //Learn a regression model on the train features
            RegressionResult res = rd.regress(train.getVariable("Target"),regressors);
            int [] cols = new int[regressors.size()];

            for(int i = 0; i < regressors.size();i++)
            {
                cols[i] = test.getColumn(test.getVariable(regressors.get(i).getName()));
            }


        //Apply it to the test set

        double[]pred = new double[test.getNumRows()];
            double [] actual = new double[test.getNumRows()];
            for(int i = 0; i < test.getNumRows();i++)
            {
                double [] x = new double[cols.length];
                for(int j = 0; j < cols.length;j++) {
                    x[j] = test.getDouble(i, cols[j]);

                }
                pred[i] = res.getPredictedValue(x);
                actual[i] = test.getDouble(i,test.getColumn(test.getVariable("Target")));


            }


            //Compute Root-mean-squared error

            return RMSE(pred,actual);

    }

    private static double RMSE(double [] pred, double[] actual)
    {
        double err = 0;
        for(int i = 0; i < pred.length;i++)
        {
            err += Math.pow(pred[i]-actual[i],2);
        }
        return Math.sqrt(err/pred.length);
    }

}
