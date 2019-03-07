package edu.pitt.csb.Pref_Div.Comparisons;

import edu.pitt.csb.Pref_Div.Gene;

import java.util.*;

import static org.apache.commons.collections4.CollectionUtils.intersection;

/**
 * Created by vinee_000 on 2/21/2019.
 * This is a wrapper class to compute cluster similarities between two sets of clusters
 */
public class ClusterSim {

    private Map<Gene,List<Gene>> m1;
    private Map<Gene,List<Gene>> m2;

    public ClusterSim(Map<Gene,List<Gene>> one, Map<Gene,List<Gene>> two)
    {
        this.m1 = one;
        this.m2 = two;
    }


    /***
     *
     * @param one First list of gene lists
     * @param two Second list of gene lists
     * @return Computes the cost of each matching between pathways in set one and set two
     */
    public static double[][] computeCostMatrix(List<List<Gene>> one, List<List<Gene>> two)
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

    /***
     *
     * @return The union between two gene lists as a list of strings
     */
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


    /***
     *
     * @return The intersection between two gene lists as a list of genes
     */
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


    /***
     *
     * @param map
     * @param one
     */
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



    /**
     *Utility method to convert a map from genes to list of genes to a list of lists
     */
    private List<List<Gene>> convertToGeneList(Map<Gene,List<Gene>> m)
    {
        List<List<Gene>> one = new ArrayList<List<Gene>>();
        Map<Gene,List<Gene>> map = m;
        addToGeneList((HashMap<Gene,List<Gene>>)map,one);
        return one;
    }

    /***
     *
     * @return Returns the cluster similarity between the two specified cluster sets
     */
    public double clusterSim()
    {
        List<List<Gene>> one = convertToGeneList(m1);

        List<List<Gene>> two = convertToGeneList(m2);

        double [][] costs = computeCostMatrix(one,two);
        return getSimilarity(costs);
    }

    private static double getSimilarity(double [][] costs)
    {
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

    public double friendlyClusterSim()
    {
        List<List<Gene>> one = convertToGeneList(m1);
        List<List<Gene>> two = convertToGeneList(m2);
        double [][] costs = computeCostMatrix(one,two);
        /***This mapping consists of merged clusters if two cluster sets like the same real pathway***/
        HashMap<Integer,List<Gene>> mergedClusters = new HashMap<Integer,List<Gene>>();
        //For each estimated cluster
        for(int i = 0; i < costs.length;i++)
        {
            int bestIndex = -1;
            double lowestCost = 1;
            //Look at the true clusters
            for(int j = 0; j< costs[i].length;j++)
            {
                if(lowestCost > costs[i][j])
                {
                    lowestCost = costs[i][j];
                    bestIndex = j;
                }
            }

            /***We should match this pathway to the j'th cluster***/
            if(mergedClusters.get(bestIndex)==null)
            {
                mergedClusters.put(bestIndex,one.get(i));
            }
            else
            {
                List<Gene> temp = mergedClusters.get(bestIndex);
                temp.addAll(one.get(i));
                mergedClusters.put(bestIndex,temp);
            }

        }

        /***Convert Int to List of Gene to Gene to List of Gene***/
        HashMap<Gene,List<Gene>> newMap = new HashMap<Gene,List<Gene>>();
        for(Integer x: mergedClusters.keySet())
        {
            List<Gene> curr = mergedClusters.get(x);
            Set<String> gNames = new HashSet<String>();
            List<Gene> finalList = new ArrayList<Gene>();
            /**Ensure there are no duplicate genes in this list***/
            for(int i = 0; i < curr.size();i++)
            {
                if(gNames.add(curr.get(i).symbol))
                {
                    finalList.add(curr.get(i));
                }
            }
            Gene first = finalList.get(0);
            finalList.remove(0);
            newMap.put(first,finalList);
        }
        return getSimilarity(computeCostMatrix(convertToGeneList(newMap),two));
    }


    /**
     *
     * @return Computes average tanimoto set similarity across a list of gene lists
     */
    public static double tanimotoSim(List<List<Gene>> allGenes)
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
