package edu.pitt.csb.mgm;

import edu.cmu.tetrad.cluster.KMeans;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradVector;

import java.util.Arrays;
import java.util.Comparator;

/**
 * Created by vinee_000 on 9/1/2016.
 */
public class ReKS {

    private DataSet data;
    private final int numNeighbors = 15;
    public ReKS(DataSet d)
    {
        data = d;
    }
    public void getTree()
    {
        TetradMatrix t = getNN(data,numNeighbors);
        //need to perform knn search somehow

    }
    private TetradMatrix getNN(DataSet d, int n)
    {
        TetradMatrix t = new TetradMatrix(d.getNumColumns(),n);
        TetradMatrix dat = d.getDoubleData();
        TetradMatrix distances = new TetradMatrix(d.getNumColumns(),d.getNumColumns());
        for(int i = 0; i < d.getNumColumns();i++)
        {
            for(int j = i; j < d.getNumColumns();j++)
            {
                distances.set(i,j,dist(dat.getColumn(i),dat.getColumn(j)));
                distances.set(j,i,distances.get(i,j));
            }
        }
        for(int i = 0; i < d.getNumColumns();i++)
        {
            t.assignRow(i,getVector(i,distances,n));
        }
//TODO
return t;
    }
    private class MyComparator implements Comparator<Integer> {

        private double [] arr2;

        public MyComparator(double [] a)
        {
            arr2 = a;
        }
        @Override
        public int compare(Integer i1, Integer i2) {
            return Double.compare(arr2[i1.intValue()],arr2[i2.intValue()]);
        }
    }
private TetradVector getVector(int ind, TetradMatrix dists, int nn)
{
    TetradVector v = new TetradVector(nn); //array to store final answer
    double [] curr = dists.getRow(ind).toArray(); //current row of distances
    Integer [] inds = new Integer[curr.length];
    for(int i = 0; i < inds.length;i++)
        inds[i] = i;
    Arrays.sort(inds,new MyComparator(curr));

    for(int i = 0; i < nn; i ++)
        v.set(i,inds[i]);
    //need to get the nn indices with smallest numbers in curr
    return v;
}
    private double dist(TetradVector v1, TetradVector v2)
    {
        if(v1.equals(v2))
            return 0;
        double sum = 0;
        for(int i = 0; i < v1.size();i++)
        {
            sum+= Math.pow(v1.get(i) - v2.get(i),2);
        }
        return Math.sqrt(sum);
    }
    //Get a matrix that is numVars x constant size
    //M(i,:) gives the top constant neighbors of variable i
    //Use clustering (Kmeans or KNN for this)
    //Anti C is another matrix that contains the actual distances between the neighbor specified in M(i,j)
    //1 - actual distance to be precise

    //Distance metric used is correlation
    //Then do 1-antiC to get C

    //affinity matrix is given by
    //e^(-1*(sin(arccos(C)).^2)

    //Make into numVars x numVars matrix with only entries in the places were distances were computed (top 15 columns for each row)
    //for i = 1:numVars
    //S_Sparse(i
    //Finally make the square matrix symmetric

}
