package edu.pitt.csb.mgm;

import edu.cmu.tetrad.graph.NodeType;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradVector;

import java.util.Arrays;
import java.util.Comparator;

/**
 * Created by vinee_000 on 9/9/2016.
 */
public class IGCI {
    public enum types{
        GAUSSIAN,UNIFORM,INTEGRAL,ENTROPY
    }
    private types rm;
    private types est;
    private TetradVector x;
    private TetradVector y;
    private double f;
    private final double GAMMA = 0.577215664901532860606512090082;
    private final double GAMMA_MINX = 1.e-12;
    private final double DIGAMMA_MINNEGX = -1250;
    private final double C_LIMIT = 49;
    private final double S_LIMIT = 1e-5;
    public IGCI(TetradVector x,TetradVector y,types refMeasure,types estimator)
    {
       rm = refMeasure;
        est = estimator;
        this.x = x;
        this.y = y;

    }

    public double getDirection()
    {
        TetradVector ytemp = y.copy();
        TetradVector xtemp = x.copy();
        double f = 0;
        if(rm==types.GAUSSIAN)
        {
            double [] mstdX = getMeanStd(xtemp);
            double [] mstdY = getMeanStd(ytemp);
            System.out.println(Arrays.toString(mstdX));
            for(int i = 0; i < xtemp.size();i++)
            {
                xtemp.set(i,(xtemp.get(i)-mstdX[0])/mstdX[1]);
                ytemp.set(i,(ytemp.get(i)-mstdY[0])/mstdY[1]);
            }

        }
        if(rm==types.UNIFORM)
        {
            double [] maxminX = getMaxMin(xtemp);
            double denomX = maxminX[0]-maxminX[1];
            for(int i = 0; i < xtemp.size();i++)
            {
                xtemp.set(i,(xtemp.get(i)-maxminX[1])/denomX);
            }
            double [] maxminY = getMaxMin(ytemp);
            double denomY = maxminY[0] - maxminY[1];
            for(int i = 0; i < ytemp.size();i++) {
                ytemp.set(i, (ytemp.get(i) - maxminY[1]) / denomY);
            }
        }
        if(est==types.ENTROPY)
        {
            double [] xtemp2 = xtemp.toArray();
            Arrays.sort(xtemp2);
            xtemp = new TetradVector(xtemp2);
            double [] ytemp2 = ytemp.toArray();
            Arrays.sort(ytemp2);
            ytemp = new TetradVector(ytemp2);
            int n1 = xtemp.size();
            double hx = 0;
            double hy = 0;
            for(int i =0; i < n1-1;i++)
            {
                if(xtemp.get(i+1)-xtemp.get(i)!= 0)
                    hx += Math.log(xtemp.get(i+1)-xtemp.get(i));
                if(ytemp.get(i+1)-ytemp.get(i)!=0)
                    hy += Math.log(ytemp.get(i+1)-ytemp.get(i));
            }
            hx = hx/(n1-1) + digamma(n1) - digamma(1);
            hy = hy/(n1-1) + digamma(n1)-digamma(1);
            f = hy-hx;
        }
        else //INTEGRAL FORMAT
        {
            double a = 0;
            double b = 0;
           /* double [] xtemp2 = xtemp.toArray();
            Arrays.sort(xtemp2);

            TetradVector x_new = new TetradVector(xtemp2);
            double [] ytemp2 = ytemp.toArray();
            Arrays.sort(ytemp2);
            TetradVector y_new = new TetradVector(ytemp2);*/
            Integer [] xInd = new Integer[xtemp.size()];
            Integer [] yInd = new Integer [ytemp.size()];
            for(int i = 0; i < xtemp.size();i++)
            {
                xInd[i] = i;
                yInd[i] = i;
            }
            Arrays.sort(xInd,new MyComparator(xtemp.toArray()));
            Arrays.sort(yInd,new MyComparator(ytemp.toArray()));
            int n1 = xtemp.size();
            for(int i = 0; i < n1-1;i++)
            {
                double x1 = xtemp.get(xInd[i]);
                double x2 = xtemp.get(xInd[i+1]);
                double y1 = ytemp.get(xInd[i]);
                double y2 = ytemp.get(xInd[i+1]);
                if(x1!=x2 && y1!=y2)
                    a+= Math.log(Math.abs((y2-y1)/(x2-x1)));
                x1 = xtemp.get(yInd[i]);
                x2 = xtemp.get(yInd[i+1]);
                y1 = ytemp.get(yInd[i]);
                y2 = ytemp.get(yInd[i+1]);
                if(x1!=x2 && y1!= y2)
                    b+= Math.log(Math.abs((x2-x1)/(y2-y1)));

            }
            System.out.println(a);
            System.out.println(b);
            f = (a-b)/n1;

        }
        return -1*f;
    }

    public double digamma(double x) {
        if (x >= 0 && x < GAMMA_MINX) {
            x = GAMMA_MINX;
        }
        if (x < DIGAMMA_MINNEGX) {
            return digamma(DIGAMMA_MINNEGX + GAMMA_MINX);
        }
        if (x > 0 && x <= S_LIMIT) {
            return -GAMMA - 1 / x;
        }

        if (x >= C_LIMIT) {
            double inv = 1 / (x * x);
            return Math.log(x) - 0.5 / x - inv
                    * ((1.0 / 12) + inv * (1.0 / 120 - inv / 252));
        }
        return digamma(x + 1) - 1 / x;
    }
    private double [] getMaxMin(TetradVector x)
    {
        double max = Double.MIN_VALUE;
        double min = Double.MAX_VALUE;

        for(int i = 0; i < x.size();i++)
        {
            if(x.get(i) > max)
                max = x.get(i);
            if(x.get(i) < min)
                min = x.get(i);
        }
        double [] result = {max,min};
        return result;
    }
    private double[] getMeanStd(TetradVector x)
    {
        double mean = 0;
        for(int i = 0; i < x.size();i++)
        {
            mean = mean + x.get(i);
        }
        mean = mean / x.size();

        double std = 0;
        for(int i = 0; i < x.size();i++)
        {
            std += Math.pow(x.get(i)-mean,2);
        }
        std = std/(x.size()-1);
        double [] result = {mean,Math.sqrt(std)};
        return result;
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
}
