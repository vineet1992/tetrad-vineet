package edu.pitt.csb.mgm;

import edu.cmu.tetrad.util.TetradVector;

/**
 * Created by vinee_000 on 9/11/2016.
 */
public class IGCI_Test {
    public static void main(String[] args) {

    double[] x = new double[100];
        double [] y = new double[100];
        for(int i = 0; i < 100; i++)
        {
            x[i] = i + 1;
            y[i] = Math.pow(x[i],3) + 5;
        }
        TetradVector xvec  = new TetradVector(x);
        TetradVector yvec = new TetradVector(y);
        IGCI i = new IGCI(xvec,yvec, IGCI.types.UNIFORM, IGCI.types.INTEGRAL);
       System.out.println(i.getDirection());
    }
}
