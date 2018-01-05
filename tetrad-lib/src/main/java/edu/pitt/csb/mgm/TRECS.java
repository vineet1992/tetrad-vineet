package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;

import java.util.ArrayList;

/**
 * Created by vinee_000 on 9/1/2016.
 */
public class TRECS {

    private DataSet data;
    private boolean featuresSet;
    private double [] selectedFeatures;
    private double fsAlpha;
    private double c1Alpha;
    private double c2Alpha;
    public TRECS(DataSet d, double alpha)
    {
        data = d;
        featuresSet = false;
        fsAlpha = alpha;
        c1Alpha = alpha;
        c2Alpha = alpha;
    }
    public TRECS (DataSet d, double a1, double a2, double a3)
    {
        data = d;
        fsAlpha = a1;
        c1Alpha = a2;
        c2Alpha = a3;

    }
    public void setInitialFeatures(double [] sf)
    {
        selectedFeatures = sf;
        featuresSet = true;
    }

    public ArrayList<double[]> run() {

        if(!featuresSet)
        {
            //Run Feature Selection Method
        }
  //      ReKS r = new ReKS(d);
//ReKSTree reks = r.getTree();
        //Run ReKS

return null;
    }

}
