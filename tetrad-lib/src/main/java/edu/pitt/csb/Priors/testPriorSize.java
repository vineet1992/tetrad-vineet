package edu.pitt.csb.Priors;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.File;
import java.util.ArrayList;

/**
 * Created by vinee_000 on 1/11/2018.
 */
public class testPriorSize {
    public static void main(String [] args) throws Exception
    {
        boolean tumors = true;
        DataSet data = null;
        if(tumors)
            data = MixedUtils.loadDataSet2("genes_with_clinical.txt");
        else
            data = MixedUtils.loadDataSet2("genes_with_clinical_normals.txt");
        //System.out.println(data);
        data = MixedUtils.completeCases(data);
        File f = new File("prior_sources");
        File [] stuff = f.listFiles();
        ArrayList<File> files = new ArrayList<File>();
        for(int i = 0; i < stuff.length;i++)
        {
            if(!stuff[i].getName().contains("PAM50"))
                files.add(stuff[i]);
        }
        TetradMatrix[] priors = new TetradMatrix[files.size()];
        for(int i = 0; i < files.size();i++) {
            priors[i] = new TetradMatrix(realDataPriorTest.loadPrior(files.get(i),data.getNumColumns()));
            int count = 0;
            TetradMatrix curr = priors[i];
            for(int j =0; j < curr.rows();j++)
            {
                for(int k = j+1; k < curr.columns();k++)
                {
                    if(curr.get(j,k)!=0)
                    {
                        count++;
                    }
                }
            }
            System.out.println(files.get(i) + "\t" + count);
        }
    }
}
