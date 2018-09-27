package edu.pitt.csb.HBA;

import edu.cmu.tetrad.data.DataSet;
import edu.pitt.csb.mgm.MixedUtils;

/**
 * Created by vinee_000 on 9/20/2018.
 */
public class readInTest {
    public static void main(String [] args) throws Exception
    {
        DataSet temp = MixedUtils.loadDataSet2("DataSet_50_300_301_CGS_M.txt",5);
        System.out.println(temp.isMixed());
        System.out.println(temp);
    }
}
