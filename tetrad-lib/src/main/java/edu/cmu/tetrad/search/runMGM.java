package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.DataSet;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;

/**
 * Created by vinee_000 on 5/27/2017.
 */
public class runMGM {
    public static void main(String [] args) throws Exception
    {
        DataSet data = MixedUtils.loadDataSet2("out.txt");
        data.removeColumn(data.getVariable("ID"));
        System.out.println(data);
        double [] l = {.35,.35,.1};
        MGM m = new MGM(data,l);
        m.learnEdges(1000);
        System.out.println(m.graphFromMGM());
    }
}
