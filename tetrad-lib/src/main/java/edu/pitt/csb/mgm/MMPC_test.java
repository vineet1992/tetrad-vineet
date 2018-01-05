package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndependenceTest;

import java.util.List;

/**
 * Created by vinee_000 on 9/2/2016.
 */
public class MMPC_test {
    public static void main(String [] args) throws Exception
    {
        DataSet d = MixedUtils.loadDataSet2("C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Fall 2016/Kauffman Project/for_Lucas/for_Lucas/tetrad_500_npn.txt");
        IndependenceTest i = new IndTestMultinomialAJ(d,.1);
        MMPC m = new MMPC(d,.05,i);
        List<Node> l = m.getPC(d.getVariable("MT"));
        System.out.println(l);
    }
}
