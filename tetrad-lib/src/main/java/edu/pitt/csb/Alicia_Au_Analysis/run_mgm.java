package edu.pitt.csb.Alicia_Au_Analysis;

import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DelimiterType;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.mgm.STEPS;

import java.io.PrintStream;

/**
 * Created by vinee_000 on 10/24/2017.
 */
public class run_mgm {
    public static void main(String [] args)throws Exception
    {
        DelimiterType d2 = DelimiterType.TAB;
        DataSet d = MixedUtils.loadDataSet2("morbidity_full_max.txt",d2);
        d = MixedUtils.completeCases(d);
     //   d.removeColumn(d.getVariable("Max.MBP"));
        d.removeColumn(d.getVariable("Any.Abnl.EEG.admission..y.1."));
        d.removeColumn(d.getVariable("Any.abnl.CT.during.admission..y.1."));
        d.removeColumn(d.getVariable("Max.CNSE"));
        d.removeColumn(d.getVariable("Max.S100B"));
        System.out.println(d);

        int ns = 20;
        double g = 0.01;
        int numLambdas = 40;
        double [] lambda = new double[numLambdas];
        for(int i = 0; i < numLambdas;i++)
        {
            lambda[i] = .1+(.7/numLambdas)*i;
        }
        IndependenceTest ii = new IndTestMultinomialAJ(d,0.1);
        /*System.out.println(ii.isIndependent(ii.getVariable("Max.CNSE"),ii.getVariable("Poor_Outcome")));
        System.out.println(ii.getPValue());
        System.out.println(ii.isIndependent(ii.getVariable("Max.MBP"),ii.getVariable("Poor_Outcome")));
        System.out.println(ii.getPValue());
        System.out.println(ii.isIndependent(ii.getVariable("Max.S100B"),ii.getVariable("Poor_Outcome")));
        System.out.println(ii.getPValue());

        System.out.println(ii.isIndependent(ii.getVariable("Max.CNSE"),ii.getVariable("Neuromorbidity.Devel.During.Admission.1.y.")));
        System.out.println(ii.getPValue());
        System.out.println(ii.isIndependent(ii.getVariable("Max.MBP"),ii.getVariable("Neuromorbidity.Devel.During.Admission.1.y.")));
        System.out.println(ii.getPValue());
        System.out.println(ii.isIndependent(ii.getVariable("Max.S100B"),ii.getVariable("Neuromorbidity.Devel.During.Admission.1.y.")));
        System.out.println(ii.getPValue());

        System.out.println(ii.isIndependent(ii.getVariable("Max.CNSE"),ii.getVariable("Developed.Neuromorbidity..1.y.")));
        System.out.println(ii.getPValue());
        System.out.println(ii.isIndependent(ii.getVariable("Max.MBP"),ii.getVariable("Developed.Neuromorbidity..1.y.")));
        System.out.println(ii.getPValue());
        System.out.println(ii.isIndependent(ii.getVariable("Max.S100B"),ii.getVariable("Developed.Neuromorbidity..1.y.")));
        System.out.println(ii.getPValue());


        System.exit(0);
*/
        STEPS s = new STEPS(d,lambda,g,d.getNumRows(),true);
        Graph g2 = s.runStepsPar();
        double [][] stab = s.stabilities;


        PrintStream out = new PrintStream("Morbidity_Full_Max_No_CNSE_S100B.txt");
        FciMaxP f = new FciMaxP(ii);
        f.setInitialGraph(g2);
        out.println(f.search());
        //System.exit(0);

        out = new PrintStream("Morbidity_Full_Max_MGM_No_CNSE_S100B.txt");
        out.println(g2);
        out.flush();
        out.close();
        out = new PrintStream("Stabilities_Full_Morbidity_Max_No_CNSE_S100B.txt");
        for(int i = 0; i < d.getNumColumns();i++)
        {
            out.print(d.getVariable(i).getName());
            if(i < d.getNumColumns()-1)
                out.print("\t");
            else
                out.println();
        }
        for(int i = 0; i < d.getNumColumns();i++)
        {
            out.print(d.getVariable(i).getName()+"\t");
            for(int j = 0; j < d.getNumColumns();j++)
            {
                out.print(stab[i][j]);
                if(j < d.getNumColumns()-1)
                    out.print("\t");
                else
                    out.println();
            }
        }
        out.flush();
        out.close();
        out = new PrintStream("Morbidity_Full_max_No_CNSE_S100B_cyto.txt");
        for(Edge e: g2.getEdges())
        {
            int i = d.getColumn(d.getVariable(e.getNode1().getName()));
            int j = d.getColumn(d.getVariable(e.getNode2().getName()));

            out.println(e.getNode1().getName() + "\tundir\t" + e.getNode2().getName() + "\t" + stab[i][j]);
        }
        out.flush();
        out.close();
    }
}
