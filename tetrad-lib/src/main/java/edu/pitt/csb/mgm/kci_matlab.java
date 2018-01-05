package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.TetradMatrix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by vinee_000 on 2/22/2016.
 */
public class kci_matlab implements IndependenceTest {
    double alpha;
    DataSet dat;
    int graphNum;
    double lastP;
    int timesCalled;
    int total;
    int changes;
    Graph tg;
    int correct;
    PrintStream out;

    public kci_matlab(DataSet data, double threshold, int gn, Graph trueGraph,PrintStream o)
    {
        total = 0;
        tg = trueGraph;
        alpha = threshold;
        dat = data;
        graphNum = gn;
        lastP = -1;
        out = o;
        out.println("Conditioning Size\tPredicted Dependent?\tActual\tX\tY");
        out.flush();
    }
    /**
     * Returns an Independence test for a subset of the variables.
     */
    public IndependenceTest indTestSubset(List<Node> vars)
    {
         throw new UnsupportedOperationException();
    }

    /**
     * Returns true if the given independence question is judged true, false if not. The independence question is of the
     * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are variables in the list returned by
     * getVariableNames().
     */
    public double getTotal()
    {
        return total;
    }
    public double getScore(){
        return 0;
    }
    public boolean isIndependent(Node x, Node y, List<Node> z)
    {
        total++;
        IndependenceTest indr = new IndTestMultinomialLogisticRegressionWald(dat,alpha,false);
        if(indr.isDependent(x,y,z)) {
            lastP = indr.getPValue();
            return false;
        }
        else {
            timesCalled++;
            String statement = "matlab -wait -nosplash -nodesktop -minimize -r \"indTestWrap(" + graphNum + ",";
            String first = x.getName();
            String second = y.getName();
            statement = statement + first.substring(1) + ",";
            statement = statement + second.substring(1) + ",[";
            for (int i = 0; i < z.size(); i++) {
                statement = statement + z.get(i).getName().substring(1);
                if (i < z.size() - 1)
                    statement = statement + " ";

            }

            statement = statement + "])\" -logfile output.log";
            try {
                Process proc = Runtime.getRuntime().exec(statement);
                proc.waitFor();
                File f = new File("output.log");
                BufferedReader b = new BufferedReader(new FileReader("output.log"));
                double pval = -1;
                while (b.ready()) {
                    String[] line = b.readLine().split(" ");
                    for (int i = 0; i < line.length; i++) {
                        try {
                            pval = Double.parseDouble(line[i].trim());
                        } catch (Exception e) {
                        }

                    }
                }
                f.delete();
                lastP = pval;
                if (pval >= alpha) {
                    Node one = tg.getNode(x.getName());
                    Node two = tg.getNode(y.getName());
                    LinkedList<Node> ll = new LinkedList<Node>();
                    for(int jj = 0; jj < z.size();jj++)
                        ll.add(tg.getNode(z.get(jj).getName()));
                    if(!tg.isDSeparatedFrom(one,two,ll))
                        out.println(z.size()+"\t"+0+"\t"+1 + "\t"+x.toString()+"\t"+y.toString());
                    else
                        out.println(z.size()+"\t"+0+"\t"+0+"\t"+x.toString()+"\t"+y.toString());
                    return true;
                }
                else {
                    changes++;

                    Node one = tg.getNode(x.getName());
                    Node two = tg.getNode(y.getName());
                    LinkedList<Node> ll = new LinkedList<Node>();
                    for(int jj = 0; jj < z.size();jj++)
                        ll.add(tg.getNode(z.get(jj).getName()));
                    if( !tg.isDSeparatedFrom(one,two,ll)) {
                        out.println(z.size()+"\t"+1+"\t"+1 +"\t" + x.toString() + "\t" + y.toString());
                        correct++;
                    }
                    else
                        out.println(z.size()+"\t"+1+"\t"+0+"\t"+x.toString()+"\t"+y.toString());
                    out.flush();
                    return false;
                }

            } catch (Exception e) {
                e.printStackTrace();
                return true;
            }
        }
    }

    /**
     * Returns true if the given independence question is judged true, false if not. The independence question is of the
     * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are variables in the list returned by
     * getVariableNames().
     */
    public boolean isIndependent(Node x, Node y, Node... z)
    {
        LinkedList<Node> thez = new LinkedList<Node>();
        for(Node s:z)
        thez.add(s);
        return isIndependent(x,y,thez);
    }

    /**
     * Returns true if the given independence question is judged false, true if not. The independence question is of the
     * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are variables in the list returned by
     * getVariableNames().
     */
    public boolean isDependent(Node x, Node y, List<Node> z)
    {
        return !isIndependent(x,y,z);
    }

    /**
     * Returns true if the given independence question is judged false, true if not. The independence question is of the
     * form x _||_ y | z, z = <z1,...,zn>, where x, y, z1,...,zn are variables in the list returned by
     * getVariableNames().
     */
    public boolean isDependent(Node x, Node y, Node... z)
    {
        LinkedList<Node> thez = new LinkedList<Node>();
        for(Node s:z)
            thez.add(s);
        return isDependent(x,y,thez);
    }

    /**
     * Returns the probability associated with the most recently executed independence test, of Double.NaN if p value is
     * not meaningful for tis test.
     */
    public double getPValue()
    {
        return lastP;
    }

    /**
     * Returns the list of variables over which this independence checker is capable of determinining independence
     * relations.
     */
   public  List<Node> getVariables()
    {
        return dat.getVariables();
    }

    /**
     * Returns the variable by the given name.
     */
   public  Node getVariable(String name)
    {
        return dat.getVariable(name);
    }

    /**
     * Returns the list of names for the variables in getNodesInEvidence.
     */
    public List<String> getVariableNames()
    {
        return dat.getVariableNames();
    }

    /**
     * Returns true if y is determined the variable in z.
     */
    public boolean determines(List<Node> z, Node y)
    {

        return false;
    }

    /**
     * Returns the significance level of the independence test.
     *
     * @throws UnsupportedOperationException if there is no significance level.
     */
    public double getAlpha()
    {
        return alpha;
    }

    /**
     * Sets the significance level.
     */
    public void setAlpha(double alpha2)
    {
        alpha = alpha2;
    }

    /**
     * '
     *
     * @return The data model for the independence test.
     */
    public DataModel getData()
    {
        return dat;
    }


    public ICovarianceMatrix getCov()
    {
        throw new UnsupportedOperationException();
    }

    public List<DataSet> getDataSets()
    {
        LinkedList<DataSet> L = new LinkedList<DataSet>();
        L.add(dat);
        return L;
    }

    public int getSampleSize()
    {
        return dat.getNumRows();
    }

    public List<TetradMatrix> getCovMatrices()
    {
        throw new UnsupportedOperationException();
    }

}
