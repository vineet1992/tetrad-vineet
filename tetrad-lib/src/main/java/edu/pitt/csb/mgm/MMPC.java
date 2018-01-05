package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.ChoiceGenerator;

import java.util.ArrayList;
import java.util.List;

/**
 * Implementation of the Max-Min Parents and Children Algorithm by Tsamardinos et al. (2003)
 * Created by vinee_000 on 9/2/2016.
 */
public class MMPC {
    private DataSet data;
    private double alpha;
    private IndependenceTest ind;
    private int depthL;
    public MMPC(DataSet d, double threshold, IndependenceTest i)
    {
        data = d;
        alpha = threshold;
        ind = i;
        depthL = -1;
    }
public void setDepth(int d)
{
    depthL = d;
}
    public ArrayList<Node> getPC(int target)
    {
        return getPC(data.getVariable(target));
    }
    public ArrayList<Node> getPC(Node target)
    {
        ArrayList<Node> currPC = new ArrayList<Node>();
        boolean varFound = true;
        System.out.println("Growing PC");
        while(varFound) //add variable that has maximum association based on its worst subset
        {
            double maxAssoc = Double.MIN_VALUE;
            Node maxN = null;
            for(int i = 0; i < data.getNumColumns();i++)
            {
                if(data.getVariable(i).getName().equals(target.getName()) || currPC.contains(data.getVariable(i)))
                    continue;
                double assoc = getMinAssoc(target,data.getVariable(i),currPC);
                //System.out.println("The min assoc: of " + target + " " + data.getVariable(i) + " is " + assoc);
                if(assoc > maxAssoc)
                {
                    maxAssoc = assoc;
                    maxN = data.getVariable(i);
                }
            }
            if(maxAssoc < 1-alpha)
            {
                System.out.println("No variable found to add: " + maxAssoc);
                varFound = false;
            }
            else
            {
                System.out.println("Adding variable " + maxN + ", with maxAssoc: " + maxAssoc);
                currPC.add(maxN);
                varFound = true;
            }
        }
        System.out.println("Shrinking PC");
        boolean removed = false;
        do {
            removed = false;
           A: for(int i = 0; i < currPC.size();i++)
            {
                System.out.println("Testing node " + i);
                Node curr = currPC.get(i);
                List<Node> temp = (List<Node>)currPC.clone();
                temp.remove(curr);
                for(int j = 0; j < temp.size();j++)
                {
                    if(j==0)
                    {
                        if(ind.isIndependent(curr,target)) {
                            System.out.println("Removing var: " + curr);
                            currPC.remove(curr);
                            removed = true;
                            break A;
                        }

                    }
                    else {
                        ChoiceGenerator cg = new ChoiceGenerator(temp.size(), j);
                        int [] choice = cg.next();
                        while(choice!=null)
                        {
                            List<Node> z = new ArrayList<Node>();
                            for(int k = 0; k < choice.length;k++)
                                z.add(temp.get(choice[k]));
                            if(ind.isIndependent(curr,target,z))
                            {
                                System.out.println("Removing variable: " + curr);
                                currPC.remove(curr);
                                removed = true;
                                break A;
                            }
                            choice = cg.next();
                        }
                    }
                }
            }

        }while(removed);
        System.out.println("No variable found to remove... returning PC: " + currPC);
        return currPC;
        //second step try to remove variables that have min association smaller than threshold based on everyone else in currPC

    }
    private double getMinAssoc(Node target, Node current, ArrayList<Node> currPC)
    {
        if(currPC.isEmpty())
        {
            ind.isIndependent(current,target);
            return 1-ind.getPValue();
        }
        else
        {
            if(depthL==-1)
            {
                double minAssoc = 1;
                for(int i = 0; i < currPC.size();i++)
                {
                    ChoiceGenerator cg = new ChoiceGenerator(currPC.size(),i);
                    int [] choice = cg.next();
                    while(choice!=null)
                    {
                        List<Node> z = new ArrayList<Node>();
                        for(int j = 0; j < choice.length;j++)
                            z.add(currPC.get(choice[j]));

                        ind.isIndependent(current, target, z);
                        if(minAssoc > (1-ind.getPValue()))
                            minAssoc = 1-ind.getPValue();
                        choice = cg.next();
                    }
                }
                return minAssoc;
            }
            else
            {
                double minAssoc = 1;
                for(int i = 0; i < depthL;i++)
                {
                    ChoiceGenerator cg = new ChoiceGenerator(currPC.size(),i);
                    int [] choice = cg.next();
                    while(choice!=null)
                    {
                        List<Node> z = new ArrayList<Node>();
                        for(int j = 0; j < choice.length;j++)
                            z.add(currPC.get(choice[j]));
                        ind.isIndependent(current,target,z);
                        if(minAssoc > (1-ind.getPValue()))
                            minAssoc = 1-ind.getPValue();
                        choice = cg.next();
                    }

                }
                return minAssoc;
            }
        }
    }
}
