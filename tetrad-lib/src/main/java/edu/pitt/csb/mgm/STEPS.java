package edu.pitt.csb.mgm;


import cern.colt.matrix.DoubleMatrix2D;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DelimiterType;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.CpcStable;
import edu.cmu.tetrad.search.FciMaxP;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.MathUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.stability.SearchWrappers;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;

import static edu.pitt.csb.stability.StabilityUtils.StabilitySearchPar;
import static edu.pitt.csb.stability.StabilityUtils.getSubSize;

/**
 * Created by vinee_000 on 10/5/2016.
 */

//TODO Incorporate iterLimit into all parts of StEPS searches for efficiency
public class STEPS {
    private DataSet d;
    private int N;
    public int b;
    private double [] lambda;
    private double gamma;
    private boolean includeZeros = true;
    private int iterLimit = 500;
    public double origLambda;
    public Graph pdGraph;
    public Graph lastGraph;
    public double[] lastLambda;
    private boolean leaveOneOut = false;
    public double [][] stabilities = null;
    private int [][] subs;
    private boolean computeStabs;
    private boolean runFull = false; //Don't monotonize instability and kill method as soon as stability threshold is met
    private boolean verbose = false;

    public STEPS(DataSet dat,double [] lam,double g,int numSub)
    {
        N = numSub;
        gamma = g;
        lambda = lam;
        d = dat;
        this.b = getSubSize(dat.getNumRows());



    }
    public STEPS(DataSet dat,double [] lam,double g,int numSub, int b)
    {
        N = numSub;
        gamma = g;
        lambda = lam;
        d = dat;
        this.b = b;
    }
    public STEPS(DataSet dat,double [] lam,double g,int [][] subsamples)
    {
        gamma = g;
        lambda = lam;
        d = dat;
        this.b = getSubSize(dat.getNumRows());

        this.subs = subsamples;


    }
    public STEPS(DataSet dat,double [] lam,double g,int numSub,boolean loo)
    {
        leaveOneOut = loo;
        N = numSub;
        gamma = g;
        lambda = lam;
        d = dat;
        this.b = getSubSize(dat.getNumRows());



    }

    public void runFull(){runFull = true;}
    public void setIterLimit(int iterLimit)
    {
        this.iterLimit = iterLimit;
    }
    public void setComputeStabs(boolean stabs){computeStabs = stabs;}
    public void setVerbose(){verbose = true;}
    public double [][] runStepsArrayPar()
    {
        double [][] result = new double[lambda.length][4];
        Arrays.sort(lambda);

        int currIndex = 0;
        double CC = -1;
        double CD = -1;
        double DD = -1;
        double CCMax = 100;
        double CCMaxI = -1;
        double CDMax = 100;
        double CDMaxI = -1;
        double DDMax = 100;
        double DDMaxI = -1;
        double oneLamb = -1;
        double allMax = 100;
        double allMaxI = -1;
        int p = 0;
        int q = 0;
        for(int i =0; i < d.getNumColumns();i++)
        {
            Node n = d.getVariable(i);
            if(n instanceof DiscreteVariable)
                q++;
            else
                p++;
        }
        // System.out.println("P:" + p);
        //   System.out.println("Q:" + q);
        //  System.out.println("b:" + b);
        A:while(true) //go until we break by having instability better than threshold
        {
            //System.out.println("Lambda: " + lambda[currIndex]);
            double [] lambdaCurr = {lambda[currIndex],lambda[currIndex],lambda[currIndex]};
            DoubleMatrix2D adjMat;
            if(subs!=null)
                adjMat = StabilityUtils.StabilitySearchPar(d,new SearchWrappers.MGMWrapper(lambdaCurr),subs);
            else if(leaveOneOut)
                adjMat = StabilityUtils.StabilitySearchPar(d,new SearchWrappers.MGMWrapper(lambdaCurr));
            else
                adjMat = StabilityUtils.StabilitySearchPar(d,new SearchWrappers.MGMWrapper(lambdaCurr),N,b);

            double ccDestable = 0;
            double cdDestable = 0;
            double ddDestable = 0;
            ////////////////// TODO Decide if this is too harsh
            double numCC = 0;
            double numCD = 0;
            double numDD = 0;
            //////////////////
            for(int j = 0; j < d.getNumColumns();j++)
            {
                for(int k = j+1; k < d.getNumColumns();k++)
                {
                    Node one = d.getVariable(j);
                    Node two = d.getVariable(k);
                    if(one instanceof DiscreteVariable && two instanceof DiscreteVariable)
                    {
                        //          System.out.println("DD: " + adjMat[j][k]);
                        //       if(adjMat[j][k]!=0)
                        numDD++;
                        ddDestable += 2*adjMat.get(j,k)*(1-adjMat.get(j,k));
                    }
                    else if(one instanceof DiscreteVariable || two instanceof DiscreteVariable)
                    {
                        //     if(adjMat[j][k]!=0)
                        numCD++;
                        //         System.out.println("CD: " + adjMat[j][k]);
                        cdDestable += 2*adjMat.get(j,k)*(1-adjMat.get(j,k));

                    }
                    else
                    {
                        //   if(adjMat[j][k]!=0)
                        numCC++;
                        //       System.out.println("CC:" + adjMat[j][k]);
                        ccDestable+=2*adjMat.get(j,k)*(1-adjMat.get(j,k));
                    }
                }
            }
            //System.out.println(adjMat);
            double allDestable = ccDestable + cdDestable + ddDestable;
            allDestable = allDestable/(numCC+numCD+numDD);
            if(allDestable <= gamma && oneLamb==-1)
                oneLamb = lambda[currIndex];
            if(allDestable <= allMax)
            {
                allMax = allDestable;
                allMaxI = lambda[currIndex];
            }
            /*
            ccDestable = ccDestable/ MathUtils.choose(p,2);
            cdDestable = cdDestable/(p*q);
            ddDestable = ddDestable/(MathUtils.choose(q,2));
            */
            ccDestable = ccDestable / numCC;
            cdDestable = cdDestable / numCD;
            ddDestable = ddDestable/numDD;
            result[currIndex][0] = ccDestable;
            result[currIndex][1] = cdDestable;
            result[currIndex][2] = ddDestable;
            result[currIndex][3] = allDestable;
         //    System.out.println("CC:" + ccDestable);
           //  System.out.println("CD:" + cdDestable);
             //System.out.println("DD:" + ddDestable);
            if(ccDestable <= gamma && CC == -1)
                CC = lambda[currIndex];
            if(cdDestable <= gamma && CD==-1)
                CD = lambda[currIndex];
            if(ddDestable <= gamma && DD == -1)
                DD = lambda[currIndex];
            if(ccDestable <= CCMax)
            {
                CCMax = ccDestable;
                CCMaxI = lambda[currIndex];
            }
            if(cdDestable <= CDMax)
            {
                CDMax = cdDestable;
                CDMaxI = lambda[currIndex];
            }
            if(ddDestable <= DDMax)
            {
                DDMax = ddDestable;
                DDMaxI = lambda[currIndex];
            }
            if(!runFull && (CC!=-1 && CD != -1 && DD !=-1 && oneLamb!=-1))
                break A;
            if(currIndex==lambda.length-1)
                break A;
            currIndex++;
        }
        if(CC==-1)
            CC = CCMaxI;
        if(CD==-1)
            CD = CDMaxI;
        if(DD==-1)
            DD = DDMaxI;

        if(runFull)
        {
            int [] inds = getBestVals(result,gamma);
            CC = lambda[inds[0]];
            CD = lambda[inds[1]];
            DD = lambda[inds[2]];
            allMaxI = lambda[inds[3]];
        }



        double [] lambda = {CC,CD,DD};
        System.out.println("Lambdas: " + Arrays.toString(lambda));
        if(oneLamb==-1 || runFull)
            origLambda = allMaxI;
        else
            origLambda = oneLamb;


        //  System.out.println("Final Params: (" + CC + "," + CD + "," + DD + ")");
        MGM m = new MGM(d,lambda);
        m.learnEdges(iterLimit);
        lastGraph = m.graphFromMGM();
       // DoubleMatrix2D stabs = StabilitySearchPar(d,new SearchWrappers.MGMWrapper(lambda));
       // this.stabilities = stabs.toArray();
        lastLambda = lambda;
        return result;
    }

    private int [] getBestVals(double [][] result, double gamma)
    {
        int [] lambda = new int[4];
        A:for(int i = 0; i < 4; i ++)
        {
            int maxIndex = -1;
            double max = -1;
            for(int j = 0; j < result.length;j++)
            {
                if(result[j][i]>max)
                {
                    maxIndex = j;
                    max = result[j][i];
                }
            }
            for(int j = maxIndex;j<result.length;j++)
            {
                if(result[j][i] < gamma)
                {
                    lambda[i] = j;
                    continue A;
                }
            }
            if(lambda[i]==0)
                lambda[i] = result.length-1;
        }
        return lambda;
    }
    //returns a matrix where rows correspond to lambda values, and columns correspond to: 0 = CC, 1 = CD, 2 = DD, 3 = ALL


    public double [][] runStepsArray()
    {
        double [][] result = new double[lambda.length][4];
        Arrays.sort(lambda);
        int currIndex = 0;
        double CC = -1;
        double CD = -1;
        double DD = -1;
        double CCMax = 100;
        double CCMaxI = -1;
        double CDMax = 100;
        double CDMaxI = -1;
        double DDMax = 100;
        double DDMaxI = -1;
        double oneLamb = -1;
        double allMax = 100;
        double allMaxI = -1;
        int p = 0;
        int q = 0;
        for(int i =0; i < d.getNumColumns();i++)
        {
            Node n = d.getVariable(i);
            if(n instanceof DiscreteVariable)
                q++;
            else
                p++;

        }
        // System.out.println("P:" + p);
        //   System.out.println("Q:" + q);
        //  System.out.println("b:" + b);
        A:while(true) //go until we break by having instability better than threshold
        {
               System.out.println("Lambda: " + lambda[currIndex]);
            double [][] adjMat = new double[d.getNumColumns()][d.getNumColumns()];
            for(int i = 0; i < N; i ++) //Learn MGM on N networks, update adjacency matrix after each run
            {
                Graph g = null;
                boolean again = true;
                while(again) {
                    //      System.out.println("here");
                    try {
                        DataSet temp = d.copy();
                        temp.permuteRows();
                        int[] removal = new int[b];
                        for (int j = 0; j < removal.length; j++)
                            removal[j] = j;
                        temp = temp.subsetRows(removal);
                        double[] l = {lambda[currIndex], lambda[currIndex], lambda[currIndex]};
                        MGM m = new MGM(temp, l);
                        m.learnEdges(iterLimit);
                        g = m.graphFromMGM();
                        again = false;
                    }
                    catch(Exception e)
                    {
                        e.printStackTrace();
                    }

                }


                for(int j = 0; j < d.getNumColumns();j++)
                {
                    for(int k = j+1; k < d.getNumColumns();k++)
                    {
                        if(g.isAdjacentTo(g.getNode(d.getVariable(j).getName()),g.getNode(d.getVariable(k).getName())))
                        {
                            adjMat[j][k] += 1/(double)N;
                        }

                    }

                }

            }
            double ccDestable = 0;
            double cdDestable = 0;
            double ddDestable = 0;
            double numCC = 0;
            double numCD = 0;
            double numDD = 0;
            for(int j = 0; j < d.getNumColumns();j++)
            {
                for(int k = j+1; k < d.getNumColumns();k++)
                {
                    Node one = d.getVariable(j);
                    Node two = d.getVariable(k);
                    if(one instanceof DiscreteVariable && two instanceof DiscreteVariable)
                    {
                        if(includeZeros) {
                            // if(adjMat[j][k]!=0)
                            numDD++;
                            //         System.out.println("CD: " + adjMat[j][k]);
                            ddDestable += 2 * adjMat[j][k] * (1 - adjMat[j][k]);
                        }
                        else
                        {
                            if(adjMat[j][k]!=0) {
                                numDD++;
                                ddDestable += 2 * adjMat[j][k] * (1 - adjMat[j][k]);
                            }
                        }
                    }
                    else if(one instanceof DiscreteVariable || two instanceof DiscreteVariable)
                    {
                        if(includeZeros) {
                            // if(adjMat[j][k]!=0)
                            numCD++;
                            //         System.out.println("CD: " + adjMat[j][k]);
                            cdDestable += 2 * adjMat[j][k] * (1 - adjMat[j][k]);
                        }
                        else
                        {
                            if(adjMat[j][k]!=0) {
                                numCD++;
                                cdDestable += 2 * adjMat[j][k] * (1 - adjMat[j][k]);
                            }
                        }

                    }
                    else
                    {
                        if(includeZeros) {
                            // if(adjMat[j][k]!=0)
                            numCC++;
                            //         System.out.println("CD: " + adjMat[j][k]);
                            ccDestable += 2 * adjMat[j][k] * (1 - adjMat[j][k]);
                        }
                        else
                        {
                            if(adjMat[j][k]!=0) {
                                numCC++;
                                ccDestable += 2 * adjMat[j][k] * (1 - adjMat[j][k]);
                            }
                        }
                    }
                }
            }
            double allDestable = ccDestable + cdDestable + ddDestable;
            allDestable = allDestable/(numCC+numCD+numDD);
            if(allDestable <= gamma && oneLamb==-1)
                oneLamb = lambda[currIndex];
            if(allDestable <= allMax)
            {
                allMax = allDestable;
                allMaxI = lambda[currIndex];
            }
            /*
            ccDestable = ccDestable/ MathUtils.choose(p,2);
            cdDestable = cdDestable/(p*q);
            ddDestable = ddDestable/(MathUtils.choose(q,2));
            */
            ccDestable = ccDestable / numCC;
            cdDestable = cdDestable / numCD;
            ddDestable = ddDestable/numDD;
            // System.out.println("CC:" + ccDestable);
            // System.out.println("CD:" + cdDestable);
            // System.out.println("DD:" + ddDestable);
            result[currIndex][0] = ccDestable;
            result[currIndex][1] = cdDestable;
            result[currIndex][2] = ddDestable;
            result[currIndex][3] = allDestable;


            if(ccDestable <= gamma && CC == -1 && includeZeros)
                CC = lambda[currIndex];
            if(cdDestable <= gamma && CD==-1 && includeZeros)
                CD = lambda[currIndex];
            if(ddDestable <= gamma && DD == -1 && includeZeros)
                DD = lambda[currIndex];
            if(ccDestable <= CCMax)
            {
                CCMax = ccDestable;
                CCMaxI = lambda[currIndex];
            }
            if(cdDestable <= CDMax)
            {
                CDMax = cdDestable;
                CDMaxI = lambda[currIndex];
            }
            if(ddDestable <= DDMax)
            {
                DDMax = ddDestable;
                DDMaxI = lambda[currIndex];
            }
           // if(CC!=-1 && CD != -1 && DD !=-1 && oneLamb!=-1)
             //   break A;
            if(currIndex==lambda.length-1)
                break A;
            currIndex++;
        }
        if(CC==-1)
            CC = CCMaxI;
        if(CD==-1)
            CD = CDMaxI;
        if(DD==-1)
            DD = DDMaxI;

        double [] lambda = {CC,CD,DD};
        if(oneLamb==-1)
            origLambda = allMaxI;
        else
            origLambda = oneLamb;
        //  System.out.println("Final Params: (" + CC + "," + CD + "," + DD + ")");
       MGM m = new MGM(d,lambda);
        m.learnEdges(iterLimit);
        lastGraph = m.graphFromMGM();
        lastLambda = lambda;
        return result;
    }

    private Graph runStepsLOO()
    {
        N = d.getNumRows();
        b = d.getNumRows()-1;
        Arrays.sort(lambda);
        int currIndex = 0;
        double CC = -1;
        double CD = -1;
        double DD = -1;
        double CCMax = 100;
        double CCMaxI = -1;
        double CDMax = 100;
        double CDMaxI = -1;
        double DDMax = 100;
        double DDMaxI = -1;
        double oneLamb = -1;
        double allMax = 100;
        double allMaxI = -1;
        int p = 0;
        int q = 0;
        for(int i =0; i < d.getNumColumns();i++)
        {
            Node n = d.getVariable(i);
            if(n instanceof DiscreteVariable)
                q++;
            else
                p++;
        }
        // System.out.println("P:" + p);
        //   System.out.println("Q:" + q);
        //  System.out.println("b:" + b);
        A:while(true) //go until we break by having instability better than threshold
        {
              // System.out.println("Lambda: " + lambda[currIndex]);
            double [][] adjMat = new double[d.getNumColumns()][d.getNumColumns()];
            B:for(int i = 0; i < d.getNumRows(); i ++) //Learn MGM on N networks, update adjacency matrix after each run
            {
                //System.out.println("Subsample: " + i);
                Graph g = null;
                    //      System.out.println("here");
                    try {

                        DataSet temp = d.copy();
                        int[] removal = new int[d.getNumRows()-1];
                        int count = 0;
                        for(int j = 0; j < d.getNumRows();j++)
                        {
                            if(j!=i) {
                                removal[count] = j;
                                count++;
                            }
                        }
                        temp = temp.subsetRows(removal);
                        double[] l = {lambda[currIndex], lambda[currIndex], lambda[currIndex]};
                        MGM m = new MGM(temp, l);
                        m.learnEdges(iterLimit);
                        g = m.graphFromMGM();
                    }
                    catch(Exception e)
                    {
                        e.printStackTrace();
                        continue B;
                    }


                for(int j = 0; j < d.getNumColumns();j++)
                {
                    for(int k = j+1; k < d.getNumColumns();k++)
                    {
                        if(g.isAdjacentTo(g.getNode(d.getVariable(j).getName()),g.getNode(d.getVariable(k).getName())))
                        {
                            adjMat[j][k] += 1/(double)N;
                        }

                    }

                }

            }
            double ccDestable = 0;
            double cdDestable = 0;
            double ddDestable = 0;
            ////////////////// TODO Decide if this is too harsh
            double numCC = 0;
            double numCD = 0;
            double numDD = 0;
            //////////////////
            for(int j = 0; j < d.getNumColumns();j++)
            {
                for(int k = j+1; k < d.getNumColumns();k++)
                {
                    Node one = d.getVariable(j);
                    Node two = d.getVariable(k);
                    if(one instanceof DiscreteVariable && two instanceof DiscreteVariable)
                    {
                        //          System.out.println("DD: " + adjMat[j][k]);
                        //       if(adjMat[j][k]!=0)
                        numDD++;
                        ddDestable += 2*adjMat[j][k]*(1-adjMat[j][k]);
                    }
                    else if(one instanceof DiscreteVariable || two instanceof DiscreteVariable)
                    {
                        //     if(adjMat[j][k]!=0)
                        numCD++;
                        //         System.out.println("CD: " + adjMat[j][k]);
                        cdDestable += 2*adjMat[j][k]*(1-adjMat[j][k]);

                    }
                    else
                    {
                        //   if(adjMat[j][k]!=0)
                        numCC++;
                        //       System.out.println("CC:" + adjMat[j][k]);
                        ccDestable+=2*adjMat[j][k]*(1-adjMat[j][k]);
                    }
                }
            }
            double allDestable = ccDestable + cdDestable + ddDestable;
            allDestable = allDestable/(numCC+numCD+numDD);
            if(allDestable <= gamma && oneLamb==-1)
                oneLamb = lambda[currIndex];
            if(allDestable <= allMax)
            {
                allMax = allDestable;
                allMaxI = lambda[currIndex];
            }
            /*
            ccDestable = ccDestable/ MathUtils.choose(p,2);
            cdDestable = cdDestable/(p*q);
            ddDestable = ddDestable/(MathUtils.choose(q,2));
            */
            ccDestable = ccDestable / numCC;
            cdDestable = cdDestable / numCD;
            ddDestable = ddDestable/numDD;
            // System.out.println("CC:" + ccDestable);
            // System.out.println("CD:" + cdDestable);
            // System.out.println("DD:" + ddDestable);
            if(ccDestable <= gamma && CC == -1)
                CC = lambda[currIndex];
            if(cdDestable <= gamma && CD==-1)
                CD = lambda[currIndex];
            if(ddDestable <= gamma && DD == -1)
                DD = lambda[currIndex];
            if(ccDestable <= CCMax)
            {
                CCMax = ccDestable;
                CCMaxI = lambda[currIndex];
            }
            if(cdDestable <= CDMax)
            {
                CDMax = cdDestable;
                CDMaxI = lambda[currIndex];
            }
            if(ddDestable <= DDMax)
            {
                DDMax = ddDestable;
                DDMaxI = lambda[currIndex];
            }
            if(CC!=-1 && CD != -1 && DD !=-1 && oneLamb!=-1)
                break A;
            if(currIndex==lambda.length-1)
                break A;
            currIndex++;
        }
        if(CC==-1)
            CC = CCMaxI;
        if(CD==-1)
            CD = CDMaxI;
        if(DD==-1)
            DD = DDMaxI;

        double [] lambda = {CC,CD,DD};
        System.out.println("Lambdas: " + Arrays.toString(lambda));
        if(oneLamb==-1)
            origLambda = allMaxI;
        else
            origLambda = oneLamb;
        //  System.out.println("Final Params: (" + CC + "," + CD + "," + DD + ")");
        double [][] stabilities = new double[d.getNumColumns()][d.getNumColumns()];
        int numSamples = d.getNumRows();
        for(int i = 0; i < d.getNumRows();i++)
        {
            try {
                DataSet temp = d.copy();
                int[] removal = new int[d.getNumRows()];
                int count = 0;
                for(int j = 0; j < d.getNumRows();j++)
                {
                    if(i!=j) {
                        removal[count] = j;
                        count++;
                    }

                }
                temp = temp.subsetRows(removal);
                MGM m = new MGM(temp, lambda);
                m.learnEdges(iterLimit);
                Graph g = m.graphFromMGM();
                for (int j = 0; j < d.getNumColumns(); j++){
                    for(int k = 0; k < d.getNumColumns();k++)
                    {
                        if(g.isAdjacentTo(g.getNode(d.getVariable(j).getName()),g.getNode(d.getVariable(k).getName())))
                        {
                            stabilities[j][k] += 1;
                        }

                    }
                }

            }
            catch(Exception e)
            {
                numSamples--;
            }
        }
        for(int i = 0; i < d.getNumColumns();i++)
        {
            for(int j = 0; j < d.getNumColumns();j++)
            {
                //stabilities[i][j] = stabilities[i][j]/numSamples;
            }
        }
        this.stabilities = stabilities;
        MGM m = new MGM(d,lambda);
        m.learnEdges(iterLimit);
        return m.graphFromMGM();
    }
    public Graph runSteps()
    {
        if(leaveOneOut)
            return runStepsLOO();
        Arrays.sort(lambda);

        int currIndex = 0;
        double CC = -1;
        double CD = -1;
        double DD = -1;
        double CCMax = 100;
        double CCMaxI = -1;
        double CDMax = 100;
        double CDMaxI = -1;
        double DDMax = 100;
        double DDMaxI = -1;
        double oneLamb = -1;
        double allMax = 100;
        double allMaxI = -1;
        int p = 0;
        int q = 0;
        for(int i =0; i < d.getNumColumns();i++)
        {
            Node n = d.getVariable(i);
            if(n instanceof DiscreteVariable)
                q++;
            else
                p++;
        }
       // System.out.println("P:" + p);
     //   System.out.println("Q:" + q);
      //  System.out.println("b:" + b);
        A:while(true) //go until we break by having instability better than threshold
        {
         //   System.out.println("Lambda: " + lambda[currIndex]);
            double [][] adjMat = new double[d.getNumColumns()][d.getNumColumns()];
            for(int i = 0; i < N; i ++) //Learn MGM on N networks, update adjacency matrix after each run
            {
                Graph g = null;
                boolean again = true;
                while(again) {
              //      System.out.println("here");
                    try {

                        DataSet temp = d.copy();
                        temp.permuteRows();
                        int[] removal = new int[b];
                        for (int j = 0; j < removal.length; j++)
                            removal[j] = j;
                        temp = temp.subsetRows(removal);
                        double[] l = {lambda[currIndex], lambda[currIndex], lambda[currIndex]};
                        MGM m = new MGM(temp, l);
                        m.learnEdges(iterLimit);
                        g = m.graphFromMGM();
                        again = false;
                    }
                    catch(Exception e)
                    {
                        e.printStackTrace();
                    }

                }


                for(int j = 0; j < d.getNumColumns();j++)
                {
                    for(int k = j+1; k < d.getNumColumns();k++)
                    {
                        if(g.isAdjacentTo(g.getNode(d.getVariable(j).getName()),g.getNode(d.getVariable(k).getName())))
                        {
                            adjMat[j][k] += 1/(double)N;
                        }

                    }

                }

            }
            double ccDestable = 0;
            double cdDestable = 0;
            double ddDestable = 0;
            ////////////////// TODO Decide if this is too harsh
            double numCC = 0;
            double numCD = 0;
            double numDD = 0;
            //////////////////
            for(int j = 0; j < d.getNumColumns();j++)
            {
                for(int k = j+1; k < d.getNumColumns();k++)
                {
                    Node one = d.getVariable(j);
                    Node two = d.getVariable(k);
                    if(one instanceof DiscreteVariable && two instanceof DiscreteVariable)
                    {
              //          System.out.println("DD: " + adjMat[j][k]);
                 //       if(adjMat[j][k]!=0)
                            numDD++;
                        ddDestable += 2*adjMat[j][k]*(1-adjMat[j][k]);
                    }
                    else if(one instanceof DiscreteVariable || two instanceof DiscreteVariable)
                    {
                   //     if(adjMat[j][k]!=0)
                            numCD++;
               //         System.out.println("CD: " + adjMat[j][k]);
                        cdDestable += 2*adjMat[j][k]*(1-adjMat[j][k]);

                    }
                    else
                    {
                     //   if(adjMat[j][k]!=0)
                            numCC++;
                 //       System.out.println("CC:" + adjMat[j][k]);
                        ccDestable+=2*adjMat[j][k]*(1-adjMat[j][k]);
                    }
                }
            }
            double allDestable = ccDestable + cdDestable + ddDestable;
            allDestable = allDestable/(numCC+numCD+numDD);
            if(allDestable <= gamma && oneLamb==-1)
                oneLamb = lambda[currIndex];
            if(allDestable <= allMax)
            {
                allMax = allDestable;
                allMaxI = lambda[currIndex];
            }
            /*
            ccDestable = ccDestable/ MathUtils.choose(p,2);
            cdDestable = cdDestable/(p*q);
            ddDestable = ddDestable/(MathUtils.choose(q,2));
            */
            ccDestable = ccDestable / numCC;
            cdDestable = cdDestable / numCD;
            ddDestable = ddDestable/numDD;
           // System.out.println("CC:" + ccDestable);
           // System.out.println("CD:" + cdDestable);
           // System.out.println("DD:" + ddDestable);
            if(ccDestable <= gamma && CC == -1)
                CC = lambda[currIndex];
            if(cdDestable <= gamma && CD==-1)
                CD = lambda[currIndex];
            if(ddDestable <= gamma && DD == -1)
                DD = lambda[currIndex];
            if(ccDestable <= CCMax)
            {
                CCMax = ccDestable;
                CCMaxI = lambda[currIndex];
            }
            if(cdDestable <= CDMax)
            {
                CDMax = cdDestable;
                CDMaxI = lambda[currIndex];
            }
            if(ddDestable <= DDMax)
            {
                DDMax = ddDestable;
                DDMaxI = lambda[currIndex];
            }
if(CC!=-1 && CD != -1 && DD !=-1 && oneLamb!=-1)
    break A;
if(currIndex==lambda.length-1)
    break A;
            currIndex++;
        }
        if(CC==-1)
            CC = CCMaxI;
        if(CD==-1)
            CD = CDMaxI;
        if(DD==-1)
            DD = DDMaxI;

double [] lambda = {CC,CD,DD};
System.out.println("Lambdas: " + Arrays.toString(lambda));
        if(oneLamb==-1)
            origLambda = allMaxI;
        else
            origLambda = oneLamb;


      //  System.out.println("Final Params: (" + CC + "," + CD + "," + DD + ")");
        MGM m = new MGM(d,lambda);
        m.learnEdges(iterLimit);
        lastLambda = lambda;
        return m.graphFromMGM();
    }


    public Graph runStepsPar()
    {
        Arrays.sort(lambda);

        int currIndex = 0;
        double CC = -1;
        double CD = -1;
        double DD = -1;
        double CCMax = 100;
        double CCMaxI = -1;
        double CDMax = 100;
        double CDMaxI = -1;
        double DDMax = 100;
        double DDMaxI = -1;
        double oneLamb = -1;
        double allMax = 100;
        double allMaxI = -1;
        int p = 0;
        int q = 0;
        for(int i =0; i < d.getNumColumns();i++)
        {
            Node n = d.getVariable(i);
            if(n instanceof DiscreteVariable)
                q++;
            else
                p++;
        }
        // System.out.println("P:" + p);
        //   System.out.println("Q:" + q);
        //  System.out.println("b:" + b);
        long time = System.nanoTime();
        A:while(true) //go until we break by having instability better than threshold
        {

                double [] lambdaCurr = {lambda[currIndex],lambda[currIndex],lambda[currIndex]};

            if(verbose)
            {
                System.out.println(currIndex*100/(double)lambda.length);
            }
                DoubleMatrix2D adjMat;
                if(subs!=null)
                    adjMat = StabilityUtils.StabilitySearchPar(d,new SearchWrappers.MGMWrapper(lambdaCurr),subs);
                else if(leaveOneOut)
                    adjMat = StabilityUtils.StabilitySearchPar(d,new SearchWrappers.MGMWrapper(lambdaCurr));
                else
                    adjMat = StabilityUtils.StabilitySearchPar(d,new SearchWrappers.MGMWrapper(lambdaCurr),N,b);

            double ccDestable = 0;
            double cdDestable = 0;
            double ddDestable = 0;
            ////////////////// TODO Decide if this is too harsh
            double numCC = 0;
            double numCD = 0;
            double numDD = 0;
            //////////////////

            //We assume here that the subsamples have the variables in the same order
            for(int j = 0; j < d.getNumColumns();j++)
            {
                for(int k = j+1; k < d.getNumColumns();k++)
                {
                    Node one;
                    Node two;
                    if(subs!=null) {
                        one = d.getVariable(j);
                        two = d.getVariable(k);
                    }
                    else
                    {
                        one = d.getVariable(j);
                        two = d.getVariable(k);
                    }
                    if(one instanceof DiscreteVariable && two instanceof DiscreteVariable)
                    {
                        //          System.out.println("DD: " + adjMat[j][k]);
                        //       if(adjMat[j][k]!=0)
                        numDD++;
                        ddDestable += 2*adjMat.get(j,k)*(1-adjMat.get(j,k));
                    }
                    else if(one instanceof DiscreteVariable || two instanceof DiscreteVariable)
                    {
                        //     if(adjMat[j][k]!=0)
                        numCD++;
                        //         System.out.println("CD: " + adjMat[j][k]);
                        cdDestable += 2*adjMat.get(j,k)*(1-adjMat.get(j,k));

                    }
                    else
                    {
                        //   if(adjMat[j][k]!=0)
                        numCC++;
                        //       System.out.println("CC:" + adjMat[j][k]);
                        ccDestable+=2*adjMat.get(j,k)*(1-adjMat.get(j,k));
                    }
                }
            }
            double allDestable = ccDestable + cdDestable + ddDestable;
            allDestable = allDestable/(numCC+numCD+numDD);
            if(allDestable <= gamma && oneLamb==-1)
                oneLamb = lambda[currIndex];
            if(allDestable <= allMax)
            {
                allMax = allDestable;
                allMaxI = lambda[currIndex];
            }
            /*
            ccDestable = ccDestable/ MathUtils.choose(p,2);
            cdDestable = cdDestable/(p*q);
            ddDestable = ddDestable/(MathUtils.choose(q,2));
            */
            ccDestable = ccDestable / numCC;
            cdDestable = cdDestable / numCD;
            ddDestable = ddDestable/numDD;
            if(ccDestable <= gamma && CC == -1)
                CC = lambda[currIndex];
            if(cdDestable <= gamma && CD==-1)
                CD = lambda[currIndex];
            if(ddDestable <= gamma && DD == -1)
                DD = lambda[currIndex];
            if(ccDestable <= CCMax)
            {
                CCMax = ccDestable;
                CCMaxI = lambda[currIndex];
            }
            if(cdDestable <= CDMax)
            {
                CDMax = cdDestable;
                CDMaxI = lambda[currIndex];
            }
            if(ddDestable <= DDMax)
            {
                DDMax = ddDestable;
                DDMaxI = lambda[currIndex];
            }
            if(CC!=-1 && CD != -1 && DD !=-1 && oneLamb!=-1)
                break A;
            if(currIndex==lambda.length-1)
                break A;
            currIndex++;
        }
        time = System.nanoTime()-time;
        double timePerIter = time / ((double)currIndex*N);
        if(CC==-1)
            CC = CCMaxI;
        if(CD==-1)
            CD = CDMaxI;
        if(DD==-1)
            DD = DDMaxI;


        if(verbose)
        {
            System.out.println(100-(timePerIter*100)/(time+timePerIter));
        }
        double [] lambda = {CC,CD,DD};
        System.out.println("Lambdas: " + Arrays.toString(lambda));
        if(oneLamb==-1)
            origLambda = allMaxI;
        else
            origLambda = oneLamb;


          //System.out.println("Final Params: (" + CC + "," + CD + "," + DD + ")");
        MGM m = new MGM(d,lambda);
        m.learnEdges(iterLimit);
        DoubleMatrix2D stabs;
        if(computeStabs) {
            if (subs != null)
                stabs = StabilityUtils.StabilitySearchPar(d, new SearchWrappers.MGMWrapper(lambda), subs);
            else if (leaveOneOut)
                stabs = StabilityUtils.StabilitySearchPar(d, new SearchWrappers.MGMWrapper(lambda));
            else
                stabs = StabilityUtils.StabilitySearchPar(d, new SearchWrappers.MGMWrapper(lambda), N, b);
            this.stabilities = stabs.toArray();
        }
        lastLambda = lambda;
        return m.graphFromMGM();
    }
    //FOR PREF-DIV USE ONLY

    public int[][] runStepsPD()
    {
        Arrays.sort(lambda);
        int currIndex = lambda.length-1;
        double CC = -1;
        double CD = -1;
        double DD = -1;
        double CCMax = 100;
        double CCMaxI = -1;
        double CDMax = 100;
        double CDMaxI = -1;
        double DDMax = 100;
        double DDMaxI = -1;
        double oneLamb = -1;
        double allMax = 100;
        double allMaxI = -1;
        int p = 0;
        int q = 0;
        for(int i =0; i < d.getNumColumns();i++)
        {
            Node n = d.getVariable(i);
            if(n instanceof DiscreteVariable)
                q++;
            else
                p++;
        }
        // System.out.println("P:" + p);
        //   System.out.println("Q:" + q);
        //  System.out.println("b:" + b);
        A:while(true) //go until we break by having instability better than threshold
        {
            //   System.out.println("Lambda: " + lambda[currIndex]);
            double [][] adjMat = new double[d.getNumColumns()][d.getNumColumns()];
            for(int i = 0; i < N; i ++) //Learn MGM on N networks, update adjacency matrix after each run
            {
                Graph g = null;
                boolean again = true;
                while(again) {
                    //      System.out.println("here");
                    try {
                        DataSet temp = d.copy();
                        temp.permuteRows();
                        int[] removal = new int[b];
                        for (int j = 0; j < removal.length; j++)
                            removal[j] = j;
                        temp = temp.subsetRows(removal);
                        double[] l = {lambda[currIndex], lambda[currIndex], lambda[currIndex]};
                        MGM m = new MGM(temp, l);
                        m.learnEdges(iterLimit);
                        g = m.graphFromMGM();
                        again = false;
                    }
                    catch(Exception e)
                    {
                        e.printStackTrace();
                    }

                }


                for(int j = 0; j < d.getNumColumns();j++)
                {
                    for(int k = j+1; k < d.getNumColumns();k++)
                    {
                        if(g.isAdjacentTo(g.getNode(d.getVariable(j).getName()),g.getNode(d.getVariable(k).getName())))
                        {
                            adjMat[j][k] += 1/(double)N;
                        }

                    }

                }

            }
            double ccDestable = 0;
            double cdDestable = 0;
            double ddDestable = 0;
            ////////////////// TODO Decide if this is too harsh
            double numCC = 0;
            double numCD = 0;
            double numDD = 0;
            //////////////////
            for(int j = 0; j < d.getNumColumns();j++)
            {
                for(int k = j+1; k < d.getNumColumns();k++)
                {
                    Node one = d.getVariable(j);
                    Node two = d.getVariable(k);
                    if(one instanceof DiscreteVariable && two instanceof DiscreteVariable)
                    {
                        //          System.out.println("DD: " + adjMat[j][k]);
                        if(adjMat[j][k]!=0)
                            numDD++;
                        ddDestable += 2*adjMat[j][k]*(1-adjMat[j][k]);
                    }
                    else if(one instanceof DiscreteVariable || two instanceof DiscreteVariable)
                    {
                        if(adjMat[j][k]!=0)
                            numCD++;
                        //         System.out.println("CD: " + adjMat[j][k]);
                        cdDestable += 2*adjMat[j][k]*(1-adjMat[j][k]);

                    }
                    else
                    {
                        if(adjMat[j][k]!=0)
                            numCC++;
                        //       System.out.println("CC:" + adjMat[j][k]);
                        ccDestable+=2*adjMat[j][k]*(1-adjMat[j][k]);
                    }
                }
            }
            double allDestable = ccDestable + cdDestable + ddDestable;
            allDestable = allDestable/(numCC+numCD+numDD);
            if(allDestable <= gamma && oneLamb==-1)
                oneLamb = lambda[currIndex];
            if(allDestable <= allMax)
            {
                allMax = allDestable;
                allMaxI = lambda[currIndex];
            }
            /*
            ccDestable = ccDestable/ MathUtils.choose(p,2);
            cdDestable = cdDestable/(p*q);
            ddDestable = ddDestable/(MathUtils.choose(q,2));
            */
            ccDestable = ccDestable / numCC;
            cdDestable = cdDestable / numCD;
            ddDestable = ddDestable/numDD;
            // System.out.println("CC:" + ccDestable);
            // System.out.println("CD:" + cdDestable);
            // System.out.println("DD:" + ddDestable);
            if(ccDestable <= gamma && CC == -1)
                CC = lambda[currIndex];
            if(cdDestable <= gamma && CD==-1)
                CD = lambda[currIndex];
            if(ddDestable <= gamma && DD == -1)
                DD = lambda[currIndex];
            if(ccDestable <= CCMax)
            {
                CCMax = ccDestable;
                CCMaxI = lambda[currIndex];
            }
            if(cdDestable <= CDMax)
            {
                CDMax = cdDestable;
                CDMaxI = lambda[currIndex];
            }
            if(ddDestable <= DDMax)
            {
                DDMax = ddDestable;
                DDMaxI = lambda[currIndex];
            }
            if(CC!=-1 && CD != -1 && DD !=-1 && oneLamb!=-1)
                break A;
            if(currIndex==0)
                break A;
            currIndex--;
        }
        if(CC==-1)
            CC = CCMaxI;
        if(CD==-1)
            CD = CDMaxI;
        if(DD==-1)
            DD = DDMaxI;

        double [] lambda = {CC,CD,DD};
        if(oneLamb==-1)
            origLambda = allMaxI;
        else
            origLambda = oneLamb;
        //  System.out.println("Final Params: (" + CC + "," + CD + "," + DD + ")");
        int[][] result = new int[d.getNumColumns()][d.getNumColumns()];
            double [][] adjMat = new double[d.getNumColumns()][d.getNumColumns()];
            for(int i = 0; i < N; i ++) //Learn MGM on N networks, update adjacency matrix after each run
            {
                Graph g = null;
                boolean again = true;
                while (again) {
                    //      System.out.println("here");
                    try {
                        DataSet temp = d.copy();
                        temp.permuteRows();
                        int[] removal = new int[b];
                        for (int j = 0; j < removal.length; j++)
                            removal[j] = j;
                        temp = temp.subsetRows(removal);
                        MGM m = new MGM(temp, lambda);
                        m.learnEdges(iterLimit);
                        g = m.graphFromMGM();
                        for(Edge e:g.getEdges())
                        {
                            int x = d.getColumn(d.getVariable(e.getNode1().getName()));
                            int y = d.getColumn(d.getVariable(e.getNode2().getName()));
                            result[x][y]++;
                            result[y][x]++;
                        }
                        again = false;
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                }

            }
            MGM m = new MGM(d,lambda);
            m.learnEdges(1000);
            pdGraph = m.graphFromMGM();
            return result;
    }

    public static void main(String [] args)throws Exception
    {

        DataSet d = MixedUtils.loadDataSet2("Vaccine_Data_Edited.txt");
        d.removeColumn(d.getVariable("IL-12_induced"));
        d.removeColumn(d.getVariable("IL-10_induced"));
        d.removeColumn(d.getVariable("IL-12_mDC"));
        d.removeColumn(d.getVariable("IL-10_mDC"));
        //d.removeColumn(d.getVariable("IL-10_ratio"));
        //d.removeColumn(d.getVariable("IL-12_ratio"));
       // System.exit(0);
        int numLambdas = 40;
        double [] initLambdas = new double[numLambdas];
        for(int i = 0; i < initLambdas.length;i++)
        {
            initLambdas[i] = .05 + .8*i/numLambdas;
        }
        STEPS s = new STEPS(d,initLambdas,.07,5);
        Graph g = s.runSteps();
        PrintStream out = new PrintStream("steps.txt");
        out.println(g);
        out.flush();
        out.close();
        out = new PrintStream("max.txt");
        IndependenceTest i = new IndTestMultinomialAJ(d,.1);
        FciMaxP f = new FciMaxP(i);
        f.setInitialGraph(g);
        out.println(f.search());
        out.flush();
        out.close();
        CpcStable c = new CpcStable(i);
        out = new PrintStream("cpc.txt");
        out.println(c.search());
        out.flush();
        out.close();



   /*
        PrintStream out = new PrintStream("STEPS_Stability_Include_Zeros.txt");
        PrintStream out2 = new PrintStream("STEPS_Stability_Ignore_Zeros.txt");
        PrintStream out3 = new PrintStream("STEPS_Comparison.txt");
        out3.println("CC_0\tCC\tCD_0\tCD\tDD_0\tDD\tALL_0\tALL\tLAMBDAS_0\tLAMBDAS");
        out.println("CC\tCD\tDD\tALL");
        out2.println("CC\tCD\tDD\tALL");
        int numRuns = 10;
        int numLambdas = 30;
        double [][][] withZeros = new double[numRuns][numLambdas][4];
        double [][][] noZeros = new double[numRuns][numLambdas][4];
        for(int runs = 0; runs < numRuns; runs++) {
            System.out.println(runs);
            MixedLeeHastieSimulation m = new MixedLeeHastieSimulation();
            Parameters p = new Parameters();
            m.simulate(p);
            DataSet d = m.getDataSet(0);
            // DataSet d = MixedUtils.loadDataSet("C:/Users/vinee_000/Downloads","Survival_0_data.txt");
            double[] initLambda = new double[numLambdas];
            double [] initLambda2 = new double[numLambdas];
            for (int i = 0; i < numLambdas; i++) {
                initLambda[i] = .01 + .99 * i / numLambdas;
            }
            for (int i = 0; i < numLambdas; i++) {
                initLambda2[i] = .1 + .7 * i / numLambdas;
            }
            double g = 0.05;
            int numSub = 10;
            STEPS s = new STEPS(d, initLambda2, g, numSub);
            double[][] temp = s.runStepsArray(false);
            Graph noZero = s.lastGraph;
            double [] nzlam = s.lastLambda;
            s = new STEPS(d,initLambda2,g,numSub);
            double [][] temp2 = s.runStepsArray(true);
            Graph yesZero = s.lastGraph;
            double [] yslam = s.lastLambda;
            out3.print(mgmPriors.getF1(noZero,m.getTrueGraph(),m.getDataSet(0),"CC") + "\t");
            out3.print(mgmPriors.getF1(yesZero,m.getTrueGraph(),m.getDataSet(0),"CC") + "\t");
            out3.print(mgmPriors.getF1(noZero,m.getTrueGraph(),m.getDataSet(0),"CD") + "\t");
            out3.print(mgmPriors.getF1(yesZero,m.getTrueGraph(),m.getDataSet(0),"CD") + "\t");
            out3.print(mgmPriors.getF1(noZero,m.getTrueGraph(),m.getDataSet(0),"DD") + "\t");
            out3.print(mgmPriors.getF1(yesZero,m.getTrueGraph(),m.getDataSet(0),"DD") + "\t");
            out3.print(mgmPriors.getF1(noZero,m.getTrueGraph(),m.getDataSet(0),"All") + "\t");
            out3.print(mgmPriors.getF1(yesZero,m.getTrueGraph(),m.getDataSet(0),"All") + "\t");
            out3.println(Arrays.toString(yslam)+"\t"+Arrays.toString(nzlam));

            System.out.println("CC\tCD\tDD\tAll");
            for (int i = 0; i < temp.length; i++) {
                for (int j = 0; j < temp[0].length; j++) {
                    System.out.print(temp[i][j] + "\t");
                }
                System.out.println();
            }
           withZeros[runs] = temp2;
           noZeros[runs] = temp;

        }
        for(int i =0; i < numLambdas;i++)
        {
            for(int j = 0; j < 4;j++)
            {
                double wz = 0;
                double nz = 0;
                int numWz = 0;
                int numNz = 0;
                for(int k = 0; k < numRuns;k++)
                {
                    if(!Double.isNaN(withZeros[k][i][j])) {
                        numWz++;
                        wz += withZeros[k][i][j];
                    }
                    if(!Double.isNaN(noZeros[k][i][j])) {
                        numNz++;
                        nz += noZeros[k][i][j];
                    }
                }
                if(numWz!=0)
                wz=wz/numWz;
                if(numNz!=0)
                nz=nz/numNz;
                out.print(wz+"\t");
                out2.print(nz+"\t");
            }
            out.println();
            out2.println();
        }
        out.flush();
        out2.flush();
        out3.flush();
        out.close();
        out2.close();
        out3.close();


*/

    }
}
