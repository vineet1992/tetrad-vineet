package edu.pitt.csb.mgm;

import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndTestMultinomialLogisticRegression;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradVector;
import edu.cmu.tetrad.util.dist.ChiSquare;
import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.linear.RealVector;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;

/**
 * Kernel Based Conditional Independence Test
 * Code Written by: Vineet Raghu
 * Test published in: Kernel-based Conditional Independence Test and Application in Causal Discovery (Zhang et al.)
 */
public class KCI_Ind implements IndependenceTest {
    HashMap<String,Integer> setSize;
    HashMap<String,String> corrRemovals; //correct changes
    HashMap<String,Integer> justNames;
    HashMap<String,Double> setP;
    HashMap<String,String> corrP;
    double alpha;
    DataSet dat;
    int graphNum;
    double lastP;
    int timesCalledCC;
    int timesCalledCCZ;
    int timesCalledCD;
    int timesCalledDD;
    int total;
    int changesCC;
    int changesCD;
    int changesDD;
    Graph tg;
    int correctCC;
    int correctCD;
    int cutoff = -1;
    int correctDD;

    int timesCalledCDZ;
    int timesCalledDDZ;
    int totalZ;
    int changesCCZ;
    int changesCDZ;
    int changesDDZ;
    int correctCCZ;
    int correctCDZ;
    int correctDDZ;
    IndTestMultinomialAJ ii;
    ArrayList<Double> chiSquareRand;

    public KCI_Ind(DataSet data,DataSet data2, double threshold)
    {
        total = 0;
        timesCalledCC = 0;
        timesCalledCD = 0;
        timesCalledDD = 0;
        changesCC = 0;
        changesCD = 0;
        changesDD = 0;
        correctCC = 0;
        correctCD = 0;
        correctDD = 0;
        timesCalledCCZ = 0;
        timesCalledCDZ = 0;
        timesCalledDDZ = 0;
         totalZ = 0;
         changesCCZ= 0;
         changesCDZ= 0;
         changesDDZ= 0;
         correctCCZ= 0;
         correctCDZ= 0;
         correctDDZ= 0;
        alpha = threshold;
        dat = data2;
        ii = new IndTestMultinomialAJ(data,alpha);
        lastP = -1;
        for(int i = 0; i < dat.getNumColumns();i++)
        {
            double [] curr = new double[dat.getNumRows()];
            for(int j = 0; j < dat.getNumRows();j++)
            {
                curr[j] = dat.getDouble(j,i);

            }
            double m = mean(curr);
            double std = std(curr,m);
            if(i==0){
                System.out.println("Mean: " + m);
                System.out.println("STD: " + std);
            }
            for(int j = 0; j < dat.getNumRows();j++)
            {
                dat.setDouble(j,i,(dat.getDouble(j,i)-m)/std);
            }
        }
        chiSquareRand = new ArrayList<Double>();
        ChiSquare cs = new ChiSquare(1);
        for(int i = 0; i < 1000*dat.getNumRows();i++)
        {
            chiSquareRand.add(cs.nextRandom());
        }
        //Normalize dataset data
        //Subtract by mean and divide by standard deviation

    }
    public KCI_Ind(DataSet data, DataSet data2, double threshold,double threshold2, Graph trueGraph,int x)
    {
        setSize = new HashMap<String,Integer>();
        corrRemovals = new HashMap<String,String>();//Comma separeted list of conditioning varaibles whenever we correctly save the variables
        justNames = new HashMap<String,Integer>(); //X1,X2 for all variables that we save the edge correctly
        setP = new HashMap<String,Double>();
        corrP = new HashMap<String,String>(); //Conditioning set when variable is incorrectly removed
        cutoff = x;
        tg = trueGraph;

        timesCalledCCZ = 0;
        timesCalledCDZ = 0;
        timesCalledDDZ = 0;
        totalZ = 0;
        changesCCZ= 0;
        changesCDZ= 0;
        changesDDZ= 0;
        correctCCZ= 0;
        correctCDZ= 0;
        correctDDZ= 0;
        total = 0;
        timesCalledCC = 0;
        timesCalledCD = 0;
        timesCalledDD = 0;
        changesCC = 0;
        changesCD = 0;
        changesDD = 0;
        correctCC = 0;
        correctCD = 0;
        correctDD = 0;
        alpha = threshold;
        dat = data2.copy();
        ii = new IndTestMultinomialAJ(data,threshold2);
        lastP = -1;
        int [] dis = MixedUtils.getDiscreteInds(data.getVariables());
       int count = 0;
        for(int i = 0; i < dat.getNumColumns();i++)
        {
           if(dis[count] == i)
           {
               if(count != dis.length-1)
                count++;
               continue;
           }
            double [] curr = new double[dat.getNumRows()];
            for(int j = 0; j < dat.getNumRows();j++)
            {
                curr[j] = dat.getDouble(j,i);

            }
            double m = mean(curr);
            double std = std(curr,m);
            if(i==0){
                System.out.println("Mean: " + m);
                System.out.println("STD: " + std);
            }
            for(int j = 0; j < dat.getNumRows();j++)
            {
                dat.setDouble(j,i,(dat.getDouble(j,i)-m)/std);
            }
        }
        chiSquareRand = new ArrayList<Double>();
        ChiSquare cs = new ChiSquare(1);
        for(int i = 0; i < 1000*dat.getNumRows();i++)
        {
            chiSquareRand.add(cs.nextRandom());
        }

    }
    public double getScore()
    {
        return 0;
    }
    public static double mean(double[] curr) {

        double sum = 0;
        int num = curr.length;
        for (int i = 0; i < num; i++)
            sum += curr[i];
        return sum / num;
    }

    public static double std(double[] curr, double mean) {
        double sum = 0;
        for (int i = 0; i < curr.length; i++) {
            sum += Math.pow((curr[i] - mean), 2);
        }
        sum = sum / curr.length;
        return Math.sqrt(sum);

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
    public boolean isIndependent(Node x, Node y, List<Node> z)
    {
        total++;
ArrayList<Node> zzzz = new ArrayList<Node>();
        for(Node tempo: z)
            zzzz.add(ii.getVariable(tempo.getName()));
        if(ii.isDependent(ii.getVariable(x.getName()),ii.getVariable(y.getName()),zzzz))
            return false;
        boolean xCont = false;
        boolean yCont = false;
        if(z.size() > cutoff)
            return false;
        try {

            ii.getData().setDouble(1, ii.getData().getColumn(ii.getVariable(x.getName())),ii.getData().getDouble(1, ii.getData().getColumn(ii.getVariable(x.getName()))));
            xCont = true;
           // System.out.println(x + ": is Continuous");
        }
        catch (Exception e){}
        try {
            ii.getData().setDouble(1, ii.getData().getColumn(ii.getVariable(y.getName())),ii.getData().getDouble(1, ii.getData().getColumn(ii.getVariable(y.getName()))));
            yCont = true;
         //System.out.println(y + ": is Continuous");
        }
        catch(Exception ee){}
       if(xCont && yCont) {
           if(z.isEmpty()) {
               timesCalledCCZ++;
           }
           else {
               timesCalledCC++;
           }
       }
        else if(!xCont && !yCont)
       {
           if(z.isEmpty()) {
               timesCalledDDZ++;
           }
           else {
               timesCalledDD++;
           }
       }

        else{
           if(z.isEmpty()) {
               timesCalledCDZ++;
           }
           else {
               timesCalledCD++;
           }
       }


        if(z.isEmpty())
        {
            boolean b =  isIndependentUncon(x, y);
            if(!b)
            {
                if(xCont && yCont)
                    changesCCZ++;

                else if(!xCont && !yCont)
                    changesDDZ++;
                else
                    changesCDZ++;
                ArrayList<Node> zz = new ArrayList<Node>();
                for(Node zzz: z)
                zz.add(tg.getNode(zzz.getName()));
                if(tg.isDConnectedTo(tg.getNode(x.getName()),tg.getNode(y.getName()),zz)) // said dependent is dependent
                 {
                    if(xCont && yCont) {

                            correctCCZ++;

                    }
                    else if(!xCont && !yCont) {

                        correctDDZ++;

                    }
                    else {

                         correctCDZ++;

                    }

                        corrRemovals.put(x.getName() + "," + y.getName(), "0");
                        //corrP.put(x.getName()+y.getName(),Double.toString(lastP));
                        justNames.put(x.getName() + "," + y.getName(), 0);
                     if(tg.isAdjacentTo(tg.getNode(x.getName()),tg.getNode(y.getName())))
                         System.out.println(x + " " + y + " are adjacent, and we saved the edge");
                }

            }
            return b;
        }
        else
        {
            boolean b = isIndependentCon(x, y, z);
            if(!b) {
                if(xCont && yCont)
                    changesCC++;
                else if(!xCont && !yCont)
                    changesDD++;
                else
                    changesCD++;
                ArrayList<Node> zz = new ArrayList<Node>();
                for(Node zzz: z)
                    zz.add(tg.getNode(zzz.getName()));
                if(tg.isDConnectedTo(tg.getNode(x.getName()),tg.getNode(y.getName()),zz)) //correctly said dependent for true answer dependent
                {
                    if(xCont && yCont)
                        correctCC++;
                    else if(!xCont && !yCont)
                        correctDD++;
                    else
                        correctCD++;
                    if(corrRemovals.get(x.getName()+ "," + y.getName())!=null) {
                        corrRemovals.put(x.getName() + "," + y.getName(), corrRemovals.get(x.getName() + "," + y.getName()) + "," + z.toString());
                        if(tg.isAdjacentTo(tg.getNode(x.getName()),tg.getNode(y.getName())))
                            System.out.println(x + " " + y + " edge was saved again with condition set " + z);
                        //corrP.put(x.getName()+y.getName(),corrP.get(x.getName()+y.getName()) + "," + Double.toString(lastP));
                    }
                    else{
                        corrRemovals.put(x.getName() + "," + y.getName(), z.toString());
                        if(tg.isAdjacentTo(tg.getNode(x.getName()),tg.getNode(y.getName())))
                            System.out.println(x + " " + y + " edge was saved first time with condition set " + z);
                        //corrP.put(x.getName()+y.getName(),Double.toString(lastP));
                        justNames.put(x.getName() + "," + y.getName(), 0);
                    }
                }

            }
            else if(justNames.get(x.getName()+ "," +y.getName())!=null)
            {
                ArrayList<Node> zz = new ArrayList<Node>();
                for(Node zzz: z)
                    zz.add(tg.getNode(zzz.getName()));
                if(tg.isDConnectedTo(tg.getNode(x.getName()),tg.getNode(y.getName()),zz)) //incorrectly removed edge previously saved
                {
                    corrP.put(x.getName() + "," + y.getName(),zz.toString());
                    if(tg.isAdjacentTo(tg.getNode(x.getName()),tg.getNode(y.getName())))
                        System.out.println(x + " " + y + " edge was removed incorrectly, despite previous save with condition set " + zz);
                }
                else
                {

                    //corrRemovals.put(x.getName()+y.getName(),zz.size()+1);
                }
            }
            return b;
        }
    }
    private boolean isIndependentCon(Node x, Node y, List<Node> z)
    {
        boolean unbiased = false;
        boolean GP = unbiased;
        boolean bootstrap = true;
        int colnumX = dat.getColumn(x);
        int colnumY = dat.getColumn(y);
        int T = dat.getNumRows();
        double [] xArr = new double[T];
        double[] yArr = new double[T];
        for(int i = 0; i < T; i ++)
        {
            xArr[i] = dat.getDouble(i,colnumX);
            yArr[i] = dat.getDouble(i,colnumY);
        }
        int num_eig = T;
        int T_BS = 5000;
        double lambda = 1E-3;
        double thres = 1E-5;
        int dim = z.size();
        double width = 0;
        if(T <=200)
            width = 1.2;
        else if (T < 1200)
            width = 0.7;
        else
            width = 0.4;
        double theta = 1/(width*width*dim);
        double [] [] H = new double[T][T];
        long time = System.currentTimeMillis();
        for(int i = 0; i < T; i ++)
        {
            for(int j = 0; j < T; j++)
            {
                if(i==j)
                    H[i][j] = 1.0 - 1.0 / T;
                else
                    H[i][j] = -1.0/T;
            }
        }
       TetradMatrix Hmat = new TetradMatrix(H);

        TetradMatrix kernArg = new TetradMatrix(T,z.size()+1);
        for(int i = 0; i < T;i++)
        {
            for(int j = 0; j < z.size()+1; j++)
            {
                if(j==0)
                  kernArg.set(i,j,xArr[i]);
                else
                    kernArg.set(i,j,dat.getDouble(i,dat.getColumn(z.get(j-1)))/2);
            }
        }
        //System.out.println("Time to setup preliminary kernel matrices: " + (System.currentTimeMillis()-time));
        double [] temp = new double[2];
        temp[0] = theta;
        temp[1] = 1;
        TetradMatrix KX = kernel(kernArg, kernArg, temp);
        kernArg = null;
        time = System.currentTimeMillis();
        KX = Hmat.times(KX.times(Hmat));
        double [] [] ky = kernel(yArr,yArr,temp);
        TetradMatrix KY = new TetradMatrix(ky);

        KY = Hmat.times(KY.times(Hmat));

        TetradMatrix KXZ = new TetradMatrix(1,1);
        TetradMatrix KYZ = new TetradMatrix(1,1);
        double sta = 0;
        double df = 0;
       // System.out.println("Time for first two Matrix Multiplications: " + (System.currentTimeMillis()-time));
        if(GP)
        {
            //TO DO optimize hyperparameters

        }
        else
        {
            time = System.currentTimeMillis();
            TetradMatrix KZ = new TetradMatrix(T,z.size());
            for(int i = 0; i < T; i ++)
            {
                for(int j = 0; j < z.size();j++)
                {
                    KZ.set(i,j,dat.getDouble(i,dat.getColumn(z.get(j))));
                }
            }
           KZ = kernel(KZ, KZ, temp);
            KZ = Hmat.times(KZ.times(Hmat));
            TetradMatrix eye = new TetradMatrix(T,T);
            for(int i = 0; i < T; i ++)
                eye.set(i,i,1);
            KZ = eye.minus(KZ.times((KZ.plus(eye.scalarMult(lambda)).inverse())));
            KXZ = KZ.times(KX.times(KZ.transpose()));
            KYZ = KZ.times(KY.times(KZ.transpose()));
            sta = KXZ.times(KYZ).trace();
            df = eye.minus(KZ).trace();
            //System.out.println("Time for setting up kernel matrices(Multiplication): " + (System.currentTimeMillis()-time));
           // System.out.println(sta);
            //System.out.println(df);
        }
        EigenDecomposition ed1;
        EigenDecomposition ed2;
        try{
            time = System.currentTimeMillis();
         ed1 = new EigenDecomposition(KXZ.plus(KXZ.transpose()).scalarMult(.5).getRealMatrix());
         ed2 = new EigenDecomposition(KYZ.plus(KYZ.transpose()).scalarMult(.5).getRealMatrix());
           // System.out.println("Time for EigenValue Decomp:" + (System.currentTimeMillis()-time));
        }
        catch(Exception e)
    {
        System.out.println("Eigenvalue didn't converge");
        return true;
    }



time = System.currentTimeMillis();
        double [] evalues1 = ed1.getRealEigenvalues();
        double [] evalues2 = ed2.getRealEigenvalues();

        //find all indices in evalues1 and evalues2 where the eigenvalue is greater than max eigenvalue * thresh
        //keep only those indices
        //get eigenvectors corresponding to those indices as well

        double max1 = 0;
        for(int i = 0 ; i < evalues1.length;i++)
        {
            if(evalues1[i] > max1)
                max1 = evalues1[i];

        }
        double max2 = 0;
        for(int i = 0; i < evalues2.length;i++)
        {
            if(evalues2[i] > max2)
                max2 = evalues2[i];
        }
        ArrayList<Integer>inds1 = new ArrayList<Integer>();
        for(int i = 0; i < evalues1.length;i++)
        {
            if(evalues1[i] > max1*thres)
                inds1.add(i);
        }
        ArrayList<Integer>inds2 = new ArrayList<Integer>();
        for(int i = 0; i < evalues2.length;i++)
        {
            if(evalues2[i] > max2*thres)
            {
                inds2.add(i);
            }
        }

        TetradMatrix eigKxz = new TetradMatrix(inds1.size(),inds1.size());
       TetradMatrix eigKyz = new TetradMatrix(inds2.size(),inds2.size());

        TetradMatrix tv1 = new TetradMatrix(ed1.getEigenvector(0).getDimension(),inds1.size());
        TetradMatrix tv2 = new TetradMatrix(ed2.getEigenvector(0).getDimension(),inds2.size());

            for (int i = 0; i < inds1.size(); i++) {

                eigKxz.set(i, i, Math.sqrt(evalues1[inds1.get(i)]));
                RealVector t = ed1.getEigenvector(inds1.get(i));


                for (int j = 0; j < t.getDimension(); j++) {
                    tv1.set(j, i, t.getEntry(j));
                }
                //columns of matrix are eigenvectors
                //rows are values of the e-vectors
            }

        for(int i = 0; i < inds2.size();i++)
        {
            eigKyz.set(i,i,Math.sqrt(evalues2[inds2.get(i)]));
            RealVector t = ed2.getEigenvector(inds2.get(i));
            for(int j = 0; j < t.getDimension();j++)
            {
                tv2.set(j,i,t.getEntry(j));
            }
        }

        //eigKXZ = diag(sqrt(evalues1))
        //tv1 = relevant eigenvectors from chosen eigenvalues

        TetradMatrix eiv_prodx = tv1.times(eigKxz.transpose());
        TetradMatrix eiv_prody = tv2.times(eigKyz.transpose());

        eigKxz = null;
        eigKyz = null;
        tv1 = null;
        tv2 = null;
        int numx = eiv_prodx.columns();
        int numy = eiv_prody.columns();
        int size_u = numx*numy;
        TetradMatrix uu = new TetradMatrix(T,size_u);
        for(int i = 0; i < numx; i++)
        {
            for(int j = 0; j < numy; j++)
            {
                for(int k = 0; k < T ;k++)
                {
                    uu.set(k,i*numy+j,eiv_prodx.get(k,i)*eiv_prody.get(k, j));
                }
            }
        }
        TetradMatrix uu_prod;
        if(size_u > T)
            uu_prod = uu.times(uu.transpose());
        else
            uu_prod = uu.transpose().times(uu);
        //System.out.println("Time to do a lot of multiplication before generating null: " + (System.currentTimeMillis()-time));


        if(bootstrap)
        {
            //TO DO
            time = System.currentTimeMillis();
            EigenDecomposition ee;
            try {
                ee = new EigenDecomposition(uu_prod.getRealMatrix());
            }
            catch(Exception e)
            {
                System.out.println("Eigenvalue Didn't converge conditional");
                return true;
            }
            int num = -1;
            if(T < size_u)
                num = T;
            else
                num = size_u;

            double [] evals = ee.getRealEigenvalues();
            double [] valsToKeep = new double[num];
            Arrays.sort(evals);
            int count = 0;
            for(int i = evals.length-1; i>=0;i--)
            {
                valsToKeep[count] = evals[i];
                count++;
            }
            double max = valsToKeep[0];
            ArrayList<Double>finalVals = new ArrayList<Double>();
            for(int i = 0; i < valsToKeep.length;i++)
            {
                if(valsToKeep[i] >= max*thres)
                    finalVals.add(valsToKeep[i]);
            }
            double cri = -1;
            double p_val = -1;

            if(finalVals.size()*T < 1E6)
            {
                double [] [] frand1 = new double[finalVals.size()][T_BS];
                for(int i = 0; i < finalVals.size();i++)
                {
                    for(int j = 0; j < T_BS;j++)
                    {
                        frand1[i][j] = chiSquareRand.get((i*T_BS+j)%chiSquareRand.size());
                    }
                }
                double [] [] eiguu = new double[1][finalVals.size()];
                for(int j = 0; j < finalVals.size();j++)
                {
                    eiguu[0][j] = finalVals.get(j);
                }
                TetradMatrix fr = new TetradMatrix(frand1);
                TetradMatrix eig_uu = new TetradMatrix(eiguu);
                TetradMatrix nullDist;
                if(unbiased)
                {
                    System.out.println("Can only return unbiased if hyperparameters are learned");
                    return false;
                }
                else
                {
                    nullDist = eig_uu.times(fr);
                }
                try{
                    PrintStream out = new PrintStream("debug.txt");
                    out.println(eig_uu);
                    out.flush();
                    out.close();
                }
                catch(Exception e)
                {
                    return false;
                }
                cri = Math.ceil((1-alpha)*T_BS);
                int sum = 0;
                for(int i = 0; i < nullDist.columns();i++)
                {
                    if(nullDist.get(0,i)>sta)
                        sum++;
                }
                lastP = sum/(double)T_BS;
               // System.out.println("Time to generate null and get result: " + (System.currentTimeMillis()-time));
                if(lastP > alpha) {

                    return true;
                }
                else {

                    return false;
                }
            }
            else
            {
                System.out.println("Unimplemented iteratively calculating null");
                return false;
            }
        }
        else
        {
            double cri_appr = -1;
            double p_appr = -1;
            double mean_appr = uu_prod.trace();
            double var_appr = 2*uu_prod.times(uu_prod).trace();
            double k_appr = mean_appr*mean_appr/var_appr;
            double theta_appr = var_appr/mean_appr;
            GammaDistribution g = new GammaDistribution(k_appr,theta_appr);
            cri_appr = g.inverseCumulativeProbability(1-alpha);
          //  System.out.println(cri_appr);
            p_appr =1- g.cumulativeProbability(sta);
            lastP = p_appr;
            if(p_appr > alpha)
                return true;
            else {

                return false;
            }
        }
    }

    private boolean isIndependentUncon(Node x, Node y)
    {
        //System.out.println(x + " " + y);
        int columnNum = dat.getColumn(x);
        int columnNumY = dat.getColumn(y);
        int T = dat.getNumRows();

        double [] xArr = new double[T];
        double [] yArr = new double[T];

        for(int i = 0; i < T; i++)
        {
            xArr[i] = dat.getDouble(i,columnNum);
            yArr[i] = dat.getDouble(i,columnNumY);
        }

        boolean approx = false;
        int num_eig = T;
        if(T>1000) {
            approx = true;
            num_eig = T/2;
        }
        int T_BS = 1000;
        double lambda = 1E-3;
        double thresh = 1E-6;
        double width = .8;
        if(T > 200)
            width = 0.5;
        if (T > 1200)
            width = 0.3;
        double theta = 1/(width*width);
        double [] [] H = new double[T][T];
        for(int i = 0; i < T; i ++)
        {
            for(int j = 0; j < T; j++)
            {
                if(i==j)
                    H[i][j] = 1.0 - 1.0 / T;
                else
                    H[i][j] = -1.0/T;
            }
        }
        double [] temp = new double[2];
        temp[0] = theta;
        temp[1] = 1;
        double [][] Kx = kernel(xArr,xArr,temp);
        double [][] Ky = kernel(yArr,yArr,temp);

        TetradMatrix kx = new TetradMatrix(Kx);

        TetradMatrix ky = new TetradMatrix(Ky);
        TetradMatrix h = new TetradMatrix(H);
        long time = System.currentTimeMillis();
        kx = h.times(kx);
        kx = kx.times(h);
        ky = h.times(ky);
        ky = ky.times(h);
        h = null;
        double sta = kx.times(ky).trace();
       // System.out.println("Time to calculate statistic: " + (System.currentTimeMillis()-time));
        //System.out.println(sta);
        double cri = -1;
        double p_val = -1;
        if(!approx) //Bootstrap
        {
           // System.out.println("Not approximating");
            try {
                PrintStream p = new PrintStream("debug.txt");
                //System.out.println(kx);
                p.println(kx.plus(kx.transpose()).scalarMult(0.5));

            }
            catch(Exception e)
            {
                return false;
            }
            EigenDecomposition ed1;
            EigenDecomposition ed2;
            try {
                time = System.currentTimeMillis();
             ed1 = new EigenDecomposition(kx.plus(kx.transpose()).scalarMult(0.5).getRealMatrix());
                ed2 = new EigenDecomposition(ky.plus(ky.transpose()).scalarMult(0.5).getRealMatrix());
               // System.out.println("Time to do EigenValue: " + (System.currentTimeMillis()-time));
            }
            catch(Exception e)
            {
                System.out.println("Eigenvalue didn't converge");
                return true;
            }
            time = System.currentTimeMillis();
            double [] ev1 = ed1.getRealEigenvalues();
            double [] ev2 = ed2.getRealEigenvalues();
            double [] eigProd = new double[ev1.length*ev2.length];
            double maxEig = Double.MIN_VALUE;
            for(int i = 0; i < ev1.length; i++)
            {
                for(int j = 0; j < ev2.length;j++)
                {
                    double curr = ev1[i]*ev2[j];
                    if(curr > maxEig)
                        maxEig = curr;
                    eigProd[i*ev2.length+j] = curr;

                }
            }
            ArrayList<Double> d = new ArrayList<Double>();
            for(int i = 0; i < eigProd.length;i++){
                if(eigProd[i] > maxEig*thresh)
                    d.add(eigProd[i]);
            }
           // System.out.println("Time to do some more multiplication: " + (System.currentTimeMillis()-time));
            if(d.size()*T < 1E6)
            {
                //THIS IS TAKING FOREVER
                //TO DO CACHE THE NULL AND REUSE THE CRAP OUT OF IT
                time = System.currentTimeMillis();
                //System.out.println(d);
                //ChiSquare chi = new ChiSquare(1);
                double[][] f_rand1 = new double[d.size()][T_BS];
                for(int i = 0; i < d.size();i++)
                {
                    for(int j = 0; j < T_BS;j++)
                    {
                            f_rand1[i][j] = chiSquareRand.get((i*T_BS+j)%chiSquareRand.size());
                    }
                }
               // System.out.println("Time to generate null: " + (System.currentTimeMillis()-time));
                time = System.currentTimeMillis();
                double [][] data = new double[1][d.size()];
                for(int i = 0; i < d.size();i++)
                    data[0][i] = d.get(i);

                TetradMatrix f_rand = new TetradMatrix(f_rand1);
                TetradMatrix ep = new TetradMatrix(data);
                ep = ep.scalarMult(1/(double)T);
                double [] [] nullDist = ep.times(f_rand).toArray();
                int sum = 0;
               // System.out.println("Time for some multiplication before p-value: " + (System.currentTimeMillis()-time));
                time = System.currentTimeMillis();
                for(int i = 0; i < nullDist.length;i++)
                {

                    for(int j = 0; j < nullDist[i].length;j++)
                    {
                        if(nullDist[i][j] > sta)
                            sum++;
                    }
                }
                //System.out.println("Got to bottom");
                double pval = (double)sum/T_BS;
                lastP = pval;
                //System.out.println("Time to n2 calculation of p-value: " + (System.currentTimeMillis()-time));
                if(pval < alpha) {

                    return false;
                }
                else
                    return true;
            }
            else
            {
                //TODO Iteratively calculate the null distribution!
               int length = 0;
                double [] nullDist = new double[T_BS];
                if(1E6/T > 100)
                    length = (int)( 1E6/T);
                else
                    length = 100;
                int itmax = (int)Math.floor(d.size()/length);
                TetradVector nd = new TetradVector(nullDist);

                for(int iter = 0; iter<itmax;iter++)
                {
                    ChiSquare cs = new ChiSquare(1);
                    double [] [] f_rand1 = new double[(int)length][T_BS];
                    for(int j = 0; j < length;j++)
                    {
                        for(int k = 0; k < T_BS;k++)
                        {
                            f_rand1[j][k] = cs.nextRandom();
                        }
                    }
                    int start = iter*length + 1;
                    int end = (iter+1)*length;
                    //for(int index = start; index <=end; index++)
                    //{

                    //}

                    //need to iteratively compute null_distr here
                }
                System.out.println("TO DO CALCULATE NULL ITERATIVELY");
                return false;
                //do final computation, and then compute p-value as above
            }

            //need to get num_eig top eigenvalues here of KX + KX' /2 and same for KY


        }
        else
        {
            time = System.currentTimeMillis();
            double mean_appr = kx.trace()*ky.trace()/T;
            double var_appr = 2*kx.times(kx).trace()*ky.times(ky).trace()/(T*T);//can optimize by not actually performing matrix multiplication
            double k_appr = mean_appr*mean_appr/var_appr;
            double theta_appr = var_appr/mean_appr;
            GammaDistribution g = new GammaDistribution(k_appr,theta_appr);
            double p_appr = 1 - g.cumulativeProbability(sta);
            lastP = p_appr;
           // System.out.println("Time to approximate p value: " + (System.currentTimeMillis()-time));
            if(p_appr < alpha) {

                return false;
            }
            else
                return true;
        }

    }
    private static double [][] dist(double [] x, double [] y)
    {
        double [][] sum = new double[x.length][y.length];
        for(int i = 0; i < x.length;i++) {
            for (int j = 0; j < y.length; j++) {
                sum[i][j] = Math.pow(x[i] - y[j], 2);

            }
        }
        return sum;
    }
    private static double dist2(double[]x,double[]y)
    {
        double sum = 0;
        for(int i = 0; i < x.length;i++)
        {
            sum+= Math.pow(x[i]-y[i],2);
        }
        return sum;
    }
    private static TetradMatrix kernel(TetradMatrix x, TetradMatrix xKern, double [] theta)
    {
        TetradMatrix result = new TetradMatrix(x.rows(),xKern.rows());
        for(int i = 0; i < x.rows();i++)
        {
            double[] currRow = x.getRow(i).toArray();

            for(int j = 0; j < xKern.rows();j++)
            {
                double [] secRow = xKern.getRow(j).toArray();
                result.set(i,j,Math.exp(-1*dist2(currRow,secRow)*theta[0]/2));
            }
        }
        return result;
    }
    private static double[][]kernel(double [] x, double [] y,double []theta )
    {
        double [][]n2 = dist(x,y);
        double wi2 = theta[0]/2;
        double [][] kx = new double[x.length][y.length];
        for(int i = 0; i < x.length;i++)
        {
            for(int j = 0; j < y.length;j++)
            {
                kx[i][j] = Math.exp(-1*n2[i][j]*wi2);
            }
        }
        return kx;

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
