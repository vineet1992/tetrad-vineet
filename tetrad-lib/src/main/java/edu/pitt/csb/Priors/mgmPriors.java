package edu.pitt.csb.Priors;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.RegressionDataset;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.MGM;
import edu.pitt.csb.mgm.MGM_Priors;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;

import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

/**
 * Created by vinee_000 on 4/25/2017.
 */
public class mgmPriors {


    Random rand = new Random();
    private final double normalEpsilon = 0.5;
    private double gEpsilon;
    private double noPrior = 0;
    private int numSubsamples;
    private double[][] lambdas;
    private DataSet data;
    private TetradMatrix[] priors;
    private boolean[][] havePriors; //did any source provide priors about edge i,j?
    private boolean[][][] sourcePrior; //did source k provide prior about edge i,j?
    private final int iterLimit = 100;
    public double oracleCC;
    public double oracleCD;
    public double oracleDD;
    public double oracleAll;
    public boolean splitLambdas = false;
    private DataSet [] subsamples;
    public Graph trueGraph;
    private PrintStream log;
    public double [] expertWeights;
    public double [] normalizedExpertWeights;
    public double [] normalizedTao;
    public double [] pValues;
    private boolean logging = false;
    private boolean normalize = true;
    private int normalizationSamples = 10000;
    public double [][] edgeScores;

    public void setLog(boolean b)
    {
        logging = b;
    }


    //If this is set to true, the method will normalize by computing a null distribution for each size of prior
    public void setNormalize(boolean b)
    {
        normalize = b;
    }
    public void setNormalizationSamples(int x)
    {
        normalizationSamples = x;
    }
    public mgmPriors(int numSubsamples, double[] initLambdas, DataSet data, TetradMatrix[] priors) {
        this.numSubsamples = numSubsamples;
        this.pValues = new double[priors.length];
        this.lambdas = constructLambdasPar(initLambdas, data);
        this.data = data;
        this.priors = priors;
        this.sourcePrior = new boolean[data.getNumColumns()][data.getNumColumns()][priors.length];
        havePriors = findPrior(priors);
        gEpsilon = 1 / (double) numSubsamples;
    }


    public mgmPriors(int numSubsamples, double[] initLambdas, DataSet data, TetradMatrix[] priors,DataSet [] subsamples) {
        this.subsamples = subsamples;
        this.pValues = new double[priors.length];
        this.numSubsamples = numSubsamples;
        if(logging) {
            try {
                this.log = new PrintStream("log_file.txt");
            } catch (Exception e) {
            }
        }
        this.lambdas = constructLambdasPar(initLambdas, data);

        this.data = data;
        this.priors = priors;
        this.sourcePrior = new boolean[data.getNumColumns()][data.getNumColumns()][priors.length];
        havePriors = findPrior(priors);

        gEpsilon = 1 / (double) numSubsamples;
    }


    public mgmPriors(int numSubsamples, double[] initLambdas, DataSet data, TetradMatrix[] priors, Graph truth,boolean splitLambdas,DataSet [] subsamples) {
        this.subsamples = subsamples;
        this.pValues = new double[priors.length];
        this.splitLambdas = splitLambdas;
        this.numSubsamples = numSubsamples;
        if(logging) {
            try {
                this.log = new PrintStream("log_file.txt");
            } catch (Exception e) {
            }
        }

        this.trueGraph = truth;
        this.lambdas = constructLambdasPar(initLambdas, data);
        this.data = data;
        this.priors = priors;
        this.sourcePrior = new boolean[data.getNumColumns()][data.getNumColumns()][priors.length];
        havePriors = findPrior(priors);

        gEpsilon = 1 / (double) numSubsamples;
    }

    public mgmPriors(int numSubsamples, double[] initLambdas, DataSet data, TetradMatrix[] priors, Graph truth,boolean splitLambdas) {
        this.splitLambdas = splitLambdas;
        this.numSubsamples = numSubsamples;
        this.pValues = new double[priors.length];
        if(logging) {
            try {
                this.log = new PrintStream("log_file.txt");
            } catch (Exception e) {
            }
        }
        this.trueGraph = truth;
        this.lambdas = constructLambdasPar(initLambdas, data);
        this.data = data;
        this.priors = priors;
        this.sourcePrior = new boolean[data.getNumColumns()][data.getNumColumns()][priors.length];
        havePriors = findPrior(priors);

        gEpsilon = 1 / (double) numSubsamples;
    }



    public Graph runPriors() {
        if(logging) {
            try {
                log.println("True Graph: " + trueGraph);
                log.flush();
            } catch (Exception e) {

            }
        }
        TetradMatrix [][] edgeCounts = null;
         edgeCounts = getEdgeProbabilitiesPar(); //matrix of numLambdas x numSubsamples that contains a matrix of edge appearences for each edge at each index
//edgeCounts = getEdgeProbabilities();

        TetradMatrix counts = getFullCounts(edgeCounts); // # of times each edge appeared
        //Across all lambda*subsamples graphs
        TetradMatrix variances = getVariances(counts); // variance of appearences of each edge
        //  System.out.println("Variances: " + variances);
        TetradMatrix[] phi = new TetradMatrix[priors.length];
        // TetradMatrix[] tao = new TetradMatrix[priors.length];
        double[] tao = new double[priors.length];
        double [] normalTao = new double[priors.length];
        double[] alpha = new double[priors.length];
        double[] weights = new double[priors.length];
        if (priors.length == 1) {

        }
        for (int tr = 0; tr < priors.length; tr++) //for each source of prior information
        {
            TetradMatrix currPrior = priors[tr];
            phi[tr] = getPhi(currPrior);
            //How many times do we expect an edge to appear based on the prior, for a given lambda?
            //phi = prior value * numSubSamples

            tao[tr] = getTao(phi[tr], counts,tr);

            if(normalize)
                normalTao[tr] = normalizeTao(tao[tr],counts,tr);
            //How confident are we about source tr?
            //tao = average deviation between our predicted probability from prior (phi) and actual counts u
        }
        if(normalize)
        {
            normalizedTao = normalTao;
            double [] normAlpha = getAlpha(normalTao);
            normalizedExpertWeights = getWeights(normAlpha);
            weights = normalizedExpertWeights;
            tao = normalizedTao;
        }
        else
        {
            alpha = getAlpha(tao);
            weights = getWeights(alpha);
            expertWeights = weights;
        }

        TetradMatrix u_mixture = getMeanMixture(phi, weights);
        TetradMatrix var_mixture = getVarMixture(weights, tao, u_mixture, phi);

        //   System.out.println("U Mixture: " + u_mixture);
        //   System.out.println("Variance Mixture: " + var_mixture);

        //use counts divided by numLambdas for the posterior calculation!
        TetradMatrix u_posterior = getMeanPosterior(counts, variances, u_mixture, var_mixture);

        //   System.out.println("U Posterior: " + u_posterior);
        //  System.out.println("Prior Available?");

        TetradMatrix var_posterior = getVarPosterior(variances, var_mixture);
        TetradMatrix[] f = getF(edgeCounts);
        //    System.out.println("F:" + f[0]);
        TetradMatrix[] g = getG(f);
        //   System.out.println("G:" + g[0]);
        TetradMatrix[] theta = getTheta(edgeCounts, u_posterior, counts, var_posterior);
        //    System.out.println("Theta: " + theta[0]);

        double[][] scoresNP = getScoresNP(theta, g); // 0th row has CC, 1st row has CD, 2nd row has DD for no prior edges
        if(logging)
        {
            try {
                log.println("Scores for Edges without prior knowledge:");
                log.println(Arrays.toString(scoresNP[2]));
            } catch (Exception e) {
            }
        }
        double[][] scoresWP = getScoresWP(theta, g); // 0th row has CC, 1st row has CD, 2nd row has DD for edges w/prior
        if(logging) {
            try {
                log.println("Scores for Edges with Prior Knoweldge:");
                log.println(Arrays.toString(scoresWP[2]));
            } catch (Exception e) {

            }
        }
        double[][] NPIndex = new double[3][2]; //CC, CD, DD rows and Index, Max Value columns
        double[][] WPIndex = new double[3][2];

        for (int i = 0; i < lambdas[0].length; i++) {
            for (int j = 0; j < 3; j++) {
                if (scoresNP[j][i] > NPIndex[j][1]) {
                    NPIndex[j][1] = scoresNP[j][i];
                    NPIndex[j][0] = i;
                }
                if (scoresWP[j][i] > WPIndex[j][1]) {
                    WPIndex[j][1] = scoresWP[j][i];
                    WPIndex[j][0] = i;
                }
            }
        }


        double[] npLambdas = {lambdas[0][(int) NPIndex[0][0]], lambdas[1][(int) NPIndex[1][0]], lambdas[2][(int) NPIndex[2][0]]};
        double[] wpLambdas = {lambdas[0][(int) WPIndex[0][0]], lambdas[1][(int) WPIndex[1][0]], lambdas[2][(int) WPIndex[2][0]]};
        System.out.println("Lambdas no prior: " + Arrays.toString(npLambdas));
        System.out.println("Lambdas with prior: " + Arrays.toString(wpLambdas));

        if(logging) {
            try {
                log.println("Lambdas no prior:");
                log.println(Arrays.toString(npLambdas));
                log.println("Lambdas with prior:");
                log.println(Arrays.toString(wpLambdas));
                log.flush();
                log.close();
            } catch (Exception e) {

            }
        }

        edgeScores = constructEdgeScores(npLambdas,wpLambdas,havePriors);
        MGM_Priors m = new MGM_Priors(data, npLambdas, wpLambdas, havePriors);
        m.learnEdges(1000);
        return m.graphFromMGM();


        //We now have the max indices for each edge type and with or without priors.
        //Now we need to figure out how to modify MGM to use the 6 lambdas instead of 3
    }

    private double[] getAlpha(double[] tao) {
        //HEURISTIC HACK TO ENSURE THAT THERE IS NO ZERO VALUE
        double min = Double.MAX_VALUE;
        for(int i = 0; i < tao.length;i++)
        {
            if(tao[i]!=0 && tao[i] < min)
                min = tao[i];
        }
        for(int i = 0; i < tao.length;i++)
        {
            if(tao[i]==0)
                tao[i] = min/2;
        }
        double taoSum = StatUtils.sum(tao);
        double[] alpha = new double[tao.length];

        for (int i = 0; i < tao.length; i++)
            alpha[i] = taoSum / tao[i];
        return alpha;
    }



    private double [][] constructEdgeScores(final double [] npLambdas, final double [] wpLambdas,final boolean [][] havePriors)
    {
        final double [][] edgeCounts = new double[data.getNumColumns()][data.getNumColumns()];
        final ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

        class StabilityAction extends RecursiveAction {
            private int chunk;
            private int from;
            private int to;

            public StabilityAction(int chunk, int from, int to){
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }


            private synchronized void updateEdgeCounts(double [][] t, Graph g)
            {
                for(int i = 0; i < data.getNumColumns();i++)
                {
                    for(int j = i+1; j < data.getNumColumns();j++)
                    {
                        if(g.getEdge(g.getNode(data.getVariable(i).getName()),g.getNode(data.getVariable(j).getName()))!=null)
                        {
                            t[i][j]++;
                            t[j][i]++;
                        }

                    }
                }
            }

            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        //temp is an int array with the samples for the current subsampling
                        DataSet current = subsamples[s];
                        MGM_Priors m = new MGM_Priors(current,npLambdas,wpLambdas,havePriors);
                        m.learnEdges(iterLimit);
                        Graph g = m.graphFromMGM();


                        updateEdgeCounts(edgeCounts,g);

                    }

                    return;
                } else {
                    List<StabilityAction> tasks = new ArrayList<>();

                    final int mid = (to + from) / 2;

                    tasks.add(new StabilityAction(chunk, from, mid));
                    tasks.add(new StabilityAction(chunk, mid, to));

                    invokeAll(tasks);

                    return;
                }
            }

        }

        final int chunk = 2;

        pool.invoke(new StabilityAction(chunk, 0, numSubsamples));

        if(logging)
            log.flush();
        return edgeCounts;
    }
    private double[] getWeights(double[] alpha) {
        double alpSum = StatUtils.sum(alpha);
        double[] weights = new double[alpha.length];
        for (int i = 0; i < alpha.length; i++)
            weights[i] = alpha[i] / alpSum;
        return weights;
    }


    private double getTao(TetradMatrix phi, TetradMatrix counts, int tr) {
        double tao = 0;
        int numPriors = 0;
        for (int i = 0; i < counts.rows(); i++) {
            for (int j = i + 1; j < counts.columns(); j++) {
                if (sourcePrior[i][j][tr]) {
                    tao += Math.abs(phi.get(i, j) - counts.get(i, j) / (double) lambdas[0].length);
                    numPriors++;
                }
            }
        }
        tao = tao / numPriors;
        return tao;
    }

    //Normalizes the deviance score by computing a null distribution with the same number of entries as the current prior
    private double normalizeTao(double tao,TetradMatrix counts, int tr)
    {
        double [] hist = new double[normalizationSamples];
       for(int i = 0; i < normalizationSamples;i++)
        {
            TetradMatrix rand = generateRandomPrior(tr);
            rand = getPhi(rand);
            double currTao = 0;
            int numEdges = 0;
            if(counts.rows()!=rand.rows()||counts.columns()!=rand.columns())
                System.exit(-1);
            for(int j = 0; j < counts.rows();j++)
            {
                for(int k = j+1; k<counts.columns();k++)
                {
                    if(rand.get(j,k)!=0)
                    {
                        currTao+=Math.abs(rand.get(j,k) -  counts.get(j,k)/(double)lambdas[0].length);
                        numEdges++;
                    }
                }
            }
            currTao=currTao/(double)numEdges;
            hist[i] = currTao;
        }
        Arrays.sort(hist);
       System.out.println(tao);
       System.out.println(Arrays.toString(hist));
        int index =0;
        while(tao > hist[index] && index < hist.length)
        {
            index++;
        }
        pValues[tr] = index/(double)normalizationSamples;
        tao = tao/StatUtils.mean(hist);

        return tao;
    }

    //Generates a random hard prior with the same number of entries as sourcePrior[][][k] (The kth expert)
    private TetradMatrix generateRandomPrior(int k)
    {
        int count = 0;
        for(int i = 0; i < data.getNumColumns();i++)
        {
            for(int j = i+1; j < data.getNumColumns();j++)
            {
                if(sourcePrior[i][j][k])
                    count++;
            }
        }
        TetradMatrix t = new TetradMatrix(data.getNumColumns(),data.getNumColumns());
        while(count >0)
        {
            int x = rand.nextInt(data.getNumColumns());
            int y = rand.nextInt(data.getNumColumns());
            while(t.get(x,y)==1 || x==y)
            {
                x = rand.nextInt(data.getNumColumns());
                y = rand.nextInt(data.getNumColumns());
            }
            t.set(x,y,1);
            t.set(y,x,1);
            count--;
        }
        return t;

    }

    private double[][] getScoresWP(TetradMatrix[] theta, TetradMatrix[] g) {
        double[][] scores = new double[3][lambdas[0].length];
        for (int i = 0; i < theta.length; i++) {
            TetradMatrix currTheta = theta[i];
            TetradMatrix currG = g[i];
            for (int j = 0; j < data.getNumColumns(); j++) {
                for (int k = j + 1; k < data.getNumColumns(); k++) {

                    if (currG.get(j, k) == 0)
                        currG.set(j, k, gEpsilon);
                    if (havePriors[j][k]) {
                        if (data.getVariable(j) instanceof ContinuousVariable && data.getVariable(k) instanceof ContinuousVariable)
                            scores[0][i] += currTheta.get(j, k) * (1 - currG.get(j, k));
                        else if (data.getVariable(j) instanceof DiscreteVariable && data.getVariable(k) instanceof DiscreteVariable)
                            scores[2][i] += currTheta.get(j, k) * (1 - currG.get(j, k));
                        else
                            scores[1][i] += currTheta.get(j, k) * (1 - currG.get(j, k));
                    }
                }
            }
        }
        return scores;
    }

    private double[][] getScoresNP(TetradMatrix[] theta, TetradMatrix[] g) {
        double[][] scores = new double[3][lambdas[0].length];
        //What to do if G = 0
        for (int i = 0; i < theta.length; i++) {
            TetradMatrix currTheta = theta[i];
            TetradMatrix currG = g[i];
            for (int j = 0; j < data.getNumColumns(); j++) {
                for (int k = j + 1; k < data.getNumColumns(); k++) {
                    if (currG.get(j, k) == 0) //TODO WHAT TO DO IF G == 0
                        currG.set(j, k, gEpsilon);
                    if (!havePriors[j][k]) {
                        if (data.getVariable(j) instanceof ContinuousVariable && data.getVariable(k) instanceof ContinuousVariable)
                            scores[0][i] += currTheta.get(j, k) * (1 - currG.get(j, k));
                        else if (data.getVariable(j) instanceof DiscreteVariable && data.getVariable(k) instanceof DiscreteVariable)
                            scores[2][i] += currTheta.get(j, k) * (1 - currG.get(j, k));
                        else
                            scores[1][i] += currTheta.get(j, k) * (1 - currG.get(j, k));
                    }
                }
            }
        }
        return scores;
    }

    private TetradMatrix[] getTheta(TetradMatrix[][] edgeCounts, TetradMatrix u_posterior, TetradMatrix counts, TetradMatrix var_posterior) {
        TetradMatrix[] theta = new TetradMatrix[lambdas[0].length];
        double q_j = numSubsamples * lambdas[0].length;

        TetradMatrix probs = counts.scalarMult(1 / q_j); //holds the probability of edge appearence (for no priors)
        if(logging)
         log.println(probs);
        for (int i = 0; i < lambdas[0].length; i++) {
            TetradMatrix temp = new TetradMatrix(data.getNumColumns(), data.getNumColumns());
            for (int ss = 0; ss < numSubsamples; ss++) {
                TetradMatrix thisCount = edgeCounts[i][ss];
                for (int j = 0; j < data.getNumColumns(); j++) {
                    for (int k = j + 1; k < data.getNumColumns(); k++) {
                        temp.set(j, k, temp.get(j, k) + thisCount.get(j, k)); //how many times did the edge appear for this particular lambda value?

                    }
                }
            }
            TetradMatrix currTheta = new TetradMatrix(data.getNumColumns(), data.getNumColumns());
            for (int ii = 0; ii < data.getNumColumns(); ii++) {
                for (int jj = ii + 1; jj < data.getNumColumns(); jj++) {
                    if (havePriors[ii][jj]) //TODO finish computing the probabilities for each edge, check which sources have prior information?
                    {
                        NormalDistribution b;
                        try {
                            b = new NormalDistribution(u_posterior.get(ii, jj), Math.sqrt(var_posterior.get(ii, jj)));
                            currTheta.set(ii, jj, b.probability(temp.get(ii, jj) - normalEpsilon, temp.get(ii, jj) + normalEpsilon));
                        } catch (Exception e) //If there is no variance, then only give probability if we match the posterior mean exactly
                        {
                            if (var_posterior.get(ii, jj) == 0) {
                                //           System.out.println("Temp: " + temp.get(ii,jj) + ", Posterior: " + u_posterior.get(ii,jj));
                                if ((int) temp.get(ii, jj) == (int) (u_posterior.get(ii, jj) / lambdas[0].length))
                                    currTheta.set(ii, jj, 1);
                                else
                                    currTheta.set(ii, jj, 0);
                            } else {
                                System.out.println(u_posterior.get(ii, jj) / (q_j) + "i: " + ii + ",j:" + jj);
                                throw e;
                            }
                        }

                    } else {
                        BinomialDistribution b = new BinomialDistribution(numSubsamples, probs.get(ii, jj));
                        currTheta.set(ii, jj, b.probability((int) temp.get(ii, jj)));
                        //if(data.getVariable(ii) instanceof DiscreteVariable && data.getVariable(jj) instanceof DiscreteVariable)
                        //System.out.println("Edge appeared " + temp.get(ii,jj) + ", Expected probability: " + probs.get(ii,jj) + ", Computed probability: " + b.probability((int)temp.get(ii,jj)));

                    }
                }
            }
            theta[i] = currTheta;
        }
        return theta;
    }

    private TetradMatrix[] getG(TetradMatrix[] f) {
        TetradMatrix[] g = new TetradMatrix[f.length];
        for (int i = 0; i < g.length; i++) {
            TetradMatrix temp = new TetradMatrix(f[0].rows(), f[0].columns());
            for (int j = 0; j < temp.rows(); j++) {
                for (int k = j + 1; k < temp.columns(); k++) {
                    temp.set(j, k, 4 * f[i].get(j, k) * (1 - f[i].get(j, k)));
                }
            }
            g[i] = temp;
        }
        if(logging) {
            for (int i = 0; i < g.length; i++) {
                log.println("Instability Matrix for Lambda: " + i);
                log.println(g[i]);
                log.flush();
            }
        }
        return g;
    }

    private TetradMatrix[] getF(TetradMatrix[][] edgeCounts) {
        TetradMatrix[] f = new TetradMatrix[lambdas[0].length];
        for (int i = 0; i < lambdas[0].length; i++) {
            TetradMatrix temp = new TetradMatrix(edgeCounts[0][0].rows(), edgeCounts[0][0].columns()); //the counts of each edge appearing for this lambda
            for (int j = 0; j < numSubsamples; j++) {
                TetradMatrix countMat = edgeCounts[i][j]; //the edge count matrix for this i,j combo
                for (int k = 0; k < countMat.rows(); k++) {
                    for (int m = k + 1; m < countMat.columns(); m++) {
                        temp.set(k, m, temp.get(k, m) + countMat.get(k, m));
                    }
                }
            }
            f[i] = temp.scalarMult(1 / (double) numSubsamples);
        }
        return f;
    }

    private TetradMatrix getFullCounts(TetradMatrix[][] edgeCounts) {
        TetradMatrix temp = new TetradMatrix(edgeCounts[0][0].rows(), edgeCounts[0][0].columns());
        for (int i = 0; i < temp.rows(); i++) {
            for (int j = i + 1; j < temp.columns(); j++) {
                int count = 0;
                for (int k = 0; k < edgeCounts.length; k++) {
                    for (int m = 0; m < edgeCounts[0].length; m++) {
                        count += edgeCounts[k][m].get(i, j);
                    }

                }
                temp.set(i, j, count);
            }
        }
        return temp;
    }

    private TetradMatrix getVarPosterior(TetradMatrix variances, TetradMatrix var_mixture) {
        TetradMatrix temp = new TetradMatrix(variances.rows(), variances.columns());
        for (int i = 0; i < temp.rows(); i++) {
            for (int j = i + 1; j < temp.columns(); j++) {
                temp.set(i, j, (variances.get(i, j) * var_mixture.get(i, j)) / (variances.get(i, j) + var_mixture.get(i, j)));
            }
        }
        return temp;
    }

    private TetradMatrix getMeanPosterior(TetradMatrix counts, TetradMatrix variances, TetradMatrix u_mixture, TetradMatrix var_mixture) {
        TetradMatrix temp = new TetradMatrix(counts.rows(), counts.columns());
        for (int i = 0; i < temp.rows(); i++) {
            for (int j = i + 1; j < temp.columns(); j++) {
                temp.set(i, j, ((counts.get(i, j) / lambdas[0].length) * var_mixture.get(i, j) + variances.get(i, j) * u_mixture.get(i, j)) / (variances.get(i, j) + var_mixture.get(i, j)));
            }
        }
        return temp;
    }

    private TetradMatrix getVarMixture(double[] weights, double[] tao, TetradMatrix u_mixture, TetradMatrix[] phi) {

        TetradMatrix temp = new TetradMatrix(phi[0].rows(), phi[0].columns());
        for (int i = 0; i < temp.rows(); i++) {
            for (int j = i + 1; j < temp.columns(); j++) {
                double sum = 0;
                for (int k = 0; k < priors.length; k++) {

                    if (sourcePrior[i][j][k])
                        sum += weights[k] * tao[k] * tao[k] + Math.pow(u_mixture.get(i, j) - phi[k].get(i, j), 2);
                }
                temp.set(i, j, sum);
            }
        }
        return temp;
    }

    private TetradMatrix getMeanMixture(TetradMatrix[] phi, double[] weights) {

        TetradMatrix temp = new TetradMatrix(phi[0].rows(), phi[0].columns());
        for (int i = 0; i < temp.rows(); i++) {
            for (int j = i + 1; j < temp.columns(); j++) {

                double sum = 0;
                for (int k = 0; k < priors.length; k++) {

                    if (sourcePrior[i][j][k]) {


                        sum += weights[k] * phi[k].get(i, j);
                    }
                }
                temp.set(i, j, sum);
            }
        }
        return temp;
    }

    private TetradMatrix abs(TetradMatrix t) {
        TetradMatrix temp = new TetradMatrix(t.rows(), t.columns());
        for (int i = 0; i < t.rows(); i++) {
            for (int j = 0; j < t.columns(); j++) {
                temp.set(i, j, Math.abs(t.get(i, j)));
            }
        }
        return temp;
    }



    private TetradMatrix getPhi(TetradMatrix priors) {
        TetradMatrix temp = new TetradMatrix(priors.rows(), priors.columns());
        for (int i = 0; i < priors.rows(); i++) {
            for (int j = i + 1; j < priors.rows(); j++) {
                temp.set(i, j, priors.get(i, j) * numSubsamples);
            }
        }
        return temp;
    }

    private TetradMatrix getVariances(TetradMatrix counts) {
        TetradMatrix temp = new TetradMatrix(counts.rows(), counts.columns());
        for (int i = 0; i < counts.rows(); i++) {
            for (int j = i + 1; j < counts.rows(); j++) {
                temp.set(i, j, counts.get(i, j) * (1 - (counts.get(i, j) / (numSubsamples * lambdas[0].length))));
            }
        }
        return temp;
    }

    private void shuffleArray(int[] ar) {
        // If running on Java 6 or older, use `new Random()` on RHS here
        for (int i = ar.length - 1; i > 0; i--) {
            int index = rand.nextInt(i + 1);
            // Simple swap
            int a = ar[index];
            ar[index] = ar[i];
            ar[i] = a;
        }
    }

    private TetradMatrix[][] getEdgeProbabilities() {
        TetradMatrix[][] edgeCounts = new TetradMatrix[lambdas[0].length][numSubsamples];
        for (int i = 0; i < lambdas[0].length; i++) {
            for (int j = 0; j < numSubsamples; j++) {
                edgeCounts[i][j] = new TetradMatrix(data.getNumColumns(), data.getNumColumns());
            }
        }
        int numSamples = data.getNumRows() / numSubsamples;
        int[] allRows = new int[data.getNumRows()];
        for (int i = 0; i < data.getNumRows(); i++) {
            allRows[i] = i;
        }
        shuffleArray(allRows);
        for (int i = 0; i < numSubsamples; i++) {
            int[] temp = new int[numSamples];
            for (int k = i * numSamples; k < (i + 1) * numSamples; k++) {
                temp[k - i * numSamples] = allRows[k];
            }
            //temp is an int array with the samples for the current subsampling
            DataSet current = data.subsetRows(temp);
            for (int j = 0; j < lambdas[0].length; j++) {
                double[] lambda = {lambdas[0][j], lambdas[1][j], lambdas[2][j]};
                MGM m = new MGM(current, lambda);
                m.learnEdges(iterLimit);
                Graph g = m.graphFromMGM();
                if(logging) {
                    log.println(i + "," + "Subsample: " + j);
                    log.println(g);
                    log.println(edgeCounts[j][i]);
                    log.flush();
                }
                updateCounts(edgeCounts[j][i], g);
            }
        }

        return edgeCounts;

    }

    public TetradMatrix[][] getEdgeProbabilitiesPar()
    {
        final TetradMatrix [][] edgeCounts = new TetradMatrix[lambdas[0].length][numSubsamples];
        for(int i = 0; i < lambdas[0].length;i++)
        {
            for(int j = 0; j< numSubsamples;j++)
            {
                edgeCounts[i][j] = new TetradMatrix(data.getNumColumns(),data.getNumColumns());
            }
        }
        final ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

        class StabilityAction extends RecursiveAction {
            private int chunk;
            private int from;
            private int to;

            public StabilityAction(int chunk, int from, int to){
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }


            private synchronized void updateEdgeCounts(TetradMatrix t, Graph g)
            {
                for(int i = 0; i < data.getNumColumns();i++)
                {
                    for(int j = i+1; j < data.getNumColumns();j++)
                    {
                        if(g.getEdge(g.getNode(data.getVariable(i).getName()),g.getNode(data.getVariable(j).getName()))!=null)
                        {
                            t.set(i,j,t.get(i,j)+1);
                        }

                    }
                }
            }
            private synchronized void debugLog(int i, int j, TetradMatrix[][] edgeCounts, Graph g)
            {
                log.println(i +  "," + "Subsample: " + j);
                log.println(g);
                log.println(edgeCounts[i][j]);
                log.flush();
            }

            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {

                        int j = s/lambdas[0].length;
                        int i = s%lambdas[0].length;
                        double [] lambda = {lambdas[0][i],lambdas[1][i],lambdas[2][i]};

                        //temp is an int array with the samples for the current subsampling
                        DataSet current = subsamples[j];
                        MGM m = new MGM(current,lambda);
                        m.learnEdges(iterLimit);
                        Graph g = m.graphFromMGM();


                        updateEdgeCounts(edgeCounts[i][j],g);
                        if(logging)
                            debugLog(i,j,edgeCounts,g);

                    }

                    return;
                } else {
                    List<StabilityAction> tasks = new ArrayList<>();

                    final int mid = (to + from) / 2;

                    tasks.add(new StabilityAction(chunk, from, mid));
                    tasks.add(new StabilityAction(chunk, mid, to));

                    invokeAll(tasks);

                    return;
                }
            }

        }

        final int chunk = 2;

        pool.invoke(new StabilityAction(chunk, 0, numSubsamples*lambdas[0].length));

        if(logging)
         log.flush();
        return edgeCounts;
    }
    private void updateCounts(TetradMatrix t, Graph g)
    {
        for(int i = 0; i < data.getNumColumns();i++)
        {
            for(int j = i+1; j < data.getNumColumns();j++)
            {
                if(g.getEdge(g.getNode(data.getVariable(i).getName()),g.getNode(data.getVariable(j).getName()))!=null)
                {
                    t.set(i,j,t.get(i,j)+1);
                }

            }
        }
    }
    public static double getF1(Graph est, Graph truth, DataSet d,String type)
    {
        double tp = 0;
        double fp = 0;
        double fn = 0;
        if(type.equals("CC"))
        {
            for(Edge e:est.getEdges())
            {
                if(d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable)
                {
                    Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()),truth.getNode(e.getNode2().getName()));
                    if(temp!=null)
                    {
                        tp++;
                    }
                    else
                        fp++;
                }
            }
            for(Edge e:truth.getEdges())
            {
                if(d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable)
                {
                    Edge temp = est.getEdge(est.getNode(e.getNode1().getName()),est.getNode(e.getNode2().getName()));
                    if(temp==null)
                        fn++;
                }
            }
        }
        else if(type.equals("CD"))
        {
            for(Edge e:est.getEdges())
            {
                if((d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable) ||  (d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable))
                {
                    Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()),truth.getNode(e.getNode2().getName()));
                    if(temp!=null)
                    {
                        tp++;
                    }
                    else
                        fp++;
                }
            }
            for(Edge e:truth.getEdges())
            {
                if((d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable) ||  (d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable))
                {
                    Edge temp = est.getEdge(est.getNode(e.getNode1().getName()),est.getNode(e.getNode2().getName()));
                    if(temp==null)
                        fn++;
                }
            }
        }
        else if(type.equals("DD"))
        {
            for(Edge e:est.getEdges())
            {
                if(d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable)
                {
                    Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()),truth.getNode(e.getNode2().getName()));

                    if(temp!=null)
                    {
                        tp++;
                    }
                    else
                        fp++;
                }
            }
            for(Edge e:truth.getEdges())
            {
                if(d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable)
                {
                    Edge temp = est.getEdge(est.getNode(e.getNode1().getName()),est.getNode(e.getNode2().getName()));
                    if(temp==null)
                        fn++;
                }
            }
        }
        else
        {
            for(Edge e:est.getEdges())
            {
                    Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()),truth.getNode(e.getNode2().getName()));
                    if(temp!=null)
                    {
                        tp++;
                    }
                    else
                        fp++;
            }
            for(Edge e:truth.getEdges())
            {
                    Edge temp = est.getEdge(est.getNode(e.getNode1().getName()),est.getNode(e.getNode2().getName()));
                    if(temp==null)
                        fn++;
            }
        }

        double prec = tp/(tp+fp);
        double rec = tp/(tp+fn);
        return (2*prec*rec)/(prec+rec);
    }


//TODO Test this method and make sure it's the same as the other constructLambdas
    public double [][] constructLambdasPar(final double[]init,final DataSet data)
    {
        final double [] oracles = new double[4];
        final int numLambdas = init.length;
        final double [][] numEdges = new double[init.length][3];
        final double [] bestF1 = new double[4];
        final ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

        class StabilityAction extends RecursiveAction {
            private int chunk;
            private int from;
            private int to;

            public StabilityAction(int chunk, int from, int to){
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }

            //could avoid using syncronized if we keep track of array of mats and add at end, but that needs lots of
            //memory
            private synchronized void addToMat(Graph g, double [][] numEdges, int s){
                numEdges[s][0] = getEdges(g,data,"CC");
                numEdges[s][1] = getEdges(g,data,"CD");
                numEdges[s][2] = g.getNumEdges()-numEdges[s][0]-numEdges[s][1];

            }
            private synchronized void updateOracle(Graph curr,double lbd,DataSet data, double [] bestF1,double [] oracles) {
                if (trueGraph != null) {
                    double F1CC = getF1(curr, trueGraph, data, "CC");
                    if (F1CC > bestF1[0] || (F1CC==bestF1[0]&&lbd < oracles[0])) {
                        bestF1[0] = F1CC;
                        oracles[0] = lbd;
                    }
                    double F1CD = getF1(curr, trueGraph, data, "CD");
                    if (F1CD > bestF1[1] || (F1CD == bestF1[1]&&lbd < oracles[1])) {
                        bestF1[1] = F1CD;
                        oracles[1] = lbd;
                    }
                    double F1DD = getF1(curr, trueGraph, data, "DD");
                    if (F1DD > bestF1[2] || (F1DD ==bestF1[2]&&lbd < oracles[2])) {

                        bestF1[2] = F1DD;
                        oracles[2] = lbd;
                    }
                    double F1All = getF1(curr, trueGraph, data, "All");
                    if (F1All > bestF1[3] || (F1All==bestF1[3]&&lbd < oracles[3])) {
                        bestF1[3] = F1All;
                        oracles[3] = lbd;
                    }
                }
            }
            private synchronized double getF1(Graph est, Graph truth, DataSet d,String type) {
                double tp = 0;
                double fp = 0;
                double fn = 0;
                if (type.equals("CC")) {
                    for (Edge e : est.getEdges()) {
                        if (d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable) {
                            Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()), truth.getNode(e.getNode2().getName()));
                            if (temp != null) {
                                tp++;
                            } else
                                fp++;
                        }
                    }
                    for (Edge e : truth.getEdges()) {
                        if (d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable) {
                            Edge temp = est.getEdge(est.getNode(e.getNode1().getName()), est.getNode(e.getNode2().getName()));
                            if (temp == null)
                                fn++;
                        }
                    }
                } else if (type.equals("CD")) {
                    for (Edge e : est.getEdges()) {
                        if ((d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable) || (d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable)) {
                            Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()), truth.getNode(e.getNode2().getName()));
                            if (temp != null) {
                                tp++;
                            } else
                                fp++;
                        }
                    }
                    for (Edge e : truth.getEdges()) {
                        if ((d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof ContinuousVariable) || (d.getVariable(e.getNode1().getName()) instanceof ContinuousVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable)) {
                            Edge temp = est.getEdge(est.getNode(e.getNode1().getName()), est.getNode(e.getNode2().getName()));
                            if (temp == null)
                                fn++;
                        }
                    }
                } else if (type.equals("DD")) {
                    for (Edge e : est.getEdges()) {
                        if (d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable) {
                            Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()), truth.getNode(e.getNode2().getName()));

                            if (temp != null) {
                                tp++;
                            } else
                                fp++;
                        }
                    }
                    for (Edge e : truth.getEdges()) {
                        if (d.getVariable(e.getNode1().getName()) instanceof DiscreteVariable && d.getVariable(e.getNode2().getName()) instanceof DiscreteVariable) {
                            Edge temp = est.getEdge(est.getNode(e.getNode1().getName()), est.getNode(e.getNode2().getName()));
                            if (temp == null)
                                fn++;
                        }
                    }
                } else {
                    for (Edge e : est.getEdges()) {
                        Edge temp = truth.getEdge(truth.getNode(e.getNode1().getName()), truth.getNode(e.getNode2().getName()));
                        if (temp != null) {
                            tp++;
                        } else
                            fp++;
                    }
                    for (Edge e : truth.getEdges()) {
                        Edge temp = est.getEdge(est.getNode(e.getNode1().getName()), est.getNode(e.getNode2().getName()));
                        if (temp == null)
                            fn++;
                    }
                }

                double prec = tp / (tp + fp);
                double rec = tp / (tp + fn);
                return (2 * prec * rec) / (prec + rec);
            }
            @Override
            protected void compute(){
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        double[] lambda = {init[s], init[s], init[s]};
                        MGM m = new MGM(data, lambda);
                        m.learnEdges(iterLimit);
                        Graph curr = m.graphFromMGM();

                        addToMat(curr, numEdges, s);
                        updateOracle(curr,init[s],data,bestF1,oracles);

                        }
                    return;
                }
                else {
                    List<StabilityAction> tasks = new ArrayList<>();

                    final int mid = (to + from) / 2;

                    tasks.add(new StabilityAction(chunk, from, mid));
                    tasks.add(new StabilityAction(chunk, mid, to));

                    invokeAll(tasks);

                    return;
                }
            }

        }

        final int chunk = 2;

        pool.invoke(new StabilityAction(chunk, 0, numLambdas));
        oracleCC = oracles[0];
        oracleCD = oracles[1];
        oracleDD = oracles[2];
        oracleAll = oracles[3];
        TetradMatrix t = new TetradMatrix(numEdges);


        if(splitLambdas)
        {
               TetradMatrix ccMat = new TetradMatrix(init.length,2);
               TetradMatrix cdMat = new TetradMatrix(init.length,2);
               TetradMatrix ddMat = new TetradMatrix(init.length,2);
               for(int i = 0; i < init.length;i++)
               {
                   ccMat.set(i,0,init[i]);
                   ccMat.set(i,1,numEdges[i][0]);
                   cdMat.set(i,0,init[i]);
                   cdMat.set(i,1,numEdges[i][1]);
                   ddMat.set(i,0,init[i]);
                   ddMat.set(i,1,numEdges[i][2]);
               }
               double [] limitsCC = getLimit(ccMat);
               double [] limitsCD = getLimit(cdMat);
               double [] limitsDD = getLimit(ddMat);
               double [][] realLimits = new double[3][init.length];
            for (int i = 0; i < init.length; i++) {
                realLimits[0][i] = map(init[i], init, limitsCC);
                realLimits[1][i] = map(init[i],init,limitsCD);
                realLimits[2][i] = map(init[i],init,limitsDD);
            }
            if(logging) {
                log.println("Lambda Range CC:" + Arrays.toString(realLimits[0]));
                log.println("Lambda Range CD:" + Arrays.toString(realLimits[1]));
                log.println("Lambda Range DD:" + Arrays.toString(realLimits[2]));
                log.flush();
            }
            return realLimits;
        }
        else {
            TetradMatrix allMat = new TetradMatrix(init.length,2);
            for (int i = 0; i < init.length; i++) {
                allMat.set(i, 0, init[i]);
                allMat.set(i, 1, numEdges[i][0] + numEdges[i][1] + numEdges[i][2]);

            }
            double[] limits = getLimit(allMat);

            /*double[] realLimits = new double[init.length];
            for (int i = 0; i < init.length; i++) {
                realLimits[i] = map(init[i], init, limits);
            }*/
            double[][] realLimits = new double[3][init.length];
            for (int i = 0; i < init.length; i++) {
                realLimits[0][i] = map(init[i], init, limits);
                realLimits[1][i] = map(init[i],init,limits);
                realLimits[2][i] = map(init[i],init,limits);
            }
            return realLimits;
        }




    }
    //This method takes the initial input lambdas and finds an appropriate range for lambda based on knee points in the graph of lambda vs number of edges in learned graph
    private double[] constructLambdas(double[]init,DataSet data)
    {
        double [][] numEdges = new double[init.length][3];

        double [] bestF1 = new double [3];
        double bestAll = 0;
        for(int i = 0; i < init.length;i++) //learn a// n MGM for each initial lambda and save the adjacency matrices
        {
            double [] lambda = {init[i],init[i],init[i]};
            MGM m =new MGM(data,lambda);
            m.learnEdges(iterLimit);
            Graph curr = m.graphFromMGM();
            if(trueGraph!=null)
            {
                double F1CC = mgmPriors.getF1(curr,trueGraph,data,"CC");
                if(F1CC > bestF1[0])
                {
                    bestF1[0] = F1CC;
                    oracleCC = init[i];
                }
                double F1CD = mgmPriors.getF1(curr,trueGraph,data,"CD");
                if(F1CD > bestF1[1])
                {

                    bestF1[1] = F1CD;
                    oracleCD = init[i];
                }
                double F1DD = mgmPriors.getF1(curr,trueGraph,data,"DD");
                if(F1DD > bestF1[2])
                {
                    bestF1[2] = F1DD;
                    oracleDD = init[i];
                }
                double F1All = mgmPriors.getF1(curr,trueGraph,data,"All");
                if(F1All > bestAll)
                {
                    bestAll = F1All;
                    oracleAll = init[i];
                }
            }

            numEdges[i][0] = getEdges(curr,data,"CC");
            numEdges[i][1] = getEdges(curr,data,"CD");
            numEdges[i][2] = curr.getNumEdges()-numEdges[i][0]-numEdges[i][1];
        }
        TetradMatrix t = new TetradMatrix(numEdges);
      /*  TetradMatrix ccMat = new TetradMatrix(init.length,2);
        TetradMatrix cdMat = new TetradMatrix(init.length,2);
        TetradMatrix ddMat = new TetradMatrix(init.length,2);*/
      TetradMatrix allMat = new TetradMatrix(init.length,2);

        for(int i = 0; i < init.length;i++)
        {
            allMat.set(i,0,init[i]);
            allMat.set(i,1,numEdges[i][0]+numEdges[i][1]+numEdges[i][2]);
           /* ccMat.set(i,1,numEdges[i][0]);
            cdMat.set(i,0,init[i]);
            cdMat.set(i,1,numEdges[i][1]);
            ddMat.set(i,0,init[i]);
            ddMat.set(i,1,numEdges[i][2]);*/
        }
        double [] limits = getLimit(allMat);
        /*double [] limitsCC = getLimit(ccMat);
        double [] limitsCD = getLimit(cdMat);
        double [] limitsDD = getLimit(ddMat);*/


        double [] realLimits =new double[init.length];
        for(int i = 0; i < init.length;i++)
        {
            realLimits[i] = map(init[i],init,limits);
        }
       /* double [][] fullResult = new double[3][init.length];
        double ogRange = init[init.length-1] - init[0];

        for(int i = 0; i < init.length;i++)
        {
            fullResult[0][i] = map(init[i],init,limitsCC);
                    fullResult[1][i] = map(init[i],init,limitsCD);
                            fullResult[2][i] = map(init[i],init,limitsDD);
        }
        TetradMatrix toPrint = new TetradMatrix(fullResult);
        System.out.println(Arrays.toString(limitsCC));
        System.out.println(Arrays.toString(limitsCD));
        System.out.println(Arrays.toString(limitsDD));
        System.out.println(Arrays.toString(init));
        System.out.println(toPrint);*/

        return realLimits;

    }

    private double map(double lambda, double [] init, double [] limit)
    {
        double ogMax = init[init.length-1];
        double ogMin = init[0];
        double ogLength = ogMax-ogMin;
        double lamLength = lambda-ogMin;
      return (lamLength/ogLength)*(limit[1]-limit[0]) + limit[0];

    }
    private double [] getLimit(TetradMatrix t) {
        ContinuousVariable v = new ContinuousVariable("Lambda");
        ContinuousVariable numEdges = new ContinuousVariable("Edges");
        List<Node> temp = new ArrayList<Node>();
        temp.add(v);
        temp.add(numEdges);
        double max = 0;
        int maxIndex = -1;
        for(int i = 1; i < t.rows()-1;i++) {
            //need to remove the points before i and regress, and the points after i and regress
            double[][] tempLow = t.toArray();
            double [][] tempHigh = t.toArray();
            double [][] realLow = new double[i][2];
            double [][] realHigh = new double[t.rows()-i-1][2];
            for(int j = 0; j < i;j++)
            {
                realLow[j][0] = tempLow[j][0];
                realLow[j][1] = tempLow[j][1];
            }
            for(int j = t.rows()-1; j > i; j--)
            {
                realHigh[j-i-1][0] = tempHigh[j][0];
                realHigh[j-i-1][1] = tempHigh[j][1];
            }

            RegressionResult result = null;
            RegressionResult result2 = null;
try {
    RegressionDataset r = new RegressionDataset(new TetradMatrix(realLow), temp);

    result = r.regress(temp.get(1), temp.get(0));

            RegressionDataset r2 = new RegressionDataset(new TetradMatrix(realHigh),temp);

          result2 = r2.regress(temp.get(1),temp.get(0));
}
catch(org.apache.commons.math3.linear.SingularMatrixException e)
{
   // System.out.println("no result for: " + i + " for original limit");
}

            if(result==null || result2 == null)
                continue;
            if(max < (result.getRSquared()+result2.getRSquared())/2)
            {
                max = (result.getRSquared()+result2.getRSquared())/2;
                maxIndex = i;
            }

            //trying to find the minimum value lambda right now
        }

        if(maxIndex==-1) {
            System.out.println("no max index");
            double [] result = new double[2];
            result[0] = t.get(0,0);
            result[1] = t.get(t.rows()-1,0);
            return result;
        }
        double[][] tempLow = t.toArray();
        double [][] tempHigh = t.toArray();
        double [][] realLow = new double[maxIndex][2];
        double [][] realHigh = new double[t.rows()-maxIndex-1][2];

        for(int j = 0; j < maxIndex;j++)
        {
            realLow[j][0] = tempLow[j][0];
            realLow[j][1] = tempLow[j][1];
        }
        for(int j = t.rows()-1; j > maxIndex; j--)
        {
            realHigh[j-maxIndex-1][0] = tempHigh[j][0];
            realHigh[j-maxIndex-1][1] = tempHigh[j][1];
        }


        double maxOne = 0;
        int indexOne = -1;
        double maxTwo = 0;
          int indexTwo = -1;
        TetradMatrix first = new TetradMatrix(realLow);
        TetradMatrix second = new TetradMatrix(realHigh);
            //test all points before it as a splitting point to see which has the best linear fit up until maxIndex
            for (int i = 1; i < first.rows(); i++) {
                double [][] tempL = first.toArray();
                double [][] tempH = first.toArray();
                double [][] realL = new double[i][2];
                double [][] realH = new double[first.rows()-i-1][2];
                for(int j = 0; j < i;j++)
                {
                    realL[j][0] = tempL[j][0];
                    realL[j][1] = tempL[j][1];
                }
                for(int j = first.rows()-1; j > i; j--)
                {
                    realH[j-i-1][0] = tempH[j][0];
                    realH[j-i-1][1] = tempH[j][1];
                }


                RegressionResult result = null;
                RegressionResult result2 = null;
                try {
                    RegressionDataset r = new RegressionDataset(new TetradMatrix(realL), temp);

                    result = r.regress(temp.get(1), temp.get(0));

                    RegressionDataset r2 = new RegressionDataset(new TetradMatrix(realH),temp);

                    result2 = r2.regress(temp.get(1),temp.get(0));
                }
                catch(org.apache.commons.math3.linear.SingularMatrixException e)
                {
                   // System.out.println("no result for: " + i + " for low limit");
                }

                if(result==null || result2 == null)
                    continue;

                if(maxOne < (result.getRSquared()+result2.getRSquared())/2)
                {
                    maxOne = (result.getRSquared()+result2.getRSquared())/2;
                    indexOne = i;
                }



        }


        if(indexOne==-1) //If the size of the data matrix is too small for one side of the partition or if no linear model fit any points
        {
            indexOne = rand.nextInt(maxIndex);
        }

        for (int i = 1; i < second.rows(); i++) {
            double[][] tempL = second.toArray();
            double[][] tempH = second.toArray();
            double[][] realL = new double[i][2];
            double[][] realH = new double[second.rows() - i - 1][2];
            for (int j = 0; j < i; j++) {
                realL[j][0] = tempL[j][0];
                realL[j][1] = tempL[j][1];
            }
            for (int j = second.rows() - 1; j > i; j--) {
                realH[j - i - 1][0] = tempH[j][0];
                realH[j - i - 1][1] = tempH[j][1];
            }

            RegressionResult result = null;
            RegressionResult result2 = null;
            try {
                RegressionDataset r = new RegressionDataset(new TetradMatrix(realL), temp);

                result = r.regress(temp.get(1), temp.get(0));

                RegressionDataset r2 = new RegressionDataset(new TetradMatrix(realH),temp);

                result2 = r2.regress(temp.get(1),temp.get(0));
            }
            catch(org.apache.commons.math3.linear.SingularMatrixException e)
            {
               // System.out.println("no result for: " + i + " for high limit");
            }

            if(result==null || result2 == null)
                continue;

            if (maxTwo < (result.getRSquared() + result2.getRSquared()) / 2) {

                maxTwo = (result.getRSquared() + result2.getRSquared()) / 2;
                indexTwo = i;
            }
        }
        if(indexTwo==-1)
        {
          indexTwo = rand.nextInt(second.rows());
        }
        double [] result = new double[2];
            result[0] = t.get(indexOne,0);
            result[1] = t.get(indexTwo + maxIndex + 1,0);
            return result;
        }
    private int getEdges(Graph g, DataSet data, String type)
    {
        int count = 0;
        for(Edge e:g.getEdges())
        {
            if(e.getNode1() instanceof DiscreteVariable  && e.getNode2() instanceof DiscreteVariable) {
                if(type.equals("DD"))
                count++;
            }
            else if (e.getNode1() instanceof ContinuousVariable && e.getNode2() instanceof ContinuousVariable) {
                if(type.equals("CC"))
                 count++;
            }
            else if(type.equals("CD"))
                count++;

        }
        return count;
    }
    private boolean[][] findPrior(TetradMatrix [] pInf)
    {
        boolean [][] temp = new boolean[pInf[0].rows()][pInf[0].columns()];
        for(int i = 0; i < pInf.length;i++)
        {

         TetradMatrix curr = pInf[i];

            for(int j = 0; j < curr.rows();j++)
            {
                for(int k = j +1; k < curr.columns();k++)
                {
                    if(curr.get(j,k)!=noPrior) {
                        temp[j][k] = true;
                        sourcePrior[j][k][i] = true;
                    }
                }
            }
        }
        return temp;
    }

    public double [][] getLambdas()
    {
        return this.lambdas;
    }
}
