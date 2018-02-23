package edu.pitt.csb.latents;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.Pref_Div.Functions;
import edu.pitt.csb.Pref_Div.Gene;
import edu.pitt.csb.mgm.*;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.SearchWrappers;
import edu.pitt.csb.stability.StabilityUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

import static edu.pitt.csb.mgm.MixedUtils.loadDataSet2;

/**
 * Created by vinee_000 on 9/11/2017.
 */
public class mirnaPrediction {
  /*  public static String runName = "0.5_TUMOR_SAMPLES";
    public static double lambdaLow = .2;
    public static int numLambdas = 30;
    public static double lambdaHigh = .7;
    public static int numSub = 10;
    public static double g = 0.05;
    public static double tao = 0.5;
    final public static double corrCutoff = 0.01;
    public static int numAlphas = 10;
    public static double alphaLow = 0.001;
    public static double alphaHigh = 0.1;
    public static double[] lastParams;
    public static boolean reuseCorrelations = false;
    public static boolean reuseFiles = false;


    public static void main(String[] args) throws Exception {
        File file = new File(runName);
        file.mkdir();
        String[] files;
        if (reuseCorrelations)
            files = new String[5];
        else
            files = new String[4];
        files[0] = "expression_cleaned.txt";
        files[1] = "Fold_Change.txt";
        files[2] = "all_gene_disease_associations.tsv";
        files[3] = "final_gene_data.txt";
        if (reuseCorrelations)
            files[4] = "correlations.txt";
        long time = System.nanoTime();
        final ArrayList<Gene> g2 = Functions.loadGeneData(".", files,null,false);
        time = System.nanoTime() - time;
        System.out.println("Time to compute correlations: " + time);
        System.out.println("Loaded Gene Data");
        //ArrayList<Gene> g2 = new ArrayList<Gene>();
        String expressionFile = "subset_expression_with_clinical.txt";
        if(runName.contains("TUMOR"))
            expressionFile = "subset_expression_tumor_samples.txt";
        else if(runName.contains("NORMAL"))
            expressionFile = "subset_expression_normal_samples.txt";
        final PrintStream temp = new PrintStream(runName + "/log.txt");

        DataSet d = MixedUtils.loadDataSet2(expressionFile);
        /*will comment this out when adding Pref-Div to the pipeline
       final ArrayList<Gene> pdSelected = new ArrayList<Gene>();
       for(Gene gen: g2)
       {
           if(d.getVariable(gen.symbol)!=null)
               pdSelected.add(gen);
       }

        d.removeColumn(d.getVariable("patient.gender"));
        final List<Node> discVars = MixedUtils.getDiscreteData(d).getVariables();
        System.out.println("Loaded Dataset");
        DoubleMatrix2D stabs;
        if (!reuseFiles) {
            stabs = runAlgorithms(d, temp);
            PrintStream out = new PrintStream(runName + "/stability_main.txt");
            for (int i = 0; i < d.getNumColumns(); i++) {
                out.print(d.getVariable(i).getName() + "\t");
            }
            out.println();
            for (int i = 0; i < stabs.rows(); i++) {
                for (int j = 0; j < stabs.columns(); j++) {
                    out.print(stabs.get(i, j) + "\t");
                }
                out.println();
            }
            out.flush();
            out.close();
        } else {
            BufferedReader b = new BufferedReader(new FileReader(runName + "/stability_main.txt"));
            String[] vars = b.readLine().split("\t");
            double[][] temp2 = new double[vars.length][vars.length];
            for (int i = 0; i < vars.length; i++) {
                String[] line = b.readLine().split("\t");
                for (int j = 0; j < vars.length; j++) {
                    temp2[i][j] = Double.parseDouble(line[j]);
                }
            }
            stabs = new DenseDoubleMatrix2D(temp2);
        }
        final ArrayList<Pair> latents = new ArrayList<Pair>();

        if (reuseFiles) {
            BufferedReader b = new BufferedReader(new FileReader(runName + "/latents_main.txt"));
            while (b.ready()) {

                String[] line = b.readLine().split("\t");
                if (stabs.get(d.getColumn(d.getVariable(line[0])), d.getColumn(d.getVariable(line[1]))) > tao)
                    latents.add(new Pair(d.getVariable(line[0]), d.getVariable(line[1])));
            }
        } else {
            PrintStream stable = new PrintStream(runName + "/latents_main.txt");
            for (int i = 0; i < d.getNumColumns(); i++) {
                for (int j = i + 1; j < d.getNumColumns(); j++) {
                    if (stabs.get(i, j) > tao) {
                        stable.println(d.getVariable(i) + "\t" + d.getVariable(j) + "\t" + stabs.get(i, j));
                        latents.add(new Pair(d.getVariable(i), d.getVariable(j)));
                        //Are we only focused on genes that are downregulated by miRNAs? There is evidence for upregulation as well
                    }
                }
            }
            stable.flush();
            stable.close();
        }
        Functions.indT.setAlpha(corrCutoff);
        System.out.println("Loaded Latents");

        final PrintStream output = new PrintStream(runName +"/validated_latents.txt");
        final ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

        class StabilityAction extends RecursiveAction {
            private int chunk;
            private int from;
            private int to;

            public StabilityAction(int chunk, int from, int to) {
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }

            private synchronized void writeGenes(PrintStream temp, Gene one, Gene two, ArrayList<Gene> genes) {
                temp.println("Related Gene Set To: " + one.symbol + "," + two.symbol + "," + genes);
            }

            private synchronized boolean getCorrelation(Gene one, Gene two) {
                return Functions.computePValue(one, two);
            }

            private synchronized DataSet getTempData(String filename, List<Node> discVars, ArrayList<Gene> genes, Gene one, Gene two) {
                writeTempExpression(filename, discVars, genes, one, two);
                try {
                    return MixedUtils.loadDataSet2("temp_expression.txt");
                } catch (Exception e) {
                    return null;
                }
            }

            @Override
            protected void compute() {
                if (to - from <= chunk) {
                    for (int s = from; s < to; s++) {
                        Pair p = latents.get(s);
                        Gene one = null;
                        Gene two = null;
                        for (Gene temp2 : g2) {
                            if (p.one.getName().equals(temp2.symbol))
                                one = temp2;
                            if (p.two.getName().equals(temp2.symbol))
                                two = temp2;

                        }
                        assert (one != null);
                        assert (two != null);
                        //  Queue q = new Queue(validationSize);
                        ArrayList<Gene> genes = new ArrayList<Gene>();
                        for (Gene gen : g2) {

                            if (one.symbol.equals(gen.symbol) || two.symbol.equals(gen.symbol))
                                continue;
                /*double corr1 = Functions.computePValue(one,gen);
                if(corr1>=corrCutoff)
                    continue;
                double corr2 = Functions.computePValue(two,gen);

                            if (Functions.computePValue(one, gen) || Functions.computePValue(two, gen)) //computePValue returns true if they are independnet
                                continue;
                            if (!Functions.getPotentialBlocker(one, two, gen))
                                continue;
                            genes.add(gen);
                        }

                        writeGenes(temp, one, two, genes);

                        if (genes.size() != 0) {
                            for(Gene x: pdSelected)
                            {
                                if(!genes.contains(x))
                                    genes.add(x);
                            }
                            DataSet d = getTempData("expression_with_clinical.txt", discVars, genes, one, two);
                            temp.println("Running algorithms to validate: " + one.symbol + "," + two.symbol);
                            B:while(true) {
                                try {
                                    DoubleMatrix2D stabs2 = runAlgorithms(d, temp);
                                    if (stabs2.get(d.getColumn(d.getVariable(one.symbol)), d.getColumn(d.getVariable(two.symbol))) > tao) {
                                        temp.println("Validated: " + one.symbol + "," + two.symbol);
                                        output.println(one.symbol + "\t" + two.symbol + "\t" + stabs2.get(d.getColumn(d.getVariable(one.symbol)), d.getColumn(d.getVariable(two.symbol))));
                                    } else {
                                        temp.println("Could not validate: " + one.symbol + "," + two.symbol);
                                    }
                                    break B;
                                } catch (Exception e) {

                                }
                            }

                        } else {
                            output.println(one.symbol + "\t" + two.symbol + "\tNo potential blocking genes found");
                        }

                        //TODO next need to determine theoretically which genes from the original study were used to identify the bidirected edge and maintain those
                        //TODO computational side: profile the correctness of each of the different orientation steps

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


        final int chunk = 1000;

        pool.invoke(new StabilityAction(chunk, 0, latents.size()));
//TODO include neighbors of latent pair in the expression dataset along with the top 200 genes from the correlation
        //TODO next need to determine theoretically which genes from the original study were used to identify the bidirected edge and maintain those
        //TODO computational side: profile the correctness of each of the different orientation steps
            /*if(false) {
                writeTempExpression("expression_with_clinical.txt", discVars, genes, one, two);
                d = MixedUtils.loadDataSet2("temp_expression.txt");
                temp.println("Running algorithms to validate: " + one.symbol + "," + two.symbol);
                stabs = runAlgorithms(d, temp);
                if (stabs.get(d.getColumn(d.getVariable(one.symbol)), d.getColumn(d.getVariable(two.symbol))) > tao) {
                    temp.println("Validated: " + one.symbol + "," + two.symbol);
                    output.println(one.symbol + "\t" + two.symbol + "\t" + stabs.get(d.getColumn(d.getVariable(one.symbol)), d.getColumn(d.getVariable(two.symbol))));
                } else {
                    temp.println("Didn't Validate: " + one.symbol + "," + two.symbol + "," + stabs.get(d.getColumn(d.getVariable(one.symbol)), d.getColumn(d.getVariable(two.symbol))));
                    output.println("Didn't Validate: " + one.symbol + "," + two.symbol + "," + stabs.get(d.getColumn(d.getVariable(one.symbol)), d.getColumn(d.getVariable(two.symbol))));
                    temp.println("Checking for correlation issues...");
                    //Check if there is high corrleation issues
                    double[] tempLambda = {lastParams[0], lastParams[1], lastParams[2]};
                    MGM m = new MGM(d, tempLambda);
                    m.learnEdges(1000);
                    IndependenceTest indT = new IndTestMultinomialAJ(d, lastParams[3]);
                    FciMaxP f = new FciMaxP(indT);
                    f.setInitialGraph(m.graphFromMGM());
                    Graph gTemp = f.search();
                    temp.println("Learned Causal Graph: " + gTemp);
                    Node varOne = gTemp.getNode(one.symbol);
                    List<Node> list = gTemp.getAdjacentNodes(varOne);
                    if (isMaximalClique(gTemp, list)) {
                        temp.println("Found Cluster for Node One: " + list);
                    }
                    Node varTwo = gTemp.getNode(two.symbol);
                    list = gTemp.getAdjacentNodes(varTwo);
                    if (isMaximalClique(gTemp, list)) {
                        temp.println("Found Cluster for Node Two: " + list);
                    }

                }
            }
        output.flush();


            /*Find top 200 correlated genes based on average absolute value of correlation to genes
            Re run procedure on the dataset from these genes only
            Determine if pair p remains bidirected connected
            If not, determine if either variable in pair is disconnected from the rest of the network
            Cluster these representatives together and repeat the procedure

        /*Run STEPS to get MGM parameter estimates
        Run STARS with MFM to get alpha estimate
        Run Stability search with the selected alpha and lambda estimates to get most stable bidirected edges
        Go to the correlation matrix computed from the Functions class to find other variables that should be used with this set of bidirected ones
        Repeat the procedure to find new bidirected edges, if the original remain stable then include these in the final

    }
    //[0,1,2]
    public synchronized static void writeTempExpression(String expFile,List<Node> discVars, List<Gene> genes, Gene one, Gene two)
    {
        try {
            PrintStream out = new PrintStream("temp_expression.txt");
            BufferedReader b = new BufferedReader(new FileReader(expFile));
            int [] indices = new int[discVars.size() + genes.size() + 2];
            String [] line = b.readLine().split("\t");
            int count = 0;
           A: for(int i = 0; i < line.length;i++)
            {
                for(Node n : discVars)
                {
                    if(line[i].equals(n.getName()))
                    {
                        out.print(n.getName() + "\t");
                        indices[count] = i;
                        count++;
                        continue A;
                    }
                }
                for(Gene g: genes)
                {
                    if(line[i].equals(g.symbol))
                    {
                        out.print(g.symbol + "\t");
                        indices[count] = i;
                        count++;
                        continue A;
                    }
                }
                if(line[i].equals(one.symbol) || line[i].equals(two.symbol))
                {
                    if(line[i].equals(one.symbol))
                        out.print(one.symbol + "\t");
                    else
                        out.print(two.symbol + "\t");
                    indices[count] = i;
                    count++;
                    continue A;
                }

            }
            assert(count == indices.length);
           out.println();
           while(b.ready())
           {
               String [] line2 = b.readLine().split("\t");
               int count2 = 0;
               for(int i = 0; i < line2.length;i++)
               {
                   if(indices[count2]==i) {
                       out.print(line2[i] + "\t");
                       count2++;
                   }
               }
               out.println();
           }
           out.flush();
           out.close();
           b.close();
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }

    }
    public static boolean isMaximalClique(Graph gTemp, List<Node> list)
    {
        for(Node x:list)
        {
            List<Node> temp = gTemp.getAdjacentNodes(x);
            List<Node> copy = new ArrayList<Node>();
            copy.addAll(list);
            for(Node m: temp) {
                if (!list.contains(m)) {
                    return false;
                }
                else
                    copy.remove(m);
            }
            copy.remove(x);
                if(!copy.isEmpty())
                    return false;

        }
        return true;
    }
    public static DoubleMatrix2D runAlgorithms(DataSet d, PrintStream temp)
    {
       double[] initLambdas = new double[numLambdas];
        for (int i = 0; i < numLambdas; i++) {
            initLambdas[i] = i * (lambdaHigh - lambdaLow) / numLambdas + lambdaLow;
        }

        STEPS s = new STEPS(d, initLambdas, g, numSub);
        s.runStepsArrayPar();
        double[] lambda = s.lastLambda;
        temp.print("STEPS Lambda:");
        temp.println(Arrays.toString(s.lastLambda));


        double[] initAlphas = new double[numAlphas];
        for (int i = 0; i < numAlphas; i++) {
            initAlphas[i] = i * (alphaHigh - alphaLow) / numAlphas + alphaLow;
        }
        STARS s2 = new STARS(d, initAlphas, g, numSub);
        s2.runSTARS();
        temp.print("STARS Alpha:");
        temp.println(s2.lastAlpha);
        double[] params = {lambda[0], lambda[1], lambda[2], s2.lastAlpha};
        lastParams = params;
        DataGraphSearch gs = new SearchWrappers.MFMWrapper(params);
        DoubleMatrix2D stabs = StabilityUtils.StabilitySearchParLatent(d, gs, 10, (int) (10 * Math.sqrt(d.getNumRows())));

        return stabs;
    }
    private static class Pair
    {
        public Node one;
        public Node two;
        public Pair(Node x, Node y)
        {
            one = x;
            two = y;
        }
    }
    private static class Queue
    {
        public ArrayList<Gene> genes;
        public ArrayList<Double> correlation;
        private int size;
        public Queue(int size)
        {
            genes = new ArrayList<Gene>(size);
            correlation = new ArrayList<Double>(size);
            this.size = size;
        }

        public boolean add(Gene temp,double corr)
        {
            if(correlation.size() < size)
            {
                A:for(int i = 0; i < correlation.size();i++)
                {
                    if(correlation.get(i) < corr)
                    {
                        correlation.add(i,corr);
                        genes.add(i,temp);
                        if(correlation.size() > size) {
                            correlation.remove(size);
                            genes.remove(size);
                        }
                        return true;
                    }
                }
                correlation.add(correlation.size(),corr);
                genes.add(genes.size(),temp);
                if(correlation.size() > size)
                {
                    correlation.remove(size);
                    genes.remove(size);
                }
                return true;
            }
            else if(correlation.get(correlation.size()-1) > corr)
                return false;
            else
            {
               B: for(int i = 0; i < correlation.size();i++)
                {
                    if(correlation.get(i) < corr)
                    {
                        correlation.add(i,corr);
                        genes.add(i,temp);
                        correlation.remove(size);
                        genes.remove(size);
                        break B;
                    }
                }
                return true;
            }
        }
    }*/
}
