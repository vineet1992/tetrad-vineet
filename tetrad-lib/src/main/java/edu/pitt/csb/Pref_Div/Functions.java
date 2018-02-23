package edu.pitt.csb.Pref_Div;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.search.IndTestPartialCorrelation;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradVector;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.io.*;
public class Functions
{
    private static final long [] CHROMOSOME_LENGTHS = {249250621,243199373,198022430,191154276,180915260,171115067,159138663,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415};//Chromosomal lengths in base pairs

  //TODO Compute this on the fly instead of offline
   private static final double MEDIAN_DISEASE_SCORE = 0.000543;
    public static IndependenceTest indT;
    private static float [] corrMat;
    private static int n;
    private static double[][] data;
    private static double maxDiseaseScore;
    private static double meanIntensity; //Used to place z-scored theory information onto the intensity distribution
    private static double stdIntensity; //Same as above
    private static PrintStream out;
    private static NumberFormat nf;
    private static int theorySources;
    private static int maxIntDiscrete = 5;

    //Computes Intensity Score for a gene g1
    //Needs parameter a which specifies tradeoff between using foldChange in Expression value versus gene-disease mapping
    //Needs a list of diseases of interest to gather a score based upon gene-disease relationship


    //Computes disease score as average over all relevant diseases
    public static double dIntensityScore(HashMap<Integer,Double> score,int [] diseaseID)
    {
        double disScore = 0;
        for(int i= 0; i < diseaseID.length;i++)
        {
            if(score.get(diseaseID[i])!=null)
                disScore+= score.get(diseaseID[i]);
        }
        disScore=disScore/diseaseID.length;
        return disScore;
    }

    //Given a gene object, and a list of diseases of interest, and a tradeoff parameter alpha this computes the intensity score for a given gene
    //Currently, the function normalizes the fold change value based on the maxFoldChange, and the maxDiseaseIntensity score
    public static double computeIntensity(Gene g1, double a,int [] disease)
    {
        double disScore = 0;
        if(g1.diseaseScores!=null)
        {
            for(int i = 0; i < disease.length;i++)
            {
                if(g1.diseaseScores.get(disease[i])!=null)
                    disScore+=g1.diseaseScores.get(disease[i]);

            }

            disScore=disScore/disease.length;
        }
        if(a < 0 || a > 1)
            return Double.MIN_VALUE;
        return a*Math.abs(g1.foldChange)+(1-a)*disScore;
    }
    //Computes the Distance between two genes g1 and g2
//Parameter a controls weight of gene expression correlation
//Parameter b controls weight of Disease similarity score
//Parameter c controls weight of Chromosomal Distance score
    public static double computeDistance(Gene g1, Gene g2, double a, double b,double c)
    {
        if((a+b+c > 1) || a < 0 || b < 0 || c < 0)
            return Double.MIN_VALUE;
        //a*correlation + b*diseaseScoreSim + (1-a-b)*chromosomal distance + gene family indicator
        double DS = computeScore(g1.diseaseScores,g2.diseaseScores);
        double CD = chromDist(g1,g2);
        double family = sameFamily(g1,g2);
        int i = 0; int j = 0;
        if(g1.ID > g2.ID) {
            j = g1.ID;
            i = g2.ID;
        }
        else
        {
            i = g1.ID;
            j = g2.ID;
        }
        int index = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
      //  return a*Math.abs(corrMat[index])+b*DS+c*CD + (1-(a+b+c))*family;
        return a*Math.abs(computeCorrelation(i,j))+b*DS+c*CD + (1-(a+b+c))*family;
    }

    //Returns the correlation between two genes, given that their information is loaded
    public static double getCorrelation(Gene g1, Gene g2)
    {
        int i = -1;
        int j = -1;
        if(g1.ID > g2.ID)
        {
            i = g1.ID;
            j = g2.ID;
        }
        else
        {
            i = g2.ID;
            j = g1.ID;
        }
        int index = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
       // return corrMat[index];
        return computeCorrelation(i,j);
    }

    //Computes array with each component of the distance score between two genes
    public static double[] computeDistance(Gene g1, Gene g2)
    {
        //a*correlation + b*diseaseScoreSim + (1-a-b)*chromosomal distance + gene family indicator
        double DS = computeScore(g1.diseaseScores,g2.diseaseScores);
        double CD = chromDist(g1,g2);
        double family = sameFamily(g1,g2);
        int i = 0; int j = 0;
        if(g1.ID > g2.ID) {
            j = g1.ID;
            i = g2.ID;
        }
        else
        {
            i = g1.ID;
            j = g2.ID;
        }
        int index = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
        double [] temp = {DS,CD,family,computeCorrelation(i,j)};
        out.println(nf.format(corrMat[index]) + "\t" + nf.format(DS) + "\t" + nf.format(CD) + "\t" + nf.format(family));
        //System.out.println("Corr: " + Math.abs(corrMat[g1.ID][g2.ID]) + ",DS:" + DS + ",CD:" + CD + ",FM:" + family);
        return temp;
    }
    //Computes Disease Similarity Score (now using cosine similarity)
    private static double computeScore(HashMap<Integer,Double> score1, HashMap<Integer,Double> score2)
    {
         if(score1==null && score2==null)
            return -1;

        ArrayList<Double> vecOne = new ArrayList<Double>();
        ArrayList<Double> vecTwo = new ArrayList<Double>();
        if(score1!=null) {
            for (Integer x : score1.keySet()) {
                if (score2!=null && score2.get(x) != null)
                {
                    vecOne.add(score1.get(x));
                    vecTwo.add(score2.get(x));
                }
                else {
                    vecOne.add(score1.get(x));
                    vecTwo.add(0.0);
                }
            }
        }
        if(score2!=null) {
            for (Integer y : score2.keySet()) {
                if (score1!=null && score1.get(y) != null) {
                    vecOne.add(score1.get(y));
                    vecTwo.add(score2.get(y));
                }
                else
                {
                    vecTwo.add(score2.get(y));
                    vecOne.add(0.0);
                }
            }
        }
        return cosineSim(vecOne,vecTwo);
    }


    private static double cosineSim(ArrayList<Double>one,ArrayList<Double>two)
    {

        if(one.size()!=two.size())
        {
            System.out.println("Cosine Similarity have different sizes");
            System.exit(-1);
        }
        double asq = 0;
        double bsq = 0;
        double prod = 0;
        for(int i = 0; i < one.size();i++)
        {
            asq+=Math.pow(one.get(i),2);
            bsq+=Math.pow(two.get(i),2);
            prod+=one.get(i)*two.get(i);
        }
        return prod/(Math.sqrt(asq)*Math.sqrt(bsq));
    }

    //Computes Chromosomal Distance score
    private static double chromDist(Gene g1, Gene g2)
    {
        if(g1.chromosome==-1 || g2.chromosome ==-1)
            return -1;//Be indifferent in the abscence of information
        if(g1.chromosome!=g2.chromosome)
            return -1;
        else
        {
            if(g1.chromosomeStart==-1 || g2.chromosomeStart==-1 || g1.chromosomeStop==-1 || g2.chromosomeStop==-1)
                return -1;
            if(((max(g1.chromosomeStop,g2.chromosomeStop)-min(g1.chromosomeStart,g2.chromosomeStart))/(CHROMOSOME_LENGTHS[g1.chromosome-1])) > 1)
            {
                //System.out.println(g1.chromosome);
            }
            return (max(g1.chromosomeStop,g2.chromosomeStop)-min(g1.chromosomeStart,g2.chromosomeStart))/CHROMOSOME_LENGTHS[g1.chromosome-1];
        }

    }

    //Max and Min Helper functions
    private static double max(long a, long b)
    {
        return (a > b) ? a : b;
    }
    private static double min(long a, long b)
    {
        return (a < b) ? a : b;
    }



    //Returns gene family score
    private static double sameFamily(Gene g1, Gene g2)
    {
        if(g1.familyID==null || g2.familyID==null)
            return -1;//We don't have family ID information for this gene so be indifferent
        for(int i = 0; i < g1.familyID.length;i++)
        {
            for(int j = 0; j < g2.familyID.length;j++)
            {
                if(g1.familyID[i]==g2.familyID[j])
                    return 0;
            }
        }
        return 1;
    }



    //Main function to load all gene data into an Arraylist and stores correlation matrix locally in the Functions class
    public static ArrayList<Gene> loadGeneData(String dir,String [] files,int[]diseaseID,boolean filter) throws Exception
    {


     /*   String expFile = dir +  "/TCGA_BRCA_FILES/expression_transposed_npn.txt";
        String fcFile = dir +  "/TCGA_BRCA_FILES/foldChange_subtype.txt";
        String gdFile = dir + "/all_gene_disease_associations.tsv";
        String giFile = dir + "/final_gene_data.txt";*/

        String expFile = dir +  "/" + files[0]; //expression file
        String fcFile = dir +  "/" + files[1]; //fold change file
        String gdFile = dir + "/" + files[2]; //gene-disease associations
        String giFile = dir + "/" + files[3]; //other genetic data (chromosome/gene family)


        /******************END AVERAGE DISEASE CODE*********************************************/
        //Constant indices of each piece of information from the Gene Information file//
        int familyInd = 2;
        int chromInd = 3;
        int symId = 0;
        /*******************************************************************************/



        /******************************LOAD IN GENE FAMILY AND CHROMOSOME LENGTH INFORMATION***************************/

        System.out.print("Loading Gene Family and Chromosome Info...");

        ArrayList<Gene> g = new ArrayList<Gene>();
        BufferedReader b = new BufferedReader(new FileReader(giFile));
        b.readLine();

        //Mapping from gene symbol to array of families
        HashMap<String,int[]> families = new HashMap<String,int[]>();
        //Mapping from gene symbol to String (Chrom_ID\tChrom_Start_Location\tChrom_End_Location)
        HashMap<String,String> chromosomeStuff = new HashMap<String,String>();
        int numSamples2 = 0;
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            if(line.length < familyInd + 1)
                continue;
            if(line[familyInd].length()<1)
                families.put(line[symId],null);
            else
            {
                String[] fam = line[familyInd].split("\\|");
                int [] famInt = new int[fam.length];
                for(int i = 0; i < fam.length;i++)
                    famInt[i] = Integer.parseInt(fam[i]);
                families.put(line[symId],famInt);
            }
            if(line.length < chromInd+1)
                continue;
            if(line[chromInd].equals("X"))
                line[chromInd] = "23";
            else if (line[chromInd].equals("Y"))
                line[chromInd] = "24";
            chromosomeStuff.put(line[symId],line[chromInd] + "\t" + line[chromInd+1] + "\t" + line[chromInd+2]);

        }
        b.close();

        System.out.println("Done");

        /*********************************************END GENE INFO LOADING**************************************/


        /********************************************READ GENE OUTCOME INTENSITY*********************************/
        //Read mutual information file between genes and specific outcome of interest
        //No normalization needed in the new scheme

        System.out.print("Loading Gene-Outcome Intensity Scores...");
        b = new BufferedReader(new FileReader(fcFile));
        ArrayList<Double> miValues = new ArrayList<Double>();
        int count = 0;
        HashMap<String,Double> fcs = new HashMap<String,Double>();
        while(b.ready())
        {
            String [] fol = b.readLine().split("\t");
            String name = fol[0];

            double mi = Double.parseDouble(fol[1]);
            fcs.put(name,mi);
            miValues.add(mi);
            count++;
        }
        double [] miArr = new double[count];
        for(int i = 0; i < count;i++)
            miArr[i] = miValues.get(i);
        meanIntensity = StatUtils.mean(miArr);
        stdIntensity = StatUtils.sd(miArr);
        b.close();
        System.out.println("Done");
        /*************************************************FINISHED READING***************************************/



        /**********************************************GENE DISEASE INTENSITY***********************************/

        System.out.print("Loading Gene Disease Intensity Scores...");
        HashMap<String,HashMap<Integer,Double>> geneDisease = new HashMap<String,HashMap<Integer,Double>>();


        //Put gene disease data into the final hashmap
        b = new BufferedReader(new FileReader(gdFile));
        boolean ready = false;
        while(b.ready())
        {
            if(ready)
            {
                String [] line = b.readLine().split("\t");
                if(geneDisease.get(line[1])==null)
                {
                    HashMap<Integer,Double> temp2 = new HashMap<Integer,Double>();
                    int disID = Integer.parseInt(line[3].split("C")[1]);
                    temp2.put(disID,Double.parseDouble(line[5]));
                    geneDisease.put(line[1],temp2);
                }
                else
                {
                    HashMap<Integer,Double> temp2 = geneDisease.get(line[1]);
                    int disID = Integer.parseInt(line[3].split("C")[1]);
                    temp2.put(disID,Double.parseDouble(line[5]));
                    geneDisease.put(line[1],temp2);
                }
            }
            if(b.readLine().startsWith("gene"))
                ready = true;


        }


        b.close();

        System.out.println("Done");

        /*******************************************************END GENE DISEASE MAPPING***********************************/



        /******************************************************GENE EXPRESSION LOADING AND PARSING*****************************/
        System.out.print("Loading Gene Expression Data...");
        b = new BufferedReader(new FileReader(expFile));
        b.readLine();
        int idCount = 0;

        //Loop through the expression file, rows are genes, columns are samples
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            numSamples2 = line.length-1;
            Gene curr = new Gene(idCount);

            //Get the current gene symbol and give it an ID
            String symbol = line[0];
            curr.symbol = symbol;

            //Give it a fold change value from our previous HashMap
            if(fcs.get(symbol)!=null)
            {
                curr.foldChange = fcs.get(symbol);
             //   fcs.remove(symbol);
            }
            else
            {
                System.out.println("No Mututal Information score for Gene:" + symbol);
                curr.foldChange = 0;
            }

            //Give it a gene Family ID from previous HashMap
            if(families.get(symbol)!=null)
            {
                curr.familyID = families.get(symbol);
                families.remove(symbol);
            }
            else
            {
                curr.familyID = null;
            }

            //Give it gene-disease mappings from HashMap
            if(geneDisease.get(symbol)!=null)
            {
                curr.diseaseScores = geneDisease.get(symbol);
                geneDisease.remove(symbol);
            }
            else
            {
                curr.diseaseScores = null;
            }

            //Give it Chromosomal Location data from previous HashMap
            if(chromosomeStuff.get(symbol)!=null)
            {
                String[] chromo = chromosomeStuff.get(symbol).split("\t");
                curr.chromosome = Integer.parseInt(chromo[0]);
                curr.chromosomeStart = Long.parseLong(chromo[1]);
                curr.chromosomeStop = Long.parseLong(chromo[2]);

                chromosomeStuff.remove(symbol);
            }
            else
            {
                curr.chromosome = -1;
            }
            idCount++;
            g.add(curr);
        }
        chromosomeStuff = null;
        geneDisease = null;
        //fcs = null;
        families = null;
        b = new BufferedReader(new FileReader(expFile));
        b.readLine();
        double [][] dat = new double[g.size()][numSamples2];
        System.out.println("Num Genes: " + g.size());
        System.out.println("Number Samples: " + numSamples2);
        int currIdx = 0;


        //Load in all expression data
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            for(int i = 1; i < line.length;i++)
            {
                dat[currIdx][i-1] = Double.parseDouble(line[i]);
            }
            currIdx++;
        }
        b.close();


        System.out.println("Done");
        System.out.print("Computing Correlations...");
        n = dat.length;

        corrMat = new float[(dat.length*(dat.length-1))/2];
       // indep = new byte[(dat.length*(dat.length-1))/2];


        //If we're using Pref-Div just to filter out genes to get a set for analysis, don't need to construct independence test object
        //This is only for latent variable prediction
        nf = new DecimalFormat("0.####");
        data = dat;
         out = new PrintStream("gene_disease_distribution.txt");
         PrintStream geneNames = new PrintStream("all_genes.txt");
         /*for(int i = 0; i < numSamples;i++)
        {
            Gene gene1 = g.get(rand.nextInt(g.size()));
            Gene gene2 = g.get(rand.nextInt(g.size()));
            computeDistance(gene1,gene2);
        }*/




        float [] dataIntensity = new float[g.size()];
        float [] dataTheory = new float[g.size()]; //TODO Generalize this to N user defined theory sources for separate theory function
        for(int i = 0; i < g.size();i++)
        {
            dataIntensity[i] = (float)g.get(i).foldChange;
            double disScore = 0;
            if(g.get(i).diseaseScores!=null)
            {
                for(int j = 0; j < diseaseID.length;j++)
                {
                    if(g.get(i).diseaseScores.get(diseaseID[j])!=null)
                        disScore+=g.get(i).diseaseScores.get(diseaseID[j]);
                }

                disScore=disScore/diseaseID.length;
            }
            else
                disScore = -1;
            dataTheory[i]  = (float)disScore;

        }
        dataIntensity = NPNIgnore(dataIntensity,-1);

        dataTheory = NPNIgnore(dataTheory,-1);

        double intenseMean = mean(dataIntensity);
        double intenseSD = sd(dataIntensity,intenseMean);

        for(int i = 0; i < g.size();i++)
        {
            if(dataTheory[i]!=-1.0)
            {
                double realTheory = dataTheory[i]*intenseSD+intenseMean;
                geneNames.println(g.get(i).symbol + "\t" + nf.format(dataIntensity[i]) + "\t" + nf.format(realTheory));

            }
            else
            {
                geneNames.println(g.get(i).symbol + "\t" + nf.format(dataIntensity[i]) + "\t" + nf.format(dataTheory[i]));
            }
            for(int j = i+1; j < g.size();j++)
            {
                computeDistance(g.get(i),g.get(j));
            }
        }

        corrMat = null;
        dataIntensity = null;
        dataTheory = null;


        System.out.println("Done");
        System.out.print("Loading STRING data...");
        b = new BufferedReader(new FileReader(dir+"/"+files[4]));
        int stringSources = 0;
        HashMap<String,int[]> str = new HashMap<String,int[]>();
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            if(fcs.get(line[0])!=null && fcs.get(line[1])!=null) {


                int[] temp = new int[line.length - 2];
                for (int i = 2; i < line.length; i++) {

                    if(line[i].equals("NA"))
                        temp[i-2] = -1;
                    else
                        temp[i-2] = Integer.parseInt(line[i]);

                }
                stringSources = line.length-2;
                str.put(line[0]+","+line[1],temp);
            }
        }
        b.close();

        System.out.println("Done");
        fcs = null;

        System.out.print("Creating all sources dataset...");
        b = new BufferedReader(new FileReader("gene_disease_distribution.txt"));
        PrintStream exp = new PrintStream("all_dissimilarity_sources.txt"); //Data Driven dissimialrity is the first column, 2 through N columns are the theory-driven measurements
        for(int i = 0; i < g.size();i++)
        {
            for(int j = i+1;j < g.size();j++)
            {
                String line = b.readLine();
                exp.print(line + "\t");
                int [] temp = str.get(g.get(i).symbol+","+g.get(j).symbol);
                if(temp==null)
                {
                    for(int k = 0; k < stringSources;k++)
                    {
                        exp.print("-1\t");
                    }
                }
                else
                {
                    for(int k = 0; k < stringSources;k++)
                        exp.print(nf.format(temp[k])+"\t");
                }
                exp.println();
            }
        }
        System.out.println("Done");
        str = null;

        geneNames.flush();
        geneNames.close();
        out.flush();
        out.close();
        return g;
    }

    public static float [] getTheoryMatrix(String theoryFile)
    {
        try{
            BufferedReader b = new BufferedReader(new FileReader(theoryFile));
            String [] line = b.readLine().split("\t");
            theorySources = line.length;
            b.close();
        }
        catch(Exception e)
        {
            System.err.println("Couldn't open file to count number of sources");
        }
        int numGenes = 8098; //TODO load numGenes using other function
        int length = (numGenes*(numGenes-1))/2;
        for(int i = 0; i < theorySources;i++)
        try {
            System.out.println("Normalizing theory source " + i);
            float [] temp = new float[length];
            boolean [] keep = new boolean[length];
            int total = 0;
            BufferedReader b = new BufferedReader(new FileReader(theoryFile));
            PrintStream out = new PrintStream("_temp_" + i + ".txt");
            for(int j = 0; j < length;j++)
            {
                String line = b.readLine();
                String x = line.split("\t")[i];
                    temp[j] = Float.parseFloat(x);
            }
            float [] npn = NPNIgnore(temp,-1);
            long time = System.nanoTime();
            System.out.print("Printing Out normalized data to file...");
            if(npn!=null)
            {
                for(int j = 0; j< npn.length;j++)
                {
                    out.println(npn[j]);
                }
            }
            else
            {
                for(int j = 0; j < temp.length;j++)
                {
                    out.println(temp[j]);
                }
            }
            time = System.nanoTime()-time;
            System.out.println("Done: " + time/Math.pow(10,9));
            out.flush();
            out.close();
        }
        catch(Exception e)
        {
            e.printStackTrace();
            System.err.println("Unable to read in theory file");
            return null;
        }

        //Assume unique values greater than 5 = continuous columns
        //Input: File name for the all theory data file,
        // non-paranormal transform each NON-DISCRETE column,
        // get mean z-score for each row,
        // map it to the normal distribution given by the distribution of correlations
        //return this vector
        return null;
    }


    private static double mean(float[] t)
    {
        double mean = 0;
        for(int i = 0; i < t.length;i++)
            mean+=t[i];
        return mean/t.length;
    }

    private static double sd(float [] t, double mean)
    {
        double sd = 0;
        for(int i = 0; i< t.length;i++)
        {
            sd += Math.pow(t[i]-mean,2);
        }
        return Math.sqrt(sd/(t.length-1));
    }

    private static float[] ranks(float[] x) {

        TreeMap<Float,Integer> temp = new TreeMap<Float,Integer>();
        for(int i = 0; i < x.length;i++)
        {
            if(temp.get(x[i])!=null)
                temp.put(x[i],temp.get(x[i])+1);
            else
                temp.put(x[i],1);
        }
        if(temp.keySet().size() <= maxIntDiscrete)
            return null;
        HashMap<Float,Float>ranks = new HashMap<Float,Float>();
        int total = 1;
        for(Float f:temp.navigableKeySet())
        {
            //put the average of total+count-1 and total
            int count = temp.get(f);
            float avg = (float)((total+count-1+total)/2.0);
            ranks.put(f,avg);
            total+=count;
        }

        /*ArrayIndexComparator c = new ArrayIndexComparator(x);
        Integer [] ranks = c.createIndexArray();
        Arrays.sort(ranks,c);*/


        for(int i = 0; i < x.length;i++) {
            x[i] = ranks.get(x[i]);

        }

        return x;
    }


    //TODO Needs debugged! Data Intensity Score looks whack
    //Runs NPN on t, ignoring those elements of t that equal val
    private static float [] NPNIgnore(float [] t, float val)
    {
        boolean [] keep = new boolean[t.length];
        int total = 0;
        for(int j = 0; j < t.length;j++)
        {
            if(t[j]!=val) {
                keep[j] = true;
                total++;
            }
        }
        float [] npn = new float[total];
        total = 0;
        for(int j =0; j < t.length;j++)
        {
            if(keep[j])
            {
                npn[total] = t[j];
                total++;
            }
        }
         npn = NPN(npn);
        if(npn!=null)
        {
            total = 0;
            for (int j = 0; j < t.length; j++) {
                if (keep[j]) {
                    t[j] = npn[total];
                    total++;
                }
            }
        }
        return t;
    }
    private static float [] NPN(float [] t)
    {

         double n = t.length;
        final double delta = 1.0 / (4.0 * Math.pow(n, 0.25) * Math.sqrt(Math.PI * Math.log(n)));
        final NormalDistribution normalDistribution = new NormalDistribution();


            System.out.print("Computing ranks...");
            long time = System.nanoTime();
            t = ranks(t);
            if(t==null)
                return null;
            time = System.nanoTime()-time;
            System.out.println("Done: " + time/Math.pow(10,9));


        System.out.print("Replacing ranks with normal distribution ICP...");
            time = System.nanoTime();
            for (int i = 0; i < t.length; i++) {
                t[i] /= n;
                if (t[i] < delta) t[i] = (float)delta;
                if (t[i] > (1. - delta)) t[i] = (float)(1. - delta);
                t[i] = (float)normalDistribution.inverseCumulativeProbability((double)t[i]);
            }


        time = System.nanoTime()-time;
        System.out.println("Done: " + time/Math.pow(10,9));

        double mu1 = mean(t);
        double std1 = sd(t,mu1);
        System.out.print("Normalizing by standard deviation...");
        time = System.nanoTime();
            for (int i = 0; i < t.length; i++) {
                t[i] /= std1;
                //   x[i] *= std1;
                //     x[i] += mu1;
            }
        time = System.nanoTime()-time;
        System.out.println("Done: " + time/Math.pow(10,9));

            return t;
    }
    //Loads in the correlation matrix from a file, and returns it as a 1D Array
    private static float [] loadCorrelations(double [][] data, String file)
    {
        try {
            BufferedReader b = new BufferedReader(new FileReader(file));
            float [] corrs = new float[(data.length*(data.length-1))/2];
            for(int i = 0; i < data.length;i++)
            {
                // System.out.println(i);
                // double [] vec1 = data[i];
                String [] line = b.readLine().split("\t");
                for(int j = i+1; j < data.length;j++)
                {
                    n = data.length;
                    int index = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
                    //       double [] vec2 = data[j];
                    //double gc = getCorrelation(vec1,vec2);
                    //     double gc = getCorrelationPValue(vec1,vec2);
                   corrs[index] = (float)Float.parseFloat(line[j-i-1]);
                }
            }
            b.close();
            return corrs;
        }
        catch(Exception e) {
            return null;
        }

    }

  //Converts a 2D array of expression data into a dataset that is stored locally within an Independence Test object
    private static void convertToDataSet(double [][]data,ArrayList<Gene> g)
    {
        List<Node> temp = new ArrayList<Node>();
        for(int i = 0; i < data.length;i++)
        {
            ContinuousVariable vTemp = new ContinuousVariable(g.get(i).symbol);
            temp.add(vTemp);
        }
        DataSet ds = new ColtDataSet(data[0].length,temp);
        for(int i = 0; i < data.length;i++) //Loop over genes
        {
            for(int j = 0; j < data[0].length;j++){ //Loop over samples
                ds.setDouble(j,i,data[i][j]); //In the dataset, columns are genes
            }
        }
        indT = new IndTestFisherZ(ds,0.05);
    }

    //Helper function to compute the correlation of two variables
    private static double getCorrelation(double[]vec1,double[]vec2)
    {
        double x = 0;
        double x2 = 0;
        double y = 0;
        double y2 = 0;
        double xy = 0;

        for(int i = 0; i < vec1.length;i++)
        {
            x+= vec1[i];
            x2 += Math.pow(vec1[i],2);
            y+= vec2[i];
            y2+= Math.pow(vec2[i],2);
            xy += vec1[i]*vec2[i];
        }
        return (vec1.length*xy-x*y)/(Math.sqrt(  (vec1.length*x2-Math.pow(x,2)) * (vec1.length*y2-Math.pow(y,2))    ));
    }


    public static ArrayList<Gene> loadGeneDataSimple(int [] diseaseID) throws Exception
    {

        String expFile = "expression_transposed.txt";
        String fcFile = "foldChanges.txt";
        String gdFile = "all_gene_disease_associations.tsv";
        String giFile = "final_gene_data.txt";


        /*****************THIS IS TO GET AN AVERAGE DISEASE SCORE FOR INDIFFERENCE*************/

        BufferedReader b2 = new BufferedReader(new FileReader(gdFile));
        HashMap<String,HashMap<Integer,Double>> allScores = new HashMap<String,HashMap<Integer,Double>>();
        int numSamples = 60000;
        while(b2.ready())
        {
            if(b2.readLine().startsWith("gene"))
                break;
        }
        int count = 0;
        while(b2.ready())
        {
            String [] line = b2.readLine().split("\t");
            String symbol = line[1];
            int disease = Integer.parseInt(line[3].split(":")[1].substring(1,line[3].split(":")[1].length()));
            double score = Double.parseDouble(line[5]);
            if(allScores.get(symbol)==null)
            {
                HashMap<Integer,Double> now = new HashMap<Integer,Double>();
                now.put(disease,score);
                allScores.put(symbol,now);
            }
            else
            {
                HashMap<Integer,Double> now = allScores.get(symbol);
                now.put(disease,score);
                allScores.put(symbol,now);
            }

        }
        b2.close();
        ArrayList<String> temp = new ArrayList<String>(allScores.keySet());
        //System.out.println(temp.size());
        Random rand = new Random();
        double mean = 0;
        for(int i = 0; i < numSamples;i++)
        {

            String gene1 = temp.get(rand.nextInt(temp.size()));
            String gene2 = temp.get(rand.nextInt(temp.size()));
            HashMap<Integer,Double> score1 = allScores.get(gene1);
            HashMap<Integer,Double> score2 = allScores.get(gene2);
            while(score1==null || score2==null)
            {
                gene1 = temp.get(rand.nextInt(temp.size()));
                gene2 = temp.get(rand.nextInt(temp.size()));
                score1 = allScores.get(gene1);
                score2 = allScores.get(gene2);
            }
            double rN = computeScore(score1,score2);
            if(maxDiseaseScore < rN)
                maxDiseaseScore = rN;
            mean+=rN;
				/*System.out.println(computeScore(score1,score2,"euclidean",false));
				System.out.println("Overlap Count: " + overlap(score1,score2));
				System.out.println("Size of set 1: " + score1.keySet().size());
				System.out.println("Size of set 2: " + score2.keySet().size());*/

        }
        mean=mean/numSamples;
        b2.close();
        allScores = null;

        /******************END AVERAGE DISEASE CODE*********************************************/
        int familyInd = 2;
        int chromInd = 3;
        int symId = 0;
        ArrayList<Gene> g = new ArrayList<Gene>();
        BufferedReader b = new BufferedReader(new FileReader(giFile));
        b.readLine();
        HashMap<String,int[]> families = new HashMap<String,int[]>();
        HashMap<String,String> chromosomeStuff = new HashMap<String,String>();
        int numSamples2 = 0;
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            if(line.length < familyInd + 1)
                continue;
            if(line[familyInd].length()<1)
                families.put(line[symId],null);
            else
            {
                String[] fam = line[familyInd].split("\\|");
                int [] famInt = new int[fam.length];
                for(int i = 0; i < fam.length;i++)
                    famInt[i] = Integer.parseInt(fam[i]);
                families.put(line[symId],famInt);
            }
            if(line.length < chromInd+1)
                continue;
            if(line[chromInd].equals("X"))
                line[chromInd] = "23";
            else if (line[chromInd].equals("Y"))
                line[chromInd] = "24";
            chromosomeStuff.put(line[symId],line[chromInd] + "\t" + line[chromInd+1] + "\t" + line[chromInd+2]);

        }
        b.close();
        b = new BufferedReader(new FileReader(fcFile));
        HashMap<String,Double> fcs = new HashMap<String,Double>();
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            if(line[1].equals("-Inf"))
                fcs.put(line[0],-1.0);
            else if(line[1].equals("Inf"))
                fcs.put(line[0],1.0);
            else if(line[1].equals("NA"))
                fcs.put(line[0],0.0);
            else
                fcs.put(line[0],Double.parseDouble(line[1]));
        }
        b.close();
        HashMap<String,HashMap<Integer,Double>> geneDisease = new HashMap<String,HashMap<Integer,Double>>();

        b = new BufferedReader(new FileReader(gdFile));
        boolean ready = false;
        while(b.ready())
        {
            if(ready)
            {
                String [] line = b.readLine().split("\t");
                if(geneDisease.get(line[1])==null)
                {
                    HashMap<Integer,Double> temp2 = new HashMap<Integer,Double>();
                    int disID = Integer.parseInt(line[3].split("C")[1]);
                    temp2.put(disID,Double.parseDouble(line[5]));
                    geneDisease.put(line[1],temp2);
                }
                else
                {
                    HashMap<Integer,Double> temp2 = geneDisease.get(line[1]);
                    int disID = Integer.parseInt(line[3].split("C")[1]);
                    temp2.put(disID,Double.parseDouble(line[5]));
                    geneDisease.put(line[1],temp2);
                }
            }
            if(b.readLine().startsWith("gene"))
                ready = true;


        }

        b.close();
        b = new BufferedReader(new FileReader(expFile));
        b.readLine();
        int idCount = 0;
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            numSamples2 = line.length-1;
            Gene curr = new Gene(idCount);

            String symbol = line[0];
            curr.symbol = symbol;
            if(fcs.get(symbol)!=null)
            {
                curr.foldChange = fcs.get(symbol);
                fcs.remove(symbol);
            }
            else
            {
                curr.foldChange = 0;
            }
            if(families.get(symbol)!=null)
            {
                curr.familyID = families.get(symbol);
                families.remove(symbol);
            }
            else

            {
                curr.familyID = null;
            }
            if(geneDisease.get(symbol)!=null)
            {
                curr.diseaseScores = geneDisease.get(symbol);
                geneDisease.remove(symbol);
            }
            else
            {
                curr.diseaseScores = null;
            }
            if(chromosomeStuff.get(symbol)!=null)
            {
                String[] chromo = chromosomeStuff.get(symbol).split("\t");
                curr.chromosome = Integer.parseInt(chromo[0]);
                curr.chromosomeStart = Long.parseLong(chromo[1]);
                curr.chromosomeStop = Long.parseLong(chromo[2]);

                chromosomeStuff.remove(symbol);
            }
            else
            {
                curr.chromosome = -1;
            }
            idCount++;
            g.add(curr);
        }
        chromosomeStuff = null;
        geneDisease = null;
        fcs = null;
        families = null;
        b = new BufferedReader(new FileReader(expFile));
        b.readLine();

        double [][] data = new double[g.size()][numSamples2];
        //System.out.println("Num Genes: " + g.size());
        //	System.out.println("Number Samples: " + numSamples2);
        int currIdx = 0;

        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            for(int i = 1; i < line.length;i++)
            {
                data[currIdx][i-1] = Double.parseDouble(line[i]);
            }
            currIdx++;
        }
        b.close();
       // corrMat = computeCorrelations(data);
        return g;
    }
    //Read in and compute correlation matrices and store this in memory
    //Have the correlation matrix only be one copy in memory but reference to it stored in each class


    //Computes and stores the correlation between gene at index i and gene at index j
    public static double computeCorrelation(int i ,int j)
    {
        if(i > j)
        {
            int swap = i;
            i = j;
            j = swap;
        }
        int index = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
        if(corrMat[index]!=0)
            return corrMat[index];
        double [] one = data[i];
        double [] two = data[j];
        corrMat[index] = (float)StatUtils.correlation(one,two);
        return corrMat[index];
    }

    //For Latent Variable Identification, try to find potential blockers of two genes
    public static boolean getPotentialBlocker(Gene one, Gene two, Gene gen)
    {
        return indT.isIndependent(indT.getVariable(one.symbol),indT.getVariable(two.symbol),indT.getVariable(gen.symbol));
    }

    //Return subset of a DataSet d based on the genes in the filename and in the list of clinical variables
    public static DataSet subsetData(DataSet d,String filename,ArrayList<String>clinical)
    {
        try {
            BufferedReader b = new BufferedReader(new FileReader(filename));
            List<Node> genes = new ArrayList<Node>();
            while (b.ready()) {
                genes.add(d.getVariable(b.readLine()));
            }
            b.close();
            for(String x:clinical)
                genes.add(d.getVariable(x));
            return d.subsetColumns(genes);
        }
        catch(Exception e)
        {
            return null;
        }
    }


}