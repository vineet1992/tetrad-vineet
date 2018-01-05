package edu.pitt.csb.Pref_Div;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.search.IndTestPartialCorrelation;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.util.StatUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.csb.mgm.IndTestMultinomialAJ;

import java.util.*;
import java.io.*;
public class Functions
{
    private static final long [] CHROMOSOME_LENGTHS = {249250621,243199373,198022430,191154276,180915260,171115067,159138663,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415};//Chromosomal lengths in base pairs

  //TODO Compute this on the fly instead of offline
   private static final double MEDIAN_DISEASE_SCORE = 0.000543;
    public static IndependenceTest indT;
    private static float [] corrMat;
    private static byte [] indep;
    private static int n;
    private static double[][] data;
    private static double averageDiseaseScore;//Average disease score to return when disease information is not present
    private static double maxDiseaseScore;
    private static double maxDiseaseIntensity;
    private static double maxFoldChange;
    private static boolean normalizeDisease = false;
    private static PrintStream out;
    private static PrintStream out2;

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
        out2.println(g1.foldChange + "\t" + disScore);
        out2.flush();
        return a*Math.abs(g1.foldChange/maxFoldChange)+(1-a)*disScore/maxDiseaseIntensity;
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
        out.println(DS + "\t" + CD + "\t" + family + "\t" + corrMat[index] );
        //System.out.println("Corr: " + Math.abs(corrMat[g1.ID][g2.ID]) + ",DS:" + DS + ",CD:" + CD + ",FM:" + family);
        return temp;
    }
    //Computes Disease Similarity Score (now using cosine similarity)
    private static double computeScore(HashMap<Integer,Double> score1, HashMap<Integer,Double> score2)
    {
        //If we have no information, report as -1 (unknown, will replace with median)TODO
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
                    vecTwo.add(MEDIAN_DISEASE_SCORE);
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
                    vecOne.add(MEDIAN_DISEASE_SCORE);
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
            return 0.5;//Be indifferent in the abscence of information
        if(g1.chromosome!=g2.chromosome)
            return 1;
        else
        {
            if(g1.chromosomeStart==-1 || g2.chromosomeStart==-1 || g1.chromosomeStop==-1 || g2.chromosomeStop==-1)
                return 0.5;
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
            return 0.5;//We don't have family ID information for this gene so be indifferent
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


        /*****************THIS IS TO GET AN AVERAGE DISEASE SCORE FOR INDIFFERENCE*************/


        //Read through Gene-Disease Association file
        BufferedReader b2 = new BufferedReader(new FileReader(gdFile));
        HashMap<String,HashMap<Integer,Double>> allScores = new HashMap<String,HashMap<Integer,Double>>();
        int numSamples = 500000; //Number of samples to obtain in order to compute distribution of gene-disease association scores
                                // (too many samples if you do genes x genes-1)
        while(b2.ready())
        {
            if(b2.readLine().startsWith("gene"))
                break;
        }
        int count = 0;

        //Read in the gene-disease file line-by-line
        while(b2.ready())
        {
            String [] line = b2.readLine().split("\t");
            String symbol = line[1];
            int disease = Integer.parseInt(line[3].split(":")[1].substring(1,line[3].split(":")[1].length()));
            double score = Double.parseDouble(line[5]);

            //Create a new hashmap entry if we haven't seen this gene previously
            if(allScores.get(symbol)==null)
            {
                HashMap<Integer,Double> now = new HashMap<Integer,Double>();
                now.put(disease,score);
                allScores.put(symbol,now);
            }
            //Just add this new disease pairing to the old hashmap
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
        //Get a random pair of genes and compute their gene-disease distance score
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
            //maxDiseaseScore to normalize the distribution
            if(maxDiseaseScore < rN)
                maxDiseaseScore = rN;
            mean+=rN;
        }
        //Mean of the distribution
        mean=mean/numSamples;
        b2.close();
        //null out the HashMap to save memory
        allScores = null;
        averageDiseaseScore = mean;

        /******************END AVERAGE DISEASE CODE*********************************************/
        //Should we normalize the disease score?
        normalizeDisease = false;
        //Constant indices of each piece of information from the Gene Information file//
        int familyInd = 2;
        int chromInd = 3;
        int symId = 0;
        /*******************************************************************************/
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
        b = new BufferedReader(new FileReader(fcFile));

        //Read fold change file between different phenotypic groups for each gene
        //First compute the maximumFoldChange value
        while(b.ready())
        {
            String fol = b.readLine().split("\t")[1];
            if(fol.contains("Inf") || fol.contains("NA"))
                continue;
            else
            if(Math.abs(Double.parseDouble(fol)) > maxFoldChange)
                maxFoldChange = Math.abs(Double.parseDouble(fol));
        }
        b.close();
        b = new BufferedReader(new FileReader(fcFile));
        HashMap<String,Double> fcs = new HashMap<String,Double>();
        //Load HashMap with all foldChange values
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

        //Find the maximum disease Intensity score for normalization purposes
        //dIntensityScore finds the value for the disease of interest
        for(String x: geneDisease.keySet())
        {
            HashMap<Integer,Double> t = geneDisease.get(x);
            if(dIntensityScore(t,diseaseID) > maxDiseaseIntensity)
                maxDiseaseIntensity = dIntensityScore(t,diseaseID);
        }


        b.close();
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
                fcs.remove(symbol);
            }
            else
            {
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
        fcs = null;
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

        System.out.println("Computing Correlations...");
        n = dat.length;
        System.out.println(n);

        corrMat = new float[(dat.length*(dat.length-1))/2];
        indep = new byte[(dat.length*(dat.length-1))/2];


        //If we're using Pref-Div just to filter out genes to get a set for analysis, don't need to construct independence test object
        //This is only for latent variable prediction
        if(!filter)
            convertToDataSet(dat,g);
        else
            data = dat;
        System.out.println("MFC: " + maxFoldChange);
        System.out.println("MDS: " + maxDiseaseIntensity);
         out = new PrintStream("gene_disease_distribution.txt");
         out2 = new PrintStream("intensity_distribution.txt");
         System.out.println(numSamples);
        for(int i = 0; i < numSamples;i++)
        {
            Gene gene1 = g.get(rand.nextInt(g.size()));
            Gene gene2 = g.get(rand.nextInt(g.size()));
            computeDistance(gene1,gene2);
        }
        out.flush();
        out.close();
        return g;
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
        averageDiseaseScore = mean;

        /******************END AVERAGE DISEASE CODE*********************************************/
        normalizeDisease =false;
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
        for(String x: geneDisease.keySet())
        {
            HashMap<Integer,Double> t = geneDisease.get(x);
            if(dIntensityScore(t,diseaseID) > maxDiseaseIntensity)
                maxDiseaseIntensity = dIntensityScore(t,diseaseID);
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

    //Return true if one and two are independnet, false otherwise
    public static boolean computePValue(Gene one, Gene two)
    {
       // System.out.println(one + "," + two);
        int ii = -1;
        int jj = -1;
        if(one.ID < two.ID)
        {
            ii =  one.ID;
            jj = two.ID;
        }
        else
        {
            ii = two.ID;
            jj = one.ID;
        }
        int index = (n*(n-1)/2) - (n-ii)*((n-ii)-1)/2 + jj - ii - 1;
        if(indep[index]==(byte)0) {
            indep[index] = indT.isIndependent(indT.getVariable(one.symbol),indT.getVariable(two.symbol)) ? (byte)1:(byte)-1;
        }
        if(indep[index]==(byte)1)
            return true;
        else
            return false;

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