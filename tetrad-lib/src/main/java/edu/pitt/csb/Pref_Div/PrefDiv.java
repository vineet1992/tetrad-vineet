package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.search.Fges;
import edu.cmu.tetrad.util.StatUtils;
import edu.pitt.csb.Pref_Div.Comparisons.ComparablePD;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.*;


public class PrefDiv implements ComparablePD {
	
    public int topK;
    public double accuracy;
    public double radius;
    public ArrayList<Gene> items;
    public ArrayList<Gene> result;
    public ArrayList<Gene> G = new ArrayList<Gene>();
    public ArrayList<Gene> B = new ArrayList<Gene>();

    public Map<Gene,List<Gene>> clusters;
    private float [] corrMat;
    private float [] theoryMat;
   private double alpha;
   private int preSampleSize = 5000;
   private DataSet data;
   private boolean clustering = false; //Are we clustering the top genes with near neighbors?
   private boolean clusterByCorrs = false;//Should we use all correlations to cluster instead of only high intensity genes?
   private boolean computeCorrs;
   private boolean corrsGiven;
   private boolean partialCorr;
   private double radiusWP;
   private boolean [] withPrior;
   private boolean usingPrior = false;
   private ArrayList<Float> sampledCorrs;
   private ArrayList<Gene> selected;



    public PrefDiv(ArrayList<Gene> items, int topK, double accuracy, double radius, float [] theory,double alpha,DataSet data,boolean computeAllCorrs,boolean partialCorr) {
    	this.topK = topK;
    	this.accuracy = accuracy;
    	this.radius = radius;
    	this.items = items;
    	this.result = new ArrayList<Gene>();
    	this.theoryMat = theory;
    	this.corrsGiven = false;
    	this.partialCorr = partialCorr;
    	if(computeAllCorrs)
    	    corrMat = Functions.computeAllCorrelations(items,data,partialCorr,false,false,1);
    	else {
            corrMat = new float[(items.size() * (items.size() - 1)) / 2];
        }
    	 this.alpha = alpha;
    	this.data = data;
    	this.computeCorrs = computeAllCorrs;
    }

    //Constructor for already computed data intensity and dissimilarity but using theory as well
    public PrefDiv(ArrayList<Gene> items, int topK, double accuracy, double radius, float[] corrs, float [] dissim, double alpha) {
        this.partialCorr = false;
        this.corrMat = corrs;
        this.computeCorrs = false;
        this.corrsGiven = true;
        this.topK = topK;
        this.accuracy = accuracy;
        this.radius = radius;
        this.items = items;
        this.result = new ArrayList<Gene>();
        this.theoryMat = dissim;
        this.alpha = alpha;

    }
    //Constructor when only using data
	public PrefDiv(ArrayList<Gene> items, int topK, double accuracy, double radius, float[] corrs) {
		this.partialCorr = false;
		this.corrMat = corrs;
		this.computeCorrs = false;
		this.corrsGiven = true;
        this.topK = topK;
    	this.accuracy = accuracy;
    	this.radius = radius;
    	this.items = items;
    	this.result = new ArrayList<Gene>();
    	this.theoryMat = null;
    	this.alpha = 1;

    }



    //Constructor when only using data
    public PrefDiv(ArrayList<Gene> items, int topK, double accuracy, double radius, double radiusWP, boolean [] wp, float[] corrs) {
        this.partialCorr = false;
        this.corrMat = corrs;
        this.computeCorrs = false;
        this.corrsGiven = true;
        this.topK = topK;
        this.accuracy = accuracy;
        this.radius = radius;
        this.items = items;
        this.result = new ArrayList<Gene>();
        this.theoryMat = null;
        this.alpha = 1;
        this.radiusWP = radiusWP;
        this.withPrior = wp;
        this.usingPrior = true;

    }

    public Map<Gene,List<Gene>> getLastCluster(){return clusters;}
    public List<Gene> getLastSelected(){return selected; }
    //============================================================
    //====================== Algorithm ===========================
    //============================================================
    

    public void setPreSampleSize(int ss)
    {
        this.preSampleSize = ss;
    }



    //Sets up so that we are clustering with Pref-Div
    //Input: cluster -> whether or not we cluster using all correlations or just using the standard Pref-Div approach
    public void setCluster(boolean cluster)
    {
        clustering = true;
        clusters = new HashMap<Gene,List<Gene>>();
        clusterByCorrs = cluster;
    }
	public static double findTopKIntensity(ArrayList<Gene> record, int topK){
		double result = 0;
		for(int i = 0; i < topK; i++)
			result+=record.get(i).intensityValue;
		return result;
	}
    
    
    public ArrayList<Gene> diverset() {

	    if(!computeCorrs && !corrsGiven)
	        sampleCorrelations();
	    long time = System.nanoTime();
      //  System.out.println("original item size = " + this.items.size());
    	
    	if(this.items.size() == topK){   		
            this.result = items;
            return result;
    	}

        ArrayList<Gene> result = new ArrayList<Gene>();
        // edge handling
        if (this.topK > this.items.size()) {
            this.result = this.items;
            return this.result;
        } else if (topK == 0){
            System.out.println("Error: topK can not be 0");
            this.result = result;
            return this.result;
        } else if (accuracy > 1) {
            System.out.println("Error: Accuracy should not be bigger than 1");
            this.result = result;
            return this.result;
        }
        
        double acc = this.accuracy;
        
        int setNum = Integer.valueOf(this.items.size()/this.topK);
        if (this.items.size() % this.topK != 0){
            setNum++;
        }

        ArrayList<Gene> temp = new ArrayList<Gene>();
        //Outside while loop? For loop here
        for (int i = 0; i < setNum; i++){

            ArrayList<Gene> S = new ArrayList<Gene>();
            //Add the k objects with the highest intensity from P to S
            for (int j = i*this.topK; j < Math.min((i+1)*this.topK, this.items.size()); j++) {
                S.add(this.items.get(j));
            }

            //Find all genes in S that aren't within radius distance of result TODO Add clustering
            temp = pairEliminator(result, S, this.radius);
            
            //Eliminate variables within temp that are within radius distance of another variable from temp (choose higher intensity)TODO Add clustering
            temp = eliminate(temp, radius);

            result.addAll(temp);


            // Use B according to accuracy
            if (i != 0)
                acc = acc/2;

        	//Sort all marked genes
        	Collections.sort(B, new Comparator<Gene>() {    
                @Override
                public int compare(Gene o1, Gene o2) {
                    return new Double(o2.intensityValue).compareTo(o1.intensityValue);
                }               
        	});
        	
            //Add enough marked genes to get to A*k for this iteration
            //TODO Right now we keep this as is, so it could be the case that a gene in a cluster is also a representative gene in the Top K
            //TODO Read the Pref-Div paper and decide if we need to do this (especially the else)
            int remain = (int)(acc*this.topK) - temp.size();
            if (remain > 0 && B.size() > 0) {	
                for (int j = 0; j < Math.min(remain,B.size()); j++){
                    result.add(B.get(0));
                    B.remove(0);
                }
            }/*else {
                if(B.size() > 0) {
                    result.add(B.get(0));
                }
            }*/
            
            if (i == 0) {
                G.addAll(B);
            }
            B.clear();
            // result bigger then required topK
            int oversize = result.size() - this.topK;
            if (oversize == 0) {
                break;
            } else if (oversize > 0){
                for (int j = result.size() - 1; j >= this.topK; j--) {
                    result.remove(j);
                }
                result.trimToSize();
                break;
            }
        }

        int addUp = this.topK - result.size();
        if (addUp > 0) {
      //  	System.out.println("addUp size = " + addUp);
            for (int i = 0; i < addUp; i++) {
                result.add(G.get(i));
            }
        }
        B.clear();
        G.clear();
        this.result.addAll(result);

        for(int x = 0; x < result.size();x++)
        {
            if(clusters.get(result.get(x))==null)
                clusters.put(result.get(x),null);
        }

        List<Gene> toRemove = new ArrayList<Gene>();
        for(Gene g: clusters.keySet())
        {
            if(!result.contains(g))
                toRemove.add(g);
        }
        for(int i = 0; i < toRemove.size();i++)
            clusters.remove(toRemove.get(i));

        if(clusterByCorrs)
        {
            addAllToClusters();
        }
        time = System.nanoTime()-time;
        selected = result;
       // System.out.println("Acutal Pref-Div Time: " + time/Math.pow(10,9));
        return result;
   
    }


    //Add all significantly correlated genes to each key in the hashmap
    private void addAllToClusters()
    {
        /***For each selected gene...***/
        for(Gene g:clusters.keySet())
        {
            List<Gene>list = new ArrayList<Gene>();

            /***For each gene in the gene list***/
            for(int i = 0; i < items.size();i++)
            {

                /***Get the index in the correlation matrix based upon the placement in the gene list***/
                int index = getIndex(g,items.get(i));
                if(g.ID==items.get(i).ID)
                    continue;


                if(usingPrior && withPrior[index] && distance(g,items.get(i)) < radiusWP)
                {
                    list.add(items.get(i));
                }
                else if((!usingPrior || !withPrior[index]) && distance(g,items.get(i))<radius)
                {
                    list.add(items.get(i));
                }
            }
            clusters.put(g,list);
        }
    }

    //Identify lower intensity genes to be removed from input, since another similar gene is already in input
    public ArrayList<Gene> eliminate(ArrayList<Gene> input, double radius){
        ArrayList<Gene> output = new ArrayList<Gene>(input);
        ArrayList<Gene> finalOutput = new ArrayList<Gene>();
        if(input.size() != 0){
            ArrayList<Gene> temp = new ArrayList<Gene>();
            Gene coreElement;
            int i = 0;
            while (i < output.size()) {
                coreElement = output.get(i);
                finalOutput.add(coreElement);
                if (output.size() != 1) {
                    for (int j = i+1; j < output.size(); j++) {

                        if(usingPrior && withPrior[getIndex(output.get(j),coreElement)] && distance(output.get(j),coreElement) > radiusWP)
                        {
                            temp.add(output.get(j));
                        }
                        else if(usingPrior && !withPrior[getIndex(output.get(j),coreElement)] && distance(output.get(j),coreElement)>radius)
                        {
                            temp.add(output.get(j));
                        }
                        else if (!usingPrior && distance(output.get(j),coreElement) > radius) {
                            temp.add(output.get(j));
                        } else {
                            if(clusters.get(coreElement)==null)
                            {
                                List<Gene> clust = new ArrayList<Gene>();
                                clust.add(output.get(j));
                                clusters.put(coreElement,clust);
                            }
                            else
                            {
                                List<Gene> clust = clusters.get(coreElement);
                                clust.add(output.get(j));
                                clusters.put(coreElement,clust);
                            }
                            B.add(output.get(j));
                        }
                    }
                    output.clear();
                    output.addAll(finalOutput);
                    output.addAll(temp);
                    temp.clear();
                }
                i++;
            }
        }
        return finalOutput;
    }

    //R, S, radius
    //Returns all genes from input (S) that are further than radius distance from the genes in finalResult (R)

    public ArrayList<Gene> pairEliminator(ArrayList<Gene> finalResult, ArrayList<Gene> input, double radius){
        ArrayList<Gene> output = new ArrayList<Gene>(input);
        //System.out.println("pairEliminator size:" + input.size());
        if (finalResult.size() != 0) {
            ArrayList<Gene> temp = new ArrayList<Gene>();
            for (int i = 0; i < finalResult.size(); i++) {
                Gene toBeCompare = finalResult.get(i);
                if (output.size() == 0) {
                    break;
                } else {
                    for (int j = 0; j < output.size(); j++) {
                        if(usingPrior && withPrior[getIndex(output.get(j),toBeCompare)] && distance(output.get(j),toBeCompare) > radiusWP)
                        {
                            temp.add(output.get(j));
                        }
                        else if(usingPrior && !withPrior[getIndex(output.get(j),toBeCompare)] && distance(output.get(j),toBeCompare)>radius)
                        {
                            temp.add(output.get(j));
                        }
                        else if (distance(output.get(j),toBeCompare) > radius) {
                            temp.add(output.get(j));
                        } else {
                            B.add(output.get(j));
                            if(clusters.get(toBeCompare)==null)
                            {
                                List<Gene> clust = new ArrayList<Gene>();
                                clust.add(output.get(j));
                                clusters.put(toBeCompare,clust);
                            }
                            else
                            {
                                List<Gene> clust = clusters.get(toBeCompare);
                                clust.add(output.get(j));
                                clusters.put(toBeCompare,clust);
                            }
                        }

                    }
                    output.clear();
                    output.addAll(temp);
                    temp.clear();
                }
            }
        }
        return output;
    }
    
    public static boolean isNumeric(String s) {  
        return s.matches("[-+]?\\d*\\.?\\d+");  
    }  
    
    public static String[] emptyReomver(String[] firstArray){
        List<String> list = new ArrayList<String>();
        for(String s : firstArray) {
           if(s != null && s.length() > 0 && !isNumeric(s)) {
              list.add(new String(s));
           }
        }
        return list.toArray(new String[list.size()]);
    }
    

    private void sampleCorrelations()
    {
        Random rand = new Random();
        sampledCorrs = new ArrayList<Float>();
        double [][] temp = data.getDoubleData().transpose().toArray();
        int [] indices = new int[preSampleSize];
        float [] preCorrs = new float[preSampleSize];
        for(int i = 0; i < this.preSampleSize;i++)
        {
            int x = 0;
            int y = 0;
            while(x==y)
            {
                x = rand.nextInt(items.size());
                y = rand.nextInt(items.size());
            }
            if(x > y)
            {
                int t = x;
                x = y;
                y = t;
            }
            int index = Functions.getIndex(x,y,items.size());
            float corr = (float)Math.abs(StatUtils.correlation(temp[data.getColumn(data.getVariable(items.get(x).symbol))],temp[data.getColumn(data.getVariable(items.get(y).symbol))]));
            indices[i] = index;
            preCorrs[i] = corr;
            sampledCorrs.add(corr);

        }
        preCorrs = Functions.NPN(preCorrs,true);
        for(int i = 0; i < preCorrs.length;i++)
        {
            corrMat[indices[i]] = preCorrs[i];
        }

        Collections.sort(sampledCorrs);
    }
    private float corrOnTheFly(int i, int j)
    {
        double [][] temp = data.getDoubleData().transpose().toArray();
        float corr = (float)Math.abs(StatUtils.correlation(temp[data.getColumn(data.getVariable(items.get(i).symbol))],temp[data.getColumn(data.getVariable(items.get(j).symbol))]));
        corr = Functions.NPNonTheFly(corr,sampledCorrs);
        corrMat[Functions.getIndex(i,j,items.size())] = corr;
        return corr;
    }
    public static boolean compareFloatArr(float[] fArr1,float[] fArr2){
    	if(fArr1.length != fArr2.length)
    		return false;
    	for(int i = 0; i < fArr1.length; i++){
    		if(fArr1[i] != fArr2[i])
    			return false;
    	}
    	return true;
    }


    public double distance(Gene gene1, Gene gene2) {

        int i = gene1.ID;
        int j = gene2.ID;
        if(gene1.ID > gene2.ID)
        {
            int temp = i;
            i = j;
            j = temp;
        }
        int index = Functions.getIndex(i,j,items.size());
        float c = corrMat[index];
        if(!computeCorrs && !corrsGiven)
        {
            if(c==0.0f)
            {
                c = corrOnTheFly(i,j);
            }
        }
        c = 1-c; //c is a probability of being connected, so 1-c can be a dissimilarity or a distance measure
        if(theoryMat==null)
            return c;
        if(theoryMat[Functions.getIndex(i,j,items.size())]==-1)
            return c;
        else
            return c*alpha + (1-theoryMat[Functions.getIndex(i,j,items.size())])*(1-alpha);
    }
    private int getIndex(Gene g1, Gene g2)
    {
        int i = g1.ID;
        int j = g2.ID;
        if(g1.ID > g2.ID)
        {
            int temp = i;
            i = j;
            j = temp;
        }
        return Functions.getIndex(i,j,items.size());
    }


 
    
    //======================================================================
    //======================== Statistic methods ===========================
    //======================================================================
    
    public int getSize() {
        return this.result.size();
    }
  
    public void printResultGeneIds(){
    	for(int i = 0; i < this.result.size(); i++)
    		System.out.println("Gene " + i + " id = " + this.result.get(i).ID);
    }
    public ArrayList<Gene> getResultGenes(){
    	return this.result;
    }
    
    
    public static void writeResulttoFile (String file, String result) {  	 
        BufferedWriter bw = null;
        try {
           bw = new BufferedWriter(new FileWriter(file, true));
	       bw.write(result);
	       bw.newLine();
	       bw.flush();
        } catch (IOException ioe) {
	       ioe.printStackTrace();
        } finally {                    
	       if (bw != null) try {
	          bw.close();
	       } catch (IOException ioe2) {
	          
	       }
        } 
   
     } 
    

    public double sumDistance() {
        double sum = 0;
        for (int i = 0; i < this.result.size() - 1; i++) {
            for (int j = i + 1; j < this.result.size(); j++) {
                sum += distance(this.result.get(i),this.result.get(j));
            }
        }
        return sum;
    }
    
    public double averageDistance() {
        double distSum = 0;
        double counter = 0;
        
        for (int i = 0; i < this.result.size() - 1; i++) {
            for (int j = i + 1; j < this.result.size(); j++ ) {
            	double tempDist = distance(this.result.get(i),this.result.get(j));
            	if(tempDist > 1.0) System.out.println("dist = " + tempDist);
                distSum += distance(this.result.get(i),this.result.get(j));
                counter++;
            }
        }
        return distSum/counter;
    }
    
    
    public double maxDistance() {
    	double max = Double.MIN_VALUE;
        double distSum = 0;
        int counter = 0;
        for (int i = 0; i < this.result.size() - 1; i++) {
            for (int j = i + 1; j < this.result.size(); j++ ) {
            	double dist = distance(this.result.get(i),this.result.get(j));
                if(dist > max)
                	max = dist;
            }
        }
        return max;
    }    

    public double sumIntensity(){
        double averageSum = 0;
        for (int i = 0; i < this.result.size(); i++) {
            averageSum += this.result.get(i).intensityValue;
        }
        return averageSum;
    }
    
    public double averageIntensity() {
    	double averageSum = 0;
        for (int i = 0; i < this.result.size(); i++) {
            averageSum = averageSum + this.result.get(i).intensityValue;
        }
        return averageSum/this.result.size();
    }


    public double coverage() {
    	double count = 0.0; 
    	for(int i = 0; i < items.size(); i++){
            for (int j = 0; j< result.size(); j++) {
                if (distance(result.get(j),items.get(i)) <= radius) {
                    count++;
                    break;
                }
            }
        }
        return count/(double)items.size();
    }
    //---------------------------------------------
}
