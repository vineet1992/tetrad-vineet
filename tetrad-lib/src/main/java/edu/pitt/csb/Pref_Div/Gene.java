package edu.pitt.csb.Pref_Div;

import java.util.*;

public class Gene
{
    public Gene(int d)
    {
        ID = d;
    }
    public String toString()
    {
        String x = "ID: " + Integer.toString(ID) + "\n";
        x+= "Fold Change: " + Double.toString(foldChange) + "\n";
        x+= "Gene Symbol: " + symbol + "\n";
        x+= "Family IDs: " + Arrays.toString(familyID) + "\n";
        x+= "Chromosome Number: " + chromosome + ", from " + chromosomeStart + " to " + chromosomeStop + "\n";
        x+= "Disease Scores: " + diseaseScores;
        return x;
    }
    int ID; //unique ID for the gene
    public double foldChange; //denotes the fold change in values from control to tumor patients
    public String symbol; //denotes the official Gene Symbol for this Gene
    int [] familyID; //denotes the gene family this gene belongs to
    int chromosome; //denotes the chromosome number upon which this gene is located
    long chromosomeStart; //denotes the starting base pair for this gene
    long chromosomeStop; //denotes the stopping base pair for this gene
    HashMap<Integer,Double> diseaseScores; //Denotes all
    double intensityValue;


    public double compareTo(Gene p) {
        return p.intensityValue - this.intensityValue;
    }

    public static Comparator IntensityComparator = new Comparator() {
        public int compare(Object o1, Object o2) {
            double geneIntensity1 = ((Gene) o1).intensityValue;
            double geneIntensity2 = ((Gene) o2).intensityValue;
            if((geneIntensity1 - geneIntensity2) > 0){
                return -1;
            }else if((geneIntensity1 - geneIntensity2) < 0){
                return 1;
            }else{
                return 0;
            }
        }
    };

}