package edu.pitt.csb.Pref_Div;
import java.util.*;
import java.io.*;

public class PrefExperiment
{
	public static void main(String [] args) throws Exception
	{
		boolean diverseB = false;
		double [] i_a = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};
		//double [] i_a = {0,0.25,0.5,0.75,1};
		PrintStream out2 = new PrintStream("intensity_values.txt");
		double intensity_a = 0.3;
		double radius = 0.25;
		double accuracy = 0.5;
		int [] diseaseID = {6142}; //Breast Cancer
		String [] files = new String[4];
		files[0] = "TCGA_BRCA_FILES/expression_transposed_npn.txt";
		files[1] = "TCGA_BRCA_FILES/foldChanges_new.txt";
		files[2] = "all_gene_disease_associations.tsv";
		files[3] = "final_gene_data.txt";
		ArrayList<Gene> g = Functions.loadGeneData(".",files,diseaseID,true);
		for(Gene x: g)
		{
			x.intensityValue = Functions.computeIntensity(x,intensity_a,diseaseID);
			out2.println(x.symbol+"\t"+x.intensityValue); //Print all intensity values out to file
		}
		int [] topK = {50,100,500};
		
		double [] constants = {0,0.25,0.5,0.75,1.0};
		
		for(int k = 0; k < topK.length;k++) //Loop over sizes of returned results (Top-K)
					{

						//Compute the intensity value of all genes, based on "alpha" parameter, and the input disease ID for disease of interest
						for(Gene x: g)
						{
							x.intensityValue = Functions.computeIntensity(x,intensity_a,diseaseID);
							out2.println(x.symbol+"\t"+x.intensityValue); //Print all intensity values out to file
						}
						//Sort by intensity
						Collections.sort(g,Gene.IntensityComparator);
						System.out.println(k);
						//Print out to the file for each size of returned results
						PrintStream out = new PrintStream("output_" + k + ".txt");
						PrintStream out3 = new PrintStream("clusters_" + k + ".txt");
						//Run Pref-Div algorithm to find a diverse and relevant gene set
						PrefDiv pd = new PrefDiv(g,topK[k],accuracy,radius,PrefDiv.findTopKIntensity(g,topK[k]));
						pd.setCluster(true);
						ArrayList<Gene> Result = pd.diverset();

						//Print out some statistics about the result
						System.out.println("PD diversity = " + pd.averageDistance());
						System.out.println("PD normalizedIntensity = " + pd.normalizedIntensity());
						System.out.println("PD Coverage = " + pd.coverage());
						//Print the output gene list to file
						for(Gene stuff: Result)
						{
							out.println(stuff.symbol);
						}
						for(Gene stuff:pd.clusters.keySet())
						{
							List<Gene> temp = pd.clusters.get(stuff);
							out3.print(stuff.symbol + "\t");
							for(int i = 0; i < temp.size();i++)
							{
								if(i==temp.size()-1)
									out3.println(temp.get(i).symbol);
								else
									out3.print(temp.get(i).symbol+"\t");
							}
						}
						out3.flush();
						out3.close();
						out.flush();
						out.close();
					}
		
		/*
		if(!diverseB)
		{
			for(int i = 0; i < i_a.length; i ++)
			{
				for(int j = 0; j < topK.length;j++)
				{
					for(Gene x: g)
					{
						x.intensityValue = Functions.computeIntensity(x,intensity_a,diseaseID);
						out2.println(x.symbol+"\t"+x.intensityValue);
					}
					Collections.sort(g,Gene.IntensityComparator);
						PrintStream out = new PrintStream("output_" + i + "_" + j + ".txt");
						PrefDiv pd = new PrefDiv(g,topK[j],i_a[i],radius,PrefDiv.findTopKIntensity(g,topK[j]));
						pd.diverset();
						ArrayList<Gene> Result = pd.getResultGenes();
						System.out.println("PD diversity = " + pd.averageDistance());
						System.out.println("PD normalizedIntensity = " + pd.normalizedIntensity());
						System.out.println("PD Coverage = " + pd.coverage());
						for(Gene stuff: Result)
						{
							out.println(stuff.symbol);
						}
						out.flush();
							out.close();
				}
			}
		}
		/*else
		{
			for(int a = 0; a < 4 ;a++) //choose a score to parametrize
			{
				for(int b = 0; b < 5; b++ )//choose a value (0.25,0.5,0.75,1.0)
				{
					double [] diverse = new double[3];
					diverse[0] = (1-constants[b])/3;
					diverse[1] = diverse[0];
					diverse[2] = diverse[1];
					if(a!=3)
					{
						diverse[a] = constants[b];
					}
					
					
					for(Gene x: g)
						x.intensityValue = Functions.computeIntensity(x,intensity_a,diseaseID);
					Collections.sort(g,Gene.IntensityComparator);
					for(int k = 0; k < topK.length;k++)
					{
						System.out.println(a + "_" + b + "_" + k);
						PrintStream out = new PrintStream("output_" + a + "_" + b + "_" + k + ".txt");
						PrefDiv pd = new PrefDiv(diverse,g,topK[k],accuracy,radius,PrefDiv.findTopKIntensity(g,topK[k]));
						pd.diverset();
						ArrayList<Gene> Result = pd.getResultGenes();
						System.out.println("PD diversity = " + pd.averageDistance());
						System.out.println("PD normalizedIntensity = " + pd.normalizedIntensity());
						System.out.println("PD Coverage = " + pd.coverage());
						for(Gene stuff: Result)
						{
							out.println(stuff.symbol);
						}
						out.flush();
							out.close();
					}
				}
			}
		}*/
	}
}