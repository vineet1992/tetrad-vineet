package edu.pitt.csb.Pref_Div;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.Fci;
import edu.pitt.csb.mgm.MixedUtils;
import org.apache.commons.math3.stat.StatUtils;

import java.util.*;
import java.io.*;

public class PrefExperiment
{


	//TODO Add if target variable is null, just use the highest variance genes as the intensity score, implemented need to debug
	//TODO if disease IDs are null, just use the data for intensity, implemented need to debug
	public static void main(String [] args) throws Exception
	{
		//PARAMS
		int numAlphas = 20;
		int ns = 20;
		double threshold = 0.1;
		double accuracy = 0;
		double radius = 0.5; //NEED to figure out how to set this in a principled way TODO
		int [] topK = {50};
		int [] diseases = {6142}; //Default is breast cancer
		boolean normalizeFiles = false;
		boolean loocv = false;
		boolean approxCorrelations = true;
		String dataFile = "";
		String intensityFile = "";
		String dissimilarityFile = "";
		String outputFile = "";
		String clusterFile = "";
		String graphFile = "";
		String target = "";
		String [] files = new String[5];
		//files[0] = "TCGA_BRCA_FILES/expression_transposed_npn.txt";
		files[0] = "MUC1_RNA_Seq.txt"; //File with only gene expression, rows are genes, columns are samples
		//files[1] = "MI_Genes.txt";
		//files[1] = "TCGA_BRCA_FILES/foldChanges_new.txt"; //TODO Need to compute this file here locally

		//Constant theory files if using the default theory matrices
		files[1] = "all_gene_disease_associations.tsv";
		files[2] = "final_gene_data.txt";
		files[3] = "Data Files/string_parsed.txt";
		/////////////////////////////////////////////////////////////


		int index = 0;
		while(index < args.length)
		{
			try {
				if (args[index].equals("-ns")) {
					ns = Integer.parseInt(args[index + 1]);
					index += 2;
				} else if (args[index].equals("-na")) {
					numAlphas = Integer.parseInt(args[index + 1]);
					index += 2;
				} else if (args[index].equals("-A")) {
					accuracy = Double.parseDouble(args[index + 1]);
					index += 2;
				} else if (args[index].equals("-r")) {
					radius = Double.parseDouble(args[index + 1]);
					index += 2;
				} else if (args[index].equals("-topK")) {
					index++;
					ArrayList<Integer> temp = new ArrayList<Integer>();
					while (!args[index].startsWith("-")) {
						temp.add(Integer.parseInt(args[index]));
						index++;
					}
					Integer[] t = (Integer[]) temp.toArray();
					topK = new int[t.length];
					for (int i = 0; i < t.length; i++)
						topK[i] = t[i];
				} else if (args[index].equals("-D")) {
					index++;
					ArrayList<Integer> temp = new ArrayList<Integer>();
					while (!args[index].startsWith("-")) {
						temp.add(Integer.parseInt(args[index]));
						index++;
					}
					Integer[] t = (Integer[]) temp.toArray();
					diseases = new int[t.length];
					for (int i = 0; i < t.length; i++)
						diseases[i] = t[i];
				} else if (args[index].equals("-normal")) {
					normalizeFiles = true;
					index++;
				} else if (args[index].equals("-t")) {
					target = args[index + 1];
					index += 2;
				} else if (args[index].equals("-data")) {
					dataFile = args[index + 1];
					index += 2;
				} else if (args[index].equals("-in")) {
					intensityFile = args[index + 1];
					index += 2;
				} else if (args[index].equals("-dis")) {
					dissimilarityFile = args[index + 1];
					index += 2;
				} else if (args[index].equals("-out")) {
					outputFile = args[index + 1];
					index += 2;
				} else if (args[index].equals("-outCluster")) {
					clusterFile = args[index + 1];
					index += 2;
				}
				else if(args[index].equals("-graph"))
				{
					graphFile = args[index+1];
					index+=2;
				}
				else if(args[index].equals("-g"))
				{
					threshold = Double.parseDouble(args[index+1]);
					index+=2;
				}
				else if(args[index].equals("-loocv"))
				{
					loocv = true;
					index++;
				}
				else if(args[index].equals("-approx"))
				{
					approxCorrelations = true;
					index++;
				}
			}
			catch(ArrayIndexOutOfBoundsException e)
			{
				System.err.println("No Argument given after: " + args[index]);
				System.exit(-1);
			}
			catch(Exception e)
			{
				System.err.println("Unable to parse argument for: " + args[index]);
				System.exit(-1);
			}
		}

	/*	if(target.equals(""))
		{
			System.err.println("Can't compute intensity values without a target, usage: (-t <target>)");
			System.exit(-1);
		}*/
		if(dataFile.equals(""))
		{
			System.err.println("No data file specified, usage: (-data <filename>)");
			System.exit(-1);
		}
		if(outputFile.equals(""))
		{
			outputFile = "OUT_" + dataFile;
		}

		if(clusterFile.equals(""))
		{
			clusterFile = "CLUST_" + dataFile;
		}

		DataSet data = null;
		try {
			data = MixedUtils.loadDataSet2(dataFile);
		}
		catch(Exception e)
		{
			System.err.println("Unable to load data from: " + dataFile);
			e.printStackTrace();
			System.exit(-1);
		}
		float [] dissimilarity;
		ArrayList<Gene> g;
		files[0] = dataFile;
		if(dissimilarityFile.equals("") || intensityFile.equals(""))//Use Default Theory computations (user only specifies the data file)
		{
			System.out.print("Creating Gene Data...");
			g = Functions.loadGeneData(".",files,diseases,true);
			System.out.println("Done");
			System.out.print("Creating Theory Data...");
			dissimilarity = Functions.createTheoryMatrix("all_dissimilarity_sources.txt");
			System.out.println("Done");
			intensityFile = "all_genes.txt";
			g = Functions.loadGeneData(intensityFile,normalizeFiles);
		}
		else { //Use user specified, files
			System.out.print("Loading Gene Data...");
			g = Functions.loadGeneData(intensityFile, normalizeFiles);
			System.out.println("Done");
			System.out.print("Loading Theory Matrix...");
			dissimilarity = Functions.loadTheoryMatrix(dissimilarityFile,normalizeFiles);
			System.out.println("Done");
		}
		//ArrayList<Gene> g = Functions.loadGeneData(".",files,diseases,true);
		for(int k = 0; k < topK.length;k++) //Loop over sizes of returned results (Top-K)
					{
						System.out.println("Running Pref Div for Top-" + topK[k] + " Genes");
						//Print out to the file for each size of returned results
						PrintStream out = new PrintStream(k + "_" + outputFile);
						PrintStream out2 = new PrintStream(k + "_" + graphFile);
						PrintStream out3 = new PrintStream(k + "_" + clusterFile);
						//Run Pref-Div algorithm to find a diverse and relevant gene set
						RunPrefDiv r = new RunPrefDiv(dissimilarity,g,data,target,loocv);
						r.setNS(ns);
						r.setAccuracy(accuracy);
						r.setNumAlphas(numAlphas);
						r.setRadius(radius);
						r.setTopK(topK[k]);
						r.setThreshold(threshold);
						r.setApproxCorrelations(approxCorrelations);
						Graph causal = r.getCausalGraph(target);
						ArrayList<Gene> Result = r.getLastGeneSet();

						out2.println(causal);
						out2.flush();
						out2.close();
						for(Gene stuff: Result)
						{
							out.println(stuff.symbol);
						}
						for(Gene stuff:r.getClusters().keySet())
						{
							List<Gene> temp = r.getClusters().get(stuff);
							if(temp==null)
								out3.println(stuff.symbol);
							else
								out3.print(stuff.symbol + "\t");
							if(temp!=null) {
								for (int i = 0; i < temp.size(); i++) {
									if (i == temp.size() - 1)
										out3.println(temp.get(i).symbol);
									else
										out3.print(temp.get(i).symbol + "\t");
								}
							}
						}
						out3.flush();
						out3.close();
						out.flush();
						out.close();
					}
	}
}