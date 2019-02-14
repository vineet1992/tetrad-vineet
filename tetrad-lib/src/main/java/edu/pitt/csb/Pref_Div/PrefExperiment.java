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
	/****
	 * The purpose of this class is to create a packgeable jar file to smoothly run Prior Information Pref-Div with prior matrices
	 * The input is a data matrix where rows are samples and columns are genes to be selected.
	 * In addition, the algorithm expects a set of files in a specified directory to be prior knowledge (probability of each variable being related to one another)
	 * Lastly, user must specify the number of genes to choose, a target variable of interest, and whether or not to run causal modeling as a post-processing step
	   The output is either a reduced dataset with only the variables specified, or a graphical model of the reduced variables

	 	Author: Vineet Raghu (vineet@cs.pitt.edu)
	 */



	static boolean boot = false; //Should we use bootstrap samples for PiPrefDiv
	static boolean loocv = false; //Should we use leave-one-out CV for PiPrefDiv
	static int numSamples = 20; //Number of bootstrap/sub-sampled samples to use
	static int numParams = 30;//Number of parameters to sweep over

	static int numSelected = 75; //How many true parents of the target are there?
	static boolean targetContinuous = true; //Is the target variable continuous?
	static boolean partialCorr = false;
	static boolean pdStability = false;
	static boolean noPrior = true;

	static double rLow = 0.3;
	static double rHigh = 1.0;
	static boolean useGraph = false;
	static String dataFile = "";
	static String directory = ".";
	static String priorDir = "";

	static String target = "";

	static String outputDataset = "";
	static String outputClusters = "";
	static String outputGraph = "";

	static boolean runPiMGM = false;

	public static void main(String [] args) throws Exception
	{


		int index = 0;
		while(index < args.length)
		{
			try {
				if (args[index].equals("-ns")) {
					numSamples = Integer.parseInt(args[index + 1]);
					index += 2;
				}  else if (args[index].equals("-t")) {
					target = args[index + 1];
					index += 2;
				} else if (args[index].equals("-data")) {
					dataFile = args[index + 1];
					index += 2;
				} else if (args[index].equals("-outData")) {
					outputDataset = args[index + 1];
					index += 2;
				} else if (args[index].equals("-outCluster")) {
					outputClusters = args[index + 1];
					index += 2;
				}else if (args[index].equals("-outGraph")) {
					outputGraph = args[index + 1];
					index += 2;
				}
				else if(args[index].equals("-numParams"))
				{
					numParams = Integer.parseInt(args[index+1]);
					index+=2;
				}

				else if(args[index].equals("-numSelect"))
				{
					numSelected = Integer.parseInt(args[index+1]);
					index+=2;
				}
				else if(args[index].equals("-loocv"))
				{
					loocv = true;
					index++;
				}
				else if(args[index].equals("-paramRange"))
				{
					rLow = Integer.parseInt(args[index+1]);
					rHigh = Integer.parseInt(args[index+2]);
					index+=3;
				}
				else if(args[index].equals("-priors"))
				{
					priorDir = args[index+1];
					noPrior = false;
					index+=2;
				}
				else if(args[index].equals("-dir"))
				{
					directory = args[index+1];
					index+=2;
				}
				else if(args[index].equals("-boot"))
				{
					boot = true;
					index++;
				}
				else if(args[index].equals("-useCausalGraph"))
				{
					useGraph = true;
					if(args.length==index+1)
					{

					}
					else if(!args[index+1].startsWith("-"))
					{
						if(args[index+1].equals("piMGM"))
						{
							runPiMGM = true;

						}
						index+=2;
					}
					else
						index++;
				}
				else if(args[index].equals("-partialCorr"))
				{
					partialCorr = true;
					index++;
				}
				else if(args[index].equals("-pdStability"))
				{
					pdStability = true;
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

		if(dataFile.equals(""))
		{
			System.err.println("No data file specified, usage: (-data <filename>)");
			System.exit(-1);
		}
		if(priorDir.equals(""))
		{
			System.out.println("No prior knowledge directory specified, running PiPrefDiv no Prior version");
		}

		if(outputDataset.equals(""))
		{
			outputDataset = "OUT_" + dataFile;
		}

		if(outputClusters.equals(""))
		{
			outputClusters = "CLUST_" + dataFile;
		}

		if(outputGraph.equals("") && useGraph)
		{
			outputGraph = "GRAPH_" + dataFile;
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
		PiPrefDiv4 ppd;
		System.out.print("Initializing radii values in the range: (" + rLow + "," + rHigh + "}...");
		double [] initRadii = new double[numParams];
		for(int i = 0; i < numParams;i++)
		{
			initRadii[i] = rLow + i*(rHigh-rLow)/numParams;
		}
		System.out.println("Done");


			ppd = new PiPrefDiv4(data,target,numSelected,initRadii);


		if(boot)
			numSamples = data.getNumRows();
		ppd.setLOOCV(loocv);
		ppd.setVerbose();

		ArrayList<Gene> output;
		Map<Gene,List<Gene>> map;
		DataSet summarized;
		if(noPrior)
		{
			System.out.println("Running PiPref-Div without prior information to select genes...");

			output = ppd.selectGenes(boot,numSamples);

		}else
		{
			System.out.println("Running PiPref-Div to select genes...");

			File dDir = new File(directory + "/" + priorDir);
			File [] files = dDir.listFiles();

			String [] diss = new String[files.length];
			for(int i = 0; i < files.length;i++)
			{
				diss[i] = files[i].getAbsolutePath();
			}
			output = ppd.selectGenes(boot,numSamples,diss);
		}

		map = ppd.getLastCluster();
		summarized = ppd.getSummarizedData();

		if(useGraph)
		{
			File temp = new File("temp.txt");
			PrintStream out = new PrintStream(temp);
			out.println(summarized);
			out.flush();
			out.close();
			if(runPiMGM)
			{
				//TODO figure out how to summarize the prior information sources to match the summarized dataset (same options as PiPrefDiv)
			}
			else
			{
				System.out.print("Running StEPS to learn a causal model...");

				String cmd = "java -jar causalDiscovery.jar -d temp.txt -mgm -steps -o " + outputGraph + " -maxCat 3";
				Runtime.getRuntime().exec(cmd);

				temp.deleteOnExit();
				System.out.println("Done");
			}

		}

		System.out.println(map);

		printOutput(output,map,directory);
		PrintStream out = new PrintStream(outputDataset);
		out.println(summarized);
		out.flush();
		out.close();
	}


	public static void printOutput(ArrayList<Gene> top, Map<Gene,List<Gene>> map, String currDir)
	{
		try {
			PrintStream out = new PrintStream(currDir + "/selected_genes.txt");
			for (int i = 0; i < top.size(); i++) {
				out.print(top.get(i) + "\t");
				List<Gene> list = map.get(top.get(i));
				if(list!=null) {
					for (int j = 0; j < list.size(); j++) {
						out.print(list.get(j) + "\t");

					}
				}
				out.println();

			}
			out.flush();
			out.close();
		}catch(Exception e)
		{
			System.err.println("Couldn't write output to file...printing error message");
			e.printStackTrace();
		}
	}
}