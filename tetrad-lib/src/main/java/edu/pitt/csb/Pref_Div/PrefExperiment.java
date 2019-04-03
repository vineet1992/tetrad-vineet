package edu.pitt.csb.Pref_Div;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.Fci;
import edu.pitt.csb.Pref_Div.Comparisons.ClusterSim;
import edu.pitt.csb.Priors.runPriors;
import edu.pitt.csb.mgm.MixedUtils;
import org.apache.commons.math3.ml.clustering.Cluster;
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
	static boolean partialCorr = false; //Should we use partial correlation instead of pearson correlation (doesn't really work right now)
	static boolean pdStability = false; //Should we run Pref-Div within the subsamples to get stability?
	static boolean noPrior = true; //Should we use prior information at all or just stability?

	static double rLow = 0.3; //Low value for the radius
	static double rHigh = 1.0; //High value for the radius
	static boolean useGraph = false; //Should we use causal modeling?
	static String dataFile = ""; //Data file
	static String directory = "."; //Home directory for the run
	static String priorDir = ""; //Prior knowledge directory

	static String target = ""; //Target variable name

	static String outputDataset = ""; //Should we output the summarized dataset?
	static String outputClusters = ""; //Should we output clusters?
	static String outputGraph = ""; //Should we output the causal graph?
	static String pathwayFile = ""; //Should we attempt to map aggregations to gene lists?

	static boolean runPiMGM = false; //Should we run piMGM to incorporate prior knowledge in undirected modeling (only if aggregation type is NONE)

	static RunPrefDiv.ClusterType ctype; //Aggregation type used


	static int pathwayLimit = 10; //How many genes must be present to include a pathway/gene list?
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
					if(args.length!=index+1 && !args[index+1].startsWith("-"))
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
				else if(args[index].equals("-ctype"))
				{
					if(args[index+1].toLowerCase().equals("none"))
					{
						ctype = RunPrefDiv.ClusterType.NONE;
					}else if(args[index+1].toLowerCase().equals("pca"))
					{
						ctype = RunPrefDiv.ClusterType.PCA;

					}else if(args[index+1].toLowerCase().equals("median"))
					{
						ctype = RunPrefDiv.ClusterType.MEDIAN;

					}else if(args[index+1].toLowerCase().equals("mean"))
					{
						ctype = RunPrefDiv.ClusterType.MEAN;
					}else
					{
						System.err.println("Unrecognized aggregation option...using None\nValid types are: PCA, MEDIAN, MEAN, NONE");
					}
					index+=2;
				}
				else if(args[index].equals("-pdStability"))
				{
					pdStability = true;
					index++;
				}
				else if(args[index].equals("-pathwayFile"))
				{
					pathwayFile = args[index+1];
					index+=2;
				}
				else
				{
					System.err.println("Unrecognized command line option " + args[index] + " ... ignoring");
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

		if(data.isMixed())
		{
			System.out.print("Mixed continuous and categorical dataset detected, extracting just the continuous data...");
			int [] indices = MixedUtils.getContinuousInds(data.getVariables());
			List<Node> varsToKeep = data.subsetColumns(indices).getVariables();
			varsToKeep.add(data.getVariable(target));
			data = data.subsetColumns(varsToKeep);
			System.out.print("Done");
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
		ppd.setClusterType(ctype);

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

		if(!pathwayFile.equals(""))
		{
			System.out.print("Mapping features to pathways...");
			mapFeaturesToData(summarized,data,pathwayFile);
			System.out.println("Done");
		}

		if(useGraph)
		{
			File temp = new File("temp.txt");
			PrintStream out = new PrintStream(temp);
			if(data.isMixed())
			{
				/***Add back the discrete variables that were removed***/
				DataSet discData = MixedUtils.getDiscreteData(data);
				summarized = DataUtils.concatenate(summarized,discData);

			}
			else
			{
				/***Add a dummy variable to run causal modeling ***/
				summarized = runPriors.addDummy(summarized);

			}
			out.println(summarized);
			out.flush();
			out.close();
			String cmd = "";
			if(runPiMGM)
			{

				System.out.print("Running piMGM to learn a causal model...");

				if(ctype!= RunPrefDiv.ClusterType.NONE)
				{
					System.err.println("Using single variable priors on aggregate features.. could get weird results. If using piMGM, then use None clustertype");
				}

				String newPriorDir = directory + "/" + priorDir + "_SIF";
				File npDir = new File(newPriorDir);
				if(!npDir.isDirectory())
					npDir.mkdir();

				convertPriorsToSif(priorDir,newPriorDir,summarized,data);

				String runDirectory = directory + "/piMGM";

				cmd = "java -jar runPriors.jar -run " + runDirectory + " -priors " + newPriorDir + " -makeScores -fullCounts -data temp.txt -sif";
			}
			else
			{
				System.out.print("Running StEPS to learn a causal model...");
				cmd = "java -jar causalDiscovery.jar -d temp.txt -mgm -steps -o " + outputGraph + " -maxCat 3";

			}

			Runtime rt = Runtime.getRuntime();

			Process proc = rt.exec(cmd);
			BufferedReader stdInput = new BufferedReader(new
					InputStreamReader(proc.getInputStream()));
			BufferedReader stdError = new BufferedReader(new
					InputStreamReader(proc.getErrorStream()));

			// read the output from the command
			String s = null;
			while ((s = stdInput.readLine()) != null) {
				System.out.println(s);
			}

			// read any errors from the attempted command
			while ((s = stdError.readLine()) != null) {
				System.err.println(s);
			}


			temp.deleteOnExit();
			System.out.println("Done");

		}

		System.out.println(map);

		printOutput(output,map,directory);
		PrintStream out = new PrintStream(outputDataset);
		out.println(summarized);
		out.flush();
		out.close();
	}


	/***
	 *
	 * @param priorDir A directory of prior knowledge files (these will be matrices with no labels)
	 * @param newPriorDir A directory to write new prior knowledge files
	 * @param summarized The summarized dataset after selecting features with piMGM
	 * @param data The original dataset
	 */
	private static void convertPriorsToSif(String priorDir, String newPriorDir, DataSet summarized, DataSet data)
	{

		List<String> names = summarized.getVariableNames();


		/**Loop through old priorDir and convert each to SIF for the new prior directory***/
		File priors = new File(priorDir);
		File [] oldPriors = priors.listFiles();
		for(int i = 0; i < oldPriors.length;i++)
		{
			try {
				BufferedReader b = new BufferedReader(new FileReader(oldPriors[i]));
				PrintStream out = new PrintStream(newPriorDir + "/Prior_" + i + ".sif");
				out.println(oldPriors[i].getName());

				/***Loop through original prior information source***/
				J:
				for (int j = 0; j < data.getNumColumns(); j++) {
					String[] currLine = b.readLine().split("\t");
					String var1 = data.getVariable(j).getName();

					int index1 = -1;
					for (int x = 0; x < names.size(); x++) {
						if (names.get(x).split("\\|")[0].contains(var1))
							index1 = x;
					}


					if (index1 == -1)
						continue J;

					K:
					for (int k = j + 1; k < data.getNumColumns(); k++) {
						String var2 = data.getVariable(k).getName();


						int index2 = -1;
						for (int x = 0; x < names.size(); x++) {
							if (names.get(x).split("\\|")[0].contains(var2))
								index2 = x;
						}

						if (index2 == -1)
							continue K;

						/***Output any non-NULL prior info to SIF file***/
						if (Double.parseDouble(currLine[k]) > 0)
							out.println(names.get(index1) + "\t" + names.get(index2) + "\t" + Double.parseDouble(currLine[k]));


					}
				}
				b.close();
				out.flush();
				out.close();
			}
			catch(Exception e)
			{
				System.err.println("Could not generate SIF prior " + i + ", need to debug!");
			}
		}

	}


	/****This function creates a file giving the best mapping from aggregate features to pathways based on a user-specified pathway file***/

	/***Pathway file consists of pathways in the rows, with the first column containing the name of the pathway***/
	public static void mapFeaturesToData(DataSet data,DataSet fullSet, String pathwayFile)
	{
		List<String> varNames = fullSet.getVariableNames();
		ArrayList<List<Gene>> features = new ArrayList<List<Gene>>();


		/***Load features from the variable names of the data***/
		for(int i = 0; i < features.size();i++)
		{
			ArrayList<Gene> temp = new ArrayList<Gene>();
			String [] names = data.getVariable(i).getName().split("|");
			for(String s: names)
			{
				Gene x= new Gene(0);
				x.symbol = s;
				temp.add(x);
			}
			features.add(temp);
		}
		ArrayList<List<Gene>> pathways = new ArrayList<List<Gene>>();


		/***List with only pathways that we included in the analysis(at least one gene was present in the data)***/
		List<String> keptNames = new ArrayList<String>();
		try{
			BufferedReader b = new BufferedReader(new FileReader(pathwayFile));
			while(b.ready())
			{
				String [] line = b.readLine().split("\t");
				ArrayList<Gene> currPath = new ArrayList<Gene>();
				for(int i = 1; i < line.length;i++)
				{
					Gene x = new Gene(0);
					x.symbol = b.readLine();

					if(varNames.contains(x.symbol))
						currPath.add(x);

				}


				if(currPath.size()>=pathwayLimit) {
					pathways.add(currPath);
					keptNames.add(line[0]);
				}
			}
			b.close();

		}
		catch(Exception e) {
			System.err.println("Could not load pathways from file: " + pathwayFile);
			return;
		}



		/***Load features from the pathways, but only include genes that were in the original dataset***/



		double[][] costs = ClusterSim.computeCostMatrix(features,pathways);

		try {


			/***Write out top 5 pathways and costs to file for each feature***/

			PrintStream out = new PrintStream("Pathway_Information.txt");

			out.print("Feature\tNumber\tPathway_Name\tScore");
			for (int i = 0; i < costs.length; i++) //For each feature
			{

				/***Create temporary index array***/
				Integer[] inds = new Integer[pathways.size()];
				for (int j = 0; j < inds.length; j++)
					inds[j] = j;
				final double[] currCosts = costs[i];


				/***Extract top 5 cost indices***/
				Arrays.sort(inds, new Comparator<Integer>() {
					@Override
					public int compare(final Integer o1, final Integer o2) {
						return Double.compare(currCosts[o1], currCosts[o2]);
					}
				});

				for (int j = 0; j < 5; j++) {
					out.println(data.getVariable(i) + "\t" + j + "\t" + keptNames.get(j) + "\t" + costs[j]);

				}
			}
			out.flush();
			out.close();
		}catch(Exception e)
		{
			System.err.println("Unable to write pathway information to the file");
			e.printStackTrace();
		}

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