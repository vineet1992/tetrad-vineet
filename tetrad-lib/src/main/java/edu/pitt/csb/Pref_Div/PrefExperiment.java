package edu.pitt.csb.Pref_Div;
import edu.cmu.tetrad.data.ContinuousVariable;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.regression.LogisticRegression;
import edu.cmu.tetrad.regression.LogisticRegressionResult;
import edu.cmu.tetrad.regression.RegressionResult;
import edu.cmu.tetrad.search.Fci;
import edu.pitt.csb.Pref_Div.Comparisons.ClusterSim;
import edu.pitt.csb.Pref_Div.Comparisons.PrefDivComparator;
import edu.pitt.csb.Priors.runPriors;
import edu.pitt.csb.mgm.MixedUtils;
import edu.pitt.csb.stability.StabilityUtils;
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

	static boolean innerCV = false; //Should we use an inner cross validation to determine which number of features to use?
	static int numCVFolds = 3; //Number of internal cross-validation folds to tune hyperparameters
	static int [] topKs  = new int[]{5,10,25,50}; //Range of top K values to test
	static RunPrefDiv.ClusterType ctype = RunPrefDiv.ClusterType.NONE; //Aggregation type used

	static String runName = "Pref_Div";

	static int maxCat = 3; //Maximum number of categories for a discrete variable

	static int pathwayLimit = 10; //How many genes must be present to include a pathway/gene list?
	public static void main(String [] args) throws Exception
	{


		List<String> varsToRemove = new ArrayList<String>(); //Assume this is variables to remove completely from the data
		List<String> varsToInclude = new ArrayList<String>(); //These variables will be included no matter what (in addition to the top K)
		int index = 0;
		ARGS:while(index < args.length)
		{
			try {
				if (args[index].equals("-ns")) {
					numSamples = Integer.parseInt(args[index + 1]);
					index += 2;
				}

				else if(args[index].equals("-disc"))
				{
					targetContinuous = false;
					index++;
				}
				else if(args[index].equals("-cv"))
				{
					innerCV = true;
					numCVFolds = Integer.parseInt(args[index+1]);
					if((index+2)<args.length && !args[index+2].startsWith("-"))
					{
						String [] ks = args[index+2].split(",");
						topKs = new int[ks.length];
						for(int i = 0; i < ks.length;i++)
						{
							topKs[i] = Integer.parseInt(ks[i]);
						}
					}
					index+=3;

				}
				else if(args[index].equals("-maxCat"))
				{
					maxCat = Integer.parseInt(args[index+1]);
					index+=2;
				}
				else if (args[index].equals("-t")) {
					target = args[index + 1];
					index += 2;
				} else if(args[index].equals("-name"))
				{
					runName = args[index+1];
					index+=2;
				}
				else if (args[index].equals("-data")) {
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
				else if(args[index].equals("-rv"))
				{
					index++;
					while(index < args.length && !args[index].startsWith("-"))
					{
						varsToRemove.add(args[index]);
						index++;
					}
				}
				else if(args[index].equals("-keep"))
				{
					index++;
					while(index < args.length && !args[index].startsWith("-"))
					{
						varsToInclude.add(args[index]);
						index++;
					}
					if(index>=args.length)
						break ARGS;
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
					if(index==args.length)
						break ARGS;
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

		String workingDir = directory + "/" + runName;

		File wDir = new File(workingDir);
		if(!wDir.isDirectory())
		{
			System.out.println("Creating directory for this run: " + wDir.getName());
			if(!wDir.mkdir())
			{
				System.err.println("Unable to create working directory: " + wDir);
				System.exit(-1);
			}
		}
		else
		{
			System.err.println("Working directory already exists, contents will be overwritten");
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
			outputDataset = "/OUT_" + dataFile;
		}

		if(outputClusters.equals(""))
		{
			outputClusters = "/CLUST_" + dataFile;
		}

		if(outputGraph.equals("") && useGraph)
		{
			outputGraph = "GRAPH_" + dataFile;
		}
		outputDataset = workingDir + "/" + outputDataset;
		outputClusters = workingDir + "/" + outputClusters;
		outputGraph = workingDir + "/" + outputGraph;
		dataFile = directory + "/" + dataFile;

		DataSet data = null;
		try {
			data = MixedUtils.loadDataSet2(dataFile,maxCat);
		}
		catch(Exception e)
		{
			System.err.println("Unable to load data from: " + dataFile);
			e.printStackTrace();
			System.exit(-1);
		}


		/**Remove specified variables***/
		for(String s: varsToRemove)
		{
			data.removeColumn(data.getVariable(s));
		}


		DataSet ensured = data;


		/***Keep all categorical variables***/

		System.out.print("Mixed continuous and categorical dataset detected, extracting just the continuous data...");
		int [] indices = MixedUtils.getDiscreteInds(data.getVariables());

		HashSet<Integer> indList = new HashSet<Integer>();
		for(int x = 0; x < indices.length;x++)
			indList.add(indices[x]);
		/***Ensure we keep specified variables***/
		for(String s: varsToInclude)
		{
			indList.add(data.getColumn(data.getVariable(s)));
		}

		/***Don't allow the target variable to be removed***/
		if(indList.contains(data.getColumn(data.getVariable(target))))
			indList.remove(data.getColumn(data.getVariable(target)));

		indices = new int[indList.size()];
		int count  = 0;
		for(Integer x: indList)
		{
			indices[count] = x;
			count++;
		}




		List<Node> varsToKeep = data.subsetColumns(indices).getVariables();
		//varsToKeep.add(data.getVariable(target));

		/***Get these variables into its own dataset and remove them from the PD dataset***/
		ensured = data.subsetColumns(varsToKeep);
		data.removeCols(indices);
		System.out.print("Done");


		List<String> colnames = data.getVariableNames();


		DataSet predictionData = data.copy();
		/***If target is discrete, then reload the dataset with it as a continuous variable for Pref-Div***/
		if(!targetContinuous )
		{
			DiscreteVariable targetVar = (DiscreteVariable) data.getVariable(target);
			int max = targetVar.getCategories().size()-1;
			data = MixedUtils.loadDataSet2(dataFile,max);
		}


		PiPrefDiv4 ppd;
		System.out.print("Initializing radii values in the range: (" + rLow + "," + rHigh + "}...");
		double [] initRadii = new double[numParams];
		for(int i = 0; i < numParams;i++)
		{
			initRadii[i] = rLow + i*(rHigh-rLow)/numParams;
		}
		System.out.println("Done");


		if(innerCV)
		{
			System.out.print("Running internal CV to choose number of variables to select...");
			numSelected = crossValidate(data,predictionData,initRadii,workingDir);
			System.out.println("Done, choosing " + numSelected + " variables");
		}


		ppd = new PiPrefDiv4(data,target,numSelected,initRadii);


		if(boot)
			numSamples = data.getNumRows();
		ppd.setLOOCV(loocv);
		ppd.setVerbose();
		ppd.setClusterType(ctype);
		ppd.setParallel(true);


		if(noPrior)
			System.out.println("Running PiPref-Div without prior information to select genes...");
		else
			System.out.println("Running PiPref-Div to select genes...");

		ArrayList<Gene> output = runPrefDiv(ppd,workingDir);
		Map<Gene,List<Gene>> map = ppd.getLastCluster();
		DataSet summarized = ppd.getSummarizedData();

		/***Don't use continuous target variable for summarized dataset***/
		if(!targetContinuous)
			summarized = RunPrefDiv.summarizeData(predictionData,output,map,ctype,target);

		if(!pathwayFile.equals(""))
		{
			System.out.print("Mapping features to pathways...");
			mapFeaturesToData(summarized,data,pathwayFile);
			System.out.println("Done");
		}

		/***Add back the saved variables***/
		if(varsToInclude.size()>0) {


			for(int i = 0; i < ensured.getNumColumns();i++)
			{
				boolean cont = true;
				if(ensured.getVariable(i) instanceof ContinuousVariable)
				{
					summarized.addVariable(0,new ContinuousVariable(ensured.getVariable(i).getName()));

				}
				else
				{
					cont = false;
					DiscreteVariable temp = (DiscreteVariable)ensured.getVariable(i);
					summarized.addVariable(0,new DiscreteVariable(temp.getName(),temp.getCategories()));
				}

				for(int j = 0; j < ensured.getNumRows();j++)
				{
					if(cont)
					{
						summarized.setDouble(j,0,ensured.getDouble(j,i));
					}else
					{
						summarized.setInt(j,0,ensured.getInt(j,i));
					}
				}
			}
		}


		if(useGraph)
		{
			File temp = new File("temp.txt");
			PrintStream out = new PrintStream(temp);





			if(!summarized.isMixed())
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

				convertPriorsToSif(priorDir,newPriorDir,summarized,colnames);

				String runDirectory = workingDir + "/piMGM";

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

		printOutput(output,map,workingDir);
		PrintStream out = new PrintStream(outputDataset);
		out.println(summarized);
		out.flush();
		out.close();
	}




	public static void outputPriorWeights(String wd, String [] files, double [] weights, double [] tao,double[]p, double[]adjP)
	{
		try{
			PrintStream out = new PrintStream(wd + "/Prior_Weights.txt");
			out.println("Source\tNormalized_Weight\tTao\tP-Value\tAdjusted p-Value");
			for(int i = 0; i < files.length;i++)
			{
				out.println(files[i] + "\t" + weights[i] + "\t" + tao[i] + "\t" + p[i] + "\t" + adjP[i]);
			}
			out.flush();
			out.close();
		}catch (Exception e)
		{
			System.err.println("Unable to print out prior weights");

		}
	}

	private static ArrayList<Gene> runPrefDiv(PiPrefDiv4 ppd,String workingDir)
	{
		ArrayList<Gene> output = new ArrayList<Gene>();
		if(noPrior)
		{

			output = ppd.selectGenes(boot,numSamples);

		}else
		{

			File dDir = new File(directory + "/" + priorDir);
			File [] files = dDir.listFiles();

			String [] diss = new String[files.length];
			for(int i = 0; i < files.length;i++)
			{
				diss[i] = files[i].getAbsolutePath();
			}

			output = ppd.selectGenes(boot,numSamples,diss);
			double [] weights = ppd.getLastWeights();
			double [] tao = ppd.getLastTao();
			double [] pVals = ppd.getPValues();
			double [] apVals = ppd.getAdjustP();
			outputPriorWeights(workingDir,diss,weights,tao,pVals,apVals);
		}
		return output;

	}


	/***
	 * Inner k-fold cross validation to determine the number of genes to select
	 * @param data Dataset to run CV process on
	 * @param predData When target variable is discrete, predData has the real dataset (with discrete target)
	 * @param initRadii The initialRadii parameters to test
	 * @param workingDir The working directory to print results to
	 * @return The best value for the number of genes (in terms of minimizing prediction accuracy in a linear regression)
	 */
	private static int crossValidate(DataSet data,DataSet predData, double [] initRadii,String workingDir)
	{
		double loss = Double.MAX_VALUE;
		int maxIndex = -1;
		int [][] folds = StabilityUtils.generateSubsamples(numCVFolds,data.getNumRows());
		for(int j = 0; j < numCVFolds;j++)
			Arrays.sort(folds[j]);

		for(int i = 0; i < topKs.length;i++)
		{
			int topK = topKs[i];
			double predAccuracy = 0;



			for(int j = 0; j < numCVFolds;j++)
			{
				DataSet testSet = data.subsetRows(folds[j]);
				DataSet trainSet = data.copy();
				trainSet.removeRows(folds[j]);

				PiPrefDiv4 ppd = new PiPrefDiv4(trainSet,target,topK,initRadii);
				ppd.setLOOCV(loocv);
				ppd.setClusterType(ctype);
				ppd.setParallel(true);
				ppd.setQuiet();

				ArrayList<Gene> selected = runPrefDiv(ppd,workingDir);

				testSet = predData.subsetRows(folds[j]);
				trainSet = predData.copy();
				trainSet.removeRows(folds[j]);

				standardizeContinuous(trainSet);
				standardizeContinuous(testSet);

				DataSet summarized = RunPrefDiv.summarizeData(trainSet, selected, ppd.getLastCluster(), ctype,target);
				DataSet summarizedTest = RunPrefDiv.summarizeData(testSet, selected, ppd.getLastCluster(), ctype,target);

				if(targetContinuous)
				{
					PrefDivComparator pdc = new PrefDivComparator(ppd,target,null,null);
					try {
						predAccuracy += pdc.getPredictionAccuracy(summarized, summarizedTest)/(double)numCVFolds;
					}
					catch(Exception e)
					{
						System.out.println(summarized);
						predAccuracy = Double.NaN;
					}
				}
				else
				{

					LogisticRegression lr = new LogisticRegression(summarized);
					List<Node> regressors = summarized.getVariables();
					regressors.remove(summarized.getVariable(target));
					LogisticRegression.Result result = lr.regress((DiscreteVariable)summarized.getVariable(target),regressors);

					predAccuracy += logRegPredict(summarizedTest,result);

				}







			}

			System.out.println("Accuracy at: " + topKs[i] + " is " + predAccuracy);

			if(predAccuracy < loss)
			{
				loss = predAccuracy;
				maxIndex = i;
			}
		}
		return topKs[maxIndex];
	}


	/***
	 * Standardize only the continuous variables in the dataset
	 * @param d a dataset to deal with
	 */
	private static void standardizeContinuous(DataSet d)
	{
		double [][] data = d.getDoubleData().transpose().toArray();
		for(int i = 0; i < d.getNumColumns();i++)
		{
			if(d.getVariable(i)instanceof ContinuousVariable)
			{
				double mean = StatUtils.mean(data[i]);
				double sd = Math.sqrt(StatUtils.variance(data[i],mean));

				for(int j = 0; j <d.getNumRows();j++)
				{
					d.setDouble(j,i,(data[i][j]-mean)/sd);
				}
			}

		}

	}

	/***
	 *
	 * @param test Testing dataset to apply this prediction model to
	 * @param result Result from regressing on the training set
	 * @return The accuracy from these predictions
	 */
	private static double logRegPredict(DataSet test, LogisticRegression.Result result)
	{
		double [] loss = new double[test.getNumRows()];
		for(int i = 0; i < test.getNumRows();i++)
		{
			double logit = result.getIntercept();
			double [] coefs = result.getCoefs();
			List<String> names = result.getRegressorNames();
			for(int j = 0; j < names.size();j++)
			{
				int col = test.getColumn(test.getVariable(names.get(j)));
				logit += coefs[j+1]*test.getDouble(i,col);
			}
			double prob = 1/(1+Math.exp(-1*logit));

			if(prob > 1-Math.pow(10,-15))
				prob = 1-Math.pow(10,-15);
			else if(prob < Math.pow(10,-15))
				prob = Math.pow(10,-15);
			loss[i] = logLoss(prob,test.getInt(i,test.getColumn(test.getVariable(target))));

		}

		return StatUtils.mean(loss);


	}

	private static double logLoss(double p, int y)
	{
		return -1*(y*Math.log(p) +(1-y)*Math.log(1-p));
	}
	/***
	 *
	 * @param priorDir A directory of prior knowledge files (these will be matrices with no labels)
	 * @param newPriorDir A directory to write new prior knowledge files
	 * @param summarized The summarized dataset after selecting features with piMGM
	 * @param colnames Column names of the original dataset
	 */
	private static void convertPriorsToSif(String priorDir, String newPriorDir, DataSet summarized, List<String> colnames)
	{

		List<String> names = summarized.getVariableNames();


		/**Loop through old priorDir and convert each to SIF for the new prior directory***/
		File priors = new File(priorDir);
		File [] oldPriors = priors.listFiles();
		for(int i = 0; i < oldPriors.length;i++)
		{
			System.out.println("Prior: " + oldPriors[i]);
			try {
				BufferedReader b = new BufferedReader(new FileReader(oldPriors[i]));
				PrintStream out = new PrintStream(newPriorDir + "/Prior_" + i + ".sif");
				out.println(oldPriors[i].getName());

				/***Loop through original prior information source***/
				J:
				for (int j = 0; j < colnames.size(); j++) {
					String[] currLine = b.readLine().split("\t");
					String var1 = colnames.get(j);

					int index1 = -1;
					for (int x = 0; x < names.size(); x++) {
						if (names.get(x).split("\\|")[0].contains(var1))
							index1 = x;
					}


					if (index1 == -1)
						continue J;

					K:
					for (int k = j + 1; k < colnames.size(); k++) {
						String var2 = colnames.get(k);


						int index2 = -1;
						for (int x = 0; x < names.size(); x++) {
							if (names.get(x).split("\\|")[0].contains(var2))
								index2 = x;
						}

						if (index2 == -1)
							continue K;

						/***Output any non-NULL prior info to SIF file***/
						if (Double.parseDouble(currLine[k]) > 0)
						{
							out.println(names.get(index1) + "\t" + names.get(index2) + "\t" + Double.parseDouble(currLine[k]));
							System.out.println("Connecting " + names.get(index1) + ", " + names.get(index2) + " score: " + Double.parseDouble(currLine[k]));
						}


					}
				}
				b.close();
				out.flush();
				out.close();
			}
			catch(Exception e)
			{
				e.printStackTrace();
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
						if(!list.get(j).symbol.equals(top.get(i).symbol))
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