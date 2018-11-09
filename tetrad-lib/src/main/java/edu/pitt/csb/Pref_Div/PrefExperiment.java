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

	//TODO adapt this to use PiPrefDiv and those expected input files
	public static void main(String [] args) throws Exception
	{
		//PARAMS

		int numAlphas = 20;
		int ns = 20;
		double threshold = 0.1;
		int numParams = -1;
		int numRadii = 10;
		int K = 50;
		double rLow = 0.1;
		double rHigh = 0.8;
		double tLow = 0.0001;
		double tHigh = 0.2;
		int numThreshold = 10;
		boolean normalizeFiles = false;
		boolean loocv = false;
		boolean bootstrap = false;
		boolean approxCorrelations = true;
		boolean useGraph = false;
		String dataFile = "";
		String intensityFile = "";
		String directory = ".";
		String dissimilarityDir = "";


		String outputFile = "";
		String clusterFile = "";
		String graphFile = "";
		String target = "";

		int index = 0;
		while(index < args.length)
		{
			try {
				if (args[index].equals("-ns")) {
					ns = Integer.parseInt(args[index + 1]);
					index += 2;
				}  else if (args[index].equals("-t")) {
					target = args[index + 1];
					index += 2;
				} else if (args[index].equals("-data")) {
					dataFile = args[index + 1];
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
				else if(args[index].equals("-numRadii"))
				{
					numRadii = Integer.parseInt(args[index+1]);
					index+=2;
				}
				else if(args[index].equals("-numThreshold"))
				{
					numThreshold = Integer.parseInt(args[index+1]);
					index+=2;
				}
				else if(args[index].equals("-numParams"))
				{
					numParams = Integer.parseInt(args[index+1]);
					index+=2;
				}
				else if(args[index].equals("-g"))
				{
					threshold = Double.parseDouble(args[index+1]);
					index+=2;
				}
				else if(args[index].equals("-K"))
				{
					K = Integer.parseInt(args[index+1]);
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
				else if(args[index].equals("-tRange"))
				{
					tLow = Integer.parseInt(args[index+1]);
					tHigh = Integer.parseInt(args[index+2]);
					index+=3;
				}
				else if(args[index].equals("-rRange"))
				{
					rLow = Double.parseDouble(args[index+1]);
					rHigh = Double.parseDouble(args[index+2]);
					index+=3;
				}
				else if(args[index].equals("-dissim"))
				{
					dissimilarityDir = args[index+1];
					index+=2;
				}
				else if(args[index].equals("-dir"))
				{
					directory = args[index+1];
					index+=2;
				}
				else if(args[index].equals("-boot"))
				{
					bootstrap = true;
					index++;
				}
				else if(args[index].equals("-useCausalGraph"))
				{
					useGraph = true;
					index++;
				}
				else if(args[index].equals("-iFile"))
				{
					intensityFile = args[index+1];
					index+=2;
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
		if(dissimilarityDir.equals(""))
		{
			System.err.println("No dissimilarity directory specified, please use -dissim <Dissimilarity Directory>\nAll dissimilarity prior knowledge files should be included in this directory");
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
		PiPrefDiv ppd;
		double [] initRadii = new double[numRadii];
		double [] initThreshold = new double[numThreshold];
		for(int i = 0; i < numRadii;i++)
		{
			initRadii[i] = rLow + i*(rHigh-rLow)/numRadii;
		}
		for(int i = 0; i < numThreshold;i++)
		{
			initThreshold[i] = (tLow + i*(tHigh-tLow)/(double)numThreshold);
		}


			ppd = new PiPrefDiv(data,target,K,initRadii,initThreshold);


		if(bootstrap)
			ns = data.getNumRows();
		ppd.setLOOCV(loocv);
		ppd.setVerbose();

		File dDir = new File(directory + "/" + dissimilarityDir);
		File [] files = dDir.listFiles();

		String [] diss = new String[files.length];
		for(int i = 0; i < files.length;i++)
		{
			diss[i] = files[i].getAbsolutePath();
		}
		ArrayList<Gene> output = ppd.selectGenes(bootstrap,ns,intensityFile,diss,useGraph);
		Map<Gene,List<Gene>> map = ppd.getLastCluster();

		System.out.println(map);

		printOutput(output,map,directory);


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