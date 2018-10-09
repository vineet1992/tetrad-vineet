package edu.pitt.csb.Olja_Cancer_Analysis;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import edu.cmu.tetrad.data.DataSet;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by vinee_000 on 7/27/2018.
 */
public class constructPrior {
    private static String priorDir = "C:/Users/vinee_000/Documents/CS Academic Stuff/Graduate/Current Projects/Olja_Cancer_Trials/RNA-Seq-Data/piMGM Analysis/Prediction_RNA_Seq/";
    public static void main(String[] args) throws Exception {


        boolean mSig = true;
        boolean string = true;
        boolean geneList = false;
        String dataset = "All_Data.txt";

        DataSet data = MixedUtils.loadDataSet2(dataset);
        HashMap<String,Integer> map = null;
        DoubleMatrix2D prior = null;
        if(string) {
            //CONSTRUCT STRING PRIOR
            System.out.print("Reading gene protein map...");
            HashMap<String, String> gp = readGeneProtein(); //Maps from ensembl_protein_id to gene symbol
            System.out.println("Done");
            System.out.print("Getting scores from STRING...");
            map = getScores(gp);//maps from Gene,Gene to score
            System.out.println("Done");
            gp = null;


            System.out.print("Constructing prior from map and data...");
            prior = makePrior(map, data);
            System.out.println("Done");
            outputPrior(prior, data, "string_prior_Disc.txt");
        }

        //CONSTRUCT MSIGDB PRIOR
        if(mSig) {
            System.out.print("Constructing MSig HashMap...");
            HashMap<String, Double> map2 = readMSig();
            System.out.println("Done");
            System.out.print("Constructing prior for MSig...");
            prior = makePriorMSig(map2, data);
            System.out.println("Done");
            outputPrior(prior, data, "M_Sig_prior_Disc.txt");
        }
if(geneList) {
    System.out.print("Constructing prior for gene list...");
    prior = makePriorGeneList(data);
    System.out.println("Done");
    outputPrior(prior, data, "Gene_List_Disc.txt");

}

    }

    public static DoubleMatrix2D makePriorGeneList(DataSet data)throws Exception
    {
        BufferedReader b = new BufferedReader(new FileReader(priorDir + "gene_list.txt"));
        ArrayList<String> genes = new ArrayList<String>();
        while(b.ready())
            genes.add(b.readLine());
        DoubleMatrix2D prior = new SparseDoubleMatrix2D(data.getNumColumns(),data.getNumColumns());
        int col = -1;
        if(data.getVariable("IgG_Ratio")!=null) {
            col = data.getColumn(data.getVariable("IgG_Ratio"));
        }
        else
        {
            col = data.getColumn(data.getVariable("Response"));
        }
            for (int i = 0; i < data.getNumColumns(); i++) {
                if (genes.contains(data.getVariable(i).getName()))
                    prior.set(i, col, 1);
            }

        return prior;
    }
    public static DoubleMatrix2D makePriorMSig(HashMap<String,Double> map, DataSet data)
    {
        int total = 0;
        int found = 0;
        DoubleMatrix2D prior = new SparseDoubleMatrix2D(data.getNumColumns(),data.getNumColumns());
        for(int i = 0; i < data.getNumColumns();i++)
        {
            for(int j = i+1; j < data.getNumColumns();j++)
            {
                total++;
                if(map.get(data.getVariable(i).getName()+","+data.getVariable(j).getName())==null && map.get(data.getVariable(j).getName()+","+data.getVariable(i).getName())==null)
                    continue;
                else
                {
                    if(map.get(data.getVariable(i).getName()+","+data.getVariable(j).getName())!=null)
                    {
                        found++;
                        double y = map.get(data.getVariable(i).getName()+","+data.getVariable(j).getName());
                        prior.set(i,j,y);
                        prior.set(j,i,y);
                    }
                    else
                    {
                        found++;
                        double y = map.get(data.getVariable(j).getName()+","+data.getVariable(i).getName());
                        prior.set(i,j,y);
                        prior.set(j,i,y);
                    }
                }

            }
        }
        System.out.print("Found " + found + " priors out of " + total + " possible");
        return prior;
    }
    public static HashMap<String, Double> readMSig()throws Exception
    {
        HashMap<String,Double> result = new HashMap<String,Double>();
        BufferedReader b = new BufferedReader(new FileReader(priorDir + "gene_similarity_matrix_cosine.txt"));
        String [] names = b.readLine().split("\t");
        int row = 0;
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            for(int i = 3; i < line.length;i++)
            {
                if((i-3)!=row && Double.parseDouble(line[i])!=0)
                {
                    result.put(names[row] + "," + names[i-3],Double.parseDouble(line[i]));
                }
            }
            row++;
        }
        b.close();
        return result;
    }
    public static void outputPrior(DoubleMatrix2D prior,DataSet data, String file) throws Exception
    {
        PrintStream out = new PrintStream(file);
        for(int i = 0; i < data.getNumColumns();i++)
        {
            if(i==data.getNumColumns()-1)
                out.println(data.getVariable(i).getName());
            else
            out.print(data.getVariable(i).getName() + "\t");
        }
        for(int i = 0; i < prior.rows();i++)
        {
            for(int j = 0; j < prior.columns();j++)
            {
                if(j==prior.columns()-1)
                    out.println(prior.get(i,j));
                else
                    out.print(prior.get(i,j)+"\t");
            }
        }
        out.flush();
        out.close();
    }
    public static DoubleMatrix2D makePrior(Map<String,Integer> map, DataSet data)
    {
        int total = 0;
        int found = 0;
        DoubleMatrix2D prior = new SparseDoubleMatrix2D(data.getNumColumns(),data.getNumColumns());
        for(int i = 0; i < data.getNumColumns();i++)
        {
            for(int j = i+1; j < data.getNumColumns();j++)
            {
                total++;
                if(map.get(data.getVariable(i).getName()+","+data.getVariable(j).getName())==null && map.get(data.getVariable(j).getName()+","+data.getVariable(i).getName())==null)
                    continue;
                else
                {
                    if(map.get(data.getVariable(i).getName()+","+data.getVariable(j).getName())!=null)
                    {
                        found++;
                        int x = map.get(data.getVariable(i).getName()+","+data.getVariable(j).getName());
                        double y = (double)x/1000.0;
                        prior.set(i,j,y);
                        prior.set(j,i,y);
                    }
                    else
                    {
                        found++;
                        int x = map.get(data.getVariable(j).getName()+","+data.getVariable(i).getName());
                        double y = (double)x/1000.0;
                        prior.set(i,j,y);
                        prior.set(j,i,y);
                    }
                }

            }
        }
        System.out.print("Found " + found + " priors out of " + total + " possible");
        return prior;
    }
    public static HashMap<String,Integer> getScores(Map<String,String> gp)throws Exception
    {
        HashMap<String,Integer> map = new HashMap<String,Integer>();
        BufferedReader b = new BufferedReader(new FileReader(priorDir + "9606.protein.links.v10.5.txt"));
        b.readLine();
        int total = 0;
        int found = 0;
        while(b.ready())
        {
            String [] line = b.readLine().split(" ");
            line[0] = line[0].replace("9606.","");
            line[1] = line[1].replace("9606.","");
            if(gp.get(line[0])!=null && gp.get(line[1])!=null)
            {
                map.put(gp.get(line[0])+","+gp.get(line[1]),Integer.parseInt(line[2]));
                found++;
            }
            total++;
        }
        System.out.print("Found ID's for " + found + " out of " + total + " genes");
        b.close();
        return map;
    }
    public static HashMap<String,String> readGeneProtein() throws Exception
    {
        HashMap<String,String> map = new HashMap<String,String>();
        BufferedReader b = new BufferedReader(new FileReader(priorDir + "gene_protein_mapping.txt"));
        b.readLine();
        while(b.ready())
        {
            String[]line = b.readLine().split("\t");
            if(line[1].length()>2)
                map.put(line[1],line[0]);
        }
        b.close();
        return map;
    }

}
