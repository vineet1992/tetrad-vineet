package edu.pitt.csb.Priors;
import java.util.*;
import java.io.*;

public class addClinicalData
{
    public static void main(String [] args) throws Exception
    {
        String geneFileName = "";
        boolean useGenes = true;
        if(geneFileName.equals(""))
            useGenes = false;
        Set<String> genesToUse = new HashSet<String>();
        if(useGenes) {
            BufferedReader b2 = new BufferedReader(new FileReader(geneFileName));

            while (b2.ready()) {
                genesToUse.add(b2.readLine());
            }
            b2.close();
        }
        //ER_Positive -> Subtype = 0 or 1
        //ER_Negative -> Subtype = 2 or 3
        PrintStream erPos = new PrintStream("ER_Positive.txt");
        PrintStream erNeg = new PrintStream("ER_Negative.txt");
        PrintStream [] types = new PrintStream [4];
        String [] typeString = {"Luminal A","Luminal B", "Triple-Negative", "HER2"};
        for(int i = 0; i < types.length;i++)
        {
            //Luminal A, Luminal B, Triple Negative, HER2
            types[i] = new PrintStream(typeString[i] + ".txt");
        }
        PrintStream out = new PrintStream("genes_with_clinical.txt");
        PrintStream out2 = new PrintStream("genes_with_clinical_normals.txt");
        out.println("/variables");
        out2.println("/variables");
        erPos.println("/variables");
        erNeg.println("/variables");
        for(int i = 0; i < types.length;i++)
            types[i].println("/variables");
        BufferedReader b = new BufferedReader(new FileReader("expression_subset.txt"));
        String [] samples = b.readLine().split("\t");
        int geneCount = 0;
        while(b.ready())
        {
            b.readLine();
            geneCount++;
        }
        System.out.println(geneCount);
        b.close();
        b = new BufferedReader(new FileReader("expression_subset.txt"));
        b.readLine();
        System.out.println(Arrays.toString(samples));
        double [][] data;
        if(useGenes)
            data = new double[genesToUse.size()+6][samples.length];
        else
            data = new double[geneCount][samples.length];
        ArrayList<String> idk = new ArrayList<String>();
        int count = 0;
        while(b.ready())
        {
            String [] stuf = b.readLine().replace("\"","").split("\t");
            String gene = stuf[0];
            //System.out.println(Arrays.toString(stuf));

            if(!useGenes || genesToUse.contains(gene))
            {
                for(int i = 1; i < stuf.length;i++)
                {
                    data[count][i-1] = Double.parseDouble(stuf[i]);
                }
                idk.add(gene);
                out.println(gene + ":Continuous");
                out2.println(gene + ":Continuous");
                erPos.println(gene + ":Continuous");
                erNeg.println(gene + ":Continuous");
                for(int xx = 0; xx < types.length;xx++)
                    types[xx].println(gene + ":Continuous");
                count++;
            }

        }
        b.close();

        out.println("Tumor_Stage:Continuous");//need to determine this mapping from the spreadsheet i*, ii*, iii*, iv
        out.println("Vital_Status:0 1");// 0 = dead
        //out.println("Gender:0 1"); //0 = female
        out.println("Age:Continuous");
        //out.println("Tumor:0 1");
        //out.println("Race:0 1 2 3");//Need to determine this mapping as well 0 = white, 1 = black, 2 = asian, 3 = american indian
        out.println("Subtype:0 1 2 3");
        //4 = not reported
        out.println("/data");

        out2.println("Tumor_Stage:Continuous");//need to determine this mapping from the spreadsheet i*, ii*, iii*, iv
        out2.println("Vital_Status:0 1");// 0 = dead
        //out.println("Gender:0 1"); //0 = female
        out2.println("Age:Continuous");
        //out.println("Tumor:0 1");
        //out.println("Race:0 1 2 3");//Need to determine this mapping as well 0 = white, 1 = black, 2 = asian, 3 = american indian
        out2.println("Subtype:0 1 2 3");
        //4 = not reported
        out2.println("/data");

        erPos.println("Tumor_Stage:Continuous");//need to determine this mapping from the spreadsheet i*, ii*, iii*, iv
        erPos.println("Vital_Status:0 1");// 0 = dead
        //out.println("Gender:0 1"); //0 = female
        erPos.println("Age:Continuous");
        //out.println("Tumor:0 1");
        //out.println("Race:0 1 2 3");//Need to determine this mapping as well 0 = white, 1 = black, 2 = asian, 3 = american indian
        erPos.println("Subtype:0 1");
        //4 = not reported
        erPos.println("/data");

        erNeg.println("Tumor_Stage:Continuous");//need to determine this mapping from the spreadsheet i*, ii*, iii*, iv
        erNeg.println("Vital_Status:0 1");// 0 = dead
        //out.println("Gender:0 1"); //0 = female
        erNeg.println("Age:Continuous");
        //out.println("Tumor:0 1");
        //out.println("Race:0 1 2 3");//Need to determine this mapping as well 0 = white, 1 = black, 2 = asian, 3 = american indian
        erNeg.println("Subtype:2 3");
        //4 = not reported
        erNeg.println("/data");

        for(int i = 0; i < types.length;i++)
        {
            types[i].println("Tumor_Stage:Continuous");//need to determine this mapping from the spreadsheet i*, ii*, iii*, iv
            types[i].println("Vital_Status:0 1");// 0 = dead
            //out.println("Gender:0 1"); //0 = female
            types[i].println("Age:Continuous");
            types[i].println("/data");
        }

        for(int i = 0; i < idk.size();i++)
        {
            out.print(idk.get(i) + "\t");
            out2.print(idk.get(i) + "\t");
            erPos.print(idk.get(i) + "\t");
            erNeg.print(idk.get(i) + "\t");
            for(int j = 0; j < types.length;j++)
                types[j].print(idk.get(i) + "\t");
        }
        out.println("Tumor_Stage\tVital_Status\tAge\tSubtype");
        out2.println("Tumor_Stage\tVital_Status\tAge\tSubtype");

        erPos.println("Tumor_Stage\tVital_Status\tAge\tSubtype");
        erNeg.println("Tumor_Stage\tVital_Status\tAge\tSubtype");
        for(int j = 0; j < types.length;j++)
            types[j].println("Tumor_Stage\tVital_Status\tAge");
        // out.println("Tumor_Stage\tVital_Status\tGender\tAge\tTumor\tRace\tSubtype");
        BufferedReader b2 = new BufferedReader(new FileReader("BRCA.merged_only_clinical_clin_format.txt"));
        //getarrays of each barcode?
        //in other file we'll have barcode by barcode
        HashMap<String,Integer> map = new HashMap<String,Integer>();
        HashMap<Integer,HashMap<String,String>> temp = new HashMap<Integer,HashMap<String,String>>();
        while(b2.ready())
        {
            String [] line = b2.readLine().split("\t");
            if(line[0].contains("race"))
            {

                for(int i = 1; i < line.length;i++)
                {
                    HashMap<String,String> hm = temp.get(i);
                    if(line[i].contains("white"))
                    {
                        hm.put("Race","0");
                    }
                    else if(line[i].contains("black"))
                    {
                        hm.put("Race","1");
                    }
                    else if(line[i].contains("asian"))
                    {
                        hm.put("Race","2");
                    }
                    else if(line[i].contains("indian"))
                    {
                        hm.put("Race","3");
                    }
                    temp.put(i,hm);
                }
            }
            if(line[0].contains("stage_event.pathologic_stage"))
            {
                for(int i = 1; i < line.length;i++)
                {
                    HashMap<String,String> hm = temp.get(i);
                    if(line[i].contains("iii"))
                    {
                        hm.put("Stage","3");
                    }
                    else if (line[i].contains("ii"))
                    {
                        hm.put("Stage","2");
                    }
                    else if(line[i].contains("iv"))
                    {
                        hm.put("Stage","4");
                    }
                    else if(line[i].contains("i"))
                    {
                        hm.put("Stage","1");
                    }
                    temp.put(i,hm);
                }
            }
            if(line[0].contains("gender"))
            {
                for(int i = 1; i < line.length;i++)
                {
                    HashMap<String,String> hm = temp.get(i);
                    if(line[i].equals("female"))
                    {
                        hm.put("Gender","0");
                    }
                    else
                    {
                        hm.put("Gender","1");
                    }
                    temp.put(i,hm);
                }
            }
            if(line[0].contains("age_at_init"))
            {
                for(int i = 1; i < line.length;i++)
                {
                    HashMap<String,String> hm = new HashMap<String,String>();
                    hm.put("Age",line[i]);
                    temp.put(i,hm);
                }
            }
            if(line[0].equals("patient.vital_status"))
            {
                for(int i = 1; i < line.length;i++)
                {
                    HashMap<String,String> hm = temp.get(i);
                    if(line[i].equals("alive"))
                    {
                        hm.put("Vital","1");
                    }
                    else
                    {
                        hm.put("Vital","0");
                    }
                    temp.put(i,hm);
                }
            }
            if(line[0].equals("patient.bcr_patient_barcode"))
            {
                for(int i = 1; i < line.length;i++)
                {
                    String code = line[i].split("-")[2].toUpperCase();
                    map.put(code,i);
                }
            }
        }
        b2 = new BufferedReader(new FileReader("12864_2016_2911_MOESM1_ESM.txt"));
        int index = -1;
        while(b2.ready())
        {
            String [] line = b2.readLine().split("\t");

            if(index==-1)
            {
                for(int i = 0; i < line.length;i++)
                {
                    if(line[i].contains("Our Classification"))
                        index = i;
                }
            }
            else
            {
                if(map.get(line[0].split("-")[2])==null)
                    continue;
                int code = map.get(line[0].split("-")[2]);
                HashMap<String,String> hm = temp.get(code);
                if(line[index].equals("NA"))
                    continue;
                else
                {
                    if(line[index].equals("Luminal A"))
                        hm.put("Subtype","0");
                    else if(line[index].equals("Luminal B"))
                        hm.put("Subtype","1");
                    else if(line[index].equals("Basal-like"))
                        hm.put("Subtype","2");
                    else if(line[index].equals("NA"))
                        continue;
                    else if(line[index].equals("HER2"))
                        hm.put("Subtype","3");
                    else
                        hm.put("Subtype",line[index]);
                }
            }

        }
        b2.close();
        for(int i = 0; i < samples.length;i++)
        {
            String code = samples[i].split("\\.")[2];
            if(map.get(code)==null)
                continue;
            else
            {
                int x = map.get(code);
                HashMap<String,String> t = temp.get(x);
                if(t.get("Gender")==null || t.get("Race")==null || t.get("Vital")==null || t.get("Age")==null || t.get("Stage")==null || t.get("Subtype")==null)
                    continue;
                int tumor = Integer.parseInt(samples[i].substring(13,14));
                //1 is normal, 0 is tumor
                if(tumor==1)
                {
                    for(int j = 0; j < idk.size();j++)
                    {
                        out2.print(data[j][i] + "\t");
                    }
                    out2.println(t.get("Stage") + "\t" + t.get("Vital") + "\t" + t.get("Age") + "\t" + t.get("Subtype"));
                }
                else
                {

                    for (int j = 0; j < idk.size(); j++) {
                        out.print(data[j][i] + "\t");
                        types[Integer.parseInt(t.get("Subtype"))].print(data[j][i] + "\t");
                    }

                    if(Integer.parseInt(t.get("Subtype"))==0 || Integer.parseInt(t.get("Subtype"))==1)
                    {
                        for(int j = 0; j < idk.size();j++)
                        {
                            erPos.print(data[j][i] + "\t");
                        }
                        erPos.println(t.get("Stage") + "\t" + t.get("Vital") + "\t" + t.get("Age") + "\t" + t.get("Subtype"));
                    }
                    else
                    {
                        for(int j = 0; j < idk.size();j++)
                        {
                            erNeg.print(data[j][i] + "\t");
                        }
                        erNeg.println(t.get("Stage") + "\t" + t.get("Vital") + "\t" + t.get("Age") + "\t" + t.get("Subtype"));
                    }

                    types[Integer.parseInt(t.get("Subtype"))].println(t.get("Stage") + "\t" + t.get("Vital") + "\t" + t.get("Age"));
                    //out.println(t.get("Stage")+"\t"+t.get("Vital")+"\t"+t.get("Gender")+"\t"+t.get("Age")+"\t"+tumor+ "\t" + t.get("Race") + "\t" + t.get("Subtype"));
                    out.println(t.get("Stage") + "\t" + t.get("Vital") + "\t" + t.get("Age") + "\t" + t.get("Subtype"));
                }
            }
        }
        out2.flush();
        out2.close();
        out.flush();
        b.close();
        b2.close();
        out.close();
        erPos.flush();
        erPos.close();
        erNeg.flush();
        erNeg.close();
        for(int i = 0; i < types.length;i++) {
            types[i].flush();
            types[i].close();
        }
    }
}