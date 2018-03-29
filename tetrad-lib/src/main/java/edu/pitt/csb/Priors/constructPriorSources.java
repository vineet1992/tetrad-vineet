package edu.pitt.csb.Priors;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;

/**
 * Created by vinee_000 on 12/6/2017.
 */
public class constructPriorSources {

    public static void main(String [] args)throws Exception {


        boolean createPathway = false;
        boolean createPAM = true;
        String priorDirectory = "pathway_lists_large";
        String outDirectory = "prior_sources";
        String expFile = "genes_with_clinical.txt";

        BufferedReader b = new BufferedReader(new FileReader(expFile));
        b.readLine();
        ArrayList<String> vars = new ArrayList<String>();
        while (b.ready()) {
            String line = b.readLine();
            if (line.equals("/data"))
                break;
            vars.add(line.split(":")[0]);

        }



        File temp = new File(outDirectory);
        if (!temp.exists())
            temp.mkdirs();
        if(createPathway) {
            File f = new File(priorDirectory);
            for (File t : f.listFiles()) {
                String[] name = t.getName().split("_");
                PrintStream out = new PrintStream(outDirectory + "/" + name[0]);
                createPrior(t, out, vars);

            }
        }
        if(createPAM) {
            BufferedReader b2 = new BufferedReader(new FileReader("PAM50.txt"));
            PrintStream out = new PrintStream(outDirectory +"/Prior_PAM50.txt");
            ArrayList<String> pam50 = new ArrayList<String>();
            while (b2.ready()) {
                pam50.add(b2.readLine());
            }
            int ind = vars.indexOf("Subtype");
            for (int i = 0; i < vars.size(); i++) {
                for (int j = 0; j < vars.size(); j++) {
                    if ((i == ind && pam50.contains(vars.get(j))) || (j == ind && pam50.contains(vars.get(i)))) {
                        if (j == vars.size() - 1)
                            out.println(1);
                        else
                            out.print(1 + "\t");
                    } else {
                        if (j == vars.size() - 1)
                            out.println(0);
                        else
                            out.print(0 + "\t");
                    }
                }
            }
            for (int k = 1; k < 6; k++) {
                b2 = new BufferedReader(new FileReader("Irr_PAM50_" + k + ".txt"));
                out = new PrintStream(outDirectory + "/Irr_PAM50_" + (k - 1) + ".txt");
                pam50 = new ArrayList<String>();
                while (b2.ready()) {
                    pam50.add(b2.readLine());
                }
                ind = vars.indexOf("Subtype");
                for (int i = 0; i < vars.size(); i++) {
                    for (int j = 0; j < vars.size(); j++) {
                        if ((i == ind && pam50.contains(vars.get(j))) || (j == ind && pam50.contains(vars.get(i)))) {
                            if (j == vars.size() - 1)
                                out.println(1);
                            else
                                out.print(1 + "\t");
                        } else {
                            if (j == vars.size() - 1)
                                out.println(0);
                            else
                                out.print(0 + "\t");
                        }
                    }
                }
            }
            out.flush();
            out.close();
            b2.close();
        }
    }
    public static void createPrior(File pathway, PrintStream out, ArrayList<String> vars) throws Exception
    {
        out.println("\"" + pathway.getName() + "\"");
        //loop through pathway file, and add elements to a double [][], then print it all out to the file
        double [][] prior = fileToPrior(pathway, vars);
        for(int i = 0; i < prior.length;i++)
        {
            for(int j = 0; j < prior[i].length;j++)
            {
                if(j==prior[i].length-1)
                    out.println(prior[i][j]);
                else
                    out.print(prior[i][j] + "\t");
            }
        }
        out.flush();
        out.close();

    }
    public static double [][] fileToPrior(File pathway, ArrayList<String>vars) throws Exception
    {
        double [][] prior = new double[vars.size()][vars.size()];
        BufferedReader b = new BufferedReader(new FileReader(pathway.getAbsolutePath()));
        b.readLine();//eat the title
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            int i = vars.indexOf(line[0]);
            int j = vars.indexOf(line[1]);
            if(i==-1||j==-1)
                continue;
            prior[i][j] = 1;
            prior[j][i] = 1;
        }
        b.close();
        return prior;
    }
}
