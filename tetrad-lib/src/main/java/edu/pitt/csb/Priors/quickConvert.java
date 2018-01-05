package edu.pitt.csb.Priors;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintStream;

/**
 * Created by vinee_000 on 12/11/2017.
 */
public class quickConvert {
    public static void main(String [] args) throws Exception
    {
        BufferedReader b = new BufferedReader(new FileReader("Stabilities_No_Priors.txt"));
        PrintStream out = new PrintStream("Stabilities_No_Priors_Fixed.txt");
        out.println(b.readLine());
        out.print(b.readLine()+"\t");
        while(b.ready())
        {
            String [] line = b.readLine().split("\t");
            if(line.length==1)
                out.print(line[0] + "\t");
            else
            {
                out.println(line[0]);
                out.print(line[1] + "\t" + line[2] + "\t");
            }
        }
        b.close();
        out.flush();
        out.close();
    }
}
