package edu.pitt.csb.Pref_Div;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.util.StatUtils;
import edu.pitt.csb.Priors.mgmPriors;
import edu.pitt.csb.mgm.MixedUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

/**
 * Created by vinee_000 on 2/21/2018.
 */
public class npnTest {
    public static void main(String [] args) throws Exception
    {

        double [] pValues = {0.500000,0.100000,0.003000,0.050000,0.020000,0.002342,0.623400,0.234230,0.000010,0.045000};
        pValues = mgmPriors.adjustPValues(pValues);
        System.out.println(Arrays.toString(pValues));
        System.exit(0);


        double [] x = {-0.132432981149731,1.45752086658549};
        double [] y = {1.09228467598778,-0.552659389322541};
        System.out.println((float)StatUtils.correlation(x,y));
        x = new double[]{0,1};
        y = new double[]{0,0.1};
        System.out.println((float)StatUtils.correlation(x,y));


        float [] x2 = {2.5f,3.5f,1.7f,-2.32f,-6.55f};

        ArrayIndexComparator c = new ArrayIndexComparator(x2);
        Integer [] inds = c.createIndexArray();
        System.out.println(Arrays.toString(inds));
        Arrays.sort(inds,c);
        System.out.println(Arrays.toString(inds));


        Random rand = new Random();
        float [] temp = new float[10];
        for(int i = 0; i < temp.length;i++)
        {
            temp[i] = rand.nextFloat()*25 + 16;
        }

        ArrayList<Float> samps = new ArrayList<Float>();
        for(int i = 0; i < temp.length;i++)
            samps.add(temp[i]);

        Collections.sort(samps);
        System.out.println(Arrays.toString(temp));
        temp = Functions.NPN(temp,true);

        System.out.println("After NPN:" + Arrays.toString(temp));

        for(int i = 0; i < 15;i++)
        {

            float x3 = rand.nextFloat()*25+16;
            System.out.println("Samps before NPN: " + samps + ", Float to add: " + x3);
            x3 = Functions.NPNonTheFly(x3,samps);
            System.out.println("Samps After NPN: " + samps + ", NPN'd value: " + x3);
        }



        /*
        DataSet d = MixedUtils.loadDataSet2("eset_all_annotations.txt");
        double [][] data = d.getDoubleData().toArray();


        for(int i = 0; i < data.length;i++)
        {
            for(int j = 0; j < data[i].length;j++)
            {
                System.out.print(data[i][j] + "\t");
            }
            System.out.println();
        }
*/


    }
}
