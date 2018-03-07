package edu.pitt.csb.Pref_Div;

import java.util.Comparator;

/**
 * Created by vinee_000 on 2/23/2018.
 */
public class ArrayIndexComparator implements Comparator<Integer>
{
    private final float[] array;
    private final double [] array2;
    private final int [] array3;

    public ArrayIndexComparator(float[] array)
    {
        this.array = array;
        this.array2 = null;
        this.array3 = null;
    }

    public ArrayIndexComparator(double [] array)
    {
        array2 = array;
        this.array = null;
        this.array3 = null;
    }
    public ArrayIndexComparator(int [] array)
    {
        array3 = array;
        this.array = null;
        this.array2 = null;
    }
    public Integer[] createIndexArray()
    {
        if(array2!=null)
        {
            Integer[] indexes = new Integer[array2.length];
            for (int i = 0; i < array2.length; i++) {
                indexes[i] = i; // Autoboxing
            }
            return indexes;
        }
        else if(array!=null) {
            Integer[] indexes = new Integer[array.length];
            for (int i = 0; i < array.length; i++) {
                indexes[i] = i; // Autoboxing
            }
            return indexes;
        }
        else
        {
            Integer[] indexes = new Integer[array3.length];
            for (int i = 0; i < array3.length; i++) {
                indexes[i] = i; // Autoboxing
            }
            return indexes;
        }
    }

    @Override
    public int compare(Integer index1, Integer index2)
    {
        if(array2!=null)
        {
            if(array2[index1] > array2[index2])
                return 1;
            else if(array2[index1] < array2[index2])
                return -1;
            else
                return 0;
        }
        else if(array3!=null)
        {
            if(array3[index1] > array3[index2])
                return 1;
            else if(array3[index1] < array3[index2])
                return -1;
            else
                return 0;
        }
        if(array[index1] > array[index2])
            return 1;
        else if(array[index1] < array[index2])
            return -1;
        else
            return 0;
      }
}
