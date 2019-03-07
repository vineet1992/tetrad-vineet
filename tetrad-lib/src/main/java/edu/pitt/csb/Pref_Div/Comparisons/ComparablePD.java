package edu.pitt.csb.Pref_Div.Comparisons;

import edu.pitt.csb.Pref_Div.Gene;

import java.util.List;
import java.util.Map;

/**
 * Created by vinee_000 on 2/21/2019.
 * Defines PrefDiv flavors that can be compared to one another
 */
public interface ComparablePD {

    public List<Gene> getLastSelected();

    public Map<Gene,List<Gene>> getLastCluster();
}
