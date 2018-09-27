package edu.pitt.csb.mgm;

import edu.cmu.tetrad.util.TetradSerializable;
import edu.pitt.csb.stability.DataGraphSearch;
import edu.pitt.csb.stability.SearchWrappers;

import java.io.ObjectStreamException;

/**
 * Created by vinee_000 on 10/10/2017.
 */
public final class Algorithm implements TetradSerializable {
    static final long serialVersionUID = 23L;

    public static final  Algorithm FCI = new Algorithm("FCI");
    public static final Algorithm MGMFCI = new Algorithm("MGM-FCI");
    public static final Algorithm MGMFCIMAX = new Algorithm("MGM-FCI-MAX");
    public static final Algorithm PCS = new Algorithm("PC-Stable");
    public static final Algorithm CPC = new Algorithm("CPC");
    public static final Algorithm PCMAX = new Algorithm("PC-Max");
    public static final Algorithm FGS = new Algorithm("FGES");
    public static final Algorithm FCIMAX = new Algorithm("FCI-MAX");
    public static final Algorithm MGM = new Algorithm("MGM");
    public static final Algorithm NONE = new Algorithm("No type");

    /**
     * The name of this type.
     */
    private final transient String name;

    /**
     * Protected constructor for the types; this allows for extension in case
     * anyone wants to add formula types.
     */
    private Algorithm(String name) {
        this.name = name;
    }

    /**
     * Generates a simple exemplar of this class to test serialization.
     */
    public static Algorithm serializableInstance() {
        return Algorithm.PCS;
    }

    public static DataGraphSearch algToSearchWrapper(Algorithm a, double [] params)
    {
        if(a==Algorithm.CPC)
            return new SearchWrappers.CpcStableWrapper(params);
        else if(a==Algorithm.FCI)
            return new SearchWrappers.FCIWrapper(params);
        else if(a==Algorithm.FCIMAX)
            return new SearchWrappers.FCIMAXWrapper(params);
        else if(a==Algorithm.FGS)
            return new SearchWrappers.FgesWrapper(params);
        else if(a==Algorithm.MGM)
            return new SearchWrappers.MGMWrapper(params);
        else if(a==Algorithm.PCMAX)
            return new SearchWrappers.PcMaxWrapper(params);
        else if(a==Algorithm.PCS)
            return new SearchWrappers.PcStableWrapper(params);
        else
            return null;
    }
    /**
     * Prints out the name of the type.
     */
    public String toString() {
        return name;
    }

    // Declarations required for serialization.
    private static int nextOrdinal = 0;
    private final int ordinal = nextOrdinal++;
    private static final Algorithm[] TYPES = {FCI,MGMFCI,MGMFCIMAX,PCS,CPC,PCMAX,NONE};

    Object readResolve() throws ObjectStreamException {
        return TYPES[ordinal]; // Canonicalize.
    }
}
