package edu.pitt.csb.mgm;

import edu.cmu.tetrad.util.TetradSerializable;

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
