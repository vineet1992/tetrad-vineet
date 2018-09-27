package edu.cmu.tetrad.algcomparison.graph;

import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.List;

/**
 * Creates a random graph by adding forward edges.
 *
 * @author jdramsey
 */
public class RandomForward implements RandomGraph {
    static final long serialVersionUID = 23L;

    @Override
    public Graph createGraph(Parameters parameters) {
        NormalDistribution nd = new NormalDistribution(parameters.getInt("meanDegree"),1);
        return GraphUtils.randomGraphRandomForwardEdges(
                parameters.getInt("numMeasures") + parameters.getInt("numLatents"),
                parameters.getInt("numLatents"),
                (int)(nd.sample()* parameters.getInt("numMeasures") / 2),
                parameters.getInt("maxDegree"),
                parameters.getInt("maxIndegree"),
                parameters.getInt("maxOutdegree"),
                parameters.getBoolean("connected"));
    }

    @Override
    public String getDescription() {
        return "Graph constructed by adding random forward edges";
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add("numMeasures");
        parameters.add("numLatents");
        parameters.add("avgDegree");
        parameters.add("maxDegree");
        parameters.add("maxIndegree");
        parameters.add("maxOutdegree");
        parameters.add("connected");
        return parameters;
    }
}
