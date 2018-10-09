package edu.cmu.tetrad.algcomparison.graph;

import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.util.Parameters;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by vinee_000 on 10/4/2018.
 */
public class RandomGraphUniform implements RandomGraph {
    static final long serialVersionUID = 23L;

    @Override
    public Graph createGraph(Parameters parameters) {
        NormalDistribution nd = new NormalDistribution(parameters.getDouble("meanDegree"),parameters.getDouble("devDegree"));
        return GraphUtils.randomGraphUniform(
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
        return "Graph constructed by generating a graph uniformly at random";
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
