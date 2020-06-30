package edu.cmu.tetrad.algcomparison.simulation;

import edu.cmu.tetrad.algcomparison.graph.RandomGraph;
import edu.cmu.tetrad.algcomparison.utils.HasParameters;
import edu.cmu.tetrad.data.DataModel;
import edu.cmu.tetrad.util.IM;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.sem.GeneralizedSemIm;
import edu.cmu.tetrad.sem.GeneralizedSemPm;
import edu.pitt.csb.mgm.MixedUtils;
import org.apache.commons.lang3.RandomUtils;

import java.util.*;

/**
 * A version of the Lee & Hastic simulation which is guaranteed ot generate a discrete
 * data set.
 *
 * @author jdramsey
 */
public class LeeHastieSimulation implements Simulation, HasParameters {
    static final long serialVersionUID = 23L;
    private RandomGraph randomGraph;
    private List<DataSet> dataSets = new ArrayList<>();
    private List<Graph> graphs = new ArrayList<>();
    private List<GeneralizedSemIm> ims = new ArrayList<>();
    private DataType dataType;
    private List<Node> shuffledOrder;
    private Graph initGraph; //Added to include a possible initial graph from which to generate data

    public LeeHastieSimulation(RandomGraph graph) {
        this.randomGraph = graph;
        this.initGraph = null;
    }


    //Added by Vineet for separate inclusion of an initial graph
    public void setInitialGraph(Graph g)
    {
        this.initGraph = g;
    }

    @Override
    public void createData(Parameters parameters) {
        double percentDiscrete = parameters.getDouble("percentDiscrete");

        boolean discrete = parameters.getString("dataType").equals("discrete");
        boolean continuous = parameters.getString("dataType").equals("continuous");

        if (discrete && percentDiscrete != 100.0) {
            throw new IllegalArgumentException("To simulate discrete data, 'percentDiscrete' must be set to 0.0.");
        } else if (continuous && percentDiscrete != 0.0) {
            throw new IllegalArgumentException("To simulate continuoue data, 'percentDiscrete' must be set to 100.0.");
        }

        if (discrete) this.dataType = DataType.Discrete;
        if (continuous) this.dataType = DataType.Continuous;

        this.shuffledOrder = null;
        //THIS IS CODE I ADDED TO ALLOW FOR SEPARATE INCLUSION OF AN INITIAL GRAPH (Vineet)
        Graph graph;
        if(initGraph==null)
            graph = randomGraph.createGraph(parameters);
        else
            graph = initGraph;

        dataSets = new ArrayList<>();
        graphs = new ArrayList<>();

        for (int i = 0; i < parameters.getInt("numRuns"); i++) {
            //System.out.println("Simulating dataset #" + (i + 1));

            if (parameters.getBoolean("differentGraphs") && i > 0) {
                graph = randomGraph.createGraph(parameters);
            }

            graphs.add(graph);

            DataSet dataSet = simulate(graph, parameters);
            dataSet.setName("" + (i + 1));
            dataSets.add(dataSet);
        }
    }

    @Override
    public Graph getTrueGraph(int index) {
        return graphs.get(index);
    }

    @Override
    public DataModel getDataModel(int index) {
        return dataSets.get(index);
    }

    @Override
    public IM getInstantiatedModel(int index){return ims.get(index);}

    @Override
    public String getDescription() {
        return "Lee & Hastie simulation using " + randomGraph.getDescription();
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = randomGraph.getParameters();
        parameters.add("numCategories");
        parameters.add("percentDiscrete");
        parameters.add("numRuns");
        parameters.add("differentGraphs");
        parameters.add("sampleSize");
        return parameters;
    }

    @Override
    public int getNumDataModels() {
        return dataSets.size();
    }

    @Override
    public DataType getDataType() {
        return dataType;
    }

    private DataSet simulate(Graph dag, Parameters parameters) {
        HashMap<String, Integer> nd = new HashMap<>();

        List<Node> nodes = dag.getNodes();

        Collections.shuffle(nodes);

        if (this.shuffledOrder == null) {
            List<Node> shuffledNodes = new ArrayList<>(nodes);
            Collections.shuffle(shuffledNodes);
            this.shuffledOrder = shuffledNodes;
        }

        for (int i = 0; i < nodes.size(); i++) {
            if (i < nodes.size() * parameters.getDouble("percentDiscrete") * 0.01) {
                final int minNumCategories = parameters.getInt("minCategories");
                final int maxNumCategories = parameters.getInt("maxCategories");
                final int value = pickNumCategories(minNumCategories, maxNumCategories);
                nd.put(shuffledOrder.get(i).getName(), value);
            } else {
                nd.put(shuffledOrder.get(i).getName(), 0);
            }
        }

        Graph graph = MixedUtils.makeMixedGraph(dag, nd);

        final double coefLow = parameters.getDouble("coefLow");
        final double coefHigh = parameters.getDouble("coefHigh");
        final double varLow = parameters.getDouble("varLow");
        final double varHigh = parameters.getDouble("varHigh");
        GeneralizedSemPm pm = MixedUtils.GaussianCategoricalPm(graph, "Split(" + (-1*coefHigh) + "," + (-1*coefLow) + "," + coefLow + "," + coefHigh + ")",varLow,varHigh);
        GeneralizedSemIm im = MixedUtils.GaussianCategoricalIm(pm);
        ims.add(im);

        DataSet ds = im.simulateDataAvoidInfinity(parameters.getInt("sampleSize"), false);
        return MixedUtils.makeMixedData(ds, nd);
    }

    private int pickNumCategories(int min, int max) {
        return RandomUtils.nextInt(min, max + 1);
    }
}
