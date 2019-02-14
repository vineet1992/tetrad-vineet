package edu.cmu.tetrad.algcomparison.simulation;

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.sem.GeneralizedSemIm;
import edu.cmu.tetrad.sem.GeneralizedSemPm;
import edu.cmu.tetrad.util.IM;
import edu.pitt.csb.mgm.MixedUtils;
import org.apache.commons.math3.stat.inference.TTest;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

/**
 * @author jdramsey
 */
public class MixedLeeHastieSimulation implements Simulation {
    private List<DataSet> dataSets;
    private Graph graph;
    private GeneralizedSemIm instModel;
    private boolean realPathways; //Should pathways be simulated more realistically

    public void simulate(Parameters parameters) {
        this.dataSets = new ArrayList<>();
        if (graph == null) {
            this.graph = GraphUtils.randomGraphRandomForwardEdges(
                    parameters.getInt("numMeasures"), parameters.getInt("numLatents"),
                    parameters.getInt("numEdges"),
                    parameters.getInt("maxDegree"),
                    parameters.getInt("maxIndegree"),
                    parameters.getInt("maxOutdegree"),
                    parameters.getInt("connected") == 1);
        }

        for (int i = 0; i < parameters.getInt("numRuns"); i++) {
            DataSet dataSet = simulate(graph, parameters);
            dataSets.add(dataSet);
        }
    }

    public void createData(edu.cmu.tetrad.util.Parameters parameters)  {
        this.dataSets = new ArrayList<>();
        if (graph == null) {
            this.graph = GraphUtils.randomGraphRandomForwardEdges(
                    parameters.getInt("numMeasures"), parameters.getInt("numLatents"),
                    parameters.getInt("numEdges"),
                    parameters.getInt("maxDegree"),
                    parameters.getInt("maxIndegree"),
                    parameters.getInt("maxOutdegree"),
                    parameters.getInt("connected") == 1);
        }

        for (int i = 0; i < parameters.getInt("numRuns"); i++) {
            DataSet dataSet = simulate(graph, parameters);
            dataSets.add(dataSet);
        }
    }

    public void setRealPathways(boolean realPathways){this.realPathways = realPathways;}

    private DataSet simulate(Graph dag, edu.cmu.tetrad.util.Parameters parameters) {
        HashMap<String, Integer> nd = new HashMap<>();

        List<Node> nodes = dag.getNodes();

        Collections.shuffle(nodes);

        for (int i = 0; i < nodes.size(); i++) {


            if (i < nodes.size() * parameters.getDouble("percentDiscreteForMixedSimulation") * 0.01) {
                nd.put(nodes.get(i).getName(), parameters.getInt("numCategories"));
            } else {
                nd.put(nodes.get(i).getName(), 0);
            }
        }
        Graph graph = MixedUtils.makeMixedGraph(dag, nd);

        GeneralizedSemPm pm = MixedUtils.GaussianCategoricalPm(graph, "Split(-1.5,-.5,.5,1.5)");
        GeneralizedSemIm im = MixedUtils.GaussianCategoricalIm(pm);

        DataSet ds = im.simulateDataAvoidInfinity(parameters.getInt("sampleSize"), false);
        return MixedUtils.makeMixedData(ds, nd);
    }



    //Returns cluster assignments for each Node <X32,1> means that node X32 is in cluster 1
    public HashMap<String,List<Integer>> simulateClustered(Parameters parameters,boolean evenDistribution)
    {
        this.dataSets = new ArrayList<>();
        if(graph==null)
        {
            if(realPathways)
            {
                this.graph = GraphUtils.randomGraphRealPathways(parameters.getInt("numMeasures"),parameters.getInt("numComponents"),parameters.getInt("numConnectedComponents"),
                        parameters.getInt("maxDegree"),parameters.getInt("maxIndegree"),parameters.getInt("maxOutdegree"),parameters.getInt("avgPathways"),parameters.getInt("avgCauses"));
            }else
            {
                this.graph = GraphUtils.randomGraphForwardEdgesClusters(parameters.getInt("numMeasures"),parameters.getInt("numComponents"),parameters.getInt("numConnectedComponents"),
                        parameters.getInt("numEdges"),evenDistribution,true,parameters.getInt("maxDegree")
                        ,parameters.getInt("maxIndegree"),parameters.getInt("maxOutdegree"));
            }


        }
        System.out.println(this.graph);
        for (int i = 0; i < parameters.getInt("numRuns"); i++) {
            DataSet dataSet = simulateClustered(graph, parameters,parameters.getInt("targetContinuous"));
            dataSet = DataUtils.standardizeData(dataSet);
            System.out.println("Graph has: " + graph.getNumNodes() + ", Data has " + dataSet.getNumColumns() + ", supposed to have: " + (parameters.getInt("numMeasures")+1));
            if(dataSet.getNumColumns()!=parameters.getInt("numMeasures")+1)
            {
                System.err.println("Dataset did not contain all variables, ensure your simulated graph is proper");
                System.exit(-1);
            }
            dataSets.add(dataSet);
        }
        return GraphUtils.clusters;

    }
    public void setDataSet(DataSet d, int index) {
        dataSets.set(index, d);
    }
    public DataModel getDataModel(int index) {
        return dataSets.get(index);
    }
    public int getNumDataModels(){return dataSets.size();}
    public void setTrueGraph(Graph g) {
        graph = g;
    }


    public Graph getTrueGraph() {
        return graph;
    }

    public Graph getTrueGraph(int x){return graph;}


    public DataSet getDataSet(int index) {
        return dataSets.get(index);
    }

    @Override
    public String getDescription() {
        return "Lee & Hastie simulation";
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add("numMeasures");
        parameters.add("numLatents");
        parameters.add("numEdges");
        parameters.add("maxDegree");
        parameters.add("maxIndegree");
        parameters.add("maxOutdegree");
        parameters.add("numRuns");
        parameters.add("sampleSize");
        return parameters;
    }

    public int getNumDataSets() {
        return dataSets.size();
    }

    @Override
    public DataType getDataType() {
        return DataType.Mixed;
    }


    public GeneralizedSemIm getIm()
    {
        return instModel;
    }
    private DataSet simulateClustered(Graph dag, Parameters parameters,int targetContinuous) {
        HashMap<String, Integer> nd = new HashMap<>();

        List<Node> nodes = dag.getNodes();

        Collections.shuffle(nodes);

        for (int i = 0; i < nodes.size(); i++) {


            if (i < nodes.size() * parameters.getDouble("percentDiscreteForMixedSimulation") * 0.01) {
                nd.put(nodes.get(i).getName(), parameters.getInt("numCategories"));
            } else {
                nd.put(nodes.get(i).getName(), 0);
            }

            if (targetContinuous == 1 && nodes.get(i).getName().equals("Target"))
                nd.put(nodes.get(i).getName(), 0);
            else if(nodes.get(i).getName().equals("Target"))
                nd.put(nodes.get(i).getName(), parameters.getInt("numCategories"));
        }
            Graph graph = MixedUtils.makeMixedGraph(dag, nd);

            GeneralizedSemPm pm = MixedUtils.GaussianCategoricalPm(graph, "Split(-2,-.05,.05,2)");
            GeneralizedSemIm im = MixedUtils.GaussianCategoricalIm(pm);

            instModel = im;
            DataSet ds = im.simulateDataAvoidInfinity(parameters.getInt("sampleSize"), false);
            return MixedUtils.makeMixedData(ds, nd);
        }
    private DataSet simulate(Graph dag, Parameters parameters) {
        HashMap<String, Integer> nd = new HashMap<>();

        List<Node> nodes = dag.getNodes();

        Collections.shuffle(nodes);

        for (int i = 0; i < nodes.size(); i++) {


            if (i < nodes.size() * parameters.getDouble("percentDiscreteForMixedSimulation") * 0.01) {
                nd.put(nodes.get(i).getName(), parameters.getInt("numCategories"));
            } else {
                nd.put(nodes.get(i).getName(), 0);
            }
        }
        Graph graph = MixedUtils.makeMixedGraph(dag, nd);

        GeneralizedSemPm pm = MixedUtils.GaussianCategoricalPm(graph, "Split(-1.5,-.5,.5,1.5)");
        GeneralizedSemIm im = MixedUtils.GaussianCategoricalIm(pm);

        DataSet ds = im.simulateDataAvoidInfinity(parameters.getInt("sampleSize"), false);
        return MixedUtils.makeMixedData(ds, nd);
    }

    public void setInitialGraph(Graph g){throw new UnsupportedOperationException();}
    public IM getInstantiatedModel(int index){throw new UnsupportedOperationException();}

}
/*    public DataSet simulate(Graph dag, Parameters parameters)
    {

    }
}
//TODO Can definitely consider changing this, this is only for CSI experiments
            if(nodes.get(i)instanceof ContinuousVariable)
                    nd.put(nodes.get(i).getName(),0);
                    else if(nodes.get(i) instanceof DiscreteVariable)
                    nd.put(nodes.get(i).getName(),parameters.getInt("numCategories"));*/
