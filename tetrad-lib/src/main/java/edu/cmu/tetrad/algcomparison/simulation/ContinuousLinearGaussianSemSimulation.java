package edu.cmu.tetrad.algcomparison.simulation;

import edu.cmu.tetrad.data.DataType;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.PcStable;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.util.TetradMatrix;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author jdramsey
 */
public class ContinuousLinearGaussianSemSimulation implements Simulation {
    private Graph graph;
    private List<DataSet> dataSets;
private int numDeter;
public TetradMatrix cov;
    public void setruns(int i)
    {
        numDeter = i;
    }
    public void simulate(Parameters parameters, Graph g) {
        dataSets = new ArrayList<>();
        Random rand = new Random();
        for (int i = 0; i < parameters.getInt("numRuns"); i++) {
            SemPm pm = new SemPm(g);
            SemIm im = new SemIm(pm);

            int numNodes = numDeter;
            ArrayList<Node> used = new ArrayList<Node>();
                Node curr = im.getVariableNode("B");
                im.setErrVar(curr, rand.nextDouble() * 0.01);
            System.out.println(im);
            dataSets.add(im.simulateData(parameters.getInt("sampleSize"), false));
        }

    }
    @Override
    public void simulate(Parameters parameters) {
        dataSets = new ArrayList<>();
        Random rand = new Random();
        this.graph = GraphUtils.randomGraphRandomForwardEdges(
                parameters.getInt("numMeasures"),
                parameters.getInt("numLatents"),
                parameters.getInt("numEdges"),
                parameters.getInt("maxDegree"),
                parameters.getInt("maxIndegree"),
                parameters.getInt("maxOutdegree"),
                parameters.getInt("connected") == 1);

        for (int i = 0; i < parameters.getInt("numRuns"); i++) {
            SemPm pm = new SemPm(graph);
            SemIm im = new SemIm(pm);
            int numNodes = numDeter;
            ArrayList<Node> used = new ArrayList<Node>();
            Random rando = new Random();
            B:for(int j = 0; j < numNodes;j++)
            {
                if(used.size()==numNodes)
                    break;
                int node = rando.nextInt(pm.getVariableNodes().size());

                Node curr = this.graph.getNode("X" + (node+1));
                while(used.contains(curr))
                {
                    node = rando.nextInt(pm.getVariableNodes().size());
                    curr = this.graph.getNode("X" + (node+1));
                }
                used.add(curr);
               curr =  im.getVariableNode("X" + (node+1));

            }
            if(parameters.getInt("faithful")==1)
            {
                boolean redo = false;
                D:do {


                    DataSet d = im.simulateData(parameters.getInt("sampleSize"), false);
                    IndependenceTest ind = new IndTestFisherZ(d, 0.001);
                    PcStable p = new PcStable(ind);
                    Graph g = p.search();
                System.out.println("Esimated: " + g + ", True:" + graph);
                    A:
                    for (Edge e : g.getEdges()) {
                        if (graph.getEdge(graph.getNode(e.getNode1().getName()), graph.getNode(e.getNode2().getName())) == null) {
                            redo = true;
                            System.out.println("False Positive");
                            break A;
                        }

                    }
                    if (!redo) {
                        C:
                        for (Edge e : graph.getEdges()) {
                            if (g.getEdge(g.getNode(e.getNode1().getName()), g.getNode(e.getNode2().getName())) == null) {
                                System.out.println("False Negative");
                                redo = true;
                                break C;
                            }
                        }
                    }

                    if(parameters.getInt("sampleSize")>10000)
                        break;

                    if(redo) {
                        System.out.println(parameters.getInt("sampleSize"));
                        parameters.setValue("sampleSize",(int)(parameters.getInt("sampleSize")*1.5));
                        continue D;
                    }
                    else {
                        System.out.println("Survived");
                    }
                    dataSets.add(d);
                    redo = false;
                } while(redo);

                if(redo) {
                    dataSets.add(null);
                    return;
                }

            }
            else
            dataSets.add(im.simulateData(parameters.getInt("sampleSize"), false));
           cov = im.getImplCovar(false);
        }
    }

    @Override
    public DataSet getDataSet(int index) {
        return dataSets.get(index);
    }

    @Override
    public Graph getTrueGraph() {
        return graph;
    }

    @Override
    public String getDescription() {
        return "Linear, Gaussian SEM simulation";
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

    @Override
    public int getNumDataSets() {
        return dataSets.size();
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }
}
