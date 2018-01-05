package edu.pitt.csb;

/**
 * Created by vinee_000 on 9/8/2017.
 */

import edu.cmu.tetrad.graph.*;

import java.util.HashSet;
import java.util.Set;

public class CompareGraphs {
    private Graph truth;
    private Graph est;
    private int adjTp;
    private int adjFp;
    private int adjFn;
    private int adjTn;
    private int arrowsTp;
    private int arrowsFp;
    private int arrowsFn;
    private int arrowsTn;

    public CompareGraphs(Graph truth, Graph est) {
        this.truth = truth;
        this.est = est;
        adjTp = 0;
        adjFp = 0;
        adjFn = 0;
        arrowsTp = 0;
        arrowsFp = 0;
        arrowsFn = 0;

        Set<Edge> allOriented = new HashSet<>();
        allOriented.addAll(this.truth.getEdges());
        allOriented.addAll(this.est.getEdges());

        Set<Edge> allUnoriented = new HashSet<>();
        for (Edge edge : this.truth.getEdges()) {
            allUnoriented.add(Edges.undirectedEdge(edge.getNode1(), edge.getNode2()));
        }

        for (Edge edge : this.est.getEdges()) {
            allUnoriented.add(Edges.undirectedEdge(edge.getNode1(), edge.getNode2()));
        }

        for (Edge edge : allUnoriented) {
            Node est1 = est.getNode(edge.getNode1().getName());
            Node est2 = est.getNode(edge.getNode2().getName());
            Node tru1 = truth.getNode(edge.getNode1().getName());
            Node tru2 = truth.getNode(edge.getNode2().getName());
            if (this.est.isAdjacentTo(est1,est2) &&
                    !this.truth.isAdjacentTo(tru1,tru2)) {
                adjFp++;
            }

            if (this.truth.isAdjacentTo(tru1,tru2) &&
                    !this.est.isAdjacentTo(est1,est2)) {
                adjFn++;
            }

            if (this.truth.isAdjacentTo(tru1,tru2) &&
                    this.est.isAdjacentTo(est1,est2)) {
                adjTp++;
            }
        }

        int allEdges = this.truth.getNumNodes() * (this.truth.getNumNodes() - 1) / 2;

        adjTn = allEdges - adjFn;




        for (Edge edge : allOriented) {
            Node est1 = est.getNode(edge.getNode1().getName());
            Node est2 = est.getNode(edge.getNode2().getName());
            Node tru1 = truth.getNode(edge.getNode1().getName());
            Node tru2 = truth.getNode(edge.getNode2().getName());
          //  Endpoint e1Est = edge.getProximalEndpoint(est1);
          //  Endpoint e2Est = edge.getProximalEndpoint(edge.getNode2());

            Edge edgeTrue = this.truth.getEdge(tru1,tru2);
            Edge edgeEst = this.est.getEdge(est1,est2);

            Endpoint e1True = null;
            Endpoint e2True = null;
            Endpoint e1Est = null;
            Endpoint e2Est = null;

            if (edgeTrue != null) {
                e1True = edgeTrue.getProximalEndpoint(edgeTrue.getNode1());
                e2True = edgeTrue.getProximalEndpoint(edgeTrue.getNode2());
            }

            if(edgeEst!=null)
            {
                e1Est = edgeEst.getProximalEndpoint(edgeEst.getNode1());
                e2Est = edgeEst.getProximalEndpoint(edgeEst.getNode2());
            }
            if (e1Est == Endpoint.ARROW && e1True != Endpoint.ARROW) {
                arrowsFp++;
            }

            if (e2Est == Endpoint.ARROW && e2True != Endpoint.ARROW) {
                arrowsFp++;
            }

            if (e1True == Endpoint.ARROW && e1Est != Endpoint.ARROW) {
                arrowsFn++;
            }

            if (e2True == Endpoint.ARROW && e2Est != Endpoint.ARROW) {
                arrowsFn++;
            }

            if (e1True == Endpoint.ARROW && e1Est == Endpoint.ARROW) {
                arrowsTp++;
            }

            if (e2True == Endpoint.ARROW && e2Est == Endpoint.ARROW) {
                arrowsTp++;
            }
        }

        int allEdges2 = this.truth.getNumNodes() * (this.truth.getNumNodes() - 1) / 2;
        arrowsTn = allEdges2 - arrowsFn;
    }


    public double getPrec()
    {
        return adjTp/(double)(adjTp+adjFp);
    }
    public double getRec()
    {
        return adjTp/(double)(adjTp+adjFn);
    }

    public double getArrowPrec()
    {
        return arrowsTp/(double)(arrowsTp+arrowsFp);
    }
    public double getArrowRec()
    {
        return arrowsTp/(double)(arrowsTp+arrowsFn);
    }

    public int getArrowsTp() {
        return arrowsTp;
    }

    public int getArrowsFp() {
        return arrowsFp;
    }

    public int getArrowsFn() {
        return arrowsFn;
    }

    public int getArrowsTn() {
        return arrowsTn;
    }
    public int getAdjTp() {
        return adjTp;
    }

    public int getAdjFp() {
        return adjFp;
    }

    public int getAdjFn() {
        return adjFn;
    }

    public int getAdjTn() {
        return adjTn;
    }

}

