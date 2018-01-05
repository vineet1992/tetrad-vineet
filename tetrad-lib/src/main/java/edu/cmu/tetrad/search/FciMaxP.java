package edu.cmu.tetrad.search;

import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.data.KnowledgeEdge;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.TetradLogger;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;

import static edu.cmu.tetrad.search.SearchGraphUtils.isArrowpointAllowed;

/**
 * Created by vinee_000 on 8/13/2016.
 */
public class FciMaxP {
    ///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015 by Peter Spirtes, Richard Scheines, Joseph   //
// Ramsey, and Clark Glymour.                                                //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA //
///////////////////////////////////////////////////////////////////////////////

        /* The PAG being constructed.
         */
        private Graph graph;

        /**
         * The SepsetMap being constructed.
         */
        private SepsetMap sepsets;

        /**
         * The background knowledge.
         */
        private IKnowledge knowledge = new Knowledge2();


        /**
         * The variables to search over (optional)
         */
        private List<Node> variables = new ArrayList<Node>();

        private IndependenceTest independenceTest;

        /**
         * flag for complete rule set, true if should use complete rule set, false otherwise.
         */
        private boolean completeRuleSetUsed = false;

        /**
         * True iff the possible dsep search is done.
         */
        private boolean possibleDsepSearchDone = true;

        /**
         * The maximum length for any discriminating path. -1 if unlimited; otherwise, a positive integer.
         */
        private int maxPathLength = -1;

        /**
         * The depth for the fast adjacency search.
         */
        private int depth = -1;

        /**
         * Elapsed time of last search.
         */
        private long elapsedTime;

        /**
         * The logger to use.
         */
        private TetradLogger logger = TetradLogger.getInstance();

    private double factor;
    private int parallelism;
        /**
         * True iff verbose output should be printed.
         */
        private boolean verbose = false;
        public ConcurrentHashMap<String,String> whyOrient;


        private Graph trueDag;
        private ConcurrentMap<Node, Integer> hashIndices;
        private ICovarianceMatrix covarianceMatrix;
        private double penaltyDiscount = 2;
        private SepsetMap possibleDsepSepsets = new SepsetMap();
        private Graph initialGraph;
        private int possibleDsepDepth = -1;


        private ForkJoinPool pool;
    private double chunk = 5;
        //============================CONSTRUCTORS============================//

        /**
         * Constructs a new FCI search for the given independence test and background knowledge.
         */

        public FciMaxP(IndependenceTest independenceTest) {
            if (independenceTest == null || knowledge == null) {
                throw new NullPointerException();
            }
            this.factor = 1;
            this.parallelism = Runtime.getRuntime().availableProcessors();
            this.independenceTest = independenceTest;
            this.variables.addAll(independenceTest.getVariables());
            buildIndexing(independenceTest.getVariables());
            pool = new ForkJoinPool();
        }
        public FciMaxP(IndependenceTest independenceTest,double factor,int parallelism) {
            if (independenceTest == null || knowledge == null) {
                throw new NullPointerException();
            }
this.factor = factor;
            this.parallelism = parallelism;
            this.independenceTest = independenceTest;
            this.variables.addAll(independenceTest.getVariables());
            buildIndexing(independenceTest.getVariables());
            pool = new ForkJoinPool(parallelism);
        }

        /**
         * Constructs a new FCI search for the given independence test and background knowledge and a list of variables to
         * search over.
         */
        public FciMaxP(IndependenceTest independenceTest, List<Node> searchVars) {
            if (independenceTest == null || knowledge == null) {
                throw new NullPointerException();
            }

            this.independenceTest = independenceTest;
            this.variables.addAll(independenceTest.getVariables());

            Set<Node> remVars = new HashSet<Node>();
            for (Node node1 : this.variables) {
                boolean search = false;
                for (Node node2 : searchVars) {
                    if (node1.getName().equals(node2.getName())) {
                        search = true;
                    }
                }
                if (!search) {
                    remVars.add(node1);
                }
            }
            this.variables.removeAll(remVars);
        }

        //========================PUBLIC METHODS==========================//

        public int getDepth() {
            return depth;
        }

        public void setDepth(int depth) {
            if (depth < -1) {
                throw new IllegalArgumentException(
                        "Depth must be -1 (unlimited) or >= 0: " + depth);
            }

            this.depth = depth;
        }

        public long getElapsedTime() {
            return this.elapsedTime;
        }

        public Graph search() {
            return search(getIndependenceTest().getVariables());
        }

        public Graph search(List<Node> nodes) {
            FasStableConcurrent fas = new FasStableConcurrent(initialGraph, getIndependenceTest());
            fas.setVerbose(verbose);
            return search(fas);
//        return search(new Fas(getIndependenceTest()));
        }

        public void setInitialGraph(Graph initialGraph) {
            this.initialGraph = initialGraph;
        }
        public Graph search2 ()
        {
            ConditionalGaussianScore score0 = new ConditionalGaussianScore(independenceTest.getDataSets().get(0));
            Fgs2 fgs = new Fgs2(score0);
            fgs.setKnowledge(getKnowledge());
            fgs.setVerbose(verbose);
            fgs.setNumPatternsToStore(0);
            if(graph!=null)
                fgs.setInitialGraph(graph);
//        fgs.setHeuristicSpeedup(faithfulnessAssumed);
            Graph g = fgs.search();
            this.setInitialGraph(g);
            return search(getIndependenceTest().getVariables());


        }
        public Graph search(IFas fas) {
            logger.log("info", "Starting FCI algorithm.");
            logger.log("info", "Independence test = " + getIndependenceTest() + ".");
            fas.setVerbose(verbose);

            fas.setKnowledge(getKnowledge());
            fas.setDepth(depth);
            fas.setVerbose(verbose);
            long time = System.currentTimeMillis();
            this.graph = fas.search();
           // System.out.println("FAS: "  + (System.currentTimeMillis()-time));
            this.sepsets = fas.getSepsets();
           // System.out.println(this.graph);
           // System.out.println(this.sepsets);
           // System.out.println(this.sepsets.get(graph.getNode("X3"),graph.getNode("X24")));
            Map<Node,Integer> pdsepChange = new ConcurrentHashMap<>();

            whyOrient = new ConcurrentHashMap<String,String>();
            graph.reorientAllWith(Endpoint.CIRCLE);

//        SepsetProducer sp = new SepsetsPossibleDsep(graph, independenceTest, knowledge, depth, maxPathLength);
            SepsetMap mp = new SepsetMap();
            SepsetMap mp2 = new SepsetMap();
            SepsetProducer sp = new SepsetsMaxPValuePossDsep(graph, independenceTest, mp, depth, maxPathLength);
            SepsetsMaxPValue sp2 = new SepsetsMaxPValue(graph, independenceTest, mp2, depth);
//TODO cannot assume same sepsets for second collider orientation phase
            // The original FCI, with or without JiJi Zhang's orientation rules
            //        // Optional step: Possible Dsep. (Needed for correctness but very time consuming.)
            if (isPossibleDsepSearchDone()) {
//            long time1 = System.currentTimeMillis();
//            new FciOrient(new SepsetsSet(this.sepsets, independenceTest)).ruleR0(graph);
//            SepsetProducer sepsetProducer = new SepsetsSet(sepsets, independenceTest);

                time = System.currentTimeMillis();
                addColliders(graph, sp2, knowledge);

               // System.out.println("Adding colliders " + (System.currentTimeMillis() - time));
              //  System.out.println("Possible dsep add colliders done");
                time = System.currentTimeMillis();
                for (Edge edge : new ArrayList<>(graph.getEdges())) {
                    Node x = edge.getNode1();
                    Node y = edge.getNode2();

                    List<Node> sepset = sp.getSepset(x, y);

                    if (sepset != null) {
                        graph.removeEdge(x, y);
                        pdsepChange.put(x,0);
                        pdsepChange.put(y,0);
                        //Adjacency set of x or y changed so remember this when reorienting colliders
                        sepsets.set(x, y, sepset);
                       // System.out.println("Possible DSEP Removed " + x + "--- " + y + " sepset = " + sepset);
                    }
                }
             //   System.out.println("Possible DSEP: " + (System.currentTimeMillis()-time));
              // System.out.println("Possible dsep done");


//            long time2 = System.currentTimeMillis();
//            logger.log("info", "Step C: " + (time2 - time1) / 1000. + "s");
//
//            // Step FCI D.
//            long time3 = System.currentTimeMillis();
//
//            System.out.println("Starting possible dsep search");
//            PossibleDsepFci possibleDSep = new PossibleDsepFci(graph, independenceTest);
//            possibleDSep.setMaxIndegree(getPossibleDsepDepth());
//            possibleDSep.setKnowledge(getKnowledge());
//            possibleDSep.setMaxPathLength(maxPathLength);
//            this.sepsets.addAll(possibleDSep.search());
//            long time4 = System.currentTimeMillis();
//            logger.log("info", "Step D: " + (time4 - time3) / 1000. + "s");
//            System.out.println("Starting possible dsep search");

                // Reorient all edges as o-o.
                graph.reorientAllWith(Endpoint.CIRCLE);

               // System.out.println("Reoriented with circles");
            }
            time = System.currentTimeMillis();
            addColliders(graph, pdsepChange, sp2,knowledge);
            //Add colliders again but no need to do max p-value independence test search on variables that haven't had a new adjacency set
           // System.out.println("Adding Colliders 2: " + (System.currentTimeMillis()-time));
           // System.out.println("Added colliders again");

            // Step CI C (Zhang's step F3.)
            long time5 = System.currentTimeMillis();

            long time6 = System.currentTimeMillis();
            logger.log("info", "Step CI C: " + (time6 - time5) / 1000. + "s");
            time = System.currentTimeMillis();
           // System.out.println(this.sepsets);
            final FciOrient fciOrient = new FciOrient(new SepsetsSet(this.sepsets, independenceTest),whyOrient);

            fciOrient.setCompleteRuleSetUsed(completeRuleSetUsed);
            fciOrient.setMaxPathLength(maxPathLength);
            fciOrient.setKnowledge(knowledge);
            fciOrient.doFinalOrientation(graph);
           // System.out.println("Final Orientation: " + (System.currentTimeMillis()-time));
           // System.out.println("Final orientation done");

            return graph;
        }

        /**
         * Step C of PC; orients colliders using specified sepset. That is, orients x *-* y *-* z as x *-> y <-* z just in
         * case y is in Sepset({x, z}).
         */
        public void findCollidersUsingSepsets(final SepsetProducer sepsetProducer, final Graph graph, boolean verbose, IKnowledge knowledge,final Map<Triple,Double> colliders) {
            //Map<Triple, Double> colliders = new ConcurrentHashMap<>();
            final HashMap<Integer,Integer> edgeSizes = new HashMap<>();

            final List<Node> nodes = graph.getNodes();
            int sum = 0;
            for(int i = 0; i < nodes.size();i++)
            {
                edgeSizes.put(i,graph.getAdjacentNodes(nodes.get(i)).size());
                sum+= graph.getAdjacentNodes(nodes.get(i)).size();
            }
            chunk = sum/(parallelism);
            chunk = chunk * factor;
            class colliderTask extends RecursiveTask<Boolean> {
                private double chunk;
                private int from;
                private int to;
                public colliderTask(double chunk,int from, int to)
                {
                    this.chunk = chunk;
                    this.from = from;
                    this.to = to;
                }
                protected Boolean compute()
                {
                    int numEdges = 0;
                    for(int j = from; j < to; j++)
                    {
                        numEdges+= edgeSizes.get(j);
                    }
                    if(to-from <= 2 || numEdges <= chunk)
                    {
                        for(int i = from; i < to; i ++)
                        {
                           final Node b = nodes.get(i);
                            final List<Node> adjacent = graph.getAdjacentNodes(b);
                            //only compare to nodes less than it in index
                            if(adjacent.size()<2)
                                continue;
                            ChoiceGenerator cg = new ChoiceGenerator(adjacent.size(),2);
                            int [] combination;
                            while((combination = cg.next()) != null)
                            {
                              final  Node a = adjacent.get(combination[0]);
                              final  Node c = adjacent.get(combination[1]);
                                if(graph.isAdjacentTo(a,c))
                                    continue;
                               // System.out.println("getting sepset for: " + a + "\t" + b + "\t" + c);
                                List<Node> sepset = sepsetProducer.getSepset(a,c);
                     //   System.out.println(sepset + "\t" + a + "\t" + b);
                                if (sepset!=null && !sepset.contains(b)) {


                                    //IndependenceTest test2 = new IndTestDSep(trueDag);
                                    //SepsetProducer sp2 = new SepsetsMaxScore(graph, test2, null, depth);
                                //    System.out.println(a + "\t" + b + "\t" + c + "\t" + sepset + "\t" + sepsetProducer.getPValue());
                                    colliders.put(new Triple(a, b, c), sepsetProducer.getPValue());

                                }
                            }
                        }
                        return true;
                    }
                    else
                    {
                        List<colliderTask> tasks = new ArrayList<>();
                        final int mid = (to+from)/2;
                        colliderTask t1 = new colliderTask(chunk,from,mid);
                        tasks.add(t1);
                        colliderTask t2 = new colliderTask(chunk,mid,to);
                        tasks.add(t2);
                        invokeAll(tasks);
                        return true;
                    }

                }
            }
            pool.invoke(new colliderTask(chunk,0,nodes.size()));


        }
    public void findCollidersUsingSepsets(final Map<Node,Integer> pdsep, final SepsetsMaxPValue sp, final Graph graph, boolean verbose, IKnowledge knowledge,final Map<Triple,Double> colliders) {
        //Map<Triple, Double> colliders = new ConcurrentHashMap<>();
        final HashMap<Integer,Integer> edgeSizes = new HashMap<>();

        final List<Node> nodes = graph.getNodes();
        int sum = 0;
        for(int i = 0; i < nodes.size();i++)
        {
            edgeSizes.put(i,graph.getAdjacentNodes(nodes.get(i)).size());
            sum+= graph.getAdjacentNodes(nodes.get(i)).size();
        }
        chunk = sum/(parallelism);
        chunk = chunk*factor;
        class colliderTask extends RecursiveTask<Boolean> {
            private double chunk;
            private int from;
            private int to;
            public colliderTask(double chunk,int from, int to)
            {
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }
            protected Boolean compute()
            {
                int numEdges = 0;
                for(int j = from; j < to; j++)
                {
                    numEdges+= edgeSizes.get(j);
                }
                if(to-from <= 2 || numEdges <= chunk)
                {
                    for(int i = from; i < to; i ++)
                    {
                        final Node b = nodes.get(i);
                        final List<Node> adjacent = graph.getAdjacentNodes(b);
                        //only compare to nodes less than it in index
                        if(adjacent.size()<2)
                            continue;
                        ChoiceGenerator cg = new ChoiceGenerator(adjacent.size(),2);
                        int [] combination;
                        while((combination = cg.next()) != null)
                        {
                            final  Node a = adjacent.get(combination[0]);
                            final  Node c = adjacent.get(combination[1]);
                            if(graph.isAdjacentTo(a,c))
                                continue;

                            List<Node> sepset = null;
                            double p = 0;
                            if(pdsep.get(a)!=null&&pdsep.get(c)!=null)
                            {
                                sepset = sp.testSepset(a,c);
                                p = sp.getPValue();
                            }
                        else {
                                sepset = sp.getSepset(a, c);
                                p = sp.getPValue();
                            }


                            if (sepset!=null && !sepset.contains(b)) {


                                //IndependenceTest test2 = new IndTestDSep(trueDag);
                                //SepsetProducer sp2 = new SepsetsMaxScore(graph, test2, null, depth);
                                //System.out.println(sepset);
                                colliders.put(new Triple(a, b, c), p);

                            }
                        }
                    }
                    return true;
                }
                else
                {
                    List<colliderTask> tasks = new ArrayList<>();
                    final int mid = (to+from)/2;
                    colliderTask t1 = new colliderTask(chunk,from,mid);
                    tasks.add(t1);
                    colliderTask t2 = new colliderTask(chunk,mid,to);
                    tasks.add(t2);
                    invokeAll(tasks);
                    return true;
                }

            }
        }
        pool.invoke(new colliderTask(chunk,0,nodes.size()));


    }
        private void addColliders(Graph graph, final Map<Node,Integer> pdsep, final SepsetsMaxPValue sp, IKnowledge knowledge)
        {
            final Map<Triple,Double> collidersPs = new ConcurrentHashMap<>();
            findCollidersUsingSepsets(pdsep,sp, graph, verbose, knowledge, collidersPs);

            final List<Triple> colliders = new ArrayList<>(collidersPs.keySet());

//If we have multiple decisions for a particular collider, keep anyone that said that there is a collider, since this will have a higher p-value
            //  System.out.println(collidersPs);
            for (Triple collider : colliders) {
              //  System.out.println(collidersPs.get(collider));
                if (collidersPs.get(collider) < getIndependenceTest().getAlpha()) continue;

                Node a = collider.getX();
                Node b = collider.getY();
                Node c = collider.getZ();
                //System.out.println("Collider from " + a + " and " + c + " to " + b);
                //System.out.println("PVal: " + collidersPs.get(collider));
             //   System.out.println("Is " + a + " to " + b + " allowed?");
             //   System.out.println(!isArrowpointAllowed(a,b,knowledge));
             //   System.out.println("Is " + c + " to " + b + " allowed?");
             //   System.out.println(!isArrowpointAllowed(c,b,knowledge));


            if (!(isArrowpointAllowed(a, b, knowledge) && isArrowpointAllowed(c, b, knowledge))) {
                continue;
            }

//            if (!graph.getEdge(a, b).pointsTowards(a) && !graph.getEdge(b, c).pointsTowards(c)) {
//                graph.removeEdge(a, b);
//                graph.removeEdge(c, b);
//                graph.addDirectedEdge(a, b);
//                graph.addDirectedEdge(c, b);
//            }

                graph.setEndpoint(a, b, Endpoint.ARROW);
                graph.setEndpoint(c, b, Endpoint.ARROW);
                List<Node> temp = sp.getSepset(a,c);
                if(temp==null) {
                    whyOrient.put(a.getName() + "," + b.getName(), "0,null");
                    whyOrient.put(c.getName() + "," + b.getName(), "0,null");
                }
                else
                {
                    String x = "0";
                    for(Node p:temp)
                        x+=(","+p.getName());
                    whyOrient.put(a.getName()+","+b.getName(),x);
                    whyOrient.put(c.getName()+","+b.getName(),x);
                }
            }
        }

        private void addColliders(Graph graph, final SepsetProducer sepsetProducer, IKnowledge knowledge) {
            final Map<Triple,Double> collidersPs = new ConcurrentHashMap<>();

            findCollidersUsingSepsets(sepsetProducer, graph, verbose, knowledge, collidersPs);
        //    System.out.println("Found Colliders Using sepsets");
            final List<Triple> colliders = new ArrayList<>(collidersPs.keySet());

//If we have multiple decisions for a particular collider, keep anyone that said that there is a collider, since this will have a higher p-value
          //  System.out.println(collidersPs);
            for (Triple collider : colliders) {
                if (collidersPs.get(collider) < getIndependenceTest().getAlpha()) continue;

                Node a = collider.getX();
                Node b = collider.getY();
                Node c = collider.getZ();
                //System.out.println("Collider from " + a + " and " + c + " to " + b);
                //System.out.println("PVal: " + collidersPs.get(collider));
           if (!(isArrowpointAllowed(a, b, knowledge) && isArrowpointAllowed(c, b, knowledge))) {
                continue;
            }

//            if (!graph.getEdge(a, b).pointsTowards(a) && !graph.getEdge(b, c).pointsTowards(c)) {
//                graph.removeEdge(a, b);
//                graph.removeEdge(c, b);
//                graph.addDirectedEdge(a, b);
//                graph.addDirectedEdge(c, b);
//            }

                graph.setEndpoint(a, b, Endpoint.ARROW);
                graph.setEndpoint(c, b, Endpoint.ARROW);
            }
        }

        private static List<Node> union(List<Node> nodes, Node a) {
            List<Node> union = new ArrayList<>(nodes);
            union.add(a);
            return union;
        }


        public SepsetMap getSepsets() {
            return this.sepsets;
        }

        public IKnowledge getKnowledge() {
            return knowledge;
        }

        public void setKnowledge(IKnowledge knowledge) {
            if (knowledge == null) {
                throw new NullPointerException();
            }

            this.knowledge = knowledge;
        }

        /**
         * @return true if Zhang's complete rule set should be used, false if only R1-R4 (the rule set of the original FCI)
         * should be used. False by default.
         */
        public boolean isCompleteRuleSetUsed() {
            return completeRuleSetUsed;
        }

        /**
         * @param completeRuleSetUsed set to true if Zhang's complete rule set should be used, false if only R1-R4 (the rule
         *                            set of the original FCI) should be used. False by default.
         */
        public void setCompleteRuleSetUsed(boolean completeRuleSetUsed) {
            this.completeRuleSetUsed = completeRuleSetUsed;
        }

        public boolean isPossibleDsepSearchDone() {
            return possibleDsepSearchDone;
        }

        public void setPossibleDsepSearchDone(boolean possibleDsepSearchDone) {
            this.possibleDsepSearchDone = possibleDsepSearchDone;
        }

        /**
         * @return the maximum length of any discriminating path, or -1 of unlimited.
         */
        public int getMaxPathLength() {
            return maxPathLength == Integer.MAX_VALUE ? -1 : maxPathLength;
        }

        /**
         * @param maxPathLength the maximum length of any discriminating path, or -1 if unlimited.
         */
        public void setMaxPathLength(int maxPathLength) {
            if (maxPathLength < -1) {
                throw new IllegalArgumentException("Max path length must be -1 (unlimited) or >= 0: " + maxPathLength);
            }

            this.maxPathLength = maxPathLength;
        }

        /**
         * True iff verbose output should be printed.
         */
        public boolean isVerbose() {
            return verbose;
        }

        public void setVerbose(boolean verbose) {
            this.verbose = verbose;
        }

        /**
         * The independence test.
         */
        public IndependenceTest getIndependenceTest() {
            return independenceTest;
        }

        public void setTrueDag(Graph trueDag) {
            this.trueDag = trueDag;
        }

        public double getPenaltyDiscount() {
            return penaltyDiscount;
        }

        public void setPenaltyDiscount(double penaltyDiscount) {
            this.penaltyDiscount = penaltyDiscount;
        }

        //===========================PRIVATE METHODS=========================//

        private void buildIndexing(List<Node> nodes) {
            this.hashIndices = new ConcurrentHashMap<Node, Integer>();
            for (Node node : nodes) {
                this.hashIndices.put(node, variables.indexOf(node));
            }
        }

        /**
         * Orients according to background knowledge
         */
        private void fciOrientbk(IKnowledge bk, Graph graph, List<Node> variables) {
            logger.log("info", "Starting BK Orientation.");

            for (Iterator<KnowledgeEdge> it =
                 bk.forbiddenEdgesIterator(); it.hasNext(); ) {
                KnowledgeEdge edge = it.next();

                //match strings to variables in the graph.
                Node from = SearchGraphUtils.translate(edge.getFrom(), variables);
                Node to = SearchGraphUtils.translate(edge.getTo(), variables);


                if (from == null || to == null) {
                    continue;
                }

                if (graph.getEdge(from, to) == null) {
                    continue;
                }

                // Orient to*->from
                graph.setEndpoint(to, from, Endpoint.ARROW);
                graph.setEndpoint(from, to, Endpoint.CIRCLE);
                logger.log("knowledgeOrientation", SearchLogUtils.edgeOrientedMsg("Knowledge", graph.getEdge(from, to)));
            }

            for (Iterator<KnowledgeEdge> it =
                 bk.requiredEdgesIterator(); it.hasNext(); ) {
                KnowledgeEdge edge = it.next();

                //match strings to variables in this graph
                Node from = SearchGraphUtils.translate(edge.getFrom(), variables);
                Node to = SearchGraphUtils.translate(edge.getTo(), variables);

                if (from == null || to == null) {
                    continue;
                }

                if (graph.getEdge(from, to) == null) {
                    continue;
                }

                graph.setEndpoint(to, from, Endpoint.TAIL);
                graph.setEndpoint(from, to, Endpoint.ARROW);
                logger.log("knowledgeOrientation", SearchLogUtils.edgeOrientedMsg("Knowledge", graph.getEdge(from, to)));
            }

            logger.log("info", "Finishing BK Orientation.");
        }

        public int getPossibleDsepDepth() {
            return possibleDsepDepth;
        }

        public void setPossibleDsepDepth(int possibleDsepDepth) {
            this.possibleDsepDepth = possibleDsepDepth;
        }
    }
