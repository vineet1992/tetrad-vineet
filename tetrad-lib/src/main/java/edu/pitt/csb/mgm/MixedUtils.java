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

package edu.pitt.csb.mgm;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.*;//IndependenceTest;
import edu.cmu.tetrad.sem.GeneralizedSemIm;
import edu.cmu.tetrad.sem.GeneralizedSemPm;
import edu.cmu.tetrad.sem.TemplateExpander;
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.StatUtils;


import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Constructor;
import java.text.ParseException;
import java.util.*;

/**
 * Created by ajsedgewick on 7/29/15.
 */
public class MixedUtils {
    public static DataSet loadDataSet2(String filename,DelimiterType d) throws IOException {
        File file = new File(filename);
        DataReader reader = new DataReader();
        reader.setVariablesSupplied(true);
        reader.setMaxIntegralDiscrete(5);
        reader.setDelimiter(d);
        return reader.parseTabular(file);
    }
    public static DataSet loadDataSet2(String filename,DelimiterType d, int maxCat) throws IOException {
        File file = new File(filename);
        DataReader reader = new DataReader();
        reader.setVariablesSupplied(true);
        reader.setMaxIntegralDiscrete(maxCat);
        reader.setDelimiter(d);
        return reader.parseTabular(file);
    }
    public static DataSet loadDataSet2(String filename,int maxDiscrete) throws IOException {
        File file = new File(filename);
        DataReader reader = new DataReader();
        reader.setVariablesSupplied(true);
        reader.setMaxIntegralDiscrete(maxDiscrete);
        return reader.parseTabular(file);
    }

    public static DataSet loadDataSet2(String filename) throws IOException {
        File file = new File(filename);
        DataReader reader = new DataReader();
        reader.setVariablesSupplied(true);
        reader.setMaxIntegralDiscrete(5);
        return reader.parseTabular(file);
    }


    public static DataSet getNonparanormalTransform(DataSet data)
    {
        DataSet cont = DataUtils.getNonparanormalTransformed(MixedUtils.getContinousData(data));
        int ix = 0;
        for(int i = 0; i < data.getNumColumns();i++)
        {
            if(data.getVariable(i)instanceof ContinuousVariable) {
                for (int j = 0; j < data.getNumRows(); j++) {
                    data.setDouble(j,i,cont.getDouble(j,ix));
                }
                ix++;
            }
        }
        return data;
    }
    public static DataSet completeCases(DataSet d)
    {
        ArrayList<Integer> rowsToRemove = new ArrayList<Integer>();
        A:for(int i = 0; i < d.getNumRows();i++)
        {
            for(int j = 0; j < d.getNumColumns();j++)
            {
              if(d.getInt(i,j)==-99)
              {
                  rowsToRemove.add(i);
                  continue A;
              }
            }
        }
        int [] remove = new int[rowsToRemove.size()];
        for(int i = 0; i < rowsToRemove.size();i++)
            remove[i] = rowsToRemove.get(i);
        d.removeRows(remove);
        return d;

    }
    public static int[] getDiscreteInds(List<Node> nodes){
        List<Integer> indList = new ArrayList<>();
        int curInd = 0;
        for(Node n: nodes){
            if(n instanceof DiscreteVariable){
                indList.add(curInd);
            }
            curInd++;
        }

        int[] inds = new int[indList.size()];
        for(int i = 0; i < inds.length; i++){
            inds[i] = indList.get(i);
        }
        return inds;
    }

    public static int[] getContinuousInds(List<Node> nodes){
        List<Integer> indList = new ArrayList<>();
        int curInd = 0;
        for(Node n: nodes){
            if(n instanceof ContinuousVariable){
                indList.add(curInd);
            }
            curInd++;
        }

        int[] inds = new int[indList.size()];
        for(int i = 0; i < inds.length; i++){
            inds[i] = indList.get(i);
        }
        return inds;
    }

    //Converts a Dataset with both ContinuousVariables and DiscreteVariables to only ContinuousVariables
    public static DataSet makeContinuousData(DataSet dsMix) {
        ArrayList<Node> contVars = new ArrayList<Node>();
        for(Node n: dsMix.getVariables()){
            if(n instanceof DiscreteVariable){
                ContinuousVariable nc = new ContinuousVariable(n.getName());
                contVars.add(nc);
            } else {
                contVars.add(n);
            }
        }

        return ColtDataSet.makeData(contVars, dsMix.getDoubleData());
    }

    //takes DataSet of all ContinuousVariables
    //convert variables to discrete if there is an entry with <NodeName, "Disc"> in nodeDists
    public static DataSet makeMixedData(DataSet dsCont,  Map<String, String> nodeDists, int numCategories) {
        ArrayList<Node> mixVars = new ArrayList<Node>();
        for(Node n: dsCont.getVariables()){
            if(nodeDists.get(n.getName()).equals("Disc")){
                DiscreteVariable nd = new DiscreteVariable(n.getName(), numCategories);
                mixVars.add(nd);
            } else {
                mixVars.add(n);
            }
        }

        return ColtDataSet.makeData(mixVars, dsCont.getDoubleData());
    }

    //takes DataSet of all ContinuousVariables
    //convert variables to discrete if there is an entry with <NodeName, x> with x > 0, num categories set to x
    public static DataSet makeMixedData(DataSet dsCont,  Map<String, Integer> nodeDists) {
        ArrayList<Node> mixVars = new ArrayList<Node>();
        for(Node n: dsCont.getVariables()){
            int nC = nodeDists.get(n.getName());
            if(nC>0){
                DiscreteVariable nd = new DiscreteVariable(n.getName(), nC);
                mixVars.add(nd);
            } else {
                mixVars.add(n);
            }
        }

//        MixedDataBox box = new MixedDataBox(mixVars, dsCont.getNumRows());
//
//        for (int i = 0; i < box.numRows(); i++) {
//            for (int j = 0; j < box.numCols(); j++) {
//                if (mixVars.get(j) instanceof ContinuousVariable) {
//                    box.set(i, j, dsCont.getDouble(i, j));
//                } else {
//                    box.set(i, j, (int) Math.round(dsCont.getDouble(i, j)));
//                }
//            }
//        }

//        return new BoxDataSet(box, mixVars);

        return ColtDataSet.makeData(mixVars, dsCont.getDoubleData());
    }

    /**
     * Makes a deep copy of a dataset (Nodes copied as well). Useful for paralellization
     * @param ds dataset to be copied
     * @return
     */
    public static ColtDataSet deepCopy(ColtDataSet ds){
        List<Node> vars = new ArrayList<>(ds.getNumColumns());
        for(Node n : ds.getVariables()){
            if(n instanceof ContinuousVariable)
                vars.add(new ContinuousVariable((ContinuousVariable)n));
            else if (n instanceof DiscreteVariable)
                vars.add(new DiscreteVariable((DiscreteVariable) n));
            else
                throw new IllegalArgumentException("Variable type of node " + n + "could not be determined");
        }
        ColtDataSet out = ColtDataSet.makeData(vars,ds.getDoubleData());

        return out;
    }

    //Takes a mixed dataset and returns only data corresponding to ContinuousVariables in order
    public static DataSet getContinousData(DataSet ds){
        ArrayList<Node> contVars = new ArrayList<Node>();
        for(Node n: ds.getVariables()){
            if(n instanceof ContinuousVariable)
                contVars.add(n);
        }
        return ds.subsetColumns(contVars);
    }

    //Takes a mixed dataset and returns only data corresponding to DiscreteVariables in order
    public static DataSet getDiscreteData(DataSet ds){
        ArrayList<Node> discVars = new ArrayList<Node>();
        for(Node n: ds.getVariables()){
            if(n instanceof DiscreteVariable)
                discVars.add(n);
        }
        return ds.subsetColumns(discVars);
    }

    public static int[] getDiscLevels(DataSet ds){
        //ArrayList<Integer> levels = new ArrayList<Integer>[];
        DataSet discDs = getDiscreteData(ds);
        int[] levels = new int[discDs.getNumColumns()];
        int i = 0;
        for(Node n: discDs.getVariables()){
            levels[i] = ((DiscreteVariable) n).getNumCategories();
            i++;
        }
        return levels;
    }

    /**
     * return vector of the maximum of each column in m (as ints, i.e. for discrete data)
     * @param m
     * @return
     */
    public static int[] colMax(DoubleMatrix2D m){
        int[] maxVec = new int[m.columns()];
        for(int i = 0; i < m.columns(); i++){
            double curmax = -1;
            for(int j = 0; j < m.rows(); j++){
                double curval = m.getQuick(j,i);
                if(curval>curmax){
                    curmax = curval;
                }
            }
            maxVec[i] = (int) curmax;
        }
        return maxVec;
    }

    public static double vecMax(DoubleMatrix1D vec){
        double curMax = Double.NEGATIVE_INFINITY;
        for(int i = 0; i < vec.size(); i++){
            double curVal = vec.getQuick(i);
            if(curVal>curMax){
                curMax = curVal;
            }
        }
        return curMax;
    }

    public static double numVals(DoubleMatrix1D vec){
        return valSet(vec).size();
    }

    public static Set<Double> valSet(DoubleMatrix1D vec){
        Set<Double> vals = new HashSet<>();
        for(int i = 0; i < vec.size(); i++){
            vals.add(vec.getQuick(i));
        }
        return vals;
    }

    //generate PM from trueGraph for mixed Gaussian and Trinary variables
    //Don't use, buggy
    public static GeneralizedSemPm GaussianTrinaryPm(Graph trueGraph, HashMap<String, String> nodeDists, int maxSample, String paramTemplate) throws IllegalStateException{

        GeneralizedSemPm semPm = new GeneralizedSemPm(trueGraph);
        try {
            List<Node> variableNodes = semPm.getVariableNodes();
            int numVars = variableNodes.size();


            semPm.setStartsWithParametersTemplate("B", paramTemplate);
            semPm.setStartsWithParametersTemplate("D", paramTemplate);

            // empirically should give us a stddev of 1 - 2
            semPm.setStartsWithParametersTemplate("al", "U(.3,1.3)");
            semPm.setStartsWithParametersTemplate("s", "U(1,2)");

            //if we don't use NB error, we could do this instead
            String templateDisc = "DiscError(err, (TSUM(NEW(B)*$)), (TSUM(NEW(B)*$)), (TSUM(NEW(B)*$)))";
            //String templateDisc0 = "DiscError(err, 2,2,2)";
            String templateDisc0 = "DiscError(err, .001,.001,.001)";


            for (Node node : variableNodes) {

                List<Node> parents = trueGraph.getParents(node);
                System.out.println("Parents:" + parents);
                //System.out.println("nParents: " + parents.size() );
                Node eNode = semPm.getErrorNode(node);

                //normal and nb work like normal sems
                String curEx = semPm.getNodeExpressionString(node);
                String errEx = semPm.getNodeExpressionString(eNode);
                String newTemp = "";

                //System.out.println("Node: " + node + "Type: " + nodeDists.get(node));

                if(nodeDists.get(node.getName()).equals("Disc")){
                    if(parents.size() == 0){
                        newTemp = templateDisc0;
                    } else {
                        newTemp = templateDisc;
                    }
                    newTemp = newTemp.replaceAll("err", eNode.getName());
                    curEx = TemplateExpander.getInstance().expandTemplate(newTemp, semPm, node);
                    //System.out.println("Disc CurEx: " + curEx);
                    errEx = TemplateExpander.getInstance().expandTemplate("U(0,1)", semPm, eNode);
                }

                //now for every discrete parent, swap for discrete params
                newTemp = "";
                if(parents.size() != 0) {
                    for (Node parNode : parents){
                        if(nodeDists.get(parNode.getName()).equals("Disc")){
                            //String curName = trueGraph.getParents(node).get(0).toString();
                            String curName = parNode.getName();
                            String disRep = "IF(" + curName + "=0,NEW(D),IF("+curName+"=1,NEW(D),NEW(D)))";
                            newTemp = curEx.replaceAll("(B[0-9]*\\*" + curName + ")(?![0-9])", disRep);
                        }
                    }
                }

                if(newTemp.length()!=0){
                    curEx = TemplateExpander.getInstance().expandTemplate(newTemp, semPm, node);
                }

                semPm.setNodeExpression(node, curEx);
                semPm.setNodeExpression(eNode, errEx);
                System.out.println("CurExp: " + curEx);
                System.out.println("ErrExp: " + errEx);
            }
        } catch (ParseException e) {
            throw new IllegalStateException("Parse error in fixing parameters.", e);
        }

        return semPm;
    }

    //generate PM from trueGraph for mixed Gaussian and Categorical variables
    //public static GeneralizedSemPm GaussianCategoricalPm(Graph trueGraph, HashMap<String, Integer> nodeDists, String paramTemplate) throws IllegalStateException{
    public static GeneralizedSemPm GaussianCategoricalPm(Graph trueGraph, String paramTemplate) throws IllegalStateException{

        Map<String, Integer> nodeDists = getNodeDists(trueGraph);

        GeneralizedSemPm semPm = new GeneralizedSemPm(trueGraph);
        try {
            List<Node> variableNodes = semPm.getVariableNodes();
            int numVars = variableNodes.size();


            semPm.setStartsWithParametersTemplate("B", paramTemplate);
            semPm.setStartsWithParametersTemplate("C", paramTemplate);
            semPm.setStartsWithParametersTemplate("D", paramTemplate);

            // empirically should give us a stddev of 1 - 2
            semPm.setStartsWithParametersTemplate("s", "U(1,2)");

            //if we don't use NB error, we could do this instead
            //String templateDisc = "DiscError(err, (TSUM(NEW(B)*$)), (TSUM(NEW(B)*$)), (TSUM(NEW(B)*$)))";
//            String templateDisc0 = "DiscError(err, 1,1,1)";

            String templateDisc0 = "DiscError(err, ";

            for (Node node : variableNodes) {

                List<Node> parents = trueGraph.getParents(node);
                //System.out.println("nParents: " + parents.size() );
                Node eNode = semPm.getErrorNode(node);

                //normal and nb work like normal sems
                String curEx = semPm.getNodeExpressionString(node);
                String errEx = semPm.getNodeExpressionString(eNode);
                String newTemp = "";

                //System.out.println("Node: " + node + "Type: " + nodeDists.get(node));

                //dist of 0 means Gaussian
                int curDist = nodeDists.get(node.getName());
                if(curDist == 1)
                    throw new IllegalArgumentException("Dist for node " + node.getName() + " is set to one (i.e. constant) which is not supported.");


                //for each discrete node use DiscError for categorical draw
                if(curDist>0){
                    if(parents.size() == 0){
                        newTemp = "DiscError(err";
                        for(int l = 0; l < curDist; l++){
                            newTemp += ",1";
                        }
                        newTemp += ")";
//                        newTemp = templateDisc0;
                    } else {
                        newTemp = "DiscError(err";
                        for(int l = 0; l < curDist; l++){
                            newTemp += ", TSUM(NEW(C)*$)";
                        }
                        newTemp += ")";
                    }
                    newTemp = newTemp.replaceAll("err", eNode.getName());
                    curEx = TemplateExpander.getInstance().expandTemplate(newTemp, semPm, node);
                    //System.out.println("Disc CurEx: " + curEx);
                    errEx = TemplateExpander.getInstance().expandTemplate("U(0,1)", semPm, eNode);
                }

                //now for every discrete parent, swap for discrete params
                newTemp = curEx;
                if(parents.size() != 0) {
                    for (Node parNode : parents){
                        int parDist = nodeDists.get(parNode.getName());

                        if(parDist>0){
                            //String curName = trueGraph.getParents(node).get(0).toString();
                            String curName = parNode.getName();
                            String disRep = "Switch(" + curName;
                            for(int l = 0; l < parDist; l++){
                                if(curDist>0) {
                                    disRep += ",NEW(D)";
                                } else {
                                    disRep += ",NEW(C)";
                                }
                            }
                            disRep += ")";

                            //replaces BX * curName with new discrete expression
                            if(curDist > 0){
                                newTemp = newTemp.replaceAll("(C[0-9]*\\*" + curName + ")(?![0-9])", disRep);
                            } else {
                                newTemp = newTemp.replaceAll("(B[0-9]*\\*" + curName + ")(?![0-9])", disRep);
                            }
                        }
                    }
                }

                if(newTemp.length()!=0){
                    //System.out.println(newTemp);
                    curEx = TemplateExpander.getInstance().expandTemplate(newTemp, semPm, node);
                }

                semPm.setNodeExpression(node, curEx);
                semPm.setNodeExpression(eNode, errEx);
            }
        } catch (ParseException e) {
            throw new IllegalStateException("Parse error in fixing parameters.", e);
        }

        return semPm;
    }

    /**
     * Set all existing parameters that begins with sta to template and also set template for any new parameters
     *
     * @param sta
     * @param template
     * @param pm
     * @return
     */
    public static void setStartsWith(String sta, String template, GeneralizedSemPm pm){
        try {
            pm.setStartsWithParametersTemplate(sta, template);
            for (String param : pm.getParameters()) {
                if (param.startsWith(sta)) {
                    pm.setParameterExpression(param, template);
                }
            }
        } catch(Throwable t){
            t.printStackTrace();
        }
        return;
    }

    //legacy
    public static GeneralizedSemIm GaussianCategoricalIm(GeneralizedSemPm pm){
        return GaussianCategoricalIm(pm, true);
    }

    /**
     *    This method is needed to normalize edge parameters for an Instantiated Mixed Model
     *    Generates edge parameters for c-d and d-d edges from a single weight, abs(w), drawn by the normal IM constructor.
     *    Abs(w) is used for d-d edges.
     *
     *    For deterministic, c-d are evenly spaced between -w and w, and d-d are a matrix with w on the diagonal and
     *    -w/(categories-1) in the rest.
     *    For random, c-d params are uniformly drawn from 0 to 1 then transformed to have w as max value and sum to 0.
     *
     * @param pm
     * @param discParamRand true for random edge generation behavior, false for deterministic
     * @return
     */
    public static GeneralizedSemIm GaussianCategoricalIm(GeneralizedSemPm pm, boolean discParamRand){

        Map<String, Integer> nodeDists = getNodeDists(pm.getGraph());

        GeneralizedSemIm im = new GeneralizedSemIm(pm);
        //System.out.println(im);
        List<Node> nodes = pm.getVariableNodes();

        //this needs to be changed for cyclic graphs...
        for(Node n: nodes){
            Set<Node> parNodes = pm.getReferencedNodes(n);
            if(parNodes.size()==0){
                continue;
            }
            for(Node par: parNodes){
                if(par.getNodeType()==NodeType.ERROR){
                    continue;
                }
                int cL = nodeDists.get(n.getName());
                int pL = nodeDists.get(par.getName());

                // c-c edges don't need params changed
                if(cL==0 && pL==0){
                    continue;
                }

                List<String> params = getEdgeParams(n, par, pm);
                // just use the first parameter as the "weight" for the whole edge
                double w = im.getParameterValue(params.get(0));
                // double[] newWeights;

                // d-d edges use one vector and permute edges, could use different strategy
                if(cL > 0 && pL > 0) {
                    double[][] newWeights = new double[cL][pL];
                    //List<Integer> indices = new ArrayList<Integer>(pL);
                    //PermutationGenerator pg = new PermutationGenerator(pL);
                    //int[] permInd = pg.next();
                    w = Math.abs(w);
                    double bgW = w/((double) pL - 1.0);
                    double[] weightVals;

                    /*if(discParamRand)
                        weightVals = generateMixedEdgeParams(w, pL);
                    else
                        weightVals = evenSplitVector(w, pL);
                    */
                    int [] weightInds = new int[cL];
                    for(int i = 0; i < cL; i++){
                        if(i < pL)
                            weightInds[i] = i;
                        else
                            weightInds[i] = i % pL;
                    }

                    if(discParamRand)
                        weightInds = arrayPermute(weightInds);


                    for(int i = 0; i < cL; i++){
                        for(int j = 0; j < pL; j++){
                            int index = i*pL + j;
                            if(weightInds[i]==j)
                                im.setParameterValue(params.get(index), w);
                            else
                                im.setParameterValue(params.get(index), -bgW);
                        }
                    }
                    //params for c-d edges
                } else {
                    double[] newWeights;
                    int curL = (pL > 0 ? pL: cL);
                    if(discParamRand)
                        newWeights = generateMixedEdgeParams(w, curL);
                    else
                        newWeights = evenSplitVector(w, curL);

                    int count = 0;
                    for(String p : params){
                        im.setParameterValue(p, newWeights[count]);
                        count++;
                    }
                }

            }
            //pm.

            //if(p.startsWith("B")){
            //    continue;
            //} else if(p.startsWith())
        }


        return im;
    }

    //Given two node names and a parameterized model return list of parameters corresponding to edge between them
    public static List<String> getEdgeParams(String s1, String s2, GeneralizedSemPm pm){
        Node n1 = pm.getNode(s1);
        Node n2 = pm.getNode(s2);
        return getEdgeParams(n1, n2, pm);
    }

    //randomly permute an array of doubles
    public static double[] arrayPermute(double[] a){
        double[] out = new double[a.length];
        List<Double> l = new ArrayList<Double>(a.length);
        for(int i =0; i < a.length; i++){
            l.add(i, a[i]);
        }
        Collections.shuffle(l);
        for(int i =0; i < a.length; i++){
            out[i] = l.get(i);
        }
        return out;
    }

    //randomly permute array of ints
    public static int[] arrayPermute(int[] a){
        int[] out = new int[a.length];
        List<Integer> l = new ArrayList<Integer>(a.length);
        for(int i =0; i < a.length; i++){
            l.add(i, a[i]);
        }
        Collections.shuffle(l);
        for(int i =0; i < a.length; i++){
            out[i] = l.get(i);
        }
        return out;
    }

    //generates a vector of length L that starts with -w and increases with consistent steps to w
    public static double[] evenSplitVector(double w, int L){
        double[] vec = new double[L];
        double step = 2.0*w/(L-1.0);
        for(int i = 0; i < L; i++){
            vec[i] = -w + i*step;
        }
        return vec;
    }

    //Given two nodes and a parameterized model return list of parameters corresponding to edge between them
    public static List<String> getEdgeParams(Node n1, Node n2, GeneralizedSemPm pm){
        //there may be a better way to do this using recursive calls of Expression.getExpressions
        Set<String> allParams = pm.getParameters();

        Node child;
        Node parent;
        if(pm.getReferencedNodes(n1).contains(n2)){
            child = n1;
            parent = n2;
        } else if (pm.getReferencedNodes(n2).contains(n1)){
            child = n2;
            parent = n1;
        } else {
            return null;
        }

        java.util.regex.Pattern parPat;
        if(parent instanceof DiscreteVariable){
            parPat = java.util.regex.Pattern.compile("Switch\\(" + parent.getName() + ",.*?\\)");
        } else {
            parPat = java.util.regex.Pattern.compile("([BC][0-9]*\\*" + parent.getName() + ")(?![0-9])");
        }

        ArrayList<String> paramList = new ArrayList<>();
        String ex = pm.getNodeExpressionString(child);
        java.util.regex.Matcher mat = parPat.matcher(ex);
        while(mat.find()){
            String curGroup = mat.group();
            if(parent instanceof DiscreteVariable){
                curGroup = curGroup.substring(("Switch(" + parent.getName()).length()+1, curGroup.length()-1);
                String[] pars = curGroup.split(",");
                for(String p : pars){
                    //if(!allParams.contains(p))
                    //    throw exception;
                    paramList.add(p);
                }
            } else{
                String p = curGroup.split("\\*")[0];
                paramList.add(p);
            }
        }
        //ex.
        //if(child instanceof DiscreteVariable){
        //    if(parent instanceof DiscreteVariable)
        //}

        /*Expression exp = pm.getNodeExpression(child);
        List<Expression> test = exp.getExpressions();
        for(Expression t : test){
            List<Expression> test2 = t.getExpressions();
            for(Expression t2: test2) {
                System.out.println(t2.toString());
            }
        }*/

        return paramList;
    }

    //generates a vector of length L with maximum value w that sums to 0
    public static double[] generateMixedEdgeParams(double w, int L){
        double[] vec = new double[L];
        RandomUtil ru = RandomUtil.getInstance();

        for(int i=0; i < L; i++){
            vec[i] = ru.nextUniform(0, 1);
        }

        double vMean = StatUtils.mean(vec);
        double vMax = 0;
        for(int i=0; i < L; i++){
            vec[i] = vec[i] - vMean;
            if(Math.abs(vec[i])> Math.abs(vMax))
                vMax = vec[i];
        }

        double scale = w/vMax;
        //maintain sign of w;
        if(vMax<0)
            scale *=-1;

        for(int i=0; i < L; i++){
            vec[i] *= scale;
        }

        return vec;
    }

    //labels corresponding to values from allEdgeStats
    public static final String EdgeStatHeader = "TD\tTU\tFL\tFD\tFU\tFPD\tFPU\tFND\tFNU\tBidir";

    //assumes Graphs have properly assigned variable types
    public static int[][] allEdgeStats(Graph pT, Graph pE){
        HashMap<String, String> nd = new HashMap<String, String>();

        //Estimated graph more likely to have correct node types...
        for(Node n : pE.getNodes()){
            if(n instanceof DiscreteVariable){
                nd.put(n.getName(), "Disc");
            } else {
                nd.put(n.getName(), "Norm");
            }
        }
        return allEdgeStats(pT, pE, nd);
    }

    public static int [][] allEdgeStatsBioInf(Graph pT, Graph pE,DataSet data)
    {
        HashMap<String, String> nd = new HashMap<String, String>();

        //Estimated graph more likely to have correct node types...
        for(Node n : pE.getNodes()){
            if(data.getVariable(n.getName()) instanceof DiscreteVariable){
                nd.put(n.getName(), "Disc");
            } else {
                nd.put(n.getName(), "Norm");
            }
        }
        return allEdgeStatsBioInf(pT, pE, nd);
    }

    public static int [][] allEdgeStatsBioInf(Graph pT, Graph pE, HashMap<String,String> nodeDists)
    {
        int atp = 0;
        int afp = 1;
        int afn = 2;
        int dtp = 3;
        int dfp = 4;
        int dfn = 5;
        int shd = 6;
        int[][] stats = new int[4][7]; //4 Types by Adj TP, Adj FP, Adj FN, Dir TP, Dir FP, Dir FN, SHD
        for(int i=0; i<stats.length; i++){
            for(int j=0; j<stats[0].length; j++){
                stats[i][j] = 0;
            }
        }
        //enforce patterns?
        //Graph pT = SearchGraphUtils.patternFromDag(tg);
        //Graph pE = SearchGraphUtils.patternFromDag(eg);

        //check that variable names are the same...

        Set<Edge> edgesT = pT.getEdges();
        Set<Edge> edgesE = pE.getEdges();

        //differences += Math.abs(e1.size() - e2.size());

        //for (int i = 0; i < e1.size(); i++) {
        int edgeType;
        for(Edge eT: edgesT){
            Node n1 = pE.getNode(eT.getNode1().getName());
            Node n2 = pE.getNode(eT.getNode2().getName());
            if(nodeDists.get(n1.getName()).equals("Norm") && nodeDists.get(n2.getName()).equals("Norm")) {
                edgeType = 0;
            } else if(nodeDists.get(n1.getName()).equals("Disc") && nodeDists.get(n2.getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }

            Edge eE = pE.getEdge(n1, n2);
            if (eE == null) {
                if (eT.isDirected()) {
                    stats[edgeType][shd]+=2;
                    stats[edgeType][dfn]++; //False Negative Directed -- FND
                }
                else
                    stats[edgeType][shd]++;
                stats[edgeType][afn]++; //False Negative Undirected -- FNU
            } else if (eE.isDirected()){
                if (eT.isDirected() && eT.pointsTowards(eT.getNode1()) == eE.pointsTowards(n1)){
                    stats[edgeType][dtp]++; //True Directed -- TD
                } else if (eT.isDirected()){
                    stats[edgeType][dfp]++; //FLip
                    stats[edgeType][shd]++;
                } else {
                    stats[edgeType][shd]++;
                    stats[edgeType][dfp]++; //Falsely Directed -- FD
                }
                stats[edgeType][atp]++;
            } else { //so eE is undirected
                if(eT.isDirected()) {
                    stats[edgeType][shd]++;
                    stats[edgeType][dfn]++; //Falsely Undirected -- FU
                }
                stats[edgeType][atp]++;
            }
        }

        for(Edge eE: edgesE){
            Node n1 = pT.getNode(eE.getNode1().getName());
            Node n2 = pT.getNode(eE.getNode2().getName());

            if(nodeDists.get(n1.getName()).equals("Norm") && nodeDists.get(n2.getName()).equals("Norm")) {
                edgeType = 0;
            } else if(nodeDists.get(n1.getName()).equals("Disc") && nodeDists.get(n2.getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }

            Edge eT = pT.getEdge(n1, n2);
            if(eT == null) {
                if(eE.isDirected()){
                    stats[edgeType][shd]+=2;
                    stats[edgeType][dfp]++; //False Positive Directed -- FPD
                }
                else
                    stats[edgeType][shd]++;
                stats[edgeType][afp]++; //False Positive Undirected -- FUD
            }
        }
        for(int t = 0; t < stats[0].length;t++)
        {
            stats[3][t] = stats[0][t] + stats[1][t] + stats[2][t];
        }
        return stats;
    }


    // break out stats by node distributions, here only "Norm" and "Disc"
    // so three types of possible edges, cc, cd, dd, output is edge type by stat type
    // counts bidirected
    public static int[][] allEdgeStats(Graph pT, Graph pE, HashMap<String, String> nodeDists) {
        int[][] stats = new int[3][10];
        for(int i=0; i<stats.length; i++){
            for(int j=0; j<stats[0].length; j++){
                stats[i][j] = 0;
            }
        }
        //enforce patterns?
        //Graph pT = SearchGraphUtils.patternFromDag(tg);
        //Graph pE = SearchGraphUtils.patternFromDag(eg);

        //check that variable names are the same...

        Set<Edge> edgesT = pT.getEdges();
        Set<Edge> edgesE = pE.getEdges();

        //differences += Math.abs(e1.size() - e2.size());

        //for (int i = 0; i < e1.size(); i++) {
        int edgeType;
        for(Edge eT: edgesT){
            Node n1 = pE.getNode(eT.getNode1().getName());
            Node n2 = pE.getNode(eT.getNode2().getName());
            if(nodeDists.get(n1.getName()).equals("Norm") && nodeDists.get(n2.getName()).equals("Norm")) {
                edgeType = 0;
            } else if(nodeDists.get(n1.getName()).equals("Disc") && nodeDists.get(n2.getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }

            Edge eE = pE.getEdge(n1, n2);
            if (eE == null) {
                if (eT.isDirected()) {
                    stats[edgeType][7]++; //False Negative Directed -- FND
                } else {
                    stats[edgeType][8]++; //False Negative Undirected -- FNU
                }
            } else if (eE.isDirected()){
                if (eT.isDirected() && eT.pointsTowards(eT.getNode1()) == eE.pointsTowards(n1)){
                    stats[edgeType][0]++; //True Directed -- TD
                } else if (eT.isDirected()){
                    stats[edgeType][2]++; //FLip
                } else {
                    stats[edgeType][3]++; //Falsely Directed -- FD
                }
            } else { //so eE is undirected
                if(eT.isDirected()) {
                    stats[edgeType][4]++; //Falsely Undirected -- FU
                } else {
                    stats[edgeType][1]++; //True Undirected -- TU
                }
            }
        }

        for(Edge eE: edgesE){
            Node n1 = pT.getNode(eE.getNode1().getName());
            Node n2 = pT.getNode(eE.getNode2().getName());

            if(nodeDists.get(n1.getName()).equals("Norm") && nodeDists.get(n2.getName()).equals("Norm")) {
                edgeType = 0;
            } else if(nodeDists.get(n1.getName()).equals("Disc") && nodeDists.get(n2.getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }

            if(eE.getEndpoint1()== Endpoint.ARROW && eE.getEndpoint2()==Endpoint.ARROW)
                stats[edgeType][9]++; //bidirected

            Edge eT = pT.getEdge(n1, n2);
            if(eT == null) {
                if(eE.isDirected()){
                    stats[edgeType][5]++; //False Positive Directed -- FPD
                } else {
                    stats[edgeType][6]++; //False Positive Undirected -- FUD
                }
            }
        }
        return stats;
    }

    //Utils
    public static Graph makeMixedGraph(Graph g, Map<String, Integer> m){
        List<Node> nodes = g.getNodes();
        for(int i = 0; i < nodes.size(); i++){
            Node n = nodes.get(i);
            int nL = m.get(n.getName());
            if(nL > 0){
                Node nNew = new DiscreteVariable(n.getName(), nL);
                nodes.set(i, nNew);
            }
        }

        Graph outG = new EdgeListGraph(nodes);
        for(Edge e: g.getEdges()){
            Node n1 = e.getNode1();
            Node n2 = e.getNode2();
            Edge eNew = new Edge(outG.getNode(n1.getName()), outG.getNode(n2.getName()), e.getEndpoint1(), e.getEndpoint2());
            outG.addEdge(eNew);
        }

        return outG;
    }

    public static String stringFrom2dArray(int[][] arr){
        String outStr = "";
        for(int i = 0; i < arr.length; i++){
            for(int j = 0; j < arr[i].length; j++){
                outStr += Integer.toString(arr[i][j]);
                if(j != arr[i].length-1)
                    outStr += "\t";
            }
            outStr+="\n";
        }
        return outStr;
    }

    public static DataSet loadDataSet(String dir, String filename) throws IOException {
        File file = new File(dir, filename);
        DataReader reader = new DataReader();
        reader.setVariablesSupplied(true);
        reader.setMaxIntegralDiscrete(5);
        return reader.parseTabular(file);
    }

    public static DataSet loadDelim(String dir, String filename) throws IOException {
        File file = new File(dir, filename);
        DataReader reader = new DataReader();
        reader.setVariablesSupplied(false);
        return reader.parseTabular(file);
    }

    //Gives a map of number of categories of DiscreteVariables in g. ContinuousVariables are mapped to 0
    public static Map<String, Integer> getNodeDists(Graph g){
        HashMap<String, Integer> map = new HashMap<>();
        List<Node> nodes = g.getNodes();
        for(Node n: nodes){
            if(n instanceof DiscreteVariable)
                map.put(n.getName(), ((DiscreteVariable) n).getNumCategories());
            else
                map.put(n.getName(), 0);
        }
        return map;
    }

    public static DataSet loadData(String dir, String filename) throws IOException {
        File file = new File(dir, filename);
        DataReader reader = new DataReader();
        reader.setVariablesSupplied(true);
        return reader.parseTabular(file);
    }

    /**
     * Check each pair of variables to see if correlation is 1. WARNING: calculates correlation matrix, memory heavy
     * when there are lots of variables
     *
     * @param ds
     * @param verbose
     * @return
     */
    public static boolean isColinear(DataSet ds, boolean verbose) {
        List<Node> nodes = ds.getVariables();
        boolean isco = false;
        CorrelationMatrix cor = new CorrelationMatrix(makeContinuousData(ds));
        for(int i = 0; i < nodes.size(); i++){
            for(int j = i+1; j < nodes.size(); j++){
                if(cor.getValue(i,j) == 1){
                    if(verbose){
                        isco = true;
                        System.out.println("Colinearity found between: " + nodes.get(i).getName() + " and " + nodes.get(j).getName());
                    } else {
                        return true;
                    }
                }
            }
        }
        return isco;
    }

    public static DoubleMatrix2D graphToMatrix(Graph graph, double undirectedWeight, double directedWeight) {
        // initialize matrix
        int n = graph.getNumNodes();
        DoubleMatrix2D matrix = DoubleFactory2D.dense.make(n, n, 0.0);

        // map node names in order of appearance
        HashMap<Node, Integer> map = new HashMap<Node, Integer>();
        int i = 0;
        for (Node node : graph.getNodes()) {
            map.put(node, i);
            i++;
        }

        // mark edges
        for (Edge edge : graph.getEdges()) {
            // if directed find which is parent/child
            Node node1 = edge.getNode1();
            Node node2 = edge.getNode2();

            //treat bidirected as undirected...
            if (!edge.isDirected() || (edge.getEndpoint1()== Endpoint.ARROW && edge.getEndpoint2()==Endpoint.ARROW)) {
                matrix.set(map.get(node1), map.get(node2), undirectedWeight);
                matrix.set(map.get(node2), map.get(node1), undirectedWeight);
            } else {
                if (edge.pointsTowards(node1)) {
                    matrix.set(map.get(node2), map.get(node1), directedWeight);
                } else {
                    //if (edge.pointsTowards(node2)) {
                    matrix.set(map.get(node1), map.get(node2), directedWeight);
                }
            }
        }
        return matrix;
    }

    //returns undirected skeleton matrix (symmetric
    public static DoubleMatrix2D skeletonToMatrix(Graph graph, DataSet d) {
        // initialize matrix
        int n = graph.getNumNodes();
        DoubleMatrix2D matrix = DoubleFactory2D.dense.make(n, n, 0.0);

        // map node names in order of appearance
        HashMap<Node, Integer> map = new HashMap<Node, Integer>();
        for (Node node : graph.getNodes()) {
            map.put(node, d.getColumn(d.getVariable(node.getName())));
        }

        // mark edges
        for (Edge edge : graph.getEdges()) {
            // if directed find which is parent/child
            Node node1 = edge.getNode1();
            Node node2 = edge.getNode2();

            matrix.set(map.get(node1), map.get(node2), 1.0);
            matrix.set(map.get(node2), map.get(node1), 1.0);

        }

        return matrix;
    }

    public static DoubleMatrix2D graphToMatrix(Graph graph){
        return graphToMatrix(graph, 1, 1);
    }

    /**
     * Returns independence tests by name located in edu.cmu.tetrad.search and edu.pitt.csb.mgm
     * also supports shorthand for LRT ("lrt) and t-tests ("tlin" for prefer linear (fastest) or "tlog" for prefer logistic)
     * @param name
     * @return
     */
    public static IndependenceTest IndTestFromString(String name, DataSet data, double alpha) {

        IndependenceTest test = null;

        if (name.equals("lrt")) {
            test = new IndTestMixedRegressionLrt(data, alpha);
            //test = new IndTestMultinomialLogisticRegression(data, alpha);
        } else if (name.equals("tlin")) {
            test = new edu.pitt.csb.mgm.IndTestMixedMultipleTTest(data, alpha);
            ((edu.pitt.csb.mgm.IndTestMixedMultipleTTest)test).setPreferLinear(true);
            //test = new IndTestMultinomialLogisticRegressionWald(data, alpha, true);
        } else if (name.equals("tlog")){
            test = new edu.pitt.csb.mgm.IndTestMixedMultipleTTest(data, alpha);
            ((edu.pitt.csb.mgm.IndTestMixedMultipleTTest)test).setPreferLinear(false);
            //test = new IndTestMultinomialLogisticRegressionWald(data, alpha, false);
        } else {

            // This should allow the user to call any independence test found in tetrad.search or mgm
            Class cl = null;
            try {
                cl = Class.forName("edu.cmu.tetrad.search." + name);
            } catch (ClassNotFoundException e) {
                System.out.println("Not found: " + "edu.cmu.tetrad.search." + name);
            } catch (Exception e) {
                e.printStackTrace();
            }

            if(cl==null) {
                try {
                    cl = Class.forName("edu.pitt.csb.mgm." + name);
                } catch (ClassNotFoundException e) {
                    throw new IllegalArgumentException("-test argument not recognized");
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            try {
                Constructor con = cl.getConstructor(DataSet.class, double.class);
                test = (IndependenceTest) con.newInstance(data, alpha);
            } catch (NoSuchMethodException e) {
                System.err.println("Independence Test: " + name + " not found");
            } catch (Exception e) {
                System.err.println("Independence Test: " + name + " found but not constructed");
                e.printStackTrace();
            }
        }

        return test;
    }


    public static int[][] allEdgeStatsLatent(Graph pt, Graph pe)
    {
        HashMap<String,String> nd = new HashMap<String,String>();
        for(Node n: pe.getNodes())
        {
            if(n instanceof DiscreteVariable)
            {
                nd.put(n.getName(),"Disc");

            }
            else
            {
                nd.put(n.getName(),"Norm");
            }
        }
        return allEdgeStatsLatent(pt, pe, nd);
    }




    public static double [] structuralHamming(Graph pT, Graph pE, DataSet data)
    {
        HashMap<String,String> nd = new HashMap<String,String>();
        for(Node n: pE.getNodes())
        {
            if(data.getVariable(n.getName()) instanceof DiscreteVariable)
            {
                nd.put(n.getName(),"Disc");
            }
            else
            {
                nd.put(n.getName(),"Norm");
            }
        }
        double[] stats = new double[4];
        for (int i = 0; i < stats.length; i++) {
                stats[i] = 0;
        }

        //LOOP THROUGH TRUE EDGES
        for(Edge e: pT.getEdges())
        {
            int edgeType = -1;
            if (nd.get(e.getNode1().getName()).equals("Norm") && nd.get(e.getNode2().getName()).equals("Norm")) {
                edgeType = 0;
            } else if (nd.get(e.getNode1().getName()).equals("Disc") && nd.get(e.getNode2().getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }
            Edge temp = pE.getEdge(pE.getNode(e.getNode1().getName()),pE.getNode(e.getNode2().getName()));
            if(temp==null) //False Negative Undirected
            {
                stats[edgeType]+=2;
                stats[3]+=2;
            }
            else {
                Endpoint[] trueEndpoints = new Endpoint[]{e.getProximalEndpoint(e.getNode1()), e.getDistalEndpoint(e.getNode1())};
                Endpoint[] estEndpoints = new Endpoint[]{temp.getProximalEndpoint(pE.getNode(e.getNode1().getName())), temp.getDistalEndpoint(pE.getNode(e.getNode1().getName()))};
                if(trueEndpoints[0]!=estEndpoints[0])
                {
                    stats[edgeType]++;
                    stats[3]++;
                }
                if(trueEndpoints[1]!=estEndpoints[1])
                {
                    stats[edgeType]++;
                    stats[3]++;
                }

            }

        }
        for(Edge e:pE.getEdges())
        {
            int edgeType = -1;
            if (nd.get(e.getNode1().getName()).equals("Norm") && nd.get(e.getNode2().getName()).equals("Norm")) {
                edgeType = 0;
            } else if (nd.get(e.getNode1().getName()).equals("Disc") && nd.get(e.getNode2().getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }
            Edge temp = pT.getEdge(pT.getNode(e.getNode1().getName()),pT.getNode(e.getNode2().getName()));
            if(temp==null) //False Positive Undirected
            {
                stats[edgeType]+=2;
                stats[3]+=2;
            }
        }

        return stats;
    }


    //THIS IS THE CURRENT PAIRWISE TEST FOR THE MGM-FCI-MAX EXPERIMENT
    //INCLUDES STRUCTURAL HAMMING DISTANCE TODO ADD THESE SCORES
    public static double [][] newLatentScores(Graph pT, Graph pE, Graph trueDAG,DataSet data, boolean verbose)
    {
        HashMap<String,String> nd = new HashMap<String,String>();
        for(Node n: pE.getNodes())
        {
            if(data.getVariable(n.getName()) instanceof DiscreteVariable)
            {
                nd.put(n.getName(),"Disc");
            }
            else
            {
                nd.put(n.getName(),"Norm");
            }
        }
        double[][] stats = new double[4][14];
        //TPU, FPU, FNU, TNU, CorrOr, TotalOr, BiDirPrec,BiDirRec,
        int TPU = 0;
        int FPU = 1;
        int FNU = 2;
        int ETP = 3;
        int EFP = 4;
        int EFN = 5;
        int ETN = 6;
        int SHD = 7;

        for (int i = 0; i < stats.length; i++) {
            for (int j = 0; j < stats[0].length; j++) {
                stats[i][j] = 0;
            }
        }

        //LOOP THROUGH TRUE EDGES
        for(Edge e: pT.getEdges())
        {
            int edgeType = -1;
            if (nd.get(e.getNode1().getName()).equals("Norm") && nd.get(e.getNode2().getName()).equals("Norm")) {
                edgeType = 0;
            } else if (nd.get(e.getNode1().getName()).equals("Disc") && nd.get(e.getNode2().getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }
            Edge temp = pE.getEdge(pE.getNode(e.getNode1().getName()),pE.getNode(e.getNode2().getName()));
            if(temp==null) //False Negative Undirected
            {
                stats[edgeType][SHD]+=2;
                stats[3][SHD]+=2;
                stats[edgeType][FNU]++;
                stats[3][FNU]++;
            }
            else
            {
                stats[edgeType][TPU]++;
                stats[3][TPU]++;
                Endpoint [] trueEndpoints = new Endpoint[]{e.getProximalEndpoint(e.getNode1()),e.getDistalEndpoint(e.getNode1())};
                Endpoint[] estEndpoints = new Endpoint[]{temp.getProximalEndpoint(pE.getNode(e.getNode1().getName())),temp.getDistalEndpoint(pE.getNode(e.getNode1().getName()))};
                int [] types = getTypes(trueEndpoints,estEndpoints,trueDAG,new String []{e.getNode1().getName(),e.getNode2().getName()});
                if(types[0]<7) {
                    stats[edgeType][types[0]]++;
                    stats[3][types[0]]++;
                }
                else if(types[0]==7)
                {
                    stats[edgeType][EFN]+=0.5;
                    stats[3][EFN]+=0.5;
                }
                else
                {
                    stats[edgeType][EFP]+=0.5;
                    stats[3][EFP]+=0.5;
                }
                if(types[1]<7) {
                    stats[edgeType][types[0]]++;
                    stats[3][types[0]]++;
                }
                else if(types[1]==7)
                {
                    stats[edgeType][EFN]+=0.5;
                    stats[3][EFN]+=0.5;
                }
                else
                {
                    stats[edgeType][EFP]+=0.5;
                    stats[3][EFP]+=0.5;
                }

                if(types[0]!=ETP)
                {
                    stats[edgeType][SHD]++;
                    stats[3][SHD]++;
                }
                if(types[1]!=ETP)
                {
                    stats[edgeType][SHD]++;
                    stats[3][SHD]++;
                }
            }

        }
        for(Edge e:pE.getEdges())
        {
            int edgeType = -1;
            if (nd.get(e.getNode1().getName()).equals("Norm") && nd.get(e.getNode2().getName()).equals("Norm")) {
                edgeType = 0;
            } else if (nd.get(e.getNode1().getName()).equals("Disc") && nd.get(e.getNode2().getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }
            Edge temp = pT.getEdge(pT.getNode(e.getNode1().getName()),pT.getNode(e.getNode2().getName()));
            if(temp==null) //False Positive Undirected
            {
                stats[edgeType][FPU]++;
                stats[3][FPU]++;
                stats[edgeType][SHD]+=2;
                stats[3][SHD]+=2;
            }
        }

        return stats;
    }
    private static int [] getTypes(Endpoint [] trueEnds, Endpoint[] estEnds, Graph dag, String [] nodes)
    {
        int TPU = 0;
        int FPU = 1;
        int FNU = 2;
        int ETP = 3;
        int EFP = 4;
        int EFN = 5;
        int ETN = 6;
        int HFN = 7;
        int HFP = 8;
        int [] result = new int[trueEnds.length];
        for(int i = 0; i < trueEnds.length;i++)
        {
            if(trueEnds[i]==estEnds[i])
                result[i] = ETP;
            else
            {
                if(trueEnds[i]==Endpoint.CIRCLE)
                {
                    if(!dag.existsDirectedPathFromTo(dag.getNode(nodes[0]),dag.getNode(nodes[1])))
                        result[i]=ETP;
                    else if(estEnds[i]==Endpoint.ARROW)
                        result[i]=EFN;
                    else
                        result[i]=EFP;
                }
                else if(trueEnds[i]==Endpoint.ARROW)
                {
                    if(estEnds[i]==Endpoint.CIRCLE)
                        result[i] = HFP;
                    else
                        result[i] = EFP;
                }
                else
                {
                    if(estEnds[i]==Endpoint.CIRCLE)
                        result[i]=HFN;
                    else
                        result[i] = EFN;
                }
            }
        }
        return result;
    }
    public static double [][] allEdgeStatsLatentNew(Graph pT,Graph pE,Graph truth, List<Node>latents, DataSet data)
    {
        //TODO handle bidirected edges here,
        HashMap<String,String> nd = new HashMap<String,String>();
        for(Node n: pE.getNodes())
        {
            if(data.getVariable(n.getName()) instanceof DiscreteVariable)
            {
                nd.put(n.getName(),"Disc");

            }
            else
            {
                nd.put(n.getName(),"Norm");
            }
        }
        double[][] stats = new double[4][14];
        //TPU, FPU, FNU, TNU, CorrOr, TotalOr, BiDirPrec,BiDirRec,
        int TPU = 0;
        int FPU = 1;
        int FNU = 2;
        int ETP = 3;
        int EFP = 4;
        int EFN = 5;
        int ETN = 6;
        int ITP =  7;
        int IFP = 8;
        int IFN = 9;
        int ITN = 10;
        int LTP = 11;
        int LFP = 12;
        int LFN = 13;

        for (int i = 0; i < stats.length; i++) {
            for (int j = 0; j < stats[0].length; j++) {
                stats[i][j] = 0;
            }
        }
        //enforce patterns?
        //Graph pT = SearchGraphUtils.patternFromDag(tg);
        //Graph pE = SearchGraphUtils.patternFromDag(eg);

        //check that variable names are the same...

        Set<Edge> edgesT = pT.getEdges();
        Set<Edge> edgesE = pE.getEdges();

        //differences += Math.abs(e1.size() - e2.size());

        //for (int i = 0; i < e1.size(); i++) {
        int edgeType;

        boolean debug = false;

        //LOOPING THROUGH TRUE EDGES
        for (Edge eT : edgesT) {
            Node n1 = pE.getNode(eT.getNode1().getName());
            Node n2 = pE.getNode(eT.getNode2().getName());
            if (nd.get(n1.getName()).equals("Norm") && nd.get(n2.getName()).equals("Norm")) {
                edgeType = 0;
            } else if (nd.get(n1.getName()).equals("Disc") && nd.get(n2.getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }

            Edge eE = pE.getEdge(n1, n2);

            if(debug)
                System.out.println("True Edge: " + eT);

            if (eE == null) {//Estimated = nothing, True = edge
                stats[edgeType][FNU]++;
                stats[3][FNU]++;

                if(debug)
                    System.out.println("No Estimated Edge");
            }
            else
            {

                //Estimated = A <-> B
                if(eE.getEndpoint1()==Endpoint.ARROW && eE.getEndpoint2()==Endpoint.ARROW) {

                    //True = A <-> B
                    if(eT.getEndpoint1()==Endpoint.ARROW && eT.getEndpoint2()==Endpoint.ARROW)
                    {
                        stats[edgeType][LTP]++;
                        stats[3][LTP]++;
                    }
                    else
                    {
                        boolean found = false;
                       A:for(Node x:latents)
                       {
                           if(truth.isParentOf(truth.getNode(x.getName()),truth.getNode(eE.getNode1().getName())) && truth.isParentOf(truth.getNode(x.getName()),truth.getNode(eE.getNode2().getName())))
                           {
                               found = true;

                               stats[edgeType][LTP]++;
                               stats[3][LTP]++;
                               break A;
                           }
                       }
                       if(!found)
                       {
                           stats[edgeType][LFP]++;
                           stats[3][LFP]++;
                       }

                    }
                }
                else if(eT.getEndpoint1()==Endpoint.ARROW && eT.getEndpoint1()==Endpoint.ARROW)
                {
                    if(eE.getEndpoint1()==Endpoint.TAIL || eE.getEndpoint2()==Endpoint.TAIL)
                    {
                        stats[edgeType][LFN]++;
                        stats[3][LFN]++;
                    }
                    else if(eE.getEndpoint1()==Endpoint.CIRCLE || eE.getEndpoint2()==Endpoint.CIRCLE)
                    {
                        stats[edgeType][LFN]+=0.5;
                        stats[3][LFN]+=0.5;
                    }
                    else
                    {
                        System.out.println("Shouldn't reach here");
                        return null;
                    }
                }
                stats[edgeType][TPU]++;
                stats[3][TPU]++;
                Edge temp;
                if(!eE.getNode1().getName().equals(eT.getNode1().getName())) //Don't reverse edge
                    temp = new Edge(eE.getNode2(),eE.getNode1(),eE.getEndpoint2(),eE.getEndpoint1());
                else
                    temp = eE;
                if(debug)
                System.out.println("Estimated Edge: " + temp);
                    if(temp.getEndpoint1()==eT.getEndpoint1())//Correct
                    {
                        stats[edgeType][ETP]++;
                        stats[3][ETP]++;
                        if(debug)
                            System.out.println("Correct Endpoint 1");
                    }
                    else if(temp.getEndpoint1()==Endpoint.CIRCLE) //A o-* B
                    {
                        if(eT.getEndpoint1()==Endpoint.ARROW)
                        {
                            stats[edgeType][EFP]+=.5;
                            stats[3][EFP]+=.5;
                            if(debug)
                                System.out.println("Half False Positive");
                        }
                        else
                        {
                         stats[edgeType][EFN]+=0.5;
                         stats[3][EFN]+=0.5;
                         if(debug)
                             System.out.println("Half False Negative");
                        }
                    }
                    else if(temp.getEndpoint1()==Endpoint.ARROW) // A <-* B
                    {
                        if(eT.getEndpoint1()==Endpoint.CIRCLE)
                        {
                            //TODO Make sure this is right
                            if(truth.isAncestorOf(truth.getNode(temp.getNode1().getName()),truth.getNode(temp.getNode2().getName())))
                            {
                                stats[edgeType][EFN]++;
                                stats[3][EFN]++;
                                if(debug)
                                    System.out.println(temp.getNode1() + " is an ancestor of " + temp.getNode2() + ", so FN");
                            }
                            else
                            {
                                stats[edgeType][ETP]++;
                                stats[3][ETP]++;
                                if(debug)
                                    System.out.println("Correct Endpoint 1");
                            }

                        }
                        else
                        {
                            stats[edgeType][EFN]++;
                            stats[3][EFN]++;
                            if(debug)
                                System.out.println("False Negative Endpoint 1");
                        }
                    }
                    else //A -* B (A-> B)
                    {
                        if(eT.getEndpoint1()==Endpoint.CIRCLE)
                        {
                            if(truth.isAncestorOf(truth.getNode(temp.getNode1().getName()),truth.getNode(temp.getNode2().getName())))
                            {
                                stats[edgeType][ETP]++;
                                stats[3][ETP]++;
                                if(debug)
                                    System.out.println("Correct Endpoint 1");
                            }
                            else
                            {
                                stats[edgeType][EFP]++;
                                stats[3][EFP]++;
                                if(debug)
                                    System.out.println("False Positive Endpoint 1");
                            }
                        }
                        else if(eT.getEndpoint1()==Endpoint.ARROW)
                        {
                            stats[edgeType][EFP]++;
                            stats[3][EFP]++;
                            if(debug)
                                System.out.println("False Positive Endpoint 1");
                        }
                    }
                    if(temp.getEndpoint2()==eT.getEndpoint2())//Correct
                    {
                        stats[edgeType][ETP]++;
                        stats[3][ETP]++;
                        if(debug)
                            System.out.println("Correct Endpoint 2");
                    }
                    else if(temp.getEndpoint2()==Endpoint.CIRCLE) //A o-* B
                    {
                        if(eT.getEndpoint2()==Endpoint.ARROW)
                        {
                            stats[edgeType][EFP]+=.5;
                            stats[3][EFP]+=.5;
                            if(debug)
                                System.out.println("Half False Positive");
                        }
                        else
                        {
                            stats[edgeType][EFN]+=0.5;
                            stats[3][EFN]+=0.5;
                            if(debug)
                                System.out.println("Half False Negative Endpoint 2");
                        }
                    }
                    else if(temp.getEndpoint2()==Endpoint.ARROW) // A <-* B
                    {
                        if(eT.getEndpoint2()==Endpoint.CIRCLE)
                        {
                            //TODO Make sure this is right
                            if(truth.isAncestorOf(truth.getNode(temp.getNode2().getName()),truth.getNode(temp.getNode1().getName())))
                            {
                                stats[edgeType][EFN]++;
                                stats[3][EFN]++;
                                if(debug)
                                    System.out.println(temp.getNode2() + " is an ancestor of " + temp.getNode1() + ", so False Negative Endpoint 2");
                            }
                            else
                            {
                                stats[edgeType][ETP]++;
                                stats[3][ETP]++;
                                if(debug)
                                    System.out.println("Correct Endpoint 2");
                            }

                        }
                        else
                        {
                            stats[edgeType][EFN]++;
                            stats[3][EFN]++;
                            if(debug)
                                System.out.println("False Negative Endpoint 2");
                        }
                    }
                    else //A -* B (A-> B)
                    {
                        if(eT.getEndpoint2()==Endpoint.CIRCLE)
                        {
                            if(truth.isAncestorOf(truth.getNode(temp.getNode2().getName()),truth.getNode(temp.getNode1().getName())))
                            {
                                stats[edgeType][ETP]++;
                                stats[3][ETP]++;
                                if(debug)
                                    System.out.println("Correct Endpoint 2");
                            }
                            else
                            {
                                stats[edgeType][EFP]++;
                                stats[3][EFP]++;
                                if(debug)
                                    System.out.println(temp.getNode2() + " is not an ancestor of " + temp.getNode1() + ", so False Positive Endpoint 1");
                            }
                        }
                        else if(eT.getEndpoint2()==Endpoint.ARROW)
                        {
                            stats[edgeType][EFP]++;
                            stats[3][EFP]++;
                            if(debug)
                                System.out.println("False Positive Endpoint 2");
                        }
                    }

                }
            }
            for(Edge eE: edgesE)
            {
                Node n1 = pT.getNode(eE.getNode1().getName());
                Node n2 = pT.getNode(eE.getNode2().getName());
                if (nd.get(n1.getName()).equals("Norm") && nd.get(n2.getName()).equals("Norm")) {
                    edgeType = 0;
                } else if (nd.get(n1.getName()).equals("Disc") && nd.get(n2.getName()).equals("Disc")) {
                    edgeType = 2;
                } else {
                    edgeType = 1;
                }

                Edge eT = pT.getEdge(n1, n2);
                if (eT == null) {//Estimated = edge, True = nothing
                    stats[edgeType][FPU]++;
                    stats[3][FPU]++;

                }
            }
        return stats;
    }
    public static double [] [] allEdgeStatsLatentNew(Graph pT,Graph pE, Graph truth,List<Node> latents) {
        //TODO handle bidirected edges here,
        HashMap<String,String> nd = new HashMap<String,String>();
        for(Node n: pE.getNodes())
        {
            if(n instanceof DiscreteVariable)
            {
                nd.put(n.getName(),"Disc");

            }
            else
            {
                nd.put(n.getName(),"Norm");
            }
        }
        double[][] stats = new double[4][14];
        //TPU, FPU, FNU, TNU, CorrOr, TotalOr, BiDirPrec,BiDirRec,
        int TPU = 0;
        int FPU = 1;
        int FNU = 2;
        int ETP = 3;
        int EFP = 4;
        int EFN = 5;
        int ETN = 6;
        int ITP =  7;
        int IFP = 8;
        int IFN = 9;
        int ITN = 10;
        int LTP = 11;
        int LFP = 12;
        int LFN = 13;

        for (int i = 0; i < stats.length; i++) {
            for (int j = 0; j < stats[0].length; j++) {
                stats[i][j] = 0;
            }
        }
        //enforce patterns?
        //Graph pT = SearchGraphUtils.patternFromDag(tg);
        //Graph pE = SearchGraphUtils.patternFromDag(eg);

        //check that variable names are the same...

        Set<Edge> edgesT = pT.getEdges();
        Set<Edge> edgesE = pE.getEdges();

        //differences += Math.abs(e1.size() - e2.size());

        //for (int i = 0; i < e1.size(); i++) {
        int edgeType;
        for (Edge eT : edgesT) {
            Node n1 = pE.getNode(eT.getNode1().getName());
            Node n2 = pE.getNode(eT.getNode2().getName());

            if (nd.get(n1.getName()).equals("Norm") && nd.get(n2.getName()).equals("Norm")) {
                edgeType = 0;
            } else if (nd.get(n1.getName()).equals("Disc") && nd.get(n2.getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }

            Edge eE = pE.getEdge(n1, n2);
            if (eE == null) {//Estimated = nothing, True = edge
                stats[edgeType][FNU]++;
                stats[3][FNU]++;

            } else {
                String A = eE.getNode1().getName();
                String B = eE.getNode2().getName();
                boolean reverse = false;
                if(eE.getProximalEndpoint(pE.getNode(A))==Endpoint.ARROW && eE.getDistalEndpoint(pE.getNode(A))!=Endpoint.ARROW) {
                    String temp = B;
                    B = A;
                    A = temp;
                    reverse = true;
                    System.out.println("Reversed Edge");
                }
                System.out.println("Estimated Edge: " + eE);
                System.out.println("True Edge: " + eT);
                System.out.println("Endpoint1-Estimated: " + eE.getProximalEndpoint(pE.getNode(A)) + ", Endpoint1-True: " + eT.getProximalEndpoint(pT.getNode(A)));
                System.out.println("Node 1: " + eE.getDistalNode(pE.getNode(B)) + " " + eT.getDistalNode(pT.getNode(B)));
                System.out.println("Node 2: " + eE.getDistalNode(pE.getNode(A)) + " " + eT.getDistalNode(pT.getNode(A)));
                System.out.println("Endpoint2-Estimated: " + eE.getDistalEndpoint(pE.getNode(A)) + ", Endpoint2-True: " + eT.getDistalEndpoint(pT.getNode(A)));
                stats[edgeType][TPU]++;
                stats[3][TPU]++;
                if(eE.getEndpoint1()==Endpoint.CIRCLE && eE.getEndpoint2() == Endpoint.CIRCLE) {
                    if (eT.getEndpoint1() == Endpoint.CIRCLE && eT.getEndpoint2() == Endpoint.CIRCLE) {
                        stats[edgeType][ITP]++;
                    } else if (eT.getEndpoint1() == Endpoint.CIRCLE || eT.getEndpoint2() == Endpoint.CIRCLE) {
                        boolean found = false;
                        for (Node x : latents) {
                            if (truth.isParentOf(x, truth.getNode(eE.getNode1().getName())) && truth.isParentOf(x, truth.getNode(eE.getNode2().getName())) && !found) {
                                found = true;
                                stats[edgeType][ITP] += .33;
                                stats[edgeType][IFP] += .66;
                                stats[edgeType][IFN] += .66;
                                stats[edgeType][ITN] += 1.33;

                            }
                        }
                        if (found) {
                            edgesE.remove(eE);
                            continue;
                        }
                        stats[edgeType][ITP] += .33;
                        stats[edgeType][ITN] += .66;
                        stats[edgeType][IFP] += .33;
                        stats[edgeType][IFN] += .66;
                    } else if (eT.getEndpoint1() == Endpoint.TAIL || eT.getEndpoint2() == Endpoint.TAIL) {
                        stats[edgeType][ITP] += .33;
                        stats[edgeType][IFP] += .66;
                        stats[edgeType][IFN] += .66;
                        stats[edgeType][ITN] += .66;
                    } else {
                        stats[edgeType][ITP] += .33;
                        stats[edgeType][IFP] += .66;
                        stats[edgeType][ITN] += 1.33;
                        stats[edgeType][IFN] += .66;
                    }
                }
                else if (eE.getProximalEndpoint(pE.getNode(A))==Endpoint.CIRCLE && eE.getDistalEndpoint(pE.getNode(A))==Endpoint.ARROW)
                {
                    if(eT.getProximalEndpoint(pT.getNode(A))==Endpoint.CIRCLE && eT.getDistalEndpoint(pT.getNode(A))==Endpoint.CIRCLE)
                    {
                        if(!truth.isAncestorOf(truth.getNode(B),truth.getNode(A)))
                        {
                            stats[edgeType][ETN]++;
                            stats[edgeType][ITP]++;
                        }
                        else
                        {
                            stats[edgeType][EFN]++;
                            stats[edgeType][IFP]++;
                        }
                    }
                    else if(eT.getProximalEndpoint(pT.getNode(A))==Endpoint.CIRCLE)
                    {
                        stats[edgeType][ETN]++;
                        stats[edgeType][ITP]++;
                    }
                    else if(eT.getDistalEndpoint(pT.getNode(B))==Endpoint.CIRCLE)
                    {
                        if(!truth.isAncestorOf(truth.getNode(B),truth.getNode(A)))
                        {
                            stats[edgeType][ETN]++;
                            stats[edgeType][ITN] += 1.5;
                            stats[edgeType][IFP]+=0.5;
                            stats[edgeType][ITP]+=0.5;
                            stats[edgeType][IFN]+= 0.5;

                        }
                        else
                        {
                            stats[edgeType][EFN]++;
                            stats[edgeType][IFP]++;
                            stats[edgeType][IFN]++;
                            stats[edgeType][ITN]+=0.5;
                        }
                    }
                    else if(eT.getProximalEndpoint(pT.getNode(A))==Endpoint.ARROW && eT.getDistalEndpoint(pT.getNode(A))==Endpoint.TAIL)
                    {
                        stats[edgeType][IFP]++;
                        stats[edgeType][IFN]++;
                        stats[edgeType][ITN]+=0.5;
                        stats[edgeType][EFN]++;
                        stats[edgeType][EFP]+=0.5;
                        stats[edgeType][ETN]+=0.5;
                    }
                    else if (eT.getProximalEndpoint(pT.getNode(A))==Endpoint.ARROW)
                    {
                        stats[edgeType][ETN]+=1.5;
                        stats[edgeType][EFP]+=0.5;
                        stats[edgeType][ITP]+=0.5;
                        stats[edgeType][ITN]+=1.5;
                        stats[edgeType][IFP]+=0.5;
                        stats[edgeType][IFN]+=0.5;
                    }
                    else
                    {
                        stats[edgeType][ITP]+=0.5;
                        stats[edgeType][ITN]++;
                        stats[edgeType][IFN]+=0.5;
                        stats[edgeType][IFP]+=0.5;
                    }

                }
                else if(eE.getProximalEndpoint(pE.getNode(A))==Endpoint.TAIL && eE.getDistalEndpoint(pE.getNode(A))==Endpoint.ARROW) //Estinmated A-->B
                {
                    if(eT.getProximalEndpoint(pT.getNode(A))==Endpoint.CIRCLE || eT.getProximalEndpoint(pT.getNode(A))==Endpoint.CIRCLE)
                    {
                        if(truth.isAncestorOf(truth.getNode(A),truth.getNode(B)))
                        {
                            stats[edgeType][ITP]++;
                            stats[edgeType][ETP]++;
                            stats[edgeType][ITN]++;
                            stats[edgeType][ETN]++;
                        }
                        else if(truth.isAncestorOf(truth.getNode(B),truth.getNode(A)))
                        {
                            stats[edgeType][IFP]++;
                            stats[edgeType][EFP]++;
                            stats[edgeType][IFN]++;
                            stats[edgeType][EFN]++;
                        }
                        else
                        {
                            stats[edgeType][ETN]++;
                            stats[edgeType][EFP]++;
                            stats[edgeType][ITN]++;
                            stats[edgeType][IFP]++;
                            stats[edgeType][IFN]++;
                        }
                    }
                    else if(eT.getProximalEndpoint(pT.getNode(A))==Endpoint.TAIL)
                    {
                        stats[edgeType][ITP]++;
                        stats[edgeType][ETP]++;
                        stats[edgeType][ITN]++;
                        stats[edgeType][ETN]++;
                    }
                    else if(eT.getDistalEndpoint(pT.getNode(A))==Endpoint.TAIL)
                    {
                        stats[edgeType][EFP]++;
                        stats[edgeType][EFN]++;
                        stats[edgeType][IFP]++;
                        stats[edgeType][IFN]++;
                    }
                    else
                    {
                        stats[edgeType][ETN]++;
                        stats[edgeType][EFP]++;
                        stats[edgeType][ITN]++;
                        stats[edgeType][IFP]++;
                        stats[edgeType][IFN]++;
                    }
                }
                else //Estimated Latent
                {
                    boolean found = false;
                    for (Node x : latents) {
                        if (truth.isParentOf(x, truth.getNode(eE.getNode1().getName())) && truth.isParentOf(x, truth.getNode(eE.getNode2().getName())) && !found) {
                            found = true;
                            stats[edgeType][ITP]++;
                            stats[edgeType][ITN]+=2;
                            stats[edgeType][ETN]+=2;

                        }
                    }
                    if (found)
                    {
                        edgesE.remove(eE);
                        continue;
                    }

                    stats[edgeType][IFN] ++;
                    stats[edgeType][ITN]++;
                    stats[edgeType][IFP]++;
                    stats[edgeType][ETN] ++;
                    stats[edgeType][EFN] ++;
                }

            }
            if (eE != null)
                edgesE.remove(eE);
        }

        for (Edge eE : edgesE) {
            Node n1 = pT.getNode(eE.getNode1().getName());
            Node n2 = pT.getNode(eE.getNode2().getName());

            if (nd.get(n1.getName()).equals("Norm") && nd.get(n2.getName()).equals("Norm")) {
                edgeType = 0;
            } else if (nd.get(n1.getName()).equals("Disc") && nd.get(n2.getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }


            Edge eT = pT.getEdge(n1, n2);
            if (eT == null) {//Estimated = edge, True = nothing
                stats[edgeType][FPU]++;
                stats[3][FPU]++;
            } else {
                System.out.println("We do get here");
            }

        }
        return stats;
    }
    public static double[][] allEdgeStatsLatentPairwise(Graph pT, Graph pE, Graph truth, List<Node> latents, DataSet data)
    {
        HashMap<String,String> nd = new HashMap<String,String>();
        for(Node n: pE.getNodes())
        {
            if(data.getVariable(n.getName()) instanceof DiscreteVariable)
            {
                nd.put(n.getName(),"Disc");
            }
            else
            {
                nd.put(n.getName(),"Norm");
            }
        }
        double[][] stats = new double[4][14];
        //TPU, FPU, FNU, TNU, CorrOr, TotalOr, BiDirPrec,BiDirRec,
        int TPU = 0;
        int FPU = 1;
        int FNU = 2;
        int ETP = 3;
        int EFP = 4;
        int EFN = 5;
        int ETN = 6;
        int ITP = 7;
        int IFP = 8;
        int IFN = 9;
        int ITN = 10;
        int LTP = 11;
        int LFP = 12;
        int LFN = 13;

        for (int i = 0; i < stats.length; i++) {
            for (int j = 0; j < stats[0].length; j++) {
                stats[i][j] = 0;
            }
        }
        List<Node> allNodes = pE.getNodes();
        for(int i= 0; i < allNodes.size();i++)
        {
            for(int j = i+1;j < allNodes.size();j++)
            {
                Node uno = allNodes.get(i);
                Node dos = allNodes.get(j);
                System.out.println("Node 1: " + uno);
                System.out.println("Node 2: " + dos);
                int edgeType = -1;
                if (nd.get(uno.getName()).equals("Norm") && nd.get(dos.getName()).equals("Norm")) {
                    edgeType = 0;
                } else if (nd.get(uno.getName()).equals("Disc") && nd.get(dos.getName()).equals("Disc")) {
                    edgeType = 2;
                } else {
                    edgeType = 1;
                }


                if(pE.existsDirectedPathFromTo(uno,dos)) //Uno is predicted to be an ancestor of dos
                {
                    System.out.println("Esimated: Ancestor");
                    if(pT.existsDirectedPathFromTo(pT.getNode(uno.getName()),pT.getNode(dos.getName())))
                    {
                        System.out.println("True: Ancestor");
                        stats[edgeType][ETP]++;
                        stats[3][ETP]++;
                    }
                    else if(pT.possibleAncestor(pT.getNode(uno.getName()),pT.getNode(dos.getName())))
                    {
                        if(truth.isAncestorOf(truth.getNode(uno.getName()),truth.getNode(dos.getName())))
                        {
                            System.out.println("True: Ancestor (Look up)");
                            stats[edgeType][ETP]++;
                            stats[3][ETP]++;
                        }
                        else
                        {
                            System.out.println("True: Non Ancestor (Look up)");
                            stats[edgeType][EFP]++;
                            stats[3][EFP]++;
                        }
                    }
                    else
                    {
                        System.out.println("True: Non Ancestor");
                        stats[edgeType][EFP]++;
                        stats[3][EFP]++;
                    }
                }
                else if(pE.possibleAncestor(uno,dos)) //IDK
                {
                    System.out.println("Estimated: Not Sure");
                    if(pT.existsDirectedPathFromTo(pT.getNode(uno.getName()),pT.getNode(dos.getName())))
                    {
                        System.out.println("True: Ancestor");
                        stats[edgeType][EFN]+= 0.5;
                        stats[3][EFN]+= 0.5;
                    }
                    else if(pT.possibleAncestor(pT.getNode(uno.getName()),pT.getNode(dos.getName())))
                    {
                        System.out.println("True: Don't know");
                        //NOTHING
                    }
                    else
                    {
                        System.out.println("True: Non Ancestor");
                        //NOTHING
                    }
                }
                else //Predicted NO
                {
                    System.out.println("Estimated Non Ancestry");
                    if(pT.existsDirectedPathFromTo(pT.getNode(uno.getName()),pT.getNode(dos.getName())))
                    {
                        System.out.println("True: Ancestor");
                        stats[edgeType][EFN]+= 1;
                        stats[3][EFN]+= 1;
                    }
                    else if(pT.possibleAncestor(pT.getNode(uno.getName()),pT.getNode(dos.getName())))
                    {
                        System.out.println("True: Don't know");
                        if(truth.isAncestorOf(truth.getNode(uno.getName()),truth.getNode(dos.getName())))
                        {
                            System.out.println("True: Ancestor (Look up)");
                            stats[edgeType][EFN]++;
                            stats[3][EFN]++;
                        }
                        else
                        {
                            System.out.println("True: Non Ancestor (Look up)");
                            stats[edgeType][ETN]++;
                            stats[3][ETN]++;
                        }
                    }
                    else
                    {
                        System.out.println("True: Non Ancestor");
                        stats[edgeType][ETN]++;
                        stats[3][ETN]++;
                    }
                }


                if(pE.existsDirectedPathFromTo(dos,uno)) //Uno is predicted to be an ancestor of dos
                {
                    if(pT.existsDirectedPathFromTo(pT.getNode(dos.getName()),pT.getNode(uno.getName())))
                    {
                        stats[edgeType][ETP]++;
                        stats[3][ETP]++;
                    }
                    else if(pT.possibleAncestor(pT.getNode(dos.getName()),pT.getNode(uno.getName())))
                    {
                        if(truth.isAncestorOf(truth.getNode(dos.getName()),truth.getNode(uno.getName())))
                        {
                            stats[edgeType][ETP]++;
                            stats[3][ETP]++;
                        }
                        else
                        {
                            stats[edgeType][EFP]++;
                            stats[3][EFP]++;
                        }
                    }
                    else
                    {
                        stats[edgeType][EFP]++;
                        stats[3][EFP]++;
                    }
                }
                else if(pE.possibleAncestor(dos,uno)) //IDK
                {
                    if(pT.existsDirectedPathFromTo(pT.getNode(dos.getName()),pT.getNode(uno.getName())))
                    {
                        stats[edgeType][EFN]+= 0.5;
                        stats[3][EFN]+= 0.5;
                    }
                    else if(pT.possibleAncestor(pT.getNode(dos.getName()),pT.getNode(uno.getName())))
                    {
                        //NOTHING
                    }
                    else
                    {
                        //NOTHING
                    }
                }
                else //Predicted NO
                {
                    if(pT.existsDirectedPathFromTo(pT.getNode(dos.getName()),pT.getNode(uno.getName())))
                    {
                        stats[edgeType][EFN]+= 1;
                        stats[3][EFN]+= 1;
                    }
                    else if(pT.possibleAncestor(pT.getNode(dos.getName()),pT.getNode(uno.getName())))
                    {
                        if(truth.isAncestorOf(truth.getNode(dos.getName()),truth.getNode(uno.getName())))
                        {
                            stats[edgeType][EFN]++;
                            stats[3][EFN]++;
                        }
                        else
                        {
                            stats[edgeType][ETN]++;
                            stats[3][ETN]++;
                        }
                    }
                    else
                    {
                        stats[edgeType][ETN]++;
                        stats[3][ETN]++;
                    }
                }
            }
        }
        Set<Edge> pE_edge = pE.getEdges();
        Set<Edge> pT_edge = pT.getEdges();
        for(Edge e: pE_edge)
        {
            Node x = e.getNode1();
            Node y = e.getNode2();
            int edgeType = -1;
            if (nd.get(x.getName()).equals("Norm") && nd.get(y.getName()).equals("Norm")) {
                edgeType = 0;
            } else if (nd.get(x.getName()).equals("Disc") && nd.get(y.getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }
            if(pT.getEdge(pT.getNode(x.getName()),pT.getNode(y.getName()))==null)
            {
                stats[edgeType][FPU]++;
                stats[3][FPU]++;
            }
            else
            {
                stats[edgeType][TPU]++;
                stats[3][TPU]++;
                if(!pT_edge.remove(pT.getEdge(pT.getNode(x.getName()),pT.getNode(y.getName()))))
                {
                    System.out.println("Bad Edge");
                    System.exit(0);
                }
            }

        }
        for(Edge e:pT_edge)
        {
            Node x = e.getNode1();
            Node y = e.getNode2();
            int edgeType = -1;
            if (nd.get(x.getName()).equals("Norm") && nd.get(y.getName()).equals("Norm")) {
                edgeType = 0;
            } else if (nd.get(x.getName()).equals("Disc") && nd.get(y.getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }
            stats[edgeType][FNU]++;
            stats[3][FNU]++;
        }
        return stats;
    }

    public static DoubleMatrix2D pagToLatentMatrix(Graph g)
    {
        int n = g.getNumNodes();
        DoubleMatrix2D matrix = DoubleFactory2D.dense.make(n, n, 0.0);

        // map node names in order of appearance
        HashMap<Node, Integer> map = new HashMap<Node, Integer>();
        int i = 0;
        for (Node node : g.getNodes()) {
            map.put(node, i);
            i++;
        }

        // mark edges
        for (Edge edge : g.getEdges()) {
            // if directed find which is parent/child
            Node node1 = edge.getNode1();
            Node node2 = edge.getNode2();

            if(edge.getEndpoint1()==Endpoint.ARROW && edge.getEndpoint2()==Endpoint.ARROW)
            {
                matrix.set(map.get(node1), map.get(node2), 1.0);
                matrix.set(map.get(node2), map.get(node1), 1.0);
            }
            else if(edge.getEndpoint1()==Endpoint.ARROW && edge.getEndpoint2()==Endpoint.CIRCLE)
            {
                matrix.set(map.get(node1), map.get(node2), 0.5);
                matrix.set(map.get(node2), map.get(node1), 0.5);
            }
            else if(edge.getEndpoint2()==Endpoint.ARROW && edge.getEndpoint1()==Endpoint.CIRCLE)
            {
                matrix.set(map.get(node1), map.get(node2), 0.5);
                matrix.set(map.get(node2), map.get(node1), 0.5);
            }

        }

        return matrix;
    }
    public static int [][] allEdgeStatsLatent(Graph pT, Graph pE, HashMap<String,String>nodeDists) {

        int[][] stats = new int[3][9];
        //TPU, FPU, FNU, TNU, CorrOr, TotalOr, BiDirPrec,BiDirRec,
        int TPU = 0;
        int FPU = 1;
        int FNU = 2;
        int TNU = 3;
        int OC = 4;
        int OW = 5;
        int BDTP = 6;
        int BDFP = 7;
        int BDFN = 8;
        for (int i = 0; i < stats.length; i++) {
            for (int j = 0; j < stats[0].length; j++) {
                stats[i][j] = 0;
            }
        }
        //enforce patterns?
        //Graph pT = SearchGraphUtils.patternFromDag(tg);
        //Graph pE = SearchGraphUtils.patternFromDag(eg);

        //check that variable names are the same...

        Set<Edge> edgesT = pT.getEdges();
        Set<Edge> edgesE = pE.getEdges();

        //differences += Math.abs(e1.size() - e2.size());

        //for (int i = 0; i < e1.size(); i++) {
        int edgeType;
        for (Edge eT : edgesT) {
            Node n1 = pE.getNode(eT.getNode1().getName());
            Node n2 = pE.getNode(eT.getNode2().getName());
            System.out.println(eT.getNode1());
            System.out.println(eT.getNode2());
            System.out.println(n1);
            System.out.println(n2);
            if (nodeDists.get(n1.getName()).equals("Norm") && nodeDists.get(n2.getName()).equals("Norm")) {
                edgeType = 0;
            } else if (nodeDists.get(n1.getName()).equals("Disc") && nodeDists.get(n2.getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }

            Edge eE = pE.getEdge(n1, n2);
            if (eE == null) {//Estimated = nothing, True = edge
                stats[edgeType][FNU]++;
                if(eT.getEndpoint1()==Endpoint.ARROW && eT.getEndpoint2()==Endpoint.ARROW)
                    stats[edgeType][BDFN]++;
            } else
            {
                stats[edgeType][TPU]++;
                if(eT.getEndpoint1()==Endpoint.ARROW && eT.getEndpoint2()==Endpoint.ARROW)
                {
                    if(eE.getEndpoint1()==Endpoint.ARROW && eE.getEndpoint2()==Endpoint.ARROW)
                        stats[edgeType][BDTP]++;
                    else
                        stats[edgeType][BDFN]++;
                }
                else
                {
                    if(eE.getEndpoint1()==Endpoint.ARROW && eE.getEndpoint2()==Endpoint.ARROW)
                        stats[edgeType][BDFP]++;
                }
                if(eT.getEndpoint1()==eE.getEndpoint1())
                {
                    stats[edgeType][OC]++;

                }
                else //estimated = circle, true = arrow
                {
                    stats[edgeType][OW]++;
                }
                if(eT.getEndpoint2()==eE.getEndpoint2()) {
                    stats[edgeType][OC]++;
                }
                else //estimated = circle, true = arrow
                {
                    stats[edgeType][OW]++;
                }
            }
            if(eE!=null)
                edgesE.remove(eE);
        }

        for(Edge eE: edgesE) {
            Node n1 = pT.getNode(eE.getNode1().getName());
            Node n2 = pT.getNode(eE.getNode2().getName());

            if (nodeDists.get(n1.getName()).equals("Norm") && nodeDists.get(n2.getName()).equals("Norm")) {
                edgeType = 0;
            } else if (nodeDists.get(n1.getName()).equals("Disc") && nodeDists.get(n2.getName()).equals("Disc")) {
                edgeType = 2;
            } else {
                edgeType = 1;
            }


            Edge eT = pT.getEdge(n1, n2);
            if (eT == null) {//Estimated = edge, True = nothing
                stats[edgeType][FPU]++;
                if(eE.getEndpoint1()==Endpoint.ARROW && eE.getEndpoint2()==Endpoint.ARROW)
                    stats[edgeType][BDFP]++;
            }
            else {
                stats[edgeType][TPU]++;
                if (eT.getEndpoint1() == eE.getEndpoint1() && eT.getEndpoint2() == eE.getEndpoint2()) {
                    stats[edgeType][OC]++;
                    if(eT.getEndpoint1()==Endpoint.ARROW && eT.getEndpoint2()==Endpoint.ARROW)
                    {
                        stats[edgeType][BDTP]++;
                    }

                }
                else
                {
                    stats[edgeType][OW]++;
                    if(eT.getEndpoint1()==Endpoint.ARROW && eT.getEndpoint2()==Endpoint.ARROW)
                    {
                        stats[edgeType][BDFN]++;
                    }
                    else if (eE.getEndpoint1()==Endpoint.ARROW && eE.getEndpoint2()==Endpoint.ARROW)
                    {
                        stats[edgeType][BDFP]++;
                    }
                }
            }

        }
        return stats;
       /* int numDiscrete = 0;
        int numContinuous = 0;
        for(Node n: pE.getNodes())
        {
            if(n instanceof DiscreteVariable)
            {
                numDiscrete++;
            }
            else
                numContinuous++;
        }

        //fill in True negatives with all edges not labeled in true graph
        int sum = 0;
        int totalEdges = (numDiscrete * (numDiscrete-1))/2;
            for(int j = 0; j <4; j++)
            sum = sum + stats[2][j];
        stats[2][TNU] = totalEdges - sum;
        sum = 0;
        totalEdges = (numContinuous * numDiscrete)/2;
        for(int j = 0; j < 4; j++)
            sum = sum + stats[1][j];
        stats[1][TNU] = totalEdges - sum;
        totalEdges = (numContinuous * (numContinuous-1))/2;
        sum = 0;
        for(int j = 0; j <4; j ++)
            sum = sum + stats[0][j];
        stats[0][TNU] = totalEdges - sum;
    String x = "Edge_Type\tTPU\tFPU\tFNU\tTNU\tTPD\tFPD\tFND\tTND\n";
        x = x + "Continuous-Continuous\t";
        for(int i = 0; i < 8; i++)
            x = x + stats[0][i] + "\t";
        x = x + "\n";
        x = x + "Continuous-Discrete\t";
        for(int i = 0; i < 8; i ++)
            x = x + stats[1][i] + "\t";
        x = x + "\n";
        x = x + "Discrete-Discrete\t";
        for(int i = 0; i < 8; i ++)
            x = x + stats[2][i] + "\t";
        return x;*/
    }
    // break out stats by node distributions, here only "Norm" and "Disc"
    // so three types of possible edges, cc, cd, dd, output is edge type by stat type
    // counts bidirected


    //main for testing
    public static void main(String[] args){
        //Graph g = GraphConverter.convert("X1-->X2,X2-->X3,X3-->X4");
        Graph g = GraphConverter.convert("X1-->X2,X2-->X3,X3-->X4, X5-->X4");
        //simple graph pm im gen example

        HashMap<String, Integer> nd = new HashMap<>();
        nd.put("X1", 0);
        nd.put("X2", 0);
        nd.put("X3", 4);
        nd.put("X4", 4);
        nd.put("X5", 0);

        g = makeMixedGraph(g, nd);

        /*Graph g = new EdgeListGraph();
        g.addNode(new ContinuousVariable("X1"));
        g.addNode(new ContinuousVariable("X2"));
        g.addNode(new DiscreteVariable("X3", 4));
        g.addNode(new DiscreteVariable("X4", 4));
        g.addNode(new ContinuousVariable("X5"));
        g.addDirectedEdge(g.getNode("X1"), g.getNode("X2"));
        g.addDirectedEdge(g.getNode("X2"), g.getNode("X3"));
        g.addDirectedEdge(g.getNode("X3"), g.getNode("X4"));
        g.addDirectedEdge(g.getNode("X4"), g.getNode("X5"));
        */
        GeneralizedSemPm pm = GaussianCategoricalPm(g, "Split(-1.5,-1,1,1.5)");
        System.out.println(pm);

        System.out.println("STARTS WITH");
        System.out.println(pm.getStartsWithParameterTemplate("C"));

        try{
            MixedUtils.setStartsWith("C", "Split(-.9,-.5,.5,.9)", pm);
        }
        catch (Throwable t) {t.printStackTrace();}

        System.out.println("STARTS WITH");
        System.out.println(pm.getStartsWithParameterTemplate("C"));


        System.out.println(pm);


        GeneralizedSemIm im = GaussianCategoricalIm(pm);
        System.out.println(im);

        int samps = 15;
        DataSet ds = im.simulateDataAvoidInfinity(samps, false);
        System.out.println(ds);

        System.out.println("num cats " + ((DiscreteVariable) g.getNode("X4")).getNumCategories());

        /*Graph trueGraph = DataGraphUtils.loadGraphTxt(new File(MixedUtils.class.getResource("test_data").getPath(), "DAG_0_graph.txt"));
        HashMap<String, Integer> nd = new HashMap<>();
        List<Node> nodes = trueGraph.getNodes();
        for(int i = 0; i < nodes.size(); i++){
            int coin = RandomUtil.getInstance().nextInt(2);
            int dist = (coin==0) ? 0 : 3; //continuous if coin == 0
            nd.put(nodes.get(i).getNode(), dist);
        }
        //System.out.println(getEdgeParams(g.getNode("X3"), g.getNode("X2"), pm).toString());
        //System.out.println(getEdgeParams(g.getNode("X4"), g.getNode("X3"), pm).toString());
        //System.out.println(getEdgeParams(g.getNode("X5"), g.getNode("X4"), pm).toString());
        //System.out.println("num cats " + ((DiscreteVariable) g.getNode("X4")).getNumCategories());
        /*
        HashMap<String, String> nd2 = new HashMap<>();
        nd2.put("X1", "Norm");
        nd2.put("X2", "Norm");
        nd2.put("X3", "Disc");
        nd2.put("X4", "Disc");
        nd2.put("X5", "Disc");
        GeneralizedSemPm pm2 = GaussianTrinaryPm(g, nd2, 10, "Split(-1.5,-.5,.5,1.5)");
        System.out.println("OLD pm:");
        System.out.print(pm2);
        */
    }
}
