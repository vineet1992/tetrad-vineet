package edu.pitt.csb.mgm;

import java.util.ArrayList;

/**
 * Created by vinee_000 on 9/1/2016.
 */
public class ReKSTree {
    public ReKSNode root;

    public ReKSTree()
    {
        root = new ReKSNode();
    }




    private class ReKSNode
    {
        public ReKSNode parent; //parent node
        public ArrayList<ReKSNode> children; //list of child nodes
        public int index; //if it's a leaf node, then the variable it refers to
        public boolean isLeaf;

        public ReKSNode()
        {

        }

    }
}
