/**
 * 
 */
package SiFit.metric;

import java.io.IOException;
import java.util.ArrayList;


import SiFit.BasicUtilityFunctions;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * @author hz22
 *
 */
public class CompareTrees {
	
    private static int numLeafs;
    private static STITree<Double> reference_tree;
    
    public CompareTrees(STITree<Double> ref) {
        reference_tree = ref;
        numLeafs = reference_tree.getLeafCount();
    }
    
    /**
     * Borrowing from the SymmetricDifference file already coded, returns the symmetric difference between two trees.
     * (Uses the Robinson-Foulds method for comparing two trees)
     * @param tree, the tree being compared to the reference tree
     * @return the score that measures the difference between the two trees.
     */
    public double calcSymmetricDist(STITree<Double> tree) {
        SymmetricDifference sd = new SymmetricDifference();
        sd.computeDifference(tree, reference_tree, true);
        return sd.getWeightedAverage();
    }
    
    
    public double compareSiFitTrees(String true_filename, String SiFit_filename, String doublet_filename) throws IOException{
    	BasicUtilityFunctions BUF = new BasicUtilityFunctions();
    	String trueTreeNewick = BUF.readNewickString(true_filename);
    	STITree<Double> trueTree = BUF.getTree(trueTreeNewick);
    	String SiFitTreeNewick = BUF.readNewickString(SiFit_filename);
    	STITree<Double> SiFitTree = BUF.getTree(SiFitTreeNewick);
    	trueTree.getRoot().setName("R");
    	for (STINode<Double> n : trueTree.getNodes()){
    		if (n.getName().equals("")){
    			n.setName("IN"+Integer.toString(n.getID()));
    		}
    	}
    	
    	SiFitTree.getRoot().setName("R");
    	for (STINode<Double> n : SiFitTree.getNodes()){
    		if (n.getName().equals("")){
    			n.setName("IN"+Integer.toString(n.getID()));
    		}
    	}
    	
    	ArrayList<String> doubletList = new ArrayList<>();
    	String doubletFlagString = BUF.readNewickString(doublet_filename);
    	for (int i = 0; i < doubletFlagString.length(); i++){
    		if (doubletFlagString.charAt(i) == '1'){
    			String doubletName = "sc" + Integer.toString(i+1);
    			doubletList.add(doubletName);
    		}
    	}
    	
    	double treeError = this.CompareTreesDoublets(trueTree, SiFitTree, doubletList);
    	return treeError;
    }
    
    public double CompareTreesDoublets(STITree<Double> trueTree,
    		STITree<Double> testTree, ArrayList<String> doubletList){
    	BasicUtilityFunctions BUF = new BasicUtilityFunctions();
    	STITree<Double> trueTreeDoubletsRmvd = this.getTreeDoubletsRemoved(trueTree, doubletList, BUF);
    	STITree<Double> testTreeDoubletsRmvd = this.getTreeDoubletsRemoved(testTree, doubletList, BUF);
    	CompareTrees CompareTreeObj = new CompareTrees(trueTreeDoubletsRmvd);
    	double dist = CompareTreeObj.calcSymmetricDist(testTreeDoubletsRmvd);
    	return dist;
    }
    
    public void getTreeSingleDoubletRemoved(STITree<Double> tree, String doubletName, BasicUtilityFunctions BUF){
    	
    	STINode<Double> node2remove = tree.getNode(doubletName);
    	STINode<Double> node2removeParent = node2remove.getParent();
    	STINode<Double> node2removeSib = BUF.getOtherChild(node2removeParent, node2remove);
    	STINode<Double> node2removeGrandParent = node2removeParent.getParent();
    	node2removeGrandParent.adoptChild(node2removeSib);
    	tree.removeNode(doubletName);
    	node2removeGrandParent.removeChild(node2removeParent, true);
//    	return tree;
    }
    
    public STITree<Double> getTreeDoubletsRemoved(STITree<Double> tree, ArrayList<String> doubletList, BasicUtilityFunctions BUF){
    	for (String doubletName : doubletList){
    		getTreeSingleDoubletRemoved(tree, doubletName, BUF);
    	}
    	return tree;
    }
    
    public double calcSymmetricDistInt(STITree<Integer> tree) {
        SymmetricDifference sd = new SymmetricDifference();
        sd.computeDifference(tree, reference_tree, true);
        return sd.getWeightedAverage();
    }
    
    public double calcFPDist(STITree<Double> tree) {
        SymmetricDifference sd = new SymmetricDifference();
        sd.computeDifference(tree, reference_tree, true);
        return sd.getWeightedFalsePositive();
    }
    
    public double calcFNDist(STITree<Double> tree) {
        SymmetricDifference sd = new SymmetricDifference();
        sd.computeDifference(tree, reference_tree, true);
        return sd.getWeightedFalseNegative();
    }
    
    public static void main(String[] args){
//    	String s = "((a,b),((c,d),e));";
    	String s = "(((sc6:9,(sc5:2,sc2:13):1):2,((sc7:11,sc1:12):15,(sc8:8,sc10:10):4):4):3,((sc9:7,sc4:5):9,sc3:8):16);";
    	BasicUtilityFunctions BUF = new BasicUtilityFunctions();
    	STITree<Double> Tree = BUF.getTree(s);
    	Tree.getRoot().setName("R");
    	// Name each internal node
    	for (STINode<Double> n : Tree.getNodes()){
    		if (n.getName().equals("")){
    			n.setName("IN"+Integer.toString(n.getID()));
    		}
    	}
    	System.out.println(Tree.toNewick());
    	CompareTrees CT = new CompareTrees(Tree);
    	ArrayList<String> doublets = new ArrayList<>();
    	doublets.add("sc6");
    	doublets.add("sc4");
    	STITree<Double> nuTree = CT.getTreeDoubletsRemoved(Tree, doublets, BUF);
    	System.out.println(nuTree.toNewick());
    	
    	String sd = "00101001";
    	for (int i = 0; i < sd.length(); i++){
    		if (sd.charAt(i) == '1'){
    			String cellName = "sc" + Integer.toString(i+1);
    			System.out.println(cellName);
    		}
    	}
    }

}
