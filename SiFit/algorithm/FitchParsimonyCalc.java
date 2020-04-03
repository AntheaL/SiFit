/**
 * 
 */
package SiFit.algorithm;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.Set;

import SiFit.BasicUtilityFunctions;
import SiFit.TopologyBranchPerturbations;
import SiFit.io.VariantMatrixReader;
import cern.colt.Arrays;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * @author hz22
 *
 */
public class FitchParsimonyCalc {

	/**
	 * @param args
	 */
	// inner classes
	/**
	 * This class holds parsimony-computation related information for each node in the tree.
	 */
	protected static class ParsimonyInfo {

		// fields
		/**
		 * The sequence at this node.  Generated for internal nodes, set in leaves.
		 */
		public Integer[] sequence;

		/**
		 * A data structure to hold the list of possible values that can be assigned to this node, when
		 * considering the parsimony score of its subtree.
		 */
		public ArrayList<Set<Integer>> setGTNode = new ArrayList<>();

		/**
		 * The node in the input tree that this node corresponds to.
		 */
		public TMutableNode node;

		// constructors
		public ParsimonyInfo(TMutableNode node, int len) {
			this.node = node;
			for (int i = 0; i < len; i++){
				setGTNode.add(new HashSet<Integer>());
			}
			
		}
	}
	
	public STITree<ParsimonyInfo> nuPtree;
	
	/**
	 * Compute the parsimony score of the tree given the specified set of
	 * genotype sequences.  The branch lengths will be set to be 
	 * the number of mutations along that edge.
	 * @param tree
	 * @param singleCellNames
	 * @param genotypeSeqs
	 * @return
	 */
	public int computeFitchParsimony(MutableTree tree, String[] singleCellNames, ArrayList<Integer[]> genotypeSeqs){
		return computeFitchParsimony(tree, singleCellNames, genotypeSeqs, null);
	}
	
	public int computeFitchParsimony(MutableTree tree, 
			String[] singleCellNames, ArrayList<Integer[]> genotypeSeqs,
			Hashtable<TNode,Integer[]> assignments){
		
		assert singleCellNames.length == genotypeSeqs.size();
		
		Hashtable<String, Integer[]> cellGTSeqMap = new Hashtable<>();
		for (int i = 0; i < singleCellNames.length; i++){
			cellGTSeqMap.put(singleCellNames[i], genotypeSeqs.get(i));
		}
		
		// make the tree that will hold the assignments
		STITree<ParsimonyInfo> ptree = new STITree<ParsimonyInfo>();
		
		int pscore = computeFitchParsimony(tree, cellGTSeqMap, genotypeSeqs.get(0).length, ptree);
		// record the assignments
		if(assignments != null) {
			for(STINode<ParsimonyInfo> node : ptree.getNodes()) {
				assignments.put(node.getData().node, node.getData().sequence);
			}
		}
		
		this.nuPtree = ptree;
		return pscore;
	}
	
	protected int computeFitchParsimony(MutableTree tree, 
			Hashtable<String, Integer[]> cellGTSeqMap,
			int seqLen, STITree<ParsimonyInfo> ptree){
		
		// copy the tree
		copyNode(ptree.getRoot(), tree.getRoot(), cellGTSeqMap, seqLen);
		
		// Perform bottom up and top down phases of Fitch
		for(int pos = 0; pos < seqLen; pos++) {
			BottomUpFitch(pos, ptree.getRoot());
			
			TopDownFitch(pos, ptree.getRoot(), ptree.getRoot().getData().setGTNode.get(pos).iterator().next());
		}
		
		// Compute scores
		return computePScore(ptree.getRoot());
	}
	
	/**
	 * Compute the parsimony score of the subtree pnode.  This method also updates the branch lengths of this subtree (and the node's incoming branch)
	 * to be the number of mutations along that branch.  This method *also* updates the branch lengths in the original
	 * tree provided as input.
	 * @param pnode
	 * @return
	 */
	protected int computePScore(STINode<ParsimonyInfo> pnode) {
		
		int pscore = 0;
		if(pnode.isRoot()) {
			pnode.setParentDistance(0);
		} 
		else {
			Integer[] ps = pnode.getParent().getData().sequence;
			Integer[] seq = pnode.getData().sequence;
			
			int num_diffs = 0;
			
			for(int i = 0; i < ps.length; i++) {
				if(ps[i] != seq[i]) {
					num_diffs++;
				}
			}

			pnode.setParentDistance(num_diffs);
		}
		
		// set the mirror node's distance
//		pnode.getData().node.setParentDistance(pnode.getParentDistance());
		
		// compute all children's pscores
		for(STINode<ParsimonyInfo> child : pnode.getChildren()) {
			pscore += computePScore(child) + child.getParentDistance();
		}

		return pscore;
		
	}
	
	/**
	 * Perform the bottom-up phase of Fitch's Algo for a single position
	 * @param pos
	 * @param pnode
	 */
	protected void BottomUpFitch(int pos, STINode<ParsimonyInfo> pnode){
		ParsimonyInfo pi = pnode.getData();
		
		// leaf node, populate the set with observed genotype
		if (pnode.isLeaf()){
			// Missing data
			if (pi.sequence[pos] == 3){
				pi.setGTNode.get(pos).add(0);
				pi.setGTNode.get(pos).add(1);
				pi.setGTNode.get(pos).add(2);
			}
			else{
				pi.setGTNode.get(pos).add(pi.sequence[pos]);
			}
		}
		// Root, populate with 0
		else if (pnode.isRoot()){
			
			pi.setGTNode.get(pos).add(0);
			for(STINode<ParsimonyInfo> child : pnode.getChildren()) {
				BottomUpFitch(pos, child);
			}
		}
		// Internal node, compute recursively for children
		else{
//			System.out.println("hello");
			ArrayList<STINode<ParsimonyInfo>> p_children = new ArrayList<>();
			for(STINode<ParsimonyInfo> child : pnode.getChildren()) {
				p_children.add(child);
				BottomUpFitch(pos, child);
			}
			// Intersection is empty, compute union
			if (Collections.disjoint(p_children.get(0).getData().setGTNode.get(pos), p_children.get(1).getData().setGTNode.get(pos)) == true){
				pi.setGTNode.get(pos).addAll(p_children.get(0).getData().setGTNode.get(pos));
				pi.setGTNode.get(pos).addAll(p_children.get(1).getData().setGTNode.get(pos));
			}
			// Compute intersection
			else{
				pi.setGTNode.get(pos).addAll(p_children.get(0).getData().setGTNode.get(pos));
				pi.setGTNode.get(pos).retainAll(p_children.get(1).getData().setGTNode.get(pos));
			}
		}
	}
	
	/**
	 * Performs the top-down phase of Fitch algorithm for a single pos
	 * @param pos
	 * @param pnode
	 * @param pAssignment
	 */
	protected void TopDownFitch(int pos, STINode<ParsimonyInfo> pnode, int pAssignment){
		ParsimonyInfo pi = pnode.getData();
		int a;
		
		if (pi.setGTNode.get(pos).size() >= 1){
			if (pi.setGTNode.get(pos).contains(pAssignment))
				a = pAssignment;
			else
				a = pi.setGTNode.get(pos).iterator().next();
		}
		else
			a = pAssignment;
		
		pi.sequence[pos] = a;
		
		// propagate down to children
		for(STINode<ParsimonyInfo> child : pnode.getChildren()) {
			TopDownFitch(pos, child, a);
		}
	}
	
	
	
	/**
	 * Copy the node and all its children.  This method is for creating a tree that is topologically identical
	 * to the input tree.  The resulting tree has additional fields for holding data used to compute the
	 * parsimony score of the tree.
	 * @param pnode
	 * @param n
	 * @param cellGTSeqMap
	 * @param seqLen
	 */
	protected void copyNode(STINode<ParsimonyInfo> pnode, TMutableNode n, Hashtable<String, Integer[]> cellGTSeqMap, int seqLen){
		ParsimonyInfo pi = new ParsimonyInfo(n, seqLen);
		pnode.setData(pi);
		if(n.isLeaf()) {
			pi.sequence = cellGTSeqMap.get(n.getName());
			if(pi.sequence == null) {
				throw new RuntimeException("No sequence provided for leaf node " + n.getName());
			}
		}
		else{
			pi.sequence = new Integer[seqLen];
			for(TMutableNode c : n.getChildren()) {
				STINode<ParsimonyInfo> pchild = pnode.createChild();
				copyNode(pchild, c, cellGTSeqMap, seqLen);
			}
		}
	}
	
	/**
	 * Find the parsimony tree from the data
	 * @param tree
	 * @param singleCellNames
	 * @param genotypeSeqs
	 * @param iter
	 * @param nCell
	 * @param TBP
	 * @return
	 */
	public STITree<Double> findParsimonyTree(
			STITree<Double> tree, String[] singleCellNames, 
			ArrayList<Integer[]> genotypeSeqs, int iter, 
			int nCell, TopologyBranchPerturbations TBP){
//		int bestPScore = (int) Double.POSITIVE_INFINITY;
		STITree<Double> bestTree = tree;
		int bestPScore = this.computeFitchParsimony(tree, singleCellNames, genotypeSeqs);
		for (int i = 0; i < iter; i++){
			STITree<Double> curTree = bestTree;
			STITree<Double> newTree = TBP.proposeTreeParsimony(curTree, nCell);
			int newTreeScore = this.computeFitchParsimony(newTree, singleCellNames, genotypeSeqs);
//			System.out.printf("score of new tree = %d\n", newTreeScore);
			if (newTreeScore < bestPScore){				
				bestPScore = newTreeScore;
				bestTree = newTree;
			}
			if (i%1000 == 0)
				System.out.println(bestPScore);
		}
		return bestTree;
	}
	
	/**
	 * Find the parsimony tree, to be directly used in phylotree search
	 * @param newickTree
	 * @param singleCellList
	 * @param genotypeSeqs
	 * @param iter
	 * @param nCell
	 * @param BUF
	 * @return
	 */
	public STITree<Double> findParsimonyTree(
			String newickTree, ArrayList<String> singleCellList,
			ArrayList<Integer[]> genotypeSeqs, int iter, 
			int nCell, BasicUtilityFunctions BUF){
		String[] singleCellNameArray = new String[nCell];
		for (int i=0; i< nCell; i++){
			singleCellNameArray[i] = singleCellList.get(i);
		}
		TopologyBranchPerturbations TBP = new TopologyBranchPerturbations();
		STITree<Double> randTree1 = BUF.getTree(newickTree);
		
		STITree<Double> randTree = this.findParsimonyTree(randTree1, singleCellNameArray, genotypeSeqs, iter, nCell, TBP);
		return randTree;
	}
	
	public void annotateMutations(String treeFile, String varMatFile) throws IOException{
		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		String intreeNewick = BUF.readNewickString(treeFile);
		STITree<Double> intree = BUF.getTree(intreeNewick);
		intree.getRoot().setName("R");
		int i = 1;
		for (STINode<Double> node: intree.getNodes()){			
			if (node.getName() == ""){
				node.setName("in"+ Integer.toString(i));
				i++;
			}
		}
		System.out.println(intree.toNewick());
		
		int nCell = intree.getLeafCount();
		VariantMatrixReader vr = new VariantMatrixReader(varMatFile);
		vr.populateVarGTMatrix(varMatFile, nCell);
		ArrayList<Integer[]> obsGenotypeMatrix = vr.obsVarGTMatrix;
		int nMut = vr.n_mutation;
		ArrayList<Integer[]> seqGT = BUF.getCellGenomes(obsGenotypeMatrix, nMut);
		
	}
	
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		
		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		String varFileName = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/beatSCITE/testFitch.txt";
		String treefile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/beatSCITE/testFtree.newick";
		
		

		
		FitchParsimonyCalc FPC = new FitchParsimonyCalc();
		
		FPC.annotateMutations(treefile, varFileName);
		

	}

}

