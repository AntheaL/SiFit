package SiFit;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.distribution.PoissonDistribution;

import SiFit.objects.GenotypeObservation;
import SiFit.objects.NodeDiploidDNA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

public class BasicUtilityFunctions {

	private static Random _rng;
	
	/**
	 * Constructor with PRNG.
	 */
	public BasicUtilityFunctions (Random rng) {
		this._rng = rng;
	}	
	
	/**
	 * Default Constructor:
	 */
	public BasicUtilityFunctions () {
		this._rng = new Random(0);
	}
	
	/*
	 * Functions required for data simulation
	 */
	
	/**
	 * Read the reference DNA from Fasta File
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	public String readFastaDNA(String filename) throws IOException{		
		String DNA = "";
		try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
		    String line;
		    while ((line = br.readLine()) != null) {
		    	if (line.startsWith(">")) {continue;}
		    	else{
		    		DNA = DNA.concat(line);
		    	}
		    }
		}
		return DNA;
	}
	
	/**
	 * Get the genotype array at the root (all 0s)
	 * @param nMut
	 * @return
	 */
	public Integer[] getRootGenotypeArr (int nMut){
		Integer[] rootGenotypeArr = new Integer[nMut];
		for (int i = 0; i < nMut; i++){
			rootGenotypeArr[i] = 0;
		}
		return rootGenotypeArr;
	}
	
	public STITree<Double> generateRandomTree(ArrayList<String> scNames, double thr){
		return randomTreeHelper(scNames, scNames.size(), thr);
	}
	
	/**
	 * Helper function for generateRandomTree
	 * @param scNames - ArrayList of leaves (single cell names) 
	 * @param n_cells - Number of leaves
	 * @param thr     - Branch Length Threshold
	 * @return a STITree object representing the randomly generated tree
	 */
	public static STITree<Double> randomTreeHelper(ArrayList<String> scNames, Integer n_cells, double thr){
		
		// If there is one cell, create a tree with one node
		if (n_cells == 1){
			STITree<Double> cell = new STITree<>(scNames.get(0), true);
			double bl = _rng.nextDouble()*thr;
			cell.getRoot().setParentDistance(bl);
			return cell;
		}
		
		// Randomly split the list of leafs into two lists.
		ArrayList<String> leftCells = new ArrayList<>();
		ArrayList<String> rightCells = new ArrayList<>();
		for (int i = 0; i < n_cells; i++){
			if (_rng.nextDouble() <= 0.5){
				leftCells.add(scNames.get(i));
			}
			else
				rightCells.add(scNames.get(i));				
		}
		
		// left split is empty
		if (leftCells.size() == 0){
			return randomTreeHelper(rightCells, rightCells.size(), thr);
		}
		// right split is empty
		if (rightCells.size() == 0){
			return randomTreeHelper(leftCells, leftCells.size(), thr);
		}
		
		// Root the randomly generated subtrees together.
		STITree<Double> leftTree = randomTreeHelper(leftCells, leftCells.size(), thr);
		STITree<Double> rightTree = randomTreeHelper(rightCells, rightCells.size(), thr);
		STITree<Double> newTree = new STITree<>(true);
		newTree.getRoot().adoptChild(leftTree.getRoot());
        newTree.getRoot().adoptChild(rightTree.getRoot());
        
        // Randomly assign branch lengths between the subtrees' roots and the new tree's root.
        double lbl = _rng.nextDouble()*thr;
        if (lbl < thr/5)
        	lbl = _rng.nextDouble()*thr;
        leftTree.getRoot().setParentDistance(lbl);
        double rbl = _rng.nextDouble()*thr;
        if (rbl < thr/5)
        	rbl = _rng.nextDouble()*thr;
        rightTree.getRoot().setParentDistance(rbl);
        
        return newTree;
	}
	
	/**
	 * Change the branch length of the leaf nodes whose other sibling is also leaf node
	 * @param leafList - List of leaves
	 * @param Tree - input tree
	 * @param lb - lower bound of new branch length
	 * @param ub - upper bound of new branch length
	 */
	public void resizeLeafBranches(ArrayList<String> leafList, STITree<Double> Tree, double lb, double ub){
		for (String s : leafList){
			STINode<Double> leaf = Tree.getNode(s);
			List<TNode> siblings = leaf.getSiblings();
			if (siblings.get(0).isLeaf())
				leaf.setParentDistance(lb + _rng.nextDouble()*ub);				
		}
	}
	
	/**
	 * Resize branch lengths according to the threshold
	 * @param Tree
	 * @param thr
	 */
	public void resizeBranchLengths(STITree<Double> Tree, double thr){
		double sumBL = 0;
		for (STINode<Double> node: Tree.getNodes()){
			if (!node.isRoot()){
//			System.out.println(node.getParentDistance());
			sumBL += node.getParentDistance();
			}
		}
//		System.out.println(sumBL);
		for (STINode<Double> node: Tree.getNodes()){
			if (!node.isRoot()){
			double nuBranchlength = (node.getParentDistance()*thr)/sumBL;
			node.setParentDistance(nuBranchlength);
			}
		}
	}
	
	public STITree<Double> getTumorNormalCellTree(int normalCellCount, int tumorCellCount, 
			double normalBLThr, double tumorBLThr, double sumBLThr){
		ArrayList<String> normalNames = new ArrayList<>();
		for (int i=1; i <= normalCellCount; i++){
			normalNames.add("sc" + Integer.toString(i));
		}
		ArrayList<String> tumorNames = new ArrayList<>();
		for (int i=normalCellCount+1; i <= normalCellCount+tumorCellCount; i++){
			tumorNames.add("sc" + Integer.toString(i));
		}
		
		STITree<Double> normalSubTreeMaker = this.generateRandomTree(normalNames, normalBLThr);
		STITree<Double> normalSubTree = this.getTree(normalSubTreeMaker.toNewick());
		STINode<Double> normalSubTreeRoot = normalSubTree.getRoot();
		normalSubTreeRoot.setName("normalR");
		
		STITree<Double> tumorSubTreeMaker = this.generateRandomTree(tumorNames, tumorBLThr);
		STITree<Double> tumorSubTree = this.getTree(tumorSubTreeMaker.toNewick());
		STINode<Double> tumorSubTreeRoot = tumorSubTree.getRoot();
		tumorSubTreeRoot.setName("tumorR");
		this.resizeLeafBranches(tumorNames, tumorSubTree, tumorBLThr/100, tumorBLThr/10);
		
		STITree<Double> root = new STITree<>();
		STINode<Double> rootTree = root.getRoot();
		rootTree.setName("R");
		rootTree.adoptChild(tumorSubTreeRoot);
		rootTree.adoptChild(normalSubTreeRoot);
		tumorSubTreeRoot.setParentDistance(_rng.nextDouble()*tumorBLThr);
		normalSubTreeRoot.setParentDistance(_rng.nextDouble()*normalBLThr);
		
		STITree<Double> finalTree = this.getTree(root.toNewick());
		
		this.resizeBranchLengths(finalTree, sumBLThr);
		
		return finalTree;
	}
	
	public boolean checkArrayMember(int[] arr, int num, int len){
		for (int i=0; i<len; i++){
			if (arr[i] == num){return true;}
		}
		return false;
	}
	
	/**
	 * Given the string with vcf headers returns an array containing the names of the single cells
	 * @param vcf_header, single_cell_names
	 */
	public void getSingleCellNames(String vcf_header, ArrayList<String> single_cell_names){
		String[] vcf_header_entry_arr = vcf_header.split("\t");
		for (int i = 0; i < vcf_header_entry_arr.length - 9; i++){
			single_cell_names.add(vcf_header_entry_arr[9+i]) ;
		}
	}
	
	/**
	 * Given a list of leaves (single cell names), returns a randomly generated tree.
	 * @param scNames - ArrayList of leaves (single cell names) 
	 * @return a STITree object representing the randomly generated tree
	 */
	public STITree<Double> generateRandomTree(ArrayList<String> scNames){
		return randomTreeHelper(scNames, scNames.size());
	}
	
	/**
	 * Helper function for generateRandomTree
	 * @param scNames - ArrayList of leaves (single cell names) 
	 * @param n_cells - Number of leaves
	 * @return a STITree object representing the randomly generated tree
	 */
	public STITree<Double> randomTreeHelper(ArrayList<String> scNames, Integer n_cells){
		
		// If there is one cell, create a tree with one node
		if (n_cells == 1){
			STITree<Double> cell = new STITree<>(scNames.get(0), true);
			cell.getRoot().setParentDistance(_rng.nextDouble());
			return cell;
		}
		
		// Randomly split the list of leafs into two lists.
		ArrayList<String> leftCells = new ArrayList<>();
		ArrayList<String> rightCells = new ArrayList<>();
		for (int i = 0; i < n_cells; i++){
			if (_rng.nextDouble() <= 0.5){
				leftCells.add(scNames.get(i));
			}
			else
				rightCells.add(scNames.get(i));				
		}
		
		// left split is empty
		if (leftCells.size() == 0){
			return randomTreeHelper(rightCells, rightCells.size());
		}
		// right split is empty
		if (rightCells.size() == 0){
			return randomTreeHelper(leftCells, leftCells.size());
		}
		
		// Root the randomly generated subtrees together.
		STITree<Double> leftTree = randomTreeHelper(leftCells, leftCells.size());
		STITree<Double> rightTree = randomTreeHelper(rightCells, rightCells.size());
		STITree<Double> newTree = new STITree<>(true);
		newTree.getRoot().adoptChild(leftTree.getRoot());
        newTree.getRoot().adoptChild(rightTree.getRoot());
        
        // Randomly assign branch lengths between the subtrees' roots and the new tree's root.
        leftTree.getRoot().setParentDistance(_rng.nextDouble());
        rightTree.getRoot().setParentDistance(_rng.nextDouble());
        
        return newTree;
	}
	
	/**
	 * Returns the map of DNA and single cell names
	 * @param sc_names
	 * @param sb_DNA
	 * @return DNA_map
	 */
	public Map<String,String> getMapCellDNA(ArrayList<String> sc_names, StringBuilder[] sb_DNA){
		Map<String,String> DNA_map = new HashMap<String,String>();
		for (int i = 0; i < sb_DNA.length; i++){
			DNA_map.put(sc_names.get(i), sb_DNA[i].toString());
		}
		return DNA_map;
	}
	
	/**
	 * Returns the map of DNA and single cell names from array of diploid DNA
	 * @param sc_names
	 * @param diploid_DNA_arr
	 * @param flag
	 * @return
	 */
	public Map<String, String> getMapCellDNA(ArrayList<String> sc_names, NodeDiploidDNA[] diploid_DNA_arr, int flag){
		Map<String,String> DNA_map = new HashMap<String,String>();
		if (flag == 1){
			for (int i = 0; i < sc_names.size(); i++){
				DNA_map.put(sc_names.get(i), diploid_DNA_arr[i].DNA_copy_1);
			}
		}
		else{
			for (int i = 0; i < sc_names.size(); i++){
				DNA_map.put(sc_names.get(i), diploid_DNA_arr[i].DNA_copy_2);
			}
		}
		return DNA_map;
	}
	
	/**
	 * Computes the prior of fn from input _fnstart and observed GT matrix
	 * @param GT, observed genotype matrix
	 * @param fnPrior, fn value to start with
	 * @return
	 */
	public double getFNPrior(ArrayList<Integer[]> GT, double fnPrior){
		double nume = 0.0;
		double total = GT.size()*GT.get(0).length;
		for (int i = 0; i < GT.size(); i++){
			for (int j = 0; j < GT.get(i).length; j++){
				if (GT.get(i)[j] == 1)
					nume += 1;
			}
		}		
		return ((fnPrior/(1-fnPrior))*nume*(1/total));
	}
	
	public double[][] computeTransMat(double mu, double omega, double d){
		double[][] TransMat = new double[3][3];
		TransMat[0][0] = 1 - mu - mu*mu;
		TransMat[0][1] = mu;
		TransMat[0][2] = mu*mu;
		TransMat[1][0] = omega/2 + d;
		TransMat[1][1] = 1 - omega - d;
		TransMat[1][2] = omega/2;
		TransMat[2][0] = d*d;
		TransMat[2][1] = d;
		TransMat[2][2] = 1 - TransMat[2][0] - TransMat[2][1];
		return TransMat;
	}
	public double[][] computeTransMat(double mu, double omega, double d, double t){
		double[][] TransMat = new double[3][3];
		TransMat[0][0] = 1 - mu*t - mu*mu*t;
		TransMat[0][1] = mu*t;
		TransMat[0][2] = mu*mu*t;
		TransMat[1][0] = (omega/2 + d)*t;
		TransMat[1][1] = 1 - (omega + d)*t;
		TransMat[1][2] = t*omega/2;
		TransMat[2][0] = d*d*t;
		TransMat[2][1] = d*t;
		TransMat[2][2] = 1 - TransMat[2][0] - TransMat[2][1];
		return TransMat;
	}
	
	public Integer[] getMutatedGenotypeArrFSM(
			Integer[] parentGenotypeArr,
			HashMap<Integer, Set<Integer>> MutationTypeMap,
			double[][] JCTransMat, int nMut){
		
		Integer[] childGenotypeArr = new Integer[nMut];
		for (int i = 0; i < nMut; i++){
			if (parentGenotypeArr[i] == 0){
				double rr = _rng.nextDouble();
				// 0 -> 2
				if (rr <= JCTransMat[0][2]){
//					System.out.println("0->2 happened");
					childGenotypeArr[i] = 2;
					MutationTypeMap.get(2).add(i);
				}
				// 0 -> 1
				else if (rr <= JCTransMat[0][1] + JCTransMat[0][2]){
					childGenotypeArr[i] = 1;
					MutationTypeMap.get(1).add(i);
				}
				// 0 -> 0
				else {
					childGenotypeArr[i] = 0;
				}
			}
			else if (parentGenotypeArr[i] == 1){
				double rr = _rng.nextDouble();
				// 1 -> 2
				if (rr <= JCTransMat[1][2]){
					childGenotypeArr[i] = 2;
					MutationTypeMap.get(2).add(i);
					MutationTypeMap.get(1).remove(i);
//					System.out.printf("%s\n", "infinite sites failed, Parallel mutation (1 -> 2) happened");
				}
				// 1 -> 0
				else if (rr <= JCTransMat[1][0] + JCTransMat[1][2]){
					childGenotypeArr[i] = 0;
//					System.out.printf("%s\n", "infinite sites failed, Back mutation (1 -> 0) happened");
				}
				// 1 -> 1
				else{
					childGenotypeArr[i] = 1;
				}
			}
			else {
				double rr = _rng.nextDouble();
				// 2 -> 0
				if (rr <= JCTransMat[2][0]){
					childGenotypeArr[i] = 0;
					
				}
				// 2 -> 1
				else if (rr <= JCTransMat[2][0] + JCTransMat[2][1]){
					childGenotypeArr[i] = 1;
				}
				// 2 -> 2
				else {
					childGenotypeArr[i] = parentGenotypeArr[i];
				}
			}
		}
		return childGenotypeArr;
	}
	
	public Integer[] getMutatedGenotypeArrFSM(
			Integer[] parentGenotypeArr,
			HashMap<Integer, Set<Integer>> MutationTypeMap,
			double[][] JCTransMat, double branchLength, int nMut){
		Integer[] childGenotypeArr = new Integer[nMut];
		int PoissonParam = 1+(int) Math.round(branchLength*nMut);		
		int nuMut2introduce = new PoissonDistribution(PoissonParam).sample();
		double mutProb = nuMut2introduce/(double)nMut;
//		System.out.printf("mean=%f\tmut=%d\n", mutProb, nuMut2introduce);
		for (int i = 0; i< parentGenotypeArr.length; i++){
			double rr0 = _rng.nextDouble();
			if (rr0 <= mutProb){
				if (parentGenotypeArr[i] == 0){
					double rr = _rng.nextDouble();
					// 0 -> 2
					if (rr <= JCTransMat[0][2]){
//						System.out.println("0->2 happened");
						childGenotypeArr[i] = 2;
						MutationTypeMap.get(2).add(i);
					}
					// 0 -> 1
					else if (rr <= JCTransMat[0][1] + JCTransMat[0][2]){
						childGenotypeArr[i] = 1;
						MutationTypeMap.get(1).add(i);
					}
					// 0 -> 0
					else {
						childGenotypeArr[i] = 0;
					}
				}
				else if (parentGenotypeArr[i] == 1){
					double rr = _rng.nextDouble();
					// 1 -> 2
					if (rr <= JCTransMat[1][2]){
						childGenotypeArr[i] = 2;
						MutationTypeMap.get(2).add(i);
						MutationTypeMap.get(1).remove(i);
//						System.out.printf("%s\n", "infinite sites failed, Parallel mutation (1 -> 2) happened");
					}
					// 1 -> 0
					else if (rr <= JCTransMat[1][0] + JCTransMat[1][2]){
						childGenotypeArr[i] = 0;
//						System.out.printf("%s\n", "infinite sites failed, Back mutation (1 -> 0) happened");
					}
					// 1 -> 1
					else{
						childGenotypeArr[i] = 1;
					}
				}
				else {
					double rr = _rng.nextDouble();
					// 2 -> 0
					if (rr <= JCTransMat[2][0]){
						childGenotypeArr[i] = 0;
						
					}
					// 2 -> 1
					else if (rr <= JCTransMat[2][0] + JCTransMat[2][1]){
						childGenotypeArr[i] = 1;
					}
					// 2 -> 2
					else {
						childGenotypeArr[i] = parentGenotypeArr[i];
					}
				}
			}
			else{
				childGenotypeArr[i] = parentGenotypeArr[i];
			}
		}
		return childGenotypeArr;
	}
	
	/**
	 * Mutate parent genotype vector to child one, binary data
	 * @param child
	 * @param parentGenotypeArr
	 * @param MutationTypeMap
	 * @param mutFlagArr
	 * @param recurProb
	 * @param omega
	 * @param deletion
	 * @param branchLength
	 * @param nMut
	 * @return
	 */
	public Integer[] getMutatedGenotypeArrFSMBinary(
			STINode<Double> child,
			Integer[] parentGenotypeArr,
			HashMap<Integer, Set<Integer>> MutationTypeMap,
			Integer[] mutFlagArr,
			double recurProb, double omega,
			double deletion, double branchLength,
			int nMut){
		
		Integer[] childGenotypeArr = new Integer[nMut];
		int nuMut2introduce;
		if (!child.isLeaf()){
			int PoissonParam = 1+(int) Math.round(branchLength*nMut);
			if (PoissonParam == 0)
				PoissonParam = 1;
			nuMut2introduce = new PoissonDistribution(PoissonParam).sample();
		}
		else
			nuMut2introduce = (int) Math.round(branchLength*nMut);
		for (int i = 0; i< parentGenotypeArr.length; i++){
			if (parentGenotypeArr[i] == 1){
				double rrw = _rng.nextDouble();
				if (rrw <= omega){

					double rrw1 = _rng.nextDouble();
					childGenotypeArr[i] = 0;
				}
				else{
					childGenotypeArr[i] = parentGenotypeArr[i];
				}
			}
			else{
				childGenotypeArr[i] = parentGenotypeArr[i];
			}
		}
		double rr = _rng.nextDouble();
		// Deletion happening
		if (rr <= deletion){
			for (int i = 0; i < nuMut2introduce; i++){
				if (MutationTypeMap.get(1).size() == 0)
					continue;
				ArrayList<Integer> mutatedPosList = new ArrayList<>();
				mutatedPosList.addAll(MutationTypeMap.get(1));
				int indexInSet = _rng.nextInt(mutatedPosList.size());
				int site_index = mutatedPosList.get(indexInSet);
				if (parentGenotypeArr[site_index] == 1){
					childGenotypeArr[site_index] = 0;
//					System.out.println("deletion to 0");
				}
				else
					childGenotypeArr[site_index] = parentGenotypeArr[site_index];
			}
		}
		else{
			for (int i = 0; i < nuMut2introduce; i++){
				double rr2 = _rng.nextDouble();
				// Recurrent mutation happens
				if (rr2 <= recurProb){
					if (MutationTypeMap.get(1).size() == 0)
						continue;

					ArrayList<Integer> mutatedPosList = new ArrayList<>();
					mutatedPosList.addAll(MutationTypeMap.get(1));
					int indexInSet = _rng.nextInt(mutatedPosList.size());
					int site_index = mutatedPosList.get(indexInSet);
					if (parentGenotypeArr[site_index] == 1)
						childGenotypeArr[site_index] = 0;         // back mutation
					else if (parentGenotypeArr[site_index] == 0)  // parallel mutation
						childGenotypeArr[site_index] = 1;
				}
				else{
					ArrayList<Integer> mutatedPosList = new ArrayList<>();
					// Nu mutation happens
					if (MutationTypeMap.get(0).size() > 0){
						mutatedPosList.addAll(MutationTypeMap.get(0));

						int indexInSet = _rng.nextInt(mutatedPosList.size());
						int site_index = mutatedPosList.get(indexInSet);
						childGenotypeArr[site_index] = 1;
						MutationTypeMap.get(0).remove(site_index);
						MutationTypeMap.get(1).add(site_index);
					}
					// No more new mutation
					else{
						mutatedPosList.addAll(MutationTypeMap.get(1));
						int indexInSet = _rng.nextInt(mutatedPosList.size());
						int site_index = mutatedPosList.get(indexInSet);
						if (parentGenotypeArr[site_index] == 1)
							childGenotypeArr[site_index] = 0;         // back mutation
						else if (parentGenotypeArr[site_index] == 0)  // parallel mutation
							childGenotypeArr[site_index] = 1;
						
					}
				}
			}
		}
		return childGenotypeArr;
					
	}
	
	public Integer[] getMutatedGenotypeArrFSMTest(
			STINode<Double> child,
			Integer[] parentGenotypeArr,
			HashMap<Integer, Set<Integer>> MutationTypeMap,
			Integer[] mutFlagArr,
			double recurProb, double omega,
			double deletion, double branchLength,
			int nMut){
		Integer[] childGenotypeArr = new Integer[nMut];
		int nuMut2introduce;
		if (!child.isLeaf()){
			int PoissonParam = 1+(int) Math.round(branchLength*nMut);
			if (PoissonParam == 0)
				PoissonParam = 1;
			nuMut2introduce = new PoissonDistribution(PoissonParam).sample();
		}
		else
			nuMut2introduce = (int) Math.round(branchLength*nMut);

		for (int i = 0; i< parentGenotypeArr.length; i++){
			childGenotypeArr[i] = parentGenotypeArr[i];
		}
		for (int i = 0; i < nuMut2introduce; i++){
			double rrD = _rng.nextDouble();
			if (rrD <= deletion){
				if (MutationTypeMap.get(1).size() == 0)
					continue;
				ArrayList<Integer> mutatedPosList = new ArrayList<>();
				mutatedPosList.addAll(MutationTypeMap.get(1));
				int indexInSet = _rng.nextInt(mutatedPosList.size());
				int site_index = mutatedPosList.get(indexInSet);
				if (parentGenotypeArr[site_index] == 1){
					double rr1 = _rng.nextDouble();
					if (rr1 <= 0.5){
						childGenotypeArr[site_index] = 0;

					}
					else{
						childGenotypeArr[site_index] = 2;

					}
				}
				else if (parentGenotypeArr[site_index] == 2)
					childGenotypeArr[site_index] = 1;
				else
					childGenotypeArr[site_index] = parentGenotypeArr[site_index];
			}
			else if (rrD <= deletion + omega){
				if (MutationTypeMap.get(1).size() == 0)
					continue;
				ArrayList<Integer> mutatedPosList = new ArrayList<>();
				mutatedPosList.addAll(MutationTypeMap.get(1));
				int indexInSet = _rng.nextInt(mutatedPosList.size());
				int site_index = mutatedPosList.get(indexInSet);
				if (parentGenotypeArr[site_index] == 1){
					double rr1 = _rng.nextDouble();
					if (rr1 <= 0.5){
						childGenotypeArr[site_index] = 0;

					}
					else{
						childGenotypeArr[site_index] = 2;

					}
				}
			}
			else if (rrD <= deletion + omega + recurProb){
				if (MutationTypeMap.get(1).size() == 0)
					continue;

				ArrayList<Integer> mutatedPosList = new ArrayList<>();
				mutatedPosList.addAll(MutationTypeMap.get(1));
				int indexInSet = _rng.nextInt(mutatedPosList.size());
				int site_index = mutatedPosList.get(indexInSet);
				if (parentGenotypeArr[site_index] == 1)
					childGenotypeArr[site_index] = 0;         // back mutation
				else if (parentGenotypeArr[site_index] == 0)  // parallel mutation
					childGenotypeArr[site_index] = 1;
				else
					childGenotypeArr[site_index] = 1;
			
			}
			else{
				ArrayList<Integer> mutatedPosList = new ArrayList<>();
				// Nu mutation happens
				if (MutationTypeMap.get(0).size() > 0){
					mutatedPosList.addAll(MutationTypeMap.get(0));

					int indexInSet = _rng.nextInt(mutatedPosList.size());
					int site_index = mutatedPosList.get(indexInSet);
					childGenotypeArr[site_index] = 1;
					MutationTypeMap.get(0).remove(site_index);
					MutationTypeMap.get(1).add(site_index);
				}
				// No more new mutation
				else{
					mutatedPosList.addAll(MutationTypeMap.get(1));
					int indexInSet = _rng.nextInt(mutatedPosList.size());
					int site_index = mutatedPosList.get(indexInSet);
					if (parentGenotypeArr[site_index] == 1)
						childGenotypeArr[site_index] = 0;         // back mutation
					else if (parentGenotypeArr[site_index] == 0)  // parallel mutation
						childGenotypeArr[site_index] = 1;
					else
						childGenotypeArr[site_index] = 1;
				}
			}
		}
		return childGenotypeArr;
	}
	
	
	public Integer[] getMutatedGenotypeArrFSM(
			STINode<Double> child,
			Integer[] parentGenotypeArr,
			HashMap<Integer, Set<Integer>> MutationTypeMap,
			Integer[] mutFlagArr,
			double recurProb, double omega,
			double deletion, double branchLength,
			int nMut){
		Integer[] childGenotypeArr = new Integer[nMut];
		int nuMut2introduce;
		if (!child.isLeaf()){

			int PoissonParam = 1+(int) Math.round(branchLength*nMut);

			if (PoissonParam == 0)
				PoissonParam = 1;
			nuMut2introduce = new PoissonDistribution(PoissonParam).sample();
		}
		else
			nuMut2introduce = (int) Math.round(branchLength*nMut);

		
		// LOH are introduced
		for (int i = 0; i< parentGenotypeArr.length; i++){
			if (parentGenotypeArr[i] == 1){
				double rrw = _rng.nextDouble();
				if (rrw <= omega){

					double rrw1 = _rng.nextDouble();
					if (rrw1 <= 0.5){
					childGenotypeArr[i] = 2;
					MutationTypeMap.get(2).add(i);
					}
					else{
						childGenotypeArr[i] = 0;
					}
				}
				else
					childGenotypeArr[i] = 1;
			}
			else
				childGenotypeArr[i] = parentGenotypeArr[i];
		}
		double rr = _rng.nextDouble();
		// Deletion happening
		if (rr <= deletion){
			double rr1 = _rng.nextDouble();
			// Deletion results in homozygous reference genotype
			if (rr1 <= 0.5){
				for (int i = 0; i < nuMut2introduce; i++){
					if (MutationTypeMap.get(1).size() == 0)
						continue;
					ArrayList<Integer> mutatedPosList = new ArrayList<>();
					mutatedPosList.addAll(MutationTypeMap.get(1));
					int indexInSet = _rng.nextInt(mutatedPosList.size());
					int site_index = mutatedPosList.get(indexInSet);
					if (parentGenotypeArr[site_index] == 1){
						childGenotypeArr[site_index] = 0;

					}
					else if (parentGenotypeArr[site_index] == 2)
						childGenotypeArr[site_index] = 0;
					else
						childGenotypeArr[site_index] = parentGenotypeArr[site_index];
				}
			}
			// Deletion results in homozygous variant genotype
			else{
				for (int i = 0; i < nuMut2introduce; i++){
					if (MutationTypeMap.get(1).size() == 0)
						continue;
					ArrayList<Integer> mutatedPosList = new ArrayList<>();
					mutatedPosList.addAll(MutationTypeMap.get(1));
					int indexInSet = _rng.nextInt(mutatedPosList.size());
					int site_index = mutatedPosList.get(indexInSet);
					if (parentGenotypeArr[site_index] == 1){
						childGenotypeArr[site_index] = 2;

					}
					else if (parentGenotypeArr[site_index] == 2)
						childGenotypeArr[site_index] = 0;
					else
						childGenotypeArr[site_index] = parentGenotypeArr[site_index];
				}
			}
		}
		else{
			for (int i = 0; i < nuMut2introduce; i++){
				double rr2 = _rng.nextDouble();
				// Recurrent mutation happens
				if (rr2 <= recurProb){
					if (MutationTypeMap.get(1).size() == 0)
						continue;

					ArrayList<Integer> mutatedPosList = new ArrayList<>();
					mutatedPosList.addAll(MutationTypeMap.get(1));
					int indexInSet = _rng.nextInt(mutatedPosList.size());
					int site_index = mutatedPosList.get(indexInSet);
					if (parentGenotypeArr[site_index] == 1)
						childGenotypeArr[site_index] = 0;         // back mutation
					else if (parentGenotypeArr[site_index] == 0)  // parallel mutation
						childGenotypeArr[site_index] = 1;
					else
						childGenotypeArr[site_index] = 1;
				}
				else{
					ArrayList<Integer> mutatedPosList = new ArrayList<>();
					// Nu mutation happens
					if (MutationTypeMap.get(0).size() > 0){
						mutatedPosList.addAll(MutationTypeMap.get(0));

						int indexInSet = _rng.nextInt(mutatedPosList.size());
						int site_index = mutatedPosList.get(indexInSet);
						childGenotypeArr[site_index] = 1;
						MutationTypeMap.get(0).remove(site_index);
						MutationTypeMap.get(1).add(site_index);
					}
					// No more new mutation
					else{
						mutatedPosList.addAll(MutationTypeMap.get(1));

						int indexInSet = _rng.nextInt(mutatedPosList.size());
						int site_index = mutatedPosList.get(indexInSet);
						if (parentGenotypeArr[site_index] == 1)
							childGenotypeArr[site_index] = 0;         // back mutation
						else if (parentGenotypeArr[site_index] == 0)  // parallel mutation
							childGenotypeArr[site_index] = 1;
						else
							childGenotypeArr[site_index] = 1;
					}
				}
			}
			
		}
		return childGenotypeArr;
	}
	
	/**
	 * Returns the number of mutations from parent to child
	 * @param integers
	 * @param childGenotypeArr
	 * @return
	 */
	public int countMutParent2Child(Integer[] integers, Integer[] childGenotypeArr){
		int mutCount = 0;
		for (int i = 0; i < integers.length; i++){
			if (childGenotypeArr[i] != integers[i])
				mutCount += 1;
		}
		return mutCount;
	}
	
	/**
	 * Returns the GenotypeObservation List from single cell names and genotype matrix 
	 * @param sc_names - List of single cell names
	 * @param GTMatrix - genotype matrix
	 * @return GenotypeObsMatrix, ArrayList<GenotypeObservation>
	 */
	public ArrayList<GenotypeObservation> getGenotypeObsFromGTMatrix(ArrayList<String> sc_names, ArrayList<Integer[]> GTMatrix){
		ArrayList<GenotypeObservation> GenotypeObsMatrix = new ArrayList<GenotypeObservation>();
		for (Integer[] gtMut: GTMatrix){
			GenotypeObsMatrix.add(new GenotypeObservation(sc_names, gtMut));
		}
		return GenotypeObsMatrix;
	}
	
	public ArrayList<Integer[]> getCellGenomes(ArrayList<Integer[]> mutProfile, int nMut){
		ArrayList<Integer[]> cellGenomeList = new ArrayList<>();
		for (int i = 0; i< mutProfile.get(0).length; i++){
			Integer[] thisCellGenome = new Integer[nMut];
			for (int j = 0; j< nMut; j++){
				thisCellGenome[j] = mutProfile.get(j)[i];
			}
			cellGenomeList.add(thisCellGenome);
		}
		return cellGenomeList;
		
	}
	
	/**
	 *  Given a Phylogenetic tree in String format, create a Tree data structure.
	 */
	public STITree<Double> getTree (String treeString) {
		NewickReader nr = new NewickReader(new StringReader(treeString));
        STITree<Double> tree = new STITree<Double>(true);
        try {
            nr.readTree(tree);
        }
        catch(Exception ex) {
            ex.printStackTrace();
        }
        
        return tree;
	}
	
	/**
	 *  Given a Phylogenetic tree in String format, create a Tree data structure.
	 */
	public STITree<Integer> getTreeInt (String treeString) {
		NewickReader nr = new NewickReader(new StringReader(treeString));
        STITree<Integer> tree = new STITree<Integer>(true);
        try {
            nr.readTree(tree);
        }
        catch(Exception ex) {
            ex.printStackTrace();
        }
        
        return tree;
	}
	
	/**
	 * Reads Newick string from a file, reads just one line
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	public String readNewickString(String filename) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(filename));
    	String line = null;
    	String this_tree = null;
    	while ((line = br.readLine()) != null) {
    		this_tree = line;
    	}
    	return this_tree;
	}
	
	/**
	 * Writes newick string to a file
	 * @param filename
	 * @param newickStr
	 */
	public void writeNewickTree(String filename, String newickStr){
		try{
			PrintWriter writer = new PrintWriter(filename, "UTF-8");
			writer.println(newickStr);
			writer.close();
		} catch(IOException e){
			e.printStackTrace();
		};
	}
	
	/**
	 * Write genotype matrix to file 
	 * @param filename
	 * @param allMutatedPosList
	 * @param scNames
	 * @param nodeName2ID_Map
	 * @param nodesGenotypeArrMap
	 * @throws FileNotFoundException
	 * @throws UnsupportedEncodingException
	 */
	public void writeGenotypeMat2File(String filename, 
			ArrayList<Integer> allMutatedPosList,
			ArrayList<String> scNames, Integer[][] genotypeMat,
			HashMap<String, Integer> nodeName2ID_Map,
			HashMap<Integer, Integer[]> nodesGenotypeArrMap) throws FileNotFoundException, UnsupportedEncodingException{

		PrintWriter writer = new PrintWriter(filename, "UTF-8");
		
		for (int i = 0; i < allMutatedPosList.size(); i++){
			String this_mut_row = Integer.toString(allMutatedPosList.get(i));
			genotypeMat[i][0] = allMutatedPosList.get(i);
			for (int j = 0; j < scNames.size(); j++){
				int cellID = nodeName2ID_Map.get(scNames.get(j));
				genotypeMat[i][j+1] = nodesGenotypeArrMap.get(cellID)[allMutatedPosList.get(i)];
				this_mut_row = this_mut_row + " " + String.valueOf(genotypeMat[i][j+1]);
			}
			writer.printf("%s\n", this_mut_row);
		}
		writer.close();
	}
	
	public void writeGenotypeMat2File(String filename, 
			int nMut, ArrayList<String> scNames, Integer[][] genotypeMat,
			HashMap<String, Integer[]> nodesGenotypeArrMap) throws FileNotFoundException, UnsupportedEncodingException{
		PrintWriter writer = new PrintWriter(filename, "UTF-8");
		for (int i = 0; i < nMut; i++){
			String this_mut_row = Integer.toString(i);
			genotypeMat[i][0] = i;
			for (int j = 0; j < scNames.size(); j++){
				String cell = scNames.get(j);
				genotypeMat[i][j+1] = nodesGenotypeArrMap.get(cell)[i];
				this_mut_row = this_mut_row + " " + String.valueOf(genotypeMat[i][j+1]);
			}
			writer.printf("%s\n", this_mut_row);
		}
		writer.close();
	}
	
	/**
	 * Add missing data to a genotype matrix
	 * @param ingenotypeMat
	 * @param outgenotypeMat
	 * @param missing
	 */
	public void addMissingData2GenotypeMat(Integer[][] ingenotypeMat, Integer[][] outgenotypeMat, double missing){
		int nCell = ingenotypeMat[0].length;
		for (int i = 0; i < ingenotypeMat.length; i++){
			outgenotypeMat[i][0] = ingenotypeMat[i][0];
			for (int j = 1; j < nCell; j++){
				double rr = _rng.nextDouble();
				if (rr <= missing)
					outgenotypeMat[i][j] = 3;
				else
					outgenotypeMat[i][j] = ingenotypeMat[i][j];
			}
		}
	}
	
	/**
	 * Print genotype matrix to a file
	 * @param filename
	 * @param genotypeMat
	 * @throws FileNotFoundException
	 * @throws UnsupportedEncodingException
	 */
	public void writeGenotypeMat2File(String filename, Integer[][] genotypeMat) throws FileNotFoundException, UnsupportedEncodingException{
		PrintWriter writer = new PrintWriter(filename, "UTF-8");
		for (int i = 0; i < genotypeMat.length; i++){
			String this_mut_row = "";
			for (int j = 0; j < genotypeMat[0].length; j++){
				this_mut_row = this_mut_row + " " + String.valueOf(genotypeMat[i][j]);
			}
			writer.printf("%s\n", this_mut_row);
		}
		writer.close();
	}
	
	/**
	 * Inside STITree there is a function named selectRandomNode. 
	 * I want a version that consistently produces the same output for same input + seed.
	 */
	public STINode<Double> selectRandomNode(STITree<Double> tree, boolean include_leaves, boolean include_root, Random rng) {
		ArrayList<Integer> idList = new ArrayList<Integer> ();
		
		for (int i = 0; i < tree.getNodeCount(); i++) {
			STINode<Double> curNode = tree.getNode(i);
			if (curNode != null) {
				if (curNode.isLeaf()) {
					if (include_leaves) {
						idList.add(i);
					}
				}
				else if (curNode.isRoot()) {
					if (include_root) {
						idList.add(i);
					}
				}
				else {
					idList.add(i);
				}
			}
			else {
				System.out.println("ERROR: Random node selection: Null Node.");
				return null;
			}
		}
		
		if (idList.size() == 0) {
			System.out.println("ERROR: Random node selection: Empty Set.");
			return null;
		}		
		double randDouble = rng.nextDouble();
		int randIndex = (int) (randDouble * idList.size());
		return tree.getNode(idList.get(randIndex));
	}	
	
	
	/**
	 * Return if the tree is binary or not.
	 */
	public boolean isBinary(String treeString) {
		STITree<Double> tree = getTree(treeString);
		for (int i = 0; i < tree.getNodeCount(); i++) {
			if (!tree.getNode(i).isLeaf() && tree.getNode(i).getChildCount() != 2) {
				return false;
			}
		}
		return true;
	}
	
		
	/**
	 *  Given a Phylogenetic tree in String format, possibly without any branch lengths.
	 *  Create a Tree data structure and assign random branch lengths.
	 *  "maxV" denotes the maximum allowed branch length.
	 */	
	public STITree<Double> getTreeWithRandomBranches (String treeString, double maxV) {
		STITree<Double> tree = getTree(treeString);
		for (int i = 0; i < tree.getNodeCount(); i++) {
			STINode<Double> node = tree.getNode(i);
			
			if (!node.isRoot()) {
				node.setParentDistance(maxV * _rng.nextDouble());
			}
		}
		return tree;
	}
	
	/**
	 *  Given a Phylogenetic tree in String format, possibly without any branch lengths.
	 *  Create a Tree data structure with fixed branch lengths for all.
	 */	
	public STITree<Double> getTreeWithFixedBranches (String treeString, double branchLength) {
		STITree<Double> tree = getTree(treeString);
		for (int i = 0; i < tree.getNodeCount(); i++) {
			STINode<Double> node = tree.getNode(i);
			
			if (!node.isRoot()) {
				node.setParentDistance(branchLength);
			}
		}
		return tree;
	}	
	
	
	
	/**
	 * Given a tree node, selects one of the child randomly.
	 */
	public STINode<Double> getRandomChild(STINode<Double> node) {
		assert node.getChildCount() >= 2;
		ArrayList<STINode<Double>> a = new ArrayList<STINode<Double>>();
		for (STINode<Double> child : node.getChildren()) {
			a.add(child);
		}

		int randIndex = (int) (_rng.nextDouble() * a.size());
		return a.get(randIndex);
	}
	
	
	/**
	 * 	Given a tree node and one of it's child, returns the other child.
	 */
	public STINode<Double> getOtherChild(STINode<Double> node, STINode<Double> c) {
		assert node.getChildCount() <= 3;		
		ArrayList<STINode<Double>> a = new ArrayList<STINode<Double>>();
		
		for (STINode<Double> child : node.getChildren()) {
			if (child.getID() != c.getID()) {
				a.add(child);
			}
		}			
		
		int randIndex = (int) (_rng.nextDouble() * a.size());
		return a.get(randIndex);
	}
	
	/**
	 * 	Given a tree node and one of it's child, return with equal probability the other child or the parent.
	 */
	public STINode<Double> getOtherNeighbour(STINode<Double> node, STINode<Double> c) {
		assert !node.isRoot();
		if (_rng.nextDouble() < 0.5) {
			return node.getParent();
		}
		return getOtherChild(node, c);
	}
	
	/**
	 * Basic Statistical Functions: 
	 */
    private double getMean(double [] data) 
    {
        double sum = 0.0;
        for(double a : data)
            sum += a;
        return (sum/data.length);
    }

    private double getVariance(double [] data)
    {
        double mean = getMean(data);
        double temp = 0;
        for(double a : data){
        	temp += ((mean-a)*(mean-a));
        }
        return (temp/(data.length-1));
    }
    
    /**
     * return the beta dist parameter a from mean and sd
     * @param mean
     * @param sd
     * @return
     */
    public double getBetaPriora(double mean, double sd){
    	return ((1 - mean)*mean*mean/(sd*sd)) - mean;
    }

		/**
	 * returns the geneName list 
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	public ArrayList<String> readGeneNames(String filename) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(filename));
		ArrayList<String> geneNames = new ArrayList<>();
		String line = null;
		while ((line = br.readLine()) != null) {
			geneNames.add(line);
		}
		br.close();
		return geneNames;
	}
    
    
	
}
