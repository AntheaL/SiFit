/**
 * 
 */
package SiFit.algorithm;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;

import org.apache.commons.lang3.StringUtils;

import SiFit.BasicUtilityFunctions;
import SiFit.io.VariantMatrixReader;
import SiFit.model.ComplexEvolutionModel;
import SiFit.model.JCModelSingleCell;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * @author hz22
 *
 */
public class InferAncestralStates {
	
	protected static class AncestralInfo{
		
		// fields
		/**
		 * The sequence at this node.  Generated for internal nodes, set in leaves.
		 */
		public Integer[] sequence;
		
		/**
		 * The node in the input tree that this node corresponds to.
		 */
		public TMutableNode node;
		
		/**
		 * The list of best genotypes for each mutation for the node.
		 * These will be filled during bottom up phase of the algorithm
		 * and used during the top down phase of the algorithm.
		 */
		public ArrayList<ArrayList<Integer>> setGTNode = new ArrayList<>();
		
		/**
		 * The list of likelihoods for each genotype for each mutation for the node.
		 * These will be filled during bottom up phase of the algorithm
		 * and used during the top down phase of the algorithm.
		 */
		public ArrayList<ArrayList<Double>> likelihoodGTNode = new ArrayList<>();
		
		public AncestralInfo(TMutableNode node, int len) {
			this.node = node;
			for (int i = 0; i < len; i++){
				setGTNode.add(new ArrayList<Integer>());
				likelihoodGTNode.add(new ArrayList<Double>());
			}
			
		}
	}
	
	Tree tree;
	JCModelSingleCell model;
	ComplexEvolutionModel modelC; 
	private static double fp = 1e-5;;
    private static double fn = 0.2;
    public static double LOHProb = 0.01;
    public static double delProb = 0.05;
    public static double recurProb = 0.05;
    private static double[][] rates;
    private static double mutation_rate = 0.1;
    public static String genotypes;
    public static int dataFlag = 1;
	public static String varFileName = null;			// Filename of the input matrix
	public static String treefile = null;		// File containing the true tree in newick form
	public static String cellfile = null;		// File containing the cell names
	public static String genefile = null;
	public static String expectedMatrixFile = null;

    
    ArrayList<String> geneNamesList;			// List of genes on which mutations occur
    public STITree<AncestralInfo> nuPtree;		// Copy of the original tree but with branches that correspond to number of mutations
    
    // Constructor
    public InferAncestralStates(double fpositive, double fnegative, double mut_rate, int df, Tree geneTree, JCModelSingleCell model){
    	fp = fpositive; fn = fnegative; mutation_rate = mut_rate;
        this.tree = geneTree;
        this.model = model;
        dataFlag = df;
        
        // Binary Data
        if (dataFlag == 0){
        	rates = new double[][]
            		{{(1-fp), fp},
            		 {fn, (1-fn)}};
        	genotypes = "01";
        }
        // Ternary Data
        else {
            rates = new double[][]
                    {{(1 - fp*fn/2.0 - fp),     fp,         fp*fn/2.0},
                     {fn/2.0,                   1 - fn,     fn/2.0},
                     {0.0,                      0.0,        1.0}};
            genotypes = "012";
        }
        
    }
    
    // Constructor 
    public InferAncestralStates(double fpositive, double fnegative, 
    							double delProb, double LOHProb, double recurProb,
    							int df){
    	fp = fpositive; fn = fnegative; 
//        this.tree = geneTree;
        this.modelC = new ComplexEvolutionModel(delProb, LOHProb, recurProb);
        dataFlag = df;
        
        // Binary Data
        if (dataFlag == 0){
        	rates = new double[][]
            		{{(1-fp), fp},
            		 {fn, (1-fn)}};
        	genotypes = "01";
        }
        // Ternary Data
        else {
            rates = new double[][]
                    {{(1 - fp*fn/2.0 - fp),     fp,         fp*fn/2.0},
                     {fn/2.0,                   1 - fn,     fn/2.0},
                     {0.0,                      0.0,        1.0}};
            genotypes = "012";
        }
    }
    
    /**
     * Infer the ancestral sequences and returns the number of mutations in the tree
     * @param tree - input tree
     * @param cellGTSeqMap - 
     * @param seqLen - Number of mutations
     * @return
     */
    protected int computeAncestralStates(MutableTree tree, 
			Hashtable<String, Integer[]> cellGTSeqMap,
			int seqLen){
    	STITree<AncestralInfo> ptree = new STITree<>();
    	
    	// copy the tree
    	copyNode(ptree.getRoot(), tree.getRoot(), cellGTSeqMap, seqLen);
    	
    	for (STINode<AncestralInfo> n : ptree.getNodes()){
    		n.setName(n.getData().node.getName());
    	}

    	// Perform bottom up and top down phases
    	for(int pos = 0; pos < seqLen; pos++) {
    		performBottomUp(pos, ptree.getRoot());
    		performTopDown(pos, ptree.getRoot());
//    		performTopDown(pos, ptree.getRoot(), ptree);
    	}
    	
    	this.nuPtree = ptree;
    	
    	return computePScore(ptree.getRoot());
    }
    
    protected void performBottomUp(int pos, STINode<AncestralInfo> pnode){
    	AncestralInfo pi = pnode.getData();
    	
    	// Leaf
    	if (pnode.isLeaf()){
    		// Get the transition matrix with branch length to the parent of leaf
    		double[][] JCTransMat;
    		
    		// This is for complex model
    		JCTransMat = modelC.getTransitionMatrixBinary(pi.node.getParentDistance());
    		// This is for JC model
//    		if (dataFlag == 0)
//    			JCTransMat = model.jukesCantorBinaryMatrixAll(pi.node.getParentDistance());
//    		else
//    			JCTransMat = model.jukesCantor012AllMatrix(pi.node.getParentDistance());
    		
    		// Missing Data
    		if (pi.sequence[pos] == 3){
    			for (int i = 0; i < genotypes.length(); i++){
    				double bestL = Double.NEGATIVE_INFINITY;
    				int bestGT = 0;
    				for (int j = 0; j < genotypes.length(); j++){
//    					double ll_i = JCTransMat[i][j];
    					double ll_i = Math.log(JCTransMat[i][j]);
    					if (ll_i > bestL){
    						bestL = ll_i;
    						bestGT = j;
    					}
    				}
    				pi.likelihoodGTNode.get(pos).add(bestL);
    				pi.setGTNode.get(pos).add(bestGT);
    			}
    		}
    		// Invalid Genotype
    		else if (genotypes.indexOf(Integer.toString(pi.sequence[pos])) == -1){
    			throw new RuntimeException(pi.sequence[pos] +" is not a valid genotype");
    		}
    		// Valid Genotype
    		else{
    			for (int i = 0; i < genotypes.length(); i++){
    				double bestL = Double.NEGATIVE_INFINITY;
    				int bestGT = 0;
    				for (int j = 0; j < genotypes.length(); j++){
    					// Account for both evolution and error rates
    					double ll_i = Math.log(JCTransMat[i][j]) + Math.log(rates[j][pi.sequence[pos]]);
    					if (ll_i > bestL){
    						bestL = ll_i;
    						bestGT = j;
    					}
    				}
    				pi.likelihoodGTNode.get(pos).add(bestL);
    				pi.setGTNode.get(pos).add(bestGT);
    			}
    		}
    	}
    	else if (pnode.isRoot()){
    		pi.setGTNode.get(pos).add(0);
    		for (STINode<AncestralInfo> child : pnode.getChildren()) {
    			performBottomUp(pos, child);
    		}
    	}
    	else{
    		ArrayList<STINode<AncestralInfo>> p_children = new ArrayList<>();
    		for(STINode<AncestralInfo> child : pnode.getChildren()) {				
				performBottomUp(pos, child);
				p_children.add(child);
			}
    		double[][] JCTransMat;
    		JCTransMat = modelC.getTransitionMatrixBinary(pi.node.getParentDistance());
//    		if (dataFlag == 0)
//    			JCTransMat = model.jukesCantorBinaryMatrixAll(pi.node.getParentDistance());
//    		else
//    			JCTransMat = model.jukesCantor012AllMatrix(pi.node.getParentDistance());

    		for (int i = 0; i < genotypes.length(); i++){
    			double bestL = Double.NEGATIVE_INFINITY;
				int bestGT = 0;
				double ll_i;
				for (int j = 0; j < genotypes.length(); j++){
//					ll_i = JCTransMat[i][j] * p_children.get(0).getData().likelihoodGTNode.get(pos).get(j) * p_children.get(1).getData().likelihoodGTNode.get(pos).get(j);
					ll_i = Math.log(JCTransMat[i][j]) + p_children.get(0).getData().likelihoodGTNode.get(pos).get(j) + p_children.get(1).getData().likelihoodGTNode.get(pos).get(j);
					if (ll_i > bestL){
						bestL = ll_i;
						bestGT = j;
					}
				}
				pi.likelihoodGTNode.get(pos).add(bestL);
				pi.setGTNode.get(pos).add(bestGT);				
    		}
    	}
    }
    
    protected void performTopDown(int pos, STINode<AncestralInfo> pnode){
    	AncestralInfo pi = pnode.getData();
    	if (pnode.isRoot()){
    		pi.sequence[pos] = 0;
    	}
    	else{
    		int index = pnode.getParent().getData().sequence[pos];
    		pi.sequence[pos] = pi.setGTNode.get(pos).get(index);
    	}
    	for (STINode<AncestralInfo> child : pnode.getChildren()) {
			performTopDown(pos, child);
		}
    }
    
    protected void performTopDown(int pos, STINode<AncestralInfo> pnode, STITree<AncestralInfo> ptree){
    	AncestralInfo pi = pnode.getData();
    	if (pnode.isRoot()){
    		pi.sequence[pos] = 0;
    	}
    	else{
    		int index = pnode.getParent().getData().sequence[pos];
    		if (checkTripleMutation(pos, pnode, ptree) > 1)
    			pi.sequence[pos] = pnode.getParent().getData().sequence[pos];
    		else
    			pi.sequence[pos] = pi.setGTNode.get(pos).get(index);
    	}
    	for (STINode<AncestralInfo> child : pnode.getChildren()) {
			performTopDown(pos, child, ptree);
		}
    }
    
    protected int checkTripleMutation(int pos, STINode<AncestralInfo> pnode, STITree<AncestralInfo> ptree){
    	int countMut = 0;
    	STINode<AncestralInfo> par = pnode.getParent();
    	while (par.equals(ptree.getRoot()) ==  false){
    		STINode<AncestralInfo> grandPar = par.getParent();
    		if (grandPar.getData().sequence[pos] != par.getData().sequence[pos])
    			countMut++;
    		par = grandPar;
    	}
    	return countMut;
    }

    /**
     * Compute the number of mutations in a tree
     * @param pnode
     * @return
     */
	protected int computePScore(STINode<AncestralInfo> pnode) {
		
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
		for(STINode<AncestralInfo> child : pnode.getChildren()) {
			pscore += computePScore(child) + child.getParentDistance();
		}
		
		return pscore;		
	}
	
	/**
	 * Copy the input tree
	 * @param pnode
	 * @param n
	 * @param cellGTSeqMap
	 * @param seqLen
	 */
	protected void copyNode(STINode<AncestralInfo> pnode, TMutableNode n, Hashtable<String, Integer[]> cellGTSeqMap, int seqLen){
		AncestralInfo pi = new AncestralInfo(n, seqLen);
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
				STINode<AncestralInfo> pchild = pnode.createChild();
				copyNode(pchild, c, cellGTSeqMap, seqLen);
			}
		}
	}
	
	public void annotateMutations(String treeFile, String varMatFile, String singleCellFile, String geneNamesFile, String expectedMatrixFile) throws IOException{
		
		// Classes to use
		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		
		// Get the input tree
		String intreeNewick = BUF.readNewickString(treeFile);
		STITree<Double> intree = BUF.getTree(intreeNewick);
		
		// Name the internal nodes of tree
		intree.getRoot().setName("R");
		int i = 1;
		for (STINode<Double> node: intree.getNodes()){			
			if (node.getName() == ""){
				node.setName("in"+ Integer.toString(i));
				i++;
			}
		}

//		System.out.println(intree.toNewick());
		// Create the single cell genomes to use as input to parsimony score calculator
		int nCell = intree.getLeafCount();
		VariantMatrixReader vr = new VariantMatrixReader(varMatFile);
		vr.populateVarGTMatrix(varMatFile, nCell);
		ArrayList<Integer[]> obsGenotypeMatrix = vr.obsVarGTMatrix;

		int nMut = obsGenotypeMatrix.size();
		ArrayList<Integer[]> seqGT = BUF.getCellGenomes(obsGenotypeMatrix, nMut);
		
		// Get the single cell names
		String singleCellNameString = BUF.readNewickString(singleCellFile);
		String[] cellNameArr = singleCellNameString.split(" ");
		
		Hashtable<String, Integer[]> cellGTSeqMap = new Hashtable<>();
		for (int j = 0; j < cellNameArr.length; j++){
			cellGTSeqMap.put(cellNameArr[j], seqGT.get(j));
		}
		
		// Get the Gene Names
		ArrayList<String> geneList = BUF.readGeneNames(geneNamesFile);
		this.geneNamesList = geneList;
		this.tree = intree;
		int score = this.computeAncestralStates(intree, cellGTSeqMap, nMut);

		
		for (STINode<AncestralInfo> node : this.nuPtree.getNodes()){
			if (node.isRoot() == false){

				if (node.isLeaf())
					continue;
				ArrayList<String> mutList = this.getMutationsP2C(node.getParent(), node, geneList, nMut);
				if (mutList.size()>0){
//					this.printMutationsP2C(node.getParent().getName(), node.getName(), geneList, nMut);
					System.out.println(node.getData().node.getName());
					System.out.printf("Par = %s, child = %s\n", node.getParent().getData().node.getName(), node.getData().node.getName());
					for (String mut : mutList)					
						System.out.println(mut);
					System.out.printf("child, %s acquired %d mutations\n", node.getData().node.getName(), mutList.size());
				}
			}
		}
		HashMap<String, ArrayList<String>> cellGenomeMap = new HashMap<>();
		for (STINode<AncestralInfo> node : this.nuPtree.getNodes()){
			if (node.isLeaf()){
				ArrayList<String> cellGenome = new ArrayList<>();
				cellGenome.add(node.getData().node.getName());
				for (Integer g : node.getParent().getData().sequence)
					cellGenome.add(Integer.toString(g));
				cellGenomeMap.put(node.getData().node.getName(), cellGenome);
			}
		}
		try{
			PrintWriter writer = new PrintWriter(expectedMatrixFile, "UTF-8");
			for (String s : cellGenomeMap.keySet()){
				String cell = StringUtils.join(cellGenomeMap.get(s), '\t');
				writer.println(cell);
			}
			
			writer.close();
		} catch(IOException e){
			e.printStackTrace();
		};
		
	}
	
	public ArrayList<String> getMutationsP2C(STINode<AncestralInfo> par, 
			 STINode<AncestralInfo> child,
			 ArrayList<String> geneNames,
			 int nMut){
		ArrayList<String> mutatedGenes = new ArrayList<>();
		for (int i = 0; i < nMut; i++){
			if (par.getData().sequence[i] != child.getData().sequence[i]){
				mutatedGenes.add(geneNames.get(i));
			}
		}
		return mutatedGenes;
	}

	public void printMutationsP2C(String parname, String childname, ArrayList<String> geneNames, int nMut, String mutPosFile) throws IOException{
		ArrayList<String> mutPosList = new BasicUtilityFunctions().readGeneNames(mutPosFile);
		//ArrayList<String> synoList = new BasicUtilityFunctions().readGeneNames(synoFile);
		STINode<AncestralInfo> par = this.nuPtree.getNode(parname);
		STINode<AncestralInfo> child = this.nuPtree.getNode(childname);
		System.out.printf("%s\t%s%n", parname, childname);
		for (int i = 0; i < nMut; i++){
			if (par.getData().sequence[i] != child.getData().sequence[i]){
				System.out.printf("%s\t%s\t%d\t%d\n", mutPosList.get(i), geneNames.get(i), par.getData().sequence[i], child.getData().sequence[i]);
			}
		}
	}

	public void printMutationsP2C(String parname, String childname, ArrayList<String> geneNames, int nMut) throws IOException{
		//ArrayList<String> mutPosList = new BasicUtilityFunctions().readGeneNames(mutPosFile);
		//ArrayList<String> synoList = new BasicUtilityFunctions().readGeneNames(synoFile);
		STINode<AncestralInfo> par = this.nuPtree.getNode(parname);
		STINode<AncestralInfo> child = this.nuPtree.getNode(childname);
		System.out.printf("%s\t%s%n", parname, childname);
		for (int i = 0; i < nMut; i++){
			if (par.getData().sequence[i] != child.getData().sequence[i]){
				System.out.printf("%s\t%d\t%d\n", geneNames.get(i), par.getData().sequence[i], child.getData().sequence[i]);
			}
		}
	}
		
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		readArguments(args);
		InferAncestralStates AI = new InferAncestralStates(fp, fn, delProb, LOHProb, recurProb, dataFlag);
		AI.annotateMutations(treefile, varFileName, cellfile, genefile, expectedMatrixFile);
	}
	
	
	public static void readArguments(String[] args){
		int nPar = args.length;
		for (int i = 0; i < nPar; i = i+2){
			if (args[i].equals("-fp") == true){
				fp = Double.parseDouble(args[i+1]);
				if (fp < 0 || fp > 1)
					throw new IllegalArgumentException("Invalid input for fp rate, should be between 0 and 1");
			}
			else if (args[i].equals("-fn") == true){
				fn = Double.parseDouble(args[i+1]);
				if (fn < 0 || fn > 1)
					throw new IllegalArgumentException("Invalid input for fn rate, should be between 0 and 1");
			}
			else if (args[i].equals("-w") == true){
				LOHProb = Double.parseDouble(args[i+1]);
				if (LOHProb < 0 || LOHProb > 1)
					throw new IllegalArgumentException("Invalid input for LOH rate, should be between 0 and 1");
			}
			else if (args[i].equals("-d") == true){
				delProb = Double.parseDouble(args[i+1]);
				if (delProb < 0 || delProb > 1)
					throw new IllegalArgumentException("Invalid input for delProb rate, should be between 0 and 1");
			}
			else if (args[i].equals("-df") == true){
				dataFlag = Integer.parseInt(args[i+1]); 
				if (dataFlag > 1)
					throw new IllegalArgumentException("Invalid input for dataFlag parameter, should be 0 or 1");
			}
			else if (args[i].equals("-ipMat") == true){
				varFileName = args[i+1];
			}
			else if (args[i].equals("-tree") == true){
				treefile = args[i+1];
			}
			else if (args[i].equals("-cellNames") == true){
				cellfile = args[i+1];
			}
			else if (args[i].equals("-geneNames") == true){
				genefile = args[i+1];
			}
			else if (args[i].equals("-expectedMatrix") == true){
				expectedMatrixFile = args[i+1];
			}
		}
	}

}
