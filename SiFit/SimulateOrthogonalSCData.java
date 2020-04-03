/**
 * 
 */
package SiFit;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import SiFit.objects.AmpADOObj;
import SiFit.objects.AmpFPObj;
import SiFit.objects.FormDoublets;
//import monoGenoTree.simulation.EvolutionModelSC;

import org.apache.commons.math3.distribution.NormalDistribution;

import SiFit.objects.NodeDiploidDNA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * @author hz22
 *
 */
public class SimulateOrthogonalSCData {

	/**
	 * @param args
	 */
	public static int nCell;							// Number of cells
	public static int nMut;								// Number of SNV sites
	public static double _fp = 0.01;					// FP rate, alpha
	public static double _fnMean = 0.2;					// FN rate, beta
	public static double _doublet = 0.0;				// Doublet rate
	public static double _missing = 0.0;				// Missing Data fraction
	public static double branchLengthThr = 0.12;
	public static double deletion = 0.05;				// Deletion Rate
	public static double omega = 0.01;				    // LOH Rate
	public static double recurProb = 0.05;				// Recurrent mutation Rate
	public static int dataFlag = 1;
	public static String trueTreeFilename = null;
	public static String trueGenotypeMatFilename = null;
	public static String doubletGenotypeMatFilename = null;
	public static String noisyGenotypeMatFilename = null;
	public static String missingNoisyGenotypeMatFilename = null;
	public static String doublet_file = null;
	public static Random _rng = new Random();
	public static Integer[][] trueGenotypeMatrix;
	public static Integer[][] doubletGenotypeMatrix;
	public static Integer[][] noisyGenotypeMatrix;
	
	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
		
		readArguments(args);
		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		
		
		// The genotype array at the root. All zeros
		Integer[] rootGenotypeArr = BUF.getRootGenotypeArr(nMut);
		// Mutation flag array, entry is 1 for a position that has been mutated
		Integer[] mutFlagArr = BUF.getRootGenotypeArr(nMut);
		
		// Map containing heterozygous and homozygous mutation positions
		HashMap<Integer, Set<Integer>> Genotype_flag_array = new HashMap<Integer, Set<Integer>>(); // Genotype array with info about which pos is heterozygous and which one is homozygous.
		Set<Integer> homoRef = new HashSet<Integer>();
		Set<Integer> mutated = new HashSet<Integer>();
		Set<Integer> homoalt = new HashSet<Integer>();
		for (int i = 0; i < nMut; i++){
			homoRef.add(i);
		}
		// Initialize Genotype_flag_array
		Genotype_flag_array.put(0, homoRef);
		Genotype_flag_array.put(1, mutated);
		Genotype_flag_array.put(2, homoalt);

		// Single cell error parameters
		double fnSD = _fnMean/10;
		HashMap<String, Double> ADO_rateMap = new HashMap<>();
		NormalDistribution ADO_normalDist = new NormalDistribution(_fnMean, fnSD);

		// Simulate true tree, A random binary tree with cells as leaves
		ArrayList<String> scNames = new ArrayList<>();
		for (int i=1; i <= nCell; i++){
			scNames.add("sc" + Integer.toString(i));
		}
		STITree<Double> TreeMaker = BUF.generateRandomTree(scNames, branchLengthThr);
		String Newick_tree = TreeMaker.toNewick();
		STITree<Double> Tree = BUF.getTree(Newick_tree);
		BUF.resizeLeafBranches(scNames, Tree, branchLengthThr/100, branchLengthThr/10); // resize the leaf branches
		BUF.resizeBranchLengths(Tree, 1.1);	// Resize all branches so that the threshold number of sites are mutated
		
		/*
		 *  Process the node info from the true tree
		 */
		int num_node = Tree.getNodeCount();
		Iterable<STINode<Double>> Nodes = Tree.getNodes();
		// Set the name of the root
		STINode<Double> Root = Tree.getRoot();
		Root.setName("R");
		
		// Data Structures for storing the information of all the nodes in the tree
		HashMap<String, Integer> Node_name_ID = new HashMap<String, Integer>(); 					// <name, ID>
		HashMap<Integer, String> Node_ID_name = new HashMap<Integer, String>(); 					// <ID, name>
		HashMap<Integer, Integer[]> nodesGenotypeArrMap = new HashMap<Integer, Integer[]>();   		// <ID, Mut_Flag>

		// Put Root Genotype Array in the nodesGenotypeArrMap
		nodesGenotypeArrMap.put(Root.getID(), rootGenotypeArr);
		
		// Name each internal node
		for (STINode<Double> n : Nodes){
			if (n.getName().equals("")){
				n.setName("IN"+Integer.toString(n.getID()));
			}
			Node_name_ID.put(n.getName(), n.getID());
			Node_ID_name.put(n.getID(), n.getName());
		}
		
		// Store the info regarding all nonLeaf Nodes
		ArrayList<STINode<Double>> Non_Leaf_Nodes = new ArrayList<STINode<Double>>();
		for (int i = 0; i < num_node; i++){
			STINode<Double> node = Tree.getNode(i);
			if (node.isLeaf() == false){
				Non_Leaf_Nodes.add(node);
			}
		}
		
		// Write tree info to a file
		PrintWriter tree_writer = new PrintWriter(trueTreeFilename, "UTF-8");
		tree_writer.printf("%s\n", Tree.toNewick());
		tree_writer.close();
		
		// Obtain Mutated DNA for each node
		int totalMutCount = 0;
		for (STINode<Double> p : Non_Leaf_Nodes){
			Iterable<STINode<Double>> p_Children = p.getChildren();
			for (STINode<Double> child: p_Children){
				Integer[] childGenotypeArr;
				if (dataFlag == 1){
				// This is for ternary data
				childGenotypeArr = BUF.getMutatedGenotypeArrFSM(child, nodesGenotypeArrMap.get(p.getID()), Genotype_flag_array, mutFlagArr, recurProb, omega, deletion, child.getParentDistance(), nMut);
				}
				else{
				// This is for checking our model of evolution
//				Integer[] childGenotypeArr = BUF.getMutatedGenotypeArrFSMTest(child, nodesGenotypeArrMap.get(p.getID()), Genotype_flag_array, mutFlagArr, recurProb, omega, deletion, child.getParentDistance(), nMut);
				// This is for binary data
				childGenotypeArr = BUF.getMutatedGenotypeArrFSMBinary(child, nodesGenotypeArrMap.get(p.getID()), Genotype_flag_array, mutFlagArr, recurProb, omega, deletion, child.getParentDistance(), nMut);
				}
				if (nodesGenotypeArrMap.containsKey(child.getID()) == false){
					nodesGenotypeArrMap.put(child.getID(), childGenotypeArr);
				}
				totalMutCount += BUF.countMutParent2Child(nodesGenotypeArrMap.get(p.getID()), childGenotypeArr);
			}
		}

		ArrayList<Integer> allMutatedPosList = new ArrayList<>();
		allMutatedPosList.addAll(Genotype_flag_array.get(1));
		Collections.sort(allMutatedPosList);
		trueGenotypeMatrix = new Integer[allMutatedPosList.size()][nCell+1];
		BUF.writeGenotypeMat2File(trueGenotypeMatFilename, allMutatedPosList, scNames, trueGenotypeMatrix, Node_name_ID, nodesGenotypeArrMap);

		// Introduce Doublets
		HashMap<Integer, Integer[]> afterDoubletsCellsGenotypeArrMap = new HashMap<>();
		if (_doublet > 0.0){			
			StringBuilder doubletFlagSb = new StringBuilder(nCell);
			for (String sc_name : scNames){
				int scID = Node_name_ID.get(sc_name);
				double rr = _rng.nextDouble();
				if (rr <= _doublet){
					doubletFlagSb.append("1");
					
					// Randomly choose another cell to form the doublet
					int rand_cell_index = _rng.nextInt(nCell);
					if (scNames.get(rand_cell_index).equals(sc_name))
						rand_cell_index = _rng.nextInt(nCell);
					Integer[] doubletMateGTarr = nodesGenotypeArrMap.get(Node_name_ID.get(scNames.get(rand_cell_index)));
					Integer[] afterDoubletCellGTarr = FormDoublets.getDoubletGTarr(nodesGenotypeArrMap.get(scID), doubletMateGTarr);
					afterDoubletsCellsGenotypeArrMap.put(scID, afterDoubletCellGTarr);
					
				}
				else{
					afterDoubletsCellsGenotypeArrMap.put(scID, nodesGenotypeArrMap.get(scID));
					doubletFlagSb.append("0");
				}
			}

			String doubletFlagString = doubletFlagSb.toString();
			// Write doublet flags info to a file
			if (doublet_file != null){
				PrintWriter doublet_writer = new PrintWriter(doublet_file, "UTF-8");
				doublet_writer.printf("%s\n", doubletFlagString);
				doublet_writer.close();
			}
			// Write double genotype matrix to a file
			if (doubletGenotypeMatFilename != null){
				doubletGenotypeMatrix = new Integer[allMutatedPosList.size()][nCell+1];
				BUF.writeGenotypeMat2File(doubletGenotypeMatFilename, allMutatedPosList, scNames, doubletGenotypeMatrix, Node_name_ID, afterDoubletsCellsGenotypeArrMap);
			}
		}
		
//		 Apply ADO to get ADO affected Genotypes
		HashMap<String, AmpADOObj> ADOAffectedCellObjMap = new HashMap<>();
		HashMap<Integer, Integer[]> afterADOCellsGenotypeArrMap = new HashMap<>();
		int totalFNs = 0;
		for (String sc_name : scNames){
			int scID = Node_name_ID.get(sc_name);
			double scADO_rate = ADO_normalDist.sample();
			ADO_rateMap.put(sc_name, scADO_rate);
			AmpADOObj ADOAffectedCellObj;
			if (_doublet > 0.0){
				ADOAffectedCellObj = new AmpADOObj(afterDoubletsCellsGenotypeArrMap.get(scID), sc_name);
				ADOAffectedCellObj.getADOAffectedGenotypeArr(scADO_rate, afterDoubletsCellsGenotypeArrMap.get(scID));
			}
			else{
				ADOAffectedCellObj = new AmpADOObj(nodesGenotypeArrMap.get(scID), sc_name);			
				ADOAffectedCellObj.getADOAffectedGenotypeArr(scADO_rate, nodesGenotypeArrMap.get(scID));
			}
			
			int thisCellFN = BUF.countMutParent2Child(nodesGenotypeArrMap.get(scID),ADOAffectedCellObj.afterADOGenotypeArr);

			totalFNs += thisCellFN;
			afterADOCellsGenotypeArrMap.put(scID, ADOAffectedCellObj.afterADOGenotypeArr);
			ADOAffectedCellObjMap.put(sc_name, ADOAffectedCellObj);
		}

		
		// Apply FP to get FP affected Genotypes
		HashMap<String, AmpFPObj> FPAffectedCellObjMap = new HashMap<>();
		HashMap<Integer, Integer[]> afterFPCellsGenotypeArrMap = new HashMap<>();
		
		int FP_count = 0;
		for (String sc_name : scNames){
			int scID = Node_name_ID.get(sc_name);
			AmpFPObj FPAffectedCellObj = new AmpFPObj(afterADOCellsGenotypeArrMap.get(scID), sc_name);
			FPAffectedCellObj.getFPAffectedGenotypeArr(_fp, afterADOCellsGenotypeArrMap.get(scID), Genotype_flag_array);
			afterFPCellsGenotypeArrMap.put(scID, FPAffectedCellObj.afterFPGenotypeArr);

			FP_count += FPAffectedCellObj.FPAffectedPosList.size();
			
		}
		
		// Data Structures for writing noisy genotype matrix to a file
		Set<Integer> allMutatedPosPlusNoise = new HashSet<>();
		allMutatedPosPlusNoise.addAll(Genotype_flag_array.get(1));
		allMutatedPosPlusNoise.addAll(Genotype_flag_array.get(2));
		ArrayList<Integer> allMutatedPosPlusNoiseList = new ArrayList<>();
		allMutatedPosPlusNoiseList.addAll(allMutatedPosPlusNoise);
		Collections.sort(allMutatedPosPlusNoiseList);
		
		// Write noisy genotype matrix to a file
		noisyGenotypeMatrix = new Integer[allMutatedPosPlusNoiseList.size()][nCell+1];
		BUF.writeGenotypeMat2File(noisyGenotypeMatFilename, allMutatedPosPlusNoiseList, scNames, noisyGenotypeMatrix, Node_name_ID, afterFPCellsGenotypeArrMap);
		
		// Introduce Missing Data
		Integer[][] missingNoisyGenotypeMatrix = new Integer[allMutatedPosPlusNoiseList.size()][nCell+1];
		if (_missing > 0.0){
			BUF.addMissingData2GenotypeMat(noisyGenotypeMatrix, missingNoisyGenotypeMatrix, _missing);
			BUF.writeGenotypeMat2File(missingNoisyGenotypeMatFilename, missingNoisyGenotypeMatrix);
		}
	}
	
	/**
	 * Read the arguments from command line
	 * @param args
	 */
	public static void readArguments(String[] args){
		int nPar = args.length;
		for (int i = 0; i < nPar; i = i+2){
			if (args[i].equals("-m") == true){
				nCell = Integer.parseInt(args[i+1]); 
			}
			else if (args[i].equals("-n") == true){
				nMut = Integer.parseInt(args[i+1]); 
			}
			else if (args[i].equals("-fp") == true){
				_fp = Double.parseDouble(args[i+1]);
				if (_fp < 0 || _fp > 1)
					throw new IllegalArgumentException("Invalid input for fp rate, should be between 0 and 1");
			}
			else if (args[i].equals("-fn") == true){
				_fnMean = Double.parseDouble(args[i+1]);
				if (_fnMean < 0 || _fnMean > 1)
					throw new IllegalArgumentException("Invalid input for fn rate, should be between 0 and 1");
			}
			else if (args[i].equals("-w") == true){
				omega = Double.parseDouble(args[i+1]);
				if (omega < 0 || omega > 1)
					throw new IllegalArgumentException("Invalid input for LOH rate, should be between 0 and 1");
			}
			else if (args[i].equals("-d") == true){
				deletion = Double.parseDouble(args[i+1]);
				if (deletion < 0 || deletion > 1)
					throw new IllegalArgumentException("Invalid input for Deletion rate, should be between 0 and 1");
			}
			else if (args[i].equals("-r") == true){
				recurProb = Double.parseDouble(args[i+1]);
				if (recurProb < 0 || recurProb > 1)
					throw new IllegalArgumentException("Invalid input for recurrent mutation rate, should be between 0 and 1");
			}
			else if (args[i].equals("-df") == true){
				dataFlag = Integer.parseInt(args[i+1]); 
				if (dataFlag > 1)
					throw new IllegalArgumentException("Invalid input for dataFlag parameter, should be 0 or 1");
			}
			else if (args[i].equals("-trueTree") == true){
				trueTreeFilename = args[i+1];
			}
			else if (args[i].equals("-trueGenotype") == true){
				trueGenotypeMatFilename = args[i+1];
			}
			else if (args[i].equals("-noisyGenotype") == true){
				noisyGenotypeMatFilename = args[i+1];
			}
			else if (args[i].equals("-doublet") == true){
				_doublet = Double.parseDouble(args[i+1]);
			}
			else if (args[i].equals("-doubletGenotype") == true){
				doubletGenotypeMatFilename = args[i+1];
			}
			else if (args[i].equals("-doubletFlag") == true){
				doublet_file = args[i+1];
			}
			else if (args[i].equals("-missing") == true){
				_missing = Double.parseDouble(args[i+1]);
			}
			else if (args[i].equals("-missingGenotype") == true){
				missingNoisyGenotypeMatFilename = args[i+1];
			}
			
		}
	}

}
