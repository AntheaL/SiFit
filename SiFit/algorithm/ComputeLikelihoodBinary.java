/**
 * 
 */
package SiFit.algorithm;

import java.util.ArrayList;
import java.util.HashMap;

import SiFit.objects.GenotypeObservation;
import jeigen.DenseMatrix;
import SiFit.BasicUtilityFunctions;
import SiFit.model.ComplexEvolutionModel;
import SiFit.model.JCModelSingleCell;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * @author hz22
 *
 */
public class ComputeLikelihoodBinary {
	
	/**
     *  The parameters that remain static throughout the algorithm:
     *  fp - the false positive error rate
     *  fn - the false negative error rate
     *  rates - the error rates; the i,j entry represents the probability of observing state i as state j
     *  mutation_rate - the approximated/expected mutation rate over the entire tree (needed for Jukes-Cantor)
     *  genotypes - the possible genotypes for single cell data
     */
	public Tree tree;
	public JCModelSingleCell model;
	public ComplexEvolutionModel modelC;
	private static double fp;
    private static double fn;
    private static double mu;
    private static double[][] rates;
    public static final String genotypes = "01";
    
    public ComputeLikelihoodBinary(double fpositive, double fnegative, double mut_rate, Tree t, JCModelSingleCell model){
    	fp = fpositive; fn = fnegative; mu = mut_rate;
        rates = new double[][]
        		{{(1-fp), fp},
        		 {fn, (1-fn)}};
        this.tree = t;
        this.model = model;
    }
    
    public ComputeLikelihoodBinary(double fpositive, double fnegative, double mut_rate, Tree t, ComplexEvolutionModel model){
    	fp = fpositive; fn = fnegative; mu = mut_rate;
    	rates = new double[][]
        		{{(1-fp), fp},
        		 {fn, (1-fn)}};
        this.tree = t;
        this.modelC = model;
    }
    
    private static DenseMatrix createLeafMatrix(Integer genotype){
    	// Missing genotype
    	if (genotype == 3){
            double [][] t = {{1},{1}};
            DenseMatrix D = new DenseMatrix(t);
            return D;
    	}
    	double [][] temp = {
                {0},
                {0}                           
            };
    	if (genotypes.indexOf(Integer.toString(genotype)) == -1)
    		throw new RuntimeException(genotype +" is not a valid genotype");
    	for (int i = 0; i <= 1; i++){
    		temp[i][0] = rates[i][genotype];
    	}
    	return new DenseMatrix(temp); 
    }
    
    private DenseMatrix probabilityOfNode(TNode node, GenotypeObservation obs_genotypes, int modelFlag){
    	if (node.isLeaf()){
    		return createLeafMatrix(obs_genotypes.cellGenotypeMap.get(node.getName()));
    	}
    	else{
    		double[][] probs = new double[][]{{1.0}, {1.0}};
    		for (TNode child : node.getChildren()) {
                double dist = child.getParentDistance();
                if (Double.isNaN(dist) || Double.isInfinite(dist))
                    throw new IllegalArgumentException("Node : " + child.getName() + " has non-finite dist: "+dist);
                double[] edgeP = new double[2];
                DenseMatrix childProbs = probabilityOfNode(child, obs_genotypes, modelFlag);
                double[][] JCTransMat = new double[2][2];
                
                // Now using just the complex model of evolution
                JCTransMat = modelC.getTransitionMatrixBinary(dist);

                for (int i = 0; i < 2; i++){
                	edgeP[0] += JCTransMat[0][i] * childProbs.get(i, 0);
                	edgeP[1] += JCTransMat[1][i] * childProbs.get(i, 0);
                }
                probs[0][0] *= edgeP[0]; probs[1][0] *= edgeP[1];
    		}
            DenseMatrix probMat = new DenseMatrix(probs);
            return probMat;
    	}
    }
    
    /**
     * Function that returns the likelihood of genotype observation at given locus
     * @param arg0 - A map of genotypes of the single cells
     * @return the likelihood
     */
	public double getGenotypeProbability(GenotypeObservation arg0, int modelFlag) {
		return probabilityOfNode(tree.getRoot(), arg0, modelFlag).mul(modelC.getEquilibriumVectorBinary()).sum().s(); // returns densematrix.get(0,0);
	}

	public double getGenotypeProbability(GenotypeObservation arg0,
			int modelFlag, HashMap<String, Integer> _doubletFlagMap) {
		return probabilityOfNode(tree.getRoot(), arg0, modelFlag, _doubletFlagMap).mul(modelC.getEquilibriumVectorBinary()).sum().s();
	}

    private DenseMatrix probabilityOfNode(TNode node, GenotypeObservation obs_genotypes, int modelFlag, HashMap<String, Integer> _doubletFlagMap){
    	if (node.isLeaf()){
    		if (_doubletFlagMap.get(node.getName()) == 0)
    			return createLeafMatrix(obs_genotypes.cellGenotypeMap.get(node.getName()));
    		else
    			return createLeafMatrix(obs_genotypes.cellGenotypeMap.get(node.getName()));
    	}
    	else{
    		double[][] probs = new double[][]{{1.0}, {1.0}};
    		for (TNode child : node.getChildren()) {
                double dist = child.getParentDistance();
                if (Double.isNaN(dist) || Double.isInfinite(dist))
                    throw new IllegalArgumentException("Node : " + child.getName() + " has non-finite dist: "+dist);
                double[] edgeP = new double[2];
                DenseMatrix childProbs = probabilityOfNode(child, obs_genotypes, modelFlag);
                double[][] JCTransMat = new double[2][2];
                if (modelFlag == 0){
                	JCTransMat = model.jukesCantorBinaryMatrix01(dist);
                }
                else{
                	JCTransMat = model.jukesCantorBinaryMatrixAll(dist);
                }
                for (int i = 0; i < 2; i++){
                	edgeP[0] += JCTransMat[0][i] * childProbs.get(i, 0);
                	edgeP[1] += JCTransMat[1][i] * childProbs.get(i, 0);
                }
                probs[0][0] *= edgeP[0]; probs[1][0] *= edgeP[1];
    		}
            DenseMatrix probMat = new DenseMatrix(probs);
            return probMat;
    	}
    }
    
    public static void main(String[] args){
    	String t = "((c:0.1,b:0.1):0.1,a);";
    	BasicUtilityFunctions BUF = new BasicUtilityFunctions();
    	STITree<Double> tree = BUF.getTree(t);
    	tree.setRooted(true);
//    	tree.getRoot().adoptChild(tree.getNode("b"));
    	JCModelSingleCell model = new JCModelSingleCell(0.1);
    	ComputeLikelihoodBinary CL = new ComputeLikelihoodBinary(0.01, 0.2, 0.1, tree, model);
    	
    	ArrayList<String> cells = new ArrayList<>();
    	cells.add("a");
    	cells.add("b");
    	cells.add("c");
//    	for (String s: tree.getLeaves()){
//    		cells.add(s);
//    	}
    	System.out.println(cells);
    	Integer[] genotypes = new Integer[]{1,0,1};
    	GenotypeObservation gt = new GenotypeObservation(cells, genotypes);
    	double val = CL.getGenotypeProbability(gt, 0);
    	System.out.println(val);
    }

}
