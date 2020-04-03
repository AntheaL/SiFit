/**
 * 
 */
package SiFit.algorithm;

import java.io.IOException;
import java.util.ArrayList;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.hj.api.SuspendableException;
import SiFit.BasicUtilityFunctions;
import SiFit.io.VariantMatrixReader;
import SiFit.objects.GenotypeObservation;

/**
 * @author hz22
 *
 */
public class FindSCITETreeLikelihood {

	/**
	 * @param args
	 * @throws IOException 
	 * @throws SuspendableException 
	 */
	public static void main(String[] args) throws IOException, SuspendableException {
		// TODO Auto-generated method stub
		double fn = 0.10529;
		double fp = 0.1;
		double mu = 0.1;
		int df = 0;
		
		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		String SiFitTreeFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/RealDatasets/TreeImages/CO5.newick";
		String SiFitTreeNewick = BUF.readNewickString(SiFitTreeFile);
		STITree<Double> SiFitTree = BUF.getTree(SiFitTreeNewick);
		
		double sumBL = 0;
		for (STINode<Double> node: SiFitTree.getNodes()){
			if (!node.isRoot()){
			sumBL += node.getParentDistance();
			}
		}
//		System.out.println(sumBL);
		int nCell = SiFitTree.getLeafCount();
		double avgBL = sumBL/(2*nCell-1);
//		System.out.println(avgBL);
		
		String SCITETreeFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/CO5_SCITE_tree.txt";
		STITree<Double> SCITETree = BUF.getTree(BUF.readNewickString(SCITETreeFile));
		for (STINode<Double> node: SCITETree.getNodes()){
			if (!node.isRoot()){
			node.setParentDistance(avgBL);
			}
		}
		System.out.println(SCITETree.toNewick());
		
		ArrayList<String> singleCellNames = new ArrayList<>();
		for (int i = 1; i <= nCell; i++){
        	singleCellNames.add(Integer.toString(i));
        }
		
		String varMatFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/RealDatasets/CO5.master_matrix.sifit";
		VariantMatrixReader vr = new VariantMatrixReader(varMatFile);
        vr.populateVarGTMatrix(varMatFile, nCell);
        ArrayList<Integer[]> obsGenotypeMatrix = vr.obsVarGTMatrix;
		ArrayList<GenotypeObservation> obsGTObsMatrix = BUF.getGenotypeObsFromGTMatrix(singleCellNames, obsGenotypeMatrix);
		
		TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fn, mu, BUF, df, 1);
		double SCITEtreeLikelihood = jointSampleHelper.computeTreeLikelihood(obsGTObsMatrix, SCITETree, fn);
		System.out.println(SCITEtreeLikelihood);
	}

}
