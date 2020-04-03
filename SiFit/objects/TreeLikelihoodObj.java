/**
 * 
 */
package SiFit.objects;

import java.io.IOException;
import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Random;

import SiFit.BasicUtilityFunctions;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * Class for storing the combination of tree, 
 * error rate and the corresponding likelihood
 * @author hz22
 *
 */
public class TreeLikelihoodObj {
	
	public double bestBeta;
	public String bestTreeNewickString;
	public double bestTreeLikelihood;
	public double bestDelProb;
	public double bestLOHProb;
	
	/**
	 * Constructor that populates the three fields
	 * @param beta
	 * @param treeN
	 * @param log_likelihood
	 */
	public TreeLikelihoodObj(double beta, String treeN, double log_likelihood){
		this.bestBeta = beta;
		this.bestTreeNewickString = treeN;
		this.bestTreeLikelihood = log_likelihood;
	}
	
	public TreeLikelihoodObj(double beta, String treeN, double log_likelihood, double d, double l){
		this.bestBeta = beta;
		this.bestTreeNewickString = treeN;
		this.bestTreeLikelihood = log_likelihood;
		this.bestDelProb = d;
		this.bestLOHProb = l;
	}
	
	public static void main(String[] args) throws IOException{
		Comparator<TreeLikelihoodObj> comparator = new Comparator<TreeLikelihoodObj>() {
			
			@Override
			public int compare(TreeLikelihoodObj o1, TreeLikelihoodObj o2) {
				// TODO Auto-generated method stub
				if (o1.bestTreeLikelihood > o2.bestTreeLikelihood)
					return 1;
				else if (o1.bestTreeLikelihood == o2.bestTreeLikelihood)
					return 0;
				else 
					return -1;
			}
		};
		PriorityQueue<TreeLikelihoodObj> queue = 
	            new PriorityQueue<TreeLikelihoodObj>(5, comparator);
		
		Random _rng = new Random();
		for (int i = 0; i < 5; i++){
			double beta = _rng.nextDouble();
			double l = 0.1;
			String s = Integer.toString(i);
			TreeLikelihoodObj TL = new TreeLikelihoodObj(beta, s, l);
			queue.add(TL);
		}
		queue.add(new TreeLikelihoodObj(0.1, "as", 0.2));
//		while (queue.size() != 0)
//        {
//            System.out.println(queue.remove().bestTreeLikelihood);
//        }
		
		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		
//		String treeFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/CO5_rslts/nuModel/CO5.newick";
		String treeFile = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/newRealDatasets/nuModel/CRC.newick";
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
		for (STINode<Double> node : intree.getNode("in28").getChildren()){
			System.out.printf("child %s has parent distance = %f \n", node.getName(), node.getParentDistance());
		}
//		
		System.out.println(intree.toNewick());
//		
////		System.out.println(intree.getNode("MA_55").getParent().getName());
//		STINode<Double> in108 = intree.getNode("in108");
//		in108.getParent().adoptChild(intree.getNode("in107"));
//		intree.getNode("in95").adoptChild(in108);
//		in108.adoptChild(intree.getNode("MA_58"));
////		STINode<Double> in41 = intree.getNode("in41");
////		STINode<Double> in40 = intree.getNode("in40");
////		STINode<Double> in42 = intree.getNode("in42");
////		STINode<Double> in43 = intree.getNode("in43");
////		STINode<Double> MD_1 = intree.getNode("MD_1");
////		STINode<Double> MD_5 = intree.getNode("MD_5");
////////		MD_1.setParent(in40);
////		in40.adoptChild(MD_1);
//////		in41.removeChild(MD_1, false);
////		STINode<Double> MD_20 = intree.getNode("MD_20");
////		
//////		MD_20.setParent(in41);
////		in41.adoptChild(in42);
////		in42.adoptChild(MD_20);
////		MD_5.setParent(in43);
////		
////
////		
////		STINode<Double> in146 = intree.getNode("in146");
////		STINode<Double> in147 = intree.getNode("in147");
////		in147.adoptChild(in41);
////		in42.adoptChild(in146);
////		
////		STINode<Double> in45 = intree.getNode("in45");
////		STINode<Double> in44 = intree.getNode("in44");
////		in45.getParent().adoptChild(in44);
////		in45.adoptChild(intree.getNode("PD_30"));
////		in147.adoptChild(in45);
////		
////		
//		System.out.println(intree.toNewick());
//		for (STINode<Double> c: in108.getChildren())
//			System.out.println(c.getName());
//		System.out.println(MD_20.getName());
	}

}
