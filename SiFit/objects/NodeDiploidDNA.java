/**
 * 
 */
package SiFit.objects;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;

/**
 * Class for Diploid DNA associated with a Node in the tree
 * @author hz22
 *
 */
public class NodeDiploidDNA {
	public String DNA_copy_1;
	public String DNA_copy_2;
	public String node_name;
	public STINode<Double> node;
	
	/**
	 * Constructor
	 * @param s1
	 * @param s2
	 */
	public NodeDiploidDNA(String s1, String s2){
		this.DNA_copy_1 = s1;
		this.DNA_copy_2 = s2;
	}
	
	/**
	 * Set the name of the node
	 * @param s
	 */
	public void setNodeName(String s){
		this.node_name = s;
	}
	
	/**
	 * Set the Node
	 * @param node
	 */
	public void setNode(STINode<Double> node){
		this.node = node;
	}
	
	public String getNodeName(){
		return this.node_name;
	}
	
	/**
	 * Set DNA_copy_1
	 * @param s
	 */
	public void setDNA_1(String s){
		this.DNA_copy_1 = s;
	}
	
	/**
	 * Set DNA_copy_2
	 * @param s
	 */
	public void setDNA_2(String s){
		this.DNA_copy_2 = s;
	}
	
	/**
	 * Set DNA_copy_1 and DNA_copy_2
	 * @param s1
	 * @param s2
	 */
	public void setDiploidDNA(String s1, String s2){
		setDNA_1(s1);
		setDNA_2(s2);
	}
	
	public static void main(String args[]){
		NodeDiploidDNA test_DNA = new NodeDiploidDNA("ACCC", "GTTT");
		test_DNA.setNodeName("test");
		test_DNA.setDiploidDNA("AAAAT", "ACCTG");
		System.out.println(test_DNA.DNA_copy_1);
		System.out.println(test_DNA.DNA_copy_2);
	}


}
