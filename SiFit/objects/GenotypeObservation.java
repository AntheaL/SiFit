/**
 * 
 */
package SiFit.objects;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Class encoding the genotype observation at a locus of the genome.
 * the member cellGenotypeMap is a Map <cellname, genotype> for the given locus
 * @author hz22
 *
 */
public class GenotypeObservation {
	public HashMap<String, Integer> cellGenotypeMap = new HashMap<String, Integer>();
	
	/**
	 * Constructor that instantiates the class
	 * @param node_names - List of string (cell names)
	 * @param genotypes - Array of integer genotypes
	 */
	public GenotypeObservation(ArrayList<String> node_names, Integer[] genotypes){
		for (int i = 0; i < node_names.size(); i ++){
			cellGenotypeMap.put(node_names.get(i), genotypes[i]);
		}
	}
	
}
