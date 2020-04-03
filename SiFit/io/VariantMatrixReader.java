/**
 * 
 */
package SiFit.io;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.io.output.ThresholdingOutputStream;

/**
 * @author hz22
 *
 */
public class VariantMatrixReader {
	public int n_mutation;	
	public String filename;
	public ArrayList<Integer[]> obsVarGTMatrix = new ArrayList<>(); // List of genotypes for all mutations
	public ArrayList<Integer> mutPosList = new ArrayList<>();
	public ArrayList<String> scNames = new ArrayList<>();
	
	public VariantMatrixReader(String filename){
		this.filename = filename;
	}
	
	/**
	 * Read from the file and populate the entries of obsVarGTMatrix
	 * @param filename
	 * @param nCell
	 * @throws IOException
	 */
	public void populateVarGTMatrix(String filename, int nCell) throws IOException{
		try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = br.readLine()) != null) { 
            	String[] thisMutVector = line.split(" ");
            	Integer[] thisMutGTArr = new Integer[nCell];
            	for (int i = 1; i <= nCell; i++){
            		thisMutGTArr[i-1] = Integer.parseInt(thisMutVector[i]);
            	}
            	this.obsVarGTMatrix.add(thisMutGTArr);            	
            }
		}
	}
	
	/**
	 * Obtains a cell, genotype vector map from the input expected variant matrix
	 * @param filename
	 * @param nCell
	 * @return
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public HashMap<String, Integer[]> getCellGTVectorMap(String filename, int nMut) throws FileNotFoundException, IOException{
		HashMap<String, Integer[]> cellGTVectorMap = new HashMap<>();
		try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
			String line;
            while ((line = br.readLine()) != null) { 
            	String[] thisCellVector = line.split("\t");
//            	System.out.println(Arrays.toString(thisCellVector));
            	Integer[] thisCellGTVector = new Integer[nMut];
            	for (int i = 1; i <= nMut; i++){
            		thisCellGTVector[i-1] = Integer.parseInt(thisCellVector[i]);
            	}
            	cellGTVectorMap.put(thisCellVector[0], thisCellGTVector);
            	this.scNames.add(thisCellVector[0]);
            }
		}
		return cellGTVectorMap;
	}

}
