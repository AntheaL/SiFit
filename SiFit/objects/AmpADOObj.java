/**
 * 
 */
package SiFit.objects;

import java.util.ArrayList;
import java.util.Random;

/**
 * @author hz22
 *
 */
public class AmpADOObj {

	/**
	 * @param args
	 */
	public Random _rng = new Random();
	public int Cell_Node_ID;
	public String Cell_Name;
	public ArrayList<Integer> ADOAffectedPosList = new ArrayList<>(); // List of heterozygous positions affected by ADO
	public ArrayList<Integer> HeterozygousPosList; // List of heterozygous positions present
	public Integer[] beforeADOGenotypeArr;
	public Integer[] afterADOGenotypeArr;
	
	public AmpADOObj(String cellName){
		this.Cell_Name = cellName;
	}
	
	public AmpADOObj(Integer[] cellGenotypeArr, String cellName){
		this.beforeADOGenotypeArr = cellGenotypeArr;
		this.Cell_Name = cellName;
	}
	
	public void getADOAffectedGenotypeArr(double ADO_rate, Integer[] beforeADOGenotypeArr){
		this.afterADOGenotypeArr =  new Integer[beforeADOGenotypeArr.length];
		for (int i = 0; i < beforeADOGenotypeArr.length; i++){
			afterADOGenotypeArr[i] = beforeADOGenotypeArr[i];
			if (beforeADOGenotypeArr[i] == 1){
				double r_prob = _rng.nextDouble();
				// ADO happening here
				if (r_prob < ADO_rate){
					double refORalt = _rng.nextDouble();
					// Alt allele is dropped
					if (refORalt <= 0.5){
						afterADOGenotypeArr[i] = 0;
					}
					// Ref allele is dropped
					else{
						afterADOGenotypeArr[i] = 2;
					}
					ADOAffectedPosList.add(i);
				}
			}
		}
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
