package SiFit.objects;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;
import java.util.Set;

public class AmpFPObj {

	/**
	 * @param args
	 */
	public Random _rng = new Random();
	public int Cell_Node_ID;
	public String Cell_Name;
	public ArrayList<Integer> FPAffectedPosList = new ArrayList<>();
	public Integer[] beforeFPGenotypeArr;
	public Integer[] afterFPGenotypeArr;
	
	public AmpFPObj(Integer[] cellGenotypeArr, String cellName){
		this.beforeFPGenotypeArr = cellGenotypeArr;
		this.Cell_Name = cellName;
	}
	
	public void getFPAffectedGenotypeArr(double FP_rate, Integer[] beforeFPGenotypeArr, HashMap<Integer, Set<Integer>> Genotype_flag_array){
		this.afterFPGenotypeArr = new Integer[beforeFPGenotypeArr.length];
		for (int i = 0; i < beforeFPGenotypeArr.length; i++){
			afterFPGenotypeArr[i] = beforeFPGenotypeArr[i];
			if (beforeFPGenotypeArr[i] == 0){
				double r_prob = _rng.nextDouble();
				if (r_prob < FP_rate){
					afterFPGenotypeArr[i] = 1;
					Genotype_flag_array.get(1).add(i);
					FPAffectedPosList.add(i);
				}
			}
			else if (beforeFPGenotypeArr[i] == 1){
				double r_prob = _rng.nextDouble();
				if (r_prob < FP_rate){
					afterFPGenotypeArr[i] = 2;
					Genotype_flag_array.get(2).add(i);
					FPAffectedPosList.add(i);
				}
			}
		}
	}
	
	public void getFPAffectedGenotypeArr(double FP_rate, Integer[] beforeFPGenotypeArr){
		this.afterFPGenotypeArr = new Integer[beforeFPGenotypeArr.length];
		for (int i = 0; i < beforeFPGenotypeArr.length; i++){
			afterFPGenotypeArr[i] = beforeFPGenotypeArr[i];
			if (beforeFPGenotypeArr[i] == 0){
				double r_prob = _rng.nextDouble();
				if (r_prob < FP_rate){
					afterFPGenotypeArr[i] = 1;
					
					FPAffectedPosList.add(i);
				}
			}
//			else if (beforeFPGenotypeArr[i] == 1){
//				double r_prob = _rng.nextDouble();
//				if (r_prob < FP_rate){
//					afterFPGenotypeArr[i] = 2;
//					
//					FPAffectedPosList.add(i);
//				}
//			}
		}
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
