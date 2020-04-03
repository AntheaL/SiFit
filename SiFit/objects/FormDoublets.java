/**
 * 
 */
package SiFit.objects;

/**
 * @author hz22
 *
 */
public class FormDoublets {

	/**
	 * @param args
	 */
	public static Integer[] getDoubletGTarr(Integer[] gt1, Integer[] gt2){
		Integer[] doubletGtarr = new Integer[gt1.length];
		for (int i = 0; i < gt1.length; i++){
			doubletGtarr[i] = getDoubletGT(gt1[i], gt2[i]);
		}
		return doubletGtarr;
	}
	
	public static Integer getDoubletGT(Integer gt1, Integer gt2){
		if (gt1 == 0 && gt2 == 0)
			return 0;
		else if (gt1 == 2 && gt2 == 2)
			return 2;
		else
			return 1;
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
