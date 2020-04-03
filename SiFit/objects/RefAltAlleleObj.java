/**
 * 
 */
package SiFit.objects;

/**
 * @author hz22
 *
 */
public class RefAltAlleleObj {
	public char ref;
	public char alt;
	public int index;
	
	public RefAltAlleleObj(char ref, char alt, int i){
		this.ref = ref;
		this.alt = alt; 
		this.index = i;
	}
	public void setRef(char c){
		this.ref = c;
	}
	public void setAlt(char c){
		this.alt = c;
	}

}
