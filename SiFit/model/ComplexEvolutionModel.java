/**
 * 
 */
package SiFit.model;

import org.apache.commons.math3.distribution.NormalDistribution;

import cern.colt.Arrays;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;
import jeigen.DenseMatrix;

/**
 * @author hz22
 *
 */
public class ComplexEvolutionModel extends SubstitutionModel {
	
	DenseMatrix rateMatrix;
	public double delProb;
	public double LOHProb;
	public double recurProb;
	
	public ComplexEvolutionModel(double d, double l, double r) {
       
		this.delProb = d;
		this.LOHProb = l; 
		this.recurProb = r;
        
    }
	public ComplexEvolutionModel(double d, double l) {
	       
		this.delProb = d;
		this.LOHProb = l; 
//		this.recurProb = r;
        
    }
	
	

	/* (non-Javadoc)
	 * @see edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel#getEquilibriumVector()
	 */
	@Override
	public DenseMatrix getEquilibriumVector() {
		return new DenseMatrix(new double[][] {{1},{0},{0}});
	}
	
	public DenseMatrix getEquilibriumVectorBinary() {
		return new DenseMatrix(new double[][] {{1},{0}});
		
	}
	
	public double[][] getTransitionMatrix(double bl){
		double[][] TransMat = new double[3][3];
		
		
		// Complex Evolution model 1
		double commT = Math.pow((Math.pow(LOHProb,2) + 2*LOHProb*delProb + 2*Math.pow(delProb,2) - 2*delProb + 1), 0.5); // this is a term that appears multiple times in the eqns below.
		
		TransMat[0][0] = (Math.exp(-(bl*(LOHProb + 2*delProb + commT + 1))/2)*((Math.exp(bl*commT) - 1)*(Math.pow(LOHProb,2) + 2*LOHProb*delProb + 2*Math.pow(delProb,2) - 2*delProb + 1) - Math.exp(bl*commT) + Math.pow(delProb,2)*(2*Math.exp(bl*commT) - 2) - LOHProb*(Math.exp(bl*commT) - 1) - delProb*(Math.exp(bl*commT) - 1) + delProb*(3*Math.exp(bl*commT) + 3)*commT + 2*Math.pow(delProb,2)*Math.exp((bl*(LOHProb + 2*delProb + commT + 1))/2)*commT + LOHProb*(Math.exp(bl*commT) + 1)*commT + LOHProb*delProb*(Math.exp(bl*commT) - 1) + 2*LOHProb*delProb*Math.exp((bl*(LOHProb + 2*delProb + commT + 1))/2)*commT + 1))/(2*(LOHProb + 3*delProb + LOHProb*delProb + Math.pow(delProb,2))*commT);
		TransMat[0][1] = -(Math.exp(-(bl*(LOHProb + 2*delProb + commT + 1))/2)*(LOHProb + 2*delProb + delProb*commT - LOHProb*Math.exp(bl*commT) - 2*delProb*Math.exp(bl*commT) + Math.pow(delProb,2)*Math.exp(bl*commT) - Math.pow(delProb,2) - 2*delProb*Math.exp((bl*(LOHProb + 2*delProb + commT + 1))/2)*commT + delProb*Math.exp(bl*commT)*commT))/((LOHProb + 3*delProb + LOHProb*delProb + Math.pow(delProb,2))*commT);
		TransMat[0][2] = -(Math.exp(-(bl*(LOHProb + 2*delProb + commT + 1))/2)*(LOHProb + delProb)*(Math.exp(bl*commT) - 2*delProb - LOHProb + LOHProb*Math.exp(bl*commT) + Math.exp(bl*commT)*commT + 2*delProb*Math.exp(bl*commT) + commT - 2*Math.exp((bl*(LOHProb + 2*delProb + commT + 1))/2)*commT - 1))/(2*(LOHProb + 3*delProb + LOHProb*delProb + Math.pow(delProb,2))*commT);
		TransMat[1][0] = -(Math.exp(-(bl*(LOHProb + 2*delProb + commT + 1))/2)*(LOHProb + delProb)*(LOHProb + 2*delProb + delProb*commT - LOHProb*Math.exp(bl*commT) - 2*delProb*Math.exp(bl*commT) + Math.pow(delProb,2)*Math.exp(bl*commT) - Math.pow(delProb,2) - 2*delProb*Math.exp((bl*(LOHProb + 2*delProb + commT + 1))/2)*commT + delProb*Math.exp(bl*commT)*commT))/(2*(LOHProb + 3*delProb + LOHProb*delProb + Math.pow(delProb,2))*commT);
		TransMat[1][1] = (Math.exp(-(bl*(LOHProb + 2*delProb + commT + 1))/2)*(Math.exp(bl*commT) - (Math.exp(bl*commT) - 1)*(Math.pow(LOHProb,2) + 2*LOHProb*delProb + 2*Math.pow(delProb,2) - 2*delProb + 1) + Math.pow(delProb,3)*(2*Math.exp(bl*commT) - 2) - Math.pow(delProb,2)*(3*Math.exp(bl*commT) - 3) + LOHProb*(Math.exp(bl*commT) - 1) - LOHProb*delProb*(2*Math.exp(bl*commT) - 2) + LOHProb*Math.pow(delProb,2)*(Math.exp(bl*commT) - 1) + delProb*(Math.exp(bl*commT) + 4*Math.exp((bl*(LOHProb + 2*delProb + commT + 1))/2) + 1)*commT + Math.pow(delProb,2)*(Math.exp(bl*commT) + 1)*commT - delProb*(Math.exp(bl*commT) - 1)*(Math.pow(LOHProb,2) + 2*LOHProb*delProb + 2*Math.pow(delProb,2) - 2*delProb + 1) + LOHProb*(Math.exp(bl*commT) + 1)*commT + LOHProb*delProb*(Math.exp(bl*commT) + 1)*commT - 1))/(2*(LOHProb + 3*delProb + LOHProb*delProb + Math.pow(delProb,2))*commT);
		TransMat[1][2] = -(Math.exp(-(bl*(LOHProb + 2*delProb + commT + 1))/2)*(LOHProb + delProb)*(delProb + Math.exp(bl*commT) + LOHProb*delProb + Math.exp(bl*commT)*commT - delProb*Math.exp(bl*commT) + commT - Math.pow(delProb,2)*Math.exp(bl*commT) - 2*Math.exp((bl*(LOHProb + 2*delProb + commT + 1))/2)*commT + Math.pow(delProb,2) - LOHProb*delProb*Math.exp(bl*commT) - 1))/(2*(LOHProb + 3*delProb + LOHProb*delProb + Math.pow(delProb,2))*commT);
		TransMat[2][0] = -(delProb*Math.exp(-(bl*(LOHProb + 2*delProb + commT + 1))/2)*(LOHProb + delProb)*(Math.exp(bl*commT) - 2*delProb - LOHProb + LOHProb*Math.exp(bl*commT) + Math.exp(bl*commT)*commT + 2*delProb*Math.exp(bl*commT) + commT - 2*Math.exp((bl*(LOHProb + 2*delProb + commT + 1))/2)*commT - 1))/(2*(LOHProb + 3*delProb + LOHProb*delProb + Math.pow(delProb,2))*commT);
		TransMat[2][1] = -(Math.exp(-(bl*(LOHProb + 2*delProb + commT + 1))/2)*(delProb*commT - delProb + delProb*Math.exp(bl*commT) + LOHProb*Math.pow(delProb,2) - Math.pow(delProb,2)*Math.exp(bl*commT) - Math.pow(delProb,3)*Math.exp(bl*commT) + Math.pow(delProb,2) + Math.pow(delProb,3) - 2*delProb*Math.exp((bl*(LOHProb + 2*delProb + commT + 1))/2)*commT + delProb*Math.exp(bl*commT)*commT - LOHProb*Math.pow(delProb,2)*Math.exp(bl*commT)))/((LOHProb + 3*delProb + LOHProb*delProb + Math.pow(delProb,2))*commT);
		TransMat[2][2] = (Math.exp(-(bl*(LOHProb + 2*delProb + commT + 1))/2)*(Math.pow(delProb,2)*(Math.exp(bl*commT) - 1) - Math.pow(delProb,3)*(2*Math.exp(bl*commT) - 2) + delProb*(Math.exp(bl*commT) - 1) + 2*LOHProb*Math.exp((bl*(LOHProb + 2*delProb + commT + 1))/2)*commT - LOHProb*Math.pow(delProb,2)*(Math.exp(bl*commT) - 1) + Math.pow(delProb,2)*(Math.exp(bl*commT) + 1)*commT + delProb*(Math.exp(bl*commT) - 1)*(Math.pow(LOHProb,2) + 2*LOHProb*delProb + 2*Math.pow(delProb,2) - 2*delProb + 1) + delProb*(2*Math.exp(bl*commT) + 2*Math.exp((bl*(LOHProb + 2*delProb + commT + 1))/2) + 2)*commT + LOHProb*delProb*(Math.exp(bl*commT) - 1) + LOHProb*delProb*(Math.exp(bl*commT) + 1)*commT))/(2*(LOHProb + 3*delProb + LOHProb*delProb + Math.pow(delProb,2))*commT);
		
		

		return TransMat;
	}
	
	/**
	 * Compute the transition matrix for the binary data type
	 * @param bl, branch length representing the expected number of mutations per site
	 * @return, the transition matrix 
	 * 		new		0 		1
	 * old
	 * 0		P(0 -> 0) P(0 -> 1)			
	 * 1		P(1 -> 0) P(1 -> 1)
	 */
	public double[][] getTransitionMatrixBinary(double bl){
		double[][] TransMat = new double[2][2];
		
		//custom model 1
		double commT = (LOHProb + delProb + 2);
		
		TransMat[0][0] = (LOHProb + delProb + 2*Math.exp(-(bl*commT)/2))/commT;
		TransMat[0][1] = -(2*(Math.exp(-(bl*commT)/2) - 1))/commT; 
		TransMat[1][0] = (Math.exp(-(bl*commT)/2)*(LOHProb + delProb)*(Math.exp((bl*commT)/2) - 1))/commT;
		TransMat[1][1] = (Math.exp(-(bl*commT)/2)*(LOHProb + delProb + 2*Math.exp((bl*commT)/2)))/commT;

		
		return TransMat;
		
	}
	
	

	/* (non-Javadoc)
	 * @see edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.RateModel#getRateMatrix()
	 */
	@Override
	public DenseMatrix getRateMatrix() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel#getProbabilityMatrixIntegrated()
	 */
	@Override
	public DenseMatrix getProbabilityMatrixIntegrated() {
		// TODO Auto-generated method stub
		return null;
	}
	
	/**
	 * Propose a new value from the current value. Sampled from a normal distribution centered at the current value
	 * @param currVal
	 * @param sd
	 * @return
	 */
	public double proposeNewParamVal(double currVal, double sd){
		NormalDistribution NormalDist = new NormalDistribution(currVal, sd);
		double newVal = NormalDist.sample();
		if (newVal < 0)
			newVal = Math.abs(newVal);
		if (newVal > 1)
			newVal = newVal - 2*(newVal - 1);
		return newVal;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		ComplexEvolutionModel model = new ComplexEvolutionModel(0.05, 0.01, 0.05);
		double[][] tm = model.getTransitionMatrix(0.2);
		for (double[] arr: tm)
			System.out.println(Arrays.toString(arr));
		JCModelSingleCell jc = new JCModelSingleCell(0.1);
		System.out.println("");
		double[][] mm = jc.jukesCantor012AllMatrix(0.2);
//		double[][] mm = jc.jukesCantorBinaryMatrixAll(0.2);
		for (double[] arr: mm)
			System.out.println(Arrays.toString(arr));
		// TODO Auto-generated method stub

	}

}
