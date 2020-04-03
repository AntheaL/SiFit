/**
 * 
 */
package SiFit.model;

import jeigen.DenseMatrix;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.SubstitutionModel;

/**
 * @author hz22
 *
 */
public class JCModelSingleCell extends SubstitutionModel {
	
	DenseMatrix rateMatrix;
	double mu; // mutation rate
	
	/**
	 * Constructor 
	 * @param equilibriumFrequencies
	 * @param transitionFrequencies
	 */
	public JCModelSingleCell(double[] equilibriumFrequencies, double[] transitionFrequencies) {
        //Pass it to next constructor.
        this(transitionFrequencies[0]);
    }
	
	/**
	 * Constructor 
	 * @param mut_rate - mutation rate, mu in the formulae of JC Model
	 */
	public JCModelSingleCell(double mut_rate){
		this.mu = mut_rate;
		this.rateMatrix = createMatrix(mut_rate);
	}
	
	/**
	 * Create rateMatrix 
	 * @param transition - transition probability
	 * @return
	 */
	public DenseMatrix createMatrix(double transition){
		double[][] rateMatrixMaker = new double[3][3];
		for (int i = 0; i < rateMatrixMaker.length; i++){
			double[] row = rateMatrixMaker[i];
			for (int j = 0; j < row.length; j++){
				if (j == i)
					rateMatrixMaker[i][j] = -3.0/4 * transition;
				else
					rateMatrixMaker[i][j] = -1.0/4 * transition;
			}
		}
		
		return new DenseMatrix(rateMatrixMaker);
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
	
	/**
	 * Computes the transition probability matrix for JC
	 * All transitions are not allowed. No back substitution
	 * 		new		0			1			2
	 *  old
	 *  0		samesame	2*samediff		diffdiff
	 *  1			0		1-samediff		samediff
	 *  2			0			0				1
	 * @param distance - branch length
	 * @return
	 */
	public double[][] jukesCantor012Matrix(double distance){
		double[][] JCTransMat = new double[3][3];
		
		double sameProb = 0.25 + 0.75 * Math.exp(-1 * mu * distance);
		double diffProb = 0.25 - 0.25 * Math.exp(-1 * mu * distance);
		
		double samesame = sameProb * sameProb;
		double samediff = sameProb * diffProb;
		double diffdiff = diffProb * diffProb;
		
		JCTransMat[0][0] = samesame;
		JCTransMat[0][1] = 2 * samediff;
		JCTransMat[0][2] = diffdiff;
		JCTransMat[1][0] = 0;
		JCTransMat[1][1] = 1 - samediff;
		JCTransMat[1][2] = samediff;
		JCTransMat[2][0] = 0;
		JCTransMat[2][1] = 0;
		JCTransMat[2][2] = 1;
		
		return JCTransMat;
	}
	
	/**
	 * Computes the transition probability matrix for JC
	 * All transitions are allowed.
	 * @param distance
	 * @return
	 */
	public double[][] jukesCantor012AllMatrix(double distance){
		double[][] JCTransMat = new double[3][3];
		
		double sameProb = 0.25 + 0.75 * Math.exp(-1 * mu * distance);
		double diffProb = 0.25 - 0.25 * Math.exp(-1 * mu * distance);
		
		double samesame = sameProb * sameProb;
		double samediff = sameProb * diffProb;
		double diffdiff = diffProb * diffProb;
		
		JCTransMat[0][0] = samesame;
		JCTransMat[0][1] = 2 * samediff;
		JCTransMat[0][2] = diffdiff;
		JCTransMat[1][0] = 0.5*samediff;
		JCTransMat[1][1] = 1 - samediff;
		JCTransMat[1][2] = 0.5*samediff;
		JCTransMat[2][0] = diffdiff;
		JCTransMat[2][1] = 2 * samediff;
		JCTransMat[2][2] = samesame;
		
		return JCTransMat;
	}
	
	/**
	 * Computes transition matrix for Binary data
	 * All transitions are allowed
	 * @param distance
	 * @return
	 */
	public double[][] jukesCantorBinaryMatrixAll(double distance){
		double[][] JCTransMat = new double[2][2];
		
		double sameProb = 0.25 + 0.75 * Math.exp(-1 * mu * distance);
		double diffProb = 0.25 - 0.25 * Math.exp(-1 * mu * distance);
		double samesame = sameProb * sameProb;
		double samediff = sameProb * diffProb;
		
		JCTransMat[0][0] = samesame;
		JCTransMat[0][1] = samediff;
		JCTransMat[1][0] = samediff;
		JCTransMat[1][1] = 1-samediff;
		
		return JCTransMat;
	}
	
	/**
	 * Computes transition matrix for Binary data
	 * All transitions are not allowed, no reversal of genotypes
	 * @param distance
	 * @return
	 */
	public double[][] jukesCantorBinaryMatrix01(double distance){
		double[][] JCTransMat = new double[2][2];
		
		double sameProb = 0.25 + 0.75 * Math.exp(-1 * mu * distance);
		double diffProb = 0.25 - 0.25 * Math.exp(-1 * mu * distance);
		double samesame = sameProb * sameProb;
		double samediff = sameProb * diffProb;
		
		JCTransMat[0][0] = samesame;
		JCTransMat[0][1] = samediff;
		JCTransMat[1][0] = 0;
		JCTransMat[1][1] = 1;
		
		return JCTransMat;
	}
	public double[][] testModel(double u, double r){
		double[][] JCTransMat = new double[2][2];
		JCTransMat[0][0] = 1-u;
		JCTransMat[0][1] = u;
		JCTransMat[1][0] = r;
		JCTransMat[1][1] = 1-r;
		return JCTransMat;
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

}
