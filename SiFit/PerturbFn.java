/**
 * 
 */
package SiFit;

import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 * @author hz22
 *
 */
public class PerturbFn {

	public double fn;
	private Random _rng;
	
	public PerturbFn(double fn){
		this.fn = fn;
		this._rng = new Random();
	}
	
	/**
	 * propose new fn value from a normal dist with mean= currFn and SD
	 * @param currFn, mean of normal dist, current value of fn
	 * @param fnSD, SD of normal dist
	 * @return proposed value of fn
	 */
	public double proposeNewFn(double currFn, double fnSD){
		NormalDistribution FnNormalDist = new NormalDistribution(currFn, fnSD);
		double newFn = FnNormalDist.sample();
		if (newFn < 0)
			newFn = Math.abs(newFn);
		if (newFn > 1)
			newFn = newFn - 2*(newFn - 1);
		return newFn;
	}
	
	public double logBetaPDF(double x, double a, double b){
		BetaDistribution betaObj = new BetaDistribution(a, b);
		return Math.log(betaObj.density(x));
	}
	
	public double getProposalRatio(double currFn, double newFn, double fnSD){
		NormalDistribution forwardNormalDist = new NormalDistribution(currFn, fnSD);
		NormalDistribution backwardNormalDist = new NormalDistribution(newFn, fnSD);
		return Math.log(backwardNormalDist.density(currFn)) - Math.log(forwardNormalDist.density(newFn));
	}

}
