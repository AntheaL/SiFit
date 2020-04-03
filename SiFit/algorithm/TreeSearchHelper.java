/**
 * 
 */
package SiFit.algorithm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import SiFit.PerturbFn;
import SiFit.TopologyBranchPerturbations;
import SiFit.BasicUtilityFunctions;
import SiFit.model.ComplexEvolutionModel;
import SiFit.model.JCModelSingleCell;
import SiFit.objects.GenotypeObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.hj.api.SuspendableException;

import static edu.rice.hj.Module1.*;

/**
 * @author hz22
 *
 */
public class TreeSearchHelper {
	
	public BasicUtilityFunctions BUF; 
	public int dataFlag;
	public int modelFlag;
	public double fn;
	public double fp;
	public double mu;
	public JCModelSingleCell model;
	public ComplexEvolutionModel modelC;
	public TopologyBranchPerturbations TBP;
	public PerturbFn proposeFnObj;
	private Random _rng;
	public double newFnSameTreeLikelihood;
	public double sameFnNewTreeLikelihood;
	HashMap<String, Integer> _doubletFlagMap;
	public double newParamSameTreeLikelihood;
	
	public TreeSearchHelper(double fp, double fn, double mu, BasicUtilityFunctions BUF, int df, int mf){
		this.fp = fp;
		this.fn = fn;
		this.mu = mu;
		this.BUF = BUF;
		this.dataFlag = df;
		this.modelFlag = mf;
		this.model = new JCModelSingleCell(mu);
		
//		this.model = new ComplexEvolutionModel(0.05, 0.01, 0.05);
		this.TBP = new TopologyBranchPerturbations();
		this.proposeFnObj = new PerturbFn(fn);
		this._rng = new Random();
	}
	public TreeSearchHelper(double fp, double fn, double mu, BasicUtilityFunctions BUF, int df, int mf, HashMap<String, Integer> doubletFlagMap){
		this.fp = fp;
		this.fn = fn;
		this.mu = mu;
		this.BUF = BUF;
		this.dataFlag = df;
		this.modelFlag = mf;
		this.model = new JCModelSingleCell(mu);
//		this.model = new ComplexEvolutionModel(0.05, 0.01, 0.05);
		this.TBP = new TopologyBranchPerturbations();
		this.proposeFnObj = new PerturbFn(fn);
		this._rng = new Random();
		this._doubletFlagMap = doubletFlagMap;
	}
	
	/**
	 * This constructor is for using the complex model of evolution
	 * @param fp
	 * @param fn
	 * @param del
	 * @param LOH
	 * @param BUF
	 * @param df
	 * @param mf
	 */
	public TreeSearchHelper(double fp, double fn, double del, double LOH, BasicUtilityFunctions BUF, int df, int mf){
		this.fp = fp;
		this.fn = fn;
		this.BUF = BUF;
		this.dataFlag = df;
		this.modelFlag = mf;
		this.modelC = new ComplexEvolutionModel(del, LOH);
		this.TBP = new TopologyBranchPerturbations();
		this.proposeFnObj = new PerturbFn(fn);
		this._rng = new Random();
	}
	
	public TreeSearchHelper(double fp, double fn, double del, double LOH, double recur, BasicUtilityFunctions BUF, int df, int mf){
		this.fp = fp;
		this.fn = fn;
		this.BUF = BUF;
		this.dataFlag = df;
		this.modelFlag = mf;
		this.modelC = new ComplexEvolutionModel(del, LOH, recur);
		this.TBP = new TopologyBranchPerturbations();
		this.proposeFnObj = new PerturbFn(fn);
		this._rng = new Random();
	}
	
	/**
	 * Computes the likelihood of a Tree and Error rate
	 * @param GTMatrix, list of GenotypeObservation, input data
	 * @param tree, candidate tree
	 * @param newFn, candidate error rate
	 * @return the log-likelihood of configuration (tree, error) 
	 * @throws SuspendableException 
	 */
//	public double computeTreeLikelihood(
//			ArrayList<GenotypeObservation> GTMatrix,
//			STITree<Double> tree, double newFn) {
//		
//		this.fn = newFn;
//		double sum = 0;
//		// Binary Data
//		if (this.dataFlag == 0){
//			for (GenotypeObservation genotypeObs : GTMatrix){
//				double val = new ComputeLikelihoodBinary(fp, fn, mu, tree, model).getGenotypeProbability(genotypeObs, modelFlag);
//				double log_val;
//				if (val == 0)
//					log_val = -10000;
//				else
//					log_val = Math.log(val);
//				sum += log_val;
//			}
//		}
//		// Ternary Data
//		else{
//			for (GenotypeObservation genotypeObs : GTMatrix){
//				double val = new ComputeLikelihood(fp, fn, mu, tree, model).getGenotypeProbability(genotypeObs, modelFlag);
//				double log_val;
//				if (val == 0)
//					log_val = -10000;
//				else
//					log_val = Math.log(val);
//				sum += log_val;
//			}
//		}
//				
//		return sum;
//	}
	
	public double computeTreeLikelihood(
			ArrayList<GenotypeObservation> GTMatrix,
			STITree<Double> tree, double newFn) throws SuspendableException {
		this.fn = newFn;
		double sum = 0;
		int nMut = GTMatrix.size(); 
		double[] val = new double[nMut];
		
		if (this.dataFlag == 0){
			assert(false);
//			launchHabaneroApp(() -> {
//			ComputeLikelihood CL = new ComputeLikelihood(fp, newFn, mu, tree, model); // This is for JC model
			ComputeLikelihoodBinary CL = new ComputeLikelihoodBinary(fp, newFn, mu, tree, modelC); // This is for Complex model
			forallChunked(0, nMut -1, (i) -> {
				val[i] = CL.getGenotypeProbability(GTMatrix.get(i), modelFlag);
			});
//			});
			for (double d: val){
//				System.out.println(d);
				double log_d;
				if (d == 0)
					log_d = -10000;
				else
					log_d = Math.log(d);
				sum += log_d;
			}
			
		}
		else{
//			ComputeLikelihood CL = new ComputeLikelihood(fp, newFn, mu, tree, model);
			ComputeLikelihood CL = new ComputeLikelihood(fp, newFn, mu, tree, modelC);
//			launchHabaneroApp(() -> {
//				System.out.println(numWorkerThreads());
			forallChunked(0, nMut -1, (i) -> {
				val[i] = CL.getGenotypeProbability(GTMatrix.get(i), modelFlag);
			});
//			});
			for (double d: val){
				double log_d;
				if (d == 0)
					log_d = -10000;
				else
					log_d = Math.log(d);
				sum += log_d;
			}
		}
//		System.out.println(sum);
		return sum;
	}

	/**
	 * Returns likelihood ratio for two fn values
	 * @param obsGTObsMatrix
	 * @param currFnLogLikelihood
	 * @param newFn
	 * @param currTree
	 * @return
	 * @throws SuspendableException 
	 */
	public double getLogLikelihoodRatio(
			ArrayList<GenotypeObservation> obsGTObsMatrix,
			double currFnLogLikelihood, double newFn,
			STITree<Double> currTree) throws SuspendableException {
		
		double loglikelihood_newFn = this.computeTreeLikelihood(obsGTObsMatrix, currTree, newFn);
		this.newFnSameTreeLikelihood = loglikelihood_newFn;
		return loglikelihood_newFn - currFnLogLikelihood;
	}
	
	/**
	 * Get the loglikelihood ratio for a new value of model of evolution parameter.
	 * @param obsGTObsMatrix
	 * @param currDelProbLogLikelihood
	 * @param newdelProb
	 * @param currTree
	 * @param modelParamFlag
	 * @return
	 * @throws SuspendableException
	 */
	public double getLogLikelihoodRatio(
			ArrayList<GenotypeObservation> obsGTObsMatrix,
			double currDelProbLogLikelihood, double newdelProb,
			STITree<Double> currTree, int modelParamFlag) throws SuspendableException {
		ComplexEvolutionModel modelcopy = this.modelC;
		// delProb is being used
		if (modelParamFlag == 0){
			modelcopy.delProb = newdelProb;
		}
		// LOHProb is being used
		else if (modelParamFlag == 1){
			modelcopy.LOHProb = newdelProb;
		}
		else {
			modelcopy.recurProb = newdelProb;
		}
		// Calculate the new likelihood
		double newParamLogLikelihood = this.computeTreeLikelihood(obsGTObsMatrix, currTree, modelcopy);
		this.newParamSameTreeLikelihood = newParamLogLikelihood;
		return newParamLogLikelihood - currDelProbLogLikelihood;
	}
	
	/**
	 * Computes the likelihood of a model. tree and Fn are fixed.
	 * @param GTMatrix
	 * @param tree
	 * @param modelcopy
	 * @return
	 * @throws SuspendableException
	 */
	public double computeTreeLikelihood(ArrayList<GenotypeObservation> GTMatrix, STITree<Double> tree,
			ComplexEvolutionModel modelcopy) throws SuspendableException {
		
		double sum = 0;
		int nMut = GTMatrix.size();
		double[] val = new double[nMut];
		
		if (this.dataFlag == 0){
			assert(false);
			ComputeLikelihoodBinary CL = new ComputeLikelihoodBinary(fp, fn, mu, tree, modelcopy); // This is for Complex model
			forallChunked(0, nMut -1, (i) -> {
				val[i] = CL.getGenotypeProbability(GTMatrix.get(i), modelFlag);
			});
			for (double d: val){
				double log_d;
				if (d == 0)
					log_d = -10000;
				else
					log_d = Math.log(d);
				sum += log_d;
			}
			
		}
		else{
			ComputeLikelihood CL = new ComputeLikelihood(fp, fn, mu, tree, modelcopy);
			forallChunked(0, nMut -1, (i) -> {
				val[i] = CL.getGenotypeProbability(GTMatrix.get(i), modelFlag);
			});
			for (double d: val){
				double log_d;
				if (d == 0)
					log_d = -10000;
				else
					log_d = Math.log(d);
				sum += log_d;
			}
		}		
		return sum;
	}
	
	/**
	 * Returns likelihood ratio for two trees
	 * @param obsGTObsMatrix
	 * @param loglikelihood_currTree
	 * @param newTree
	 * @return
	 * @throws SuspendableException 
	 */
	public double getLogLikelihoodRatio(
			ArrayList<GenotypeObservation> obsGTObsMatrix,
			double loglikelihood_currTree, STITree<Double> newTree) throws SuspendableException {
		double currFn = this.fn;
		double loglikelihood_newTree = this.computeTreeLikelihood(obsGTObsMatrix, newTree, currFn);
		this.sameFnNewTreeLikelihood = loglikelihood_newTree;
		return loglikelihood_newTree - loglikelihood_currTree;
	}
	
	/**
	 * Compute Acceptance Ratio for new Fn 
	 * @param loglikelihoodRatio
	 * @param priorRatio
	 * @param proposalRatio
	 * @return
	 */
	public double computeAcceptanceRatioFn(double loglikelihoodRatio, double priorRatio, double proposalRatio){
		double r = Math.exp(loglikelihoodRatio + priorRatio + proposalRatio);
		return Math.min(1,r);
	} 
	
	/**
	 * Compute Acceptance Ratio for new tree
	 * @param loglikelihoodRatio
	 * @param proposalRatio
	 * @return
	 */
	public double computeAcceptanceRatioTree(double loglikelihoodRatio, double proposalRatio){
		double r = Math.exp(loglikelihoodRatio) * proposalRatio;
		return Math.min(1,r);
	}
	
	/**
	 * Returns flag true/false to accept/reject the sample
	 * @param acceptance_ratio
	 * @return
	 */
	public boolean getAcceptanceFlag(double acceptance_ratio){
		if (acceptance_ratio == 1){
			return true;
		}
		else {
			double u = _rng.nextDouble();
			if (u < acceptance_ratio){
				return true;
			}
			else{
				return false;
			}
		}
	}



}
