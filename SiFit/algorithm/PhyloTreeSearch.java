/**
 * 
 */
package SiFit.algorithm;

import static edu.rice.hj.Module0.launchHabaneroApp;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Random;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.hj.api.SuspendableException;
import SiFit.BasicUtilityFunctions;
import SiFit.TopologyBranchPerturbations;
import SiFit.io.VariantMatrixReader;
import SiFit.metric.CompareTrees;
import SiFit.objects.GenotypeObservation;
import SiFit.objects.TreeLikelihoodObj;

/**
 * @author hz22
 *
 */
public class PhyloTreeSearch {
	
	/**
	 * Finds the best configuration of tree and error rate that maximizes the likelihood.
	 * This is run for one restart of SiFit. 
	 * @param obsGTObsMatrix - ArrayList<GenotypeObservation> Input Genotype Matrix represented in GenotypeObservation Form
	 * @param obsGenotypeMatrix - ArrayList<Integer[]> Input Genotype Matrix needed for Parsimony calculation
	 * @param singleCellNames - arrayList of single cell names
	 * @param nCell - Number of cells
	 * @param nMut - Number of mutations
	 * @param fp - estimated FP rate
	 * @param fnADO - starting FN rate
	 * @param errorlearn - probability for choosing error rate learning move
	 * @param mutRate - mu
	 * @param n_iter - Number of iterations
	 * @param p_iter - Number of print interations
	 * @param dataFlag - Binary/Ternary data flag
	 * @param modelFlag - which Model to use (all transitions/limited transitions)
	 * @param restartIndex - index of restart, even ones are used for parsimony initiation
	 * @param BUF - Instance of utility functions
	 * @throws SuspendableException 
	 */
	public static TreeLikelihoodObj findBestTree(
			ArrayList<GenotypeObservation> obsGTObsMatrix,
			ArrayList<Integer[]> obsGenotypeMatrix,
			ArrayList<String> singleCellNames,
			int nCell, int nMut, double fp, double fnADO,
			double errorlearn, double mutRate, 
			int n_iter, int p_iter,
			int dataFlag, int modelFlag, int restartIndex,
			BasicUtilityFunctions BUF) throws SuspendableException{
		
		Random _rng = new Random();
        // fn Parameters
//        double fnPriorMean = BUF.getFNPrior(obsGenotypeMatrix, fnADO);
        double fnPriorMean = fnADO;
        double fnPriorSD = fnPriorMean*0.5;
        double betaPriora = ((1 - fnPriorMean)*fnPriorMean*fnPriorMean/(fnPriorSD*fnPriorSD)) - fnPriorMean;
        double betaPriorb = betaPriora * ((1/fnPriorMean) - 1);
        
        // Obtain matrix for running parsimony
        ArrayList<Integer[]> obsCellGenomes = BUF.getCellGenomes(obsGenotypeMatrix, nMut);
        
		// Generate random tree
		STITree<Double> randTreeMaker = BUF.generateRandomTree(singleCellNames);
		STITree<Double> randTree = BUF.getTree(randTreeMaker.toNewick());

		// For even numbered restarts, Parsimony algorithm 
		// is run to obtain a good initial tree
		if (restartIndex % 2 == 0){
			STITree<Double> randTree1 = BUF.getTree(randTreeMaker.toNewick());
			int parsimonyIter = n_iter/10;
			
			// Run Parsimony to get the initial tree
			FitchParsimonyCalc FitchObj = new FitchParsimonyCalc();
			TopologyBranchPerturbations TBP = new TopologyBranchPerturbations();
			String[] singleCellNameArray = new String[nCell];
			for (int i=0; i< nCell; i++){
				singleCellNameArray[i] = singleCellNames.get(i);
			}
			
			randTree = FitchObj.findParsimonyTree(randTree1, singleCellNameArray, obsCellGenomes, parsimonyIter, nCell, TBP);
		}

		// Data Structures to store tree, error rate and likelihood values
        ArrayList<STITree<Double>> treeList = new ArrayList<STITree<Double>>();
        double[] treeLoglikelihoodArr = new double[n_iter];
        double[] fnArr = new double[n_iter];
        treeList.add(randTree);
        fnArr[0] = fnPriorMean;
        
        TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fnPriorMean, mutRate, BUF, dataFlag, modelFlag);
        treeLoglikelihoodArr[0] = jointSampleHelper.computeTreeLikelihood(obsGTObsMatrix, randTree, fnPriorMean);
        
        double bestLikelihood = Double.NEGATIVE_INFINITY;
        STITree<Double> bestTree = randTree;
        STITree<Double> previousTree = randTree;
        int bestTreeIndex = 0;
        double bestFn = fnPriorMean;
        
        for (int i = 1; i < n_iter; i++){
//        	STITree<Double> currTree = treeList.get(i-1);
        	STITree<Double> currTree = previousTree;
        	double loglikelihood_currTree = treeLoglikelihoodArr[i-1];
        	double currFn = fnArr[i-1];
        	double rr = _rng.nextDouble();
        	
        	// Propose a new error value
        	if (rr <= errorlearn){
        		double newFn = jointSampleHelper.proposeFnObj.proposeNewFn(currFn, fnPriorSD);
        		double newFnLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newFn, betaPriora, betaPriorb);
        		double newFnPriorRatio = newFnLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currFn, betaPriora, betaPriorb);
        		double newFnProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currFn, newFn, fnPriorSD);
        		double newFnLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newFn, currTree);
        		double fnAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newFnLogLikelihoodRatio, newFnPriorRatio, newFnProposalRatio);
        		boolean fnAcceptFlag = jointSampleHelper.getAcceptanceFlag(fnAcceptanceRatio);
//        		treeList.add(currTree);
        		
        		
        		// new error value accepted
        		if(fnAcceptFlag == true){
	        		fnArr[i] = newFn;
	        		treeLoglikelihoodArr[i] = jointSampleHelper.newFnSameTreeLikelihood;
        		}
        		// new error value rejected
        		else{
        			fnArr[i] = currFn;
	        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
        		}
        	}
        	// Propose a new tree
        	else{
        		STITree<Double> newTree = jointSampleHelper.TBP.proposeTree(currTree, nCell);
        		double newTreeLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newTree);
        		double treeAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioTree(newTreeLogLikelihoodRatio, jointSampleHelper.TBP.getProposalRatio());
        		boolean treeAcceptFlag = jointSampleHelper.getAcceptanceFlag(treeAcceptanceRatio);
        		fnArr[i] = currFn;
        		
        		// new tree accepted
            	if (treeAcceptFlag == true){
//            		treeList.add(newTree);
            		previousTree = newTree;
            		treeLoglikelihoodArr[i] = jointSampleHelper.sameFnNewTreeLikelihood;
            	}
            	// new tree rejected
            	else{
//            		treeList.add(currTree);
            		
            		treeLoglikelihoodArr[i] = loglikelihood_currTree;
            	}
        	}
        	if (treeLoglikelihoodArr[i] > bestLikelihood){
        		bestLikelihood = treeLoglikelihoodArr[i];
//        		bestTree = treeList.get(i);
        		bestTree = previousTree;
        		bestTreeIndex = i;
        		bestFn = fnArr[i];
        	}
        	if (i % p_iter == 0){
        		System.out.printf("After %d iterations, likelihood of the best configuration is %f\n", i, bestLikelihood);
        	}
        	
        }
        
        
        String bestTreeNewick = bestTree.toNewick();
        TreeLikelihoodObj bestConfiguration = new TreeLikelihoodObj(bestFn, bestTreeNewick, bestLikelihood);
		return bestConfiguration;
	}
	
	/**
	 * Finds the best configuration of tree and error rate and model that maximizes the likelihood.
	 * This is run for one restart of SiFit. This one learns the parameters of the complex model of evolution.
	 * @param obsGTObsMatrix
	 * @param obsGenotypeMatrix
	 * @param singleCellNames
	 * @param nCell
	 * @param nMut
	 * @param fp
	 * @param fnADO
	 * @param errorlearn
	 * @param modelParamLearn
	 * @param delProb
	 * @param LOHProb
	 * @param n_iter
	 * @param p_iter
	 * @param dataFlag
	 * @param modelFlag
	 * @param restartIndex
	 * @param BUF
	 * @return
	 * @throws SuspendableException
	 */
	public static TreeLikelihoodObj findBestTree(
			ArrayList<GenotypeObservation> obsGTObsMatrix,
			ArrayList<Integer[]> obsGenotypeMatrix,
			ArrayList<String> singleCellNames,
			int nCell, int nMut, double fp, double fnADO,
			double errorlearn, double modelParamLearn, double mhProb,
			double delProb, double LOHProb, double recurProb,
			int n_iter, int p_iter,
			int dataFlag, int modelFlag, int restartIndex,
			BasicUtilityFunctions BUF) throws SuspendableException{
		
		Random _rng = new Random();
        // fn Parameters, learn from the input value
//        double fnPriorMean = BUF.getFNPrior(obsGenotypeMatrix, fnADO);
        double fnPriorMean = fnADO;
        double fnPriorSD = fnPriorMean*0.1;
        double betaPriora = ((1 - fnPriorMean)*fnPriorMean*fnPriorMean/(fnPriorSD*fnPriorSD)) - fnPriorMean;
        double betaPriorb = betaPriora * ((1/fnPriorMean) - 1);
        
        // Parameterize prior distribution for delProb
        double delProbPriora = BUF.getBetaPriora(delProb, delProb*0.1);
        double delProbPriorb = delProbPriora * ((1/delProb) - 1);
        
        // Parameterize prior distribution for LOHProb
        double LOHProbPriora = BUF.getBetaPriora(LOHProb, LOHProb*0.1);
        double LOHProbPriorb = LOHProbPriora * ((1/LOHProb) - 1);

        // Parameterize prior distribution for recurProb
        double recurProbPriora = BUF.getBetaPriora(recurProb, recurProb*0.3); // More variance for recurProb
        double recurProbPriorb = recurProbPriora * ((1/recurProb) - 1);
        

        
        // Obtain matrix for running parsimony
        ArrayList<Integer[]> obsCellGenomes = BUF.getCellGenomes(obsGenotypeMatrix, nMut);
        
		// Generate random tree
		STITree<Double> randTreeMaker = BUF.generateRandomTree(singleCellNames);
		STITree<Double> randTree = BUF.getTree(randTreeMaker.toNewick());
		

		// For even numbered restarts, Parsimony algorithm 
		// is run to obtain a good initial tree
		FitchParsimonyCalc FitchObj = new FitchParsimonyCalc();		
		if (restartIndex % 2 == 0){
//			int parsimonyIter = n_iter/10;
			// Now running parsimony for 10000 iterations
			int parsimonyIter = 10000;
			randTree = FitchObj.findParsimonyTree(randTreeMaker.toNewick(), singleCellNames, obsCellGenomes, parsimonyIter, nCell, BUF);			
		}
		

		// Data Structures to store tree, error rate and likelihood values
        ArrayList<STITree<Double>> treeList = new ArrayList<STITree<Double>>();
        double[] treeLoglikelihoodArr = new double[n_iter];
        double[] fnArr = new double[n_iter];
        double delProbEst = delProb;
        double LOHProbEst = LOHProb;
        double recurProbEst = recurProb;
        treeList.add(randTree);
        fnArr[0] = fnPriorMean;
        
//        TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fnPriorMean, mutRate, BUF, dataFlag, modelFlag);
        TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fnPriorMean, delProb, LOHProb, recurProb, BUF, dataFlag, modelFlag);
        treeLoglikelihoodArr[0] = jointSampleHelper.computeTreeLikelihood(obsGTObsMatrix, randTree, fnPriorMean);

        
        double bestLikelihood = Double.NEGATIVE_INFINITY;
        STITree<Double> bestTree = randTree;
        STITree<Double> previousTree = randTree;
        int bestTreeIndex = 0;
        double bestFn = fnPriorMean;
        double bestDelProb = delProbEst;
        double bestLOHProb = LOHProbEst;
        double bestRecurProb = recurProbEst;
        
        
        for (int i = 1; i < n_iter; i++){
        	STITree<Double> currTree = previousTree;
        	double loglikelihood_currTree = treeLoglikelihoodArr[i-1];
        	double currFn = fnArr[i-1];
        	double rr = _rng.nextDouble();
        	double currDelProb = delProbEst;
        	double currLOHProb = LOHProbEst;
        	double currRecurProb = recurProbEst;
        	
        	// Propose a new error value
        	if (rr <= errorlearn){
        		double newFn = jointSampleHelper.proposeFnObj.proposeNewFn(currFn, fnPriorSD);
        		double randTheta = _rng.nextDouble();
        		// Hill-climbing Step
        		if (randTheta <= 1 - mhProb){
        			double newFnLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newFn, currTree);
        			// new error value accepted
        			if (newFnLogLikelihoodRatio > 0){
        				fnArr[i] = newFn;
    	        		treeLoglikelihoodArr[i] = jointSampleHelper.newFnSameTreeLikelihood;
            		}
            		// new error value rejected
            		else{
            			fnArr[i] = currFn;
    	        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}
        		else{
	        		double newFnLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newFn, betaPriora, betaPriorb);
	        		double newFnPriorRatio = newFnLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currFn, betaPriora, betaPriorb);
	        		double newFnProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currFn, newFn, fnPriorSD);
	        		double newFnLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newFn, currTree);
	        		double fnAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newFnLogLikelihoodRatio, newFnPriorRatio, newFnProposalRatio);
	        		boolean fnAcceptFlag = jointSampleHelper.getAcceptanceFlag(fnAcceptanceRatio);
		        			        		
	        		// new error value accepted
	        		if(fnAcceptFlag == true){    
		        		fnArr[i] = newFn;
		        		treeLoglikelihoodArr[i] = jointSampleHelper.newFnSameTreeLikelihood;
	        		}
	        		// new error value rejected
	        		else{
	        			fnArr[i] = currFn;
		        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
	        		}
        		}
        	}
        	// Propose a new model parameter value
        	else if (rr > errorlearn & rr <= (errorlearn + modelParamLearn)){
        		double rr2 = _rng.nextDouble();
        		fnArr[i] = currFn;
        		// Propose a new value for delProb
        		if (rr2 <= 0.5){
        			double newDelProb = jointSampleHelper.modelC.proposeNewParamVal(currDelProb, currDelProb*0.1);
//        			double newDelProbLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newDelProb, delProbPriora, delProbPriorb);
//        			double newDelProbPriorRatio = newDelProbLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currDelProb, delProbPriora, delProbPriorb);
//        			double newDelProbProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currDelProb, newDelProb, currDelProb*0.1);
        			double newDelProbLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newDelProb, currTree, 0);
//        			double newDelProbAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newDelProbLogLikelihoodRatio, newDelProbPriorRatio, newDelProbProposalRatio);
//            		boolean AcceptFlag = jointSampleHelper.getAcceptanceFlag(newDelProbAcceptanceRatio);
            		if (newDelProbLogLikelihoodRatio > 0){
            			delProbEst = newDelProb;
            			treeLoglikelihoodArr[i] = jointSampleHelper.newParamSameTreeLikelihood;
            			jointSampleHelper.modelC.delProb = newDelProb;
            		}
            		else{
            			treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}
        		// Propose a new value for LOHProb
        		else{
//        		else if (rr2 <= 0.66){
        			double newLOHProb = jointSampleHelper.modelC.proposeNewParamVal(currLOHProb, currLOHProb*0.1);
//        			double newLOHProbLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newLOHProb, LOHProbPriora, LOHProbPriorb);
//        			double newLOHProbPriorRatio = newLOHProbLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currLOHProb, LOHProbPriora, LOHProbPriorb);
//        			double newLOHProbProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currLOHProb, newLOHProb, currLOHProb*0.1);
        			double newLOHProbLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newLOHProb, currTree, 1);
//        			double newLOHProbAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newLOHProbLogLikelihoodRatio, newLOHProbPriorRatio, newLOHProbProposalRatio);
//            		boolean AcceptFlag = jointSampleHelper.getAcceptanceFlag(newLOHProbAcceptanceRatio);
            		if (newLOHProbLogLikelihoodRatio > 0){
            			LOHProbEst = newLOHProb;
            			treeLoglikelihoodArr[i] = jointSampleHelper.newParamSameTreeLikelihood;
            			jointSampleHelper.modelC.LOHProb = newLOHProb;
            		}
            		else{
            			treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}

        	}
        	// Propose a new tree
        	else{
        		STITree<Double> newTree = jointSampleHelper.TBP.proposeTree(currTree, nCell);
        		fnArr[i] = currFn;
        		double randTheta = _rng.nextDouble();
        		// Hill-climbing Step
        		if (randTheta <= 1 - mhProb){
        			double newTreeLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newTree);
        			if (newTreeLogLikelihoodRatio > 0){
        				previousTree = newTree;
                		treeLoglikelihoodArr[i] = jointSampleHelper.sameFnNewTreeLikelihood;
        			}
        			else
        				treeLoglikelihoodArr[i] = loglikelihood_currTree;
        		}
        		else{
	        		double newTreeLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newTree);
	        		double treeAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioTree(newTreeLogLikelihoodRatio, jointSampleHelper.TBP.getProposalRatio());
	        		boolean treeAcceptFlag = jointSampleHelper.getAcceptanceFlag(treeAcceptanceRatio);
	        		
	        		
	        		
	        		// new tree accepted
	            	if (treeAcceptFlag == true){
//	            	if (newTreeLogLikelihoodRatio > 0){
	
	            		previousTree = newTree;
	            		treeLoglikelihoodArr[i] = jointSampleHelper.sameFnNewTreeLikelihood;
	            	}
	            	// new tree rejected
	            	else{
	//            		treeList.add(currTree);
	            		
	            		treeLoglikelihoodArr[i] = loglikelihood_currTree;
	            	}
        		}
        	}
        	if (treeLoglikelihoodArr[i] > bestLikelihood){
        		bestLikelihood = treeLoglikelihoodArr[i];

        		bestTree = previousTree;
        		bestTreeIndex = i;
        		bestFn = fnArr[i];
        		bestDelProb = delProbEst;
        		bestLOHProb = LOHProbEst;

        	}
        	if (i % p_iter == 0){
        		System.out.printf("After %d iterations, likelihood of the best configuration is %f \n", i, bestLikelihood);
        		
        	}
        	
        }
        
        
        String bestTreeNewick = bestTree.toNewick();
        System.out.printf("best delProb = %f\t best LOHProb = %f\n ", bestDelProb, bestLOHProb);

        TreeLikelihoodObj bestConfiguration = new TreeLikelihoodObj(bestFn, bestTreeNewick, bestLikelihood, bestDelProb, bestLOHProb);
		return bestConfiguration;
	}
	
	public static TreeLikelihoodObj findBestTreeCO5(
			ArrayList<GenotypeObservation> obsGTObsMatrix,
			ArrayList<Integer[]> obsGenotypeMatrix,
			ArrayList<String> singleCellNames,
			int nCell, int nMut, double fp, double fnADO,
			double errorlearn, double modelParamLearn, double mhProb,
			double delProb, double LOHProb, double recurProb,
			int n_iter, int p_iter,
			int dataFlag, int modelFlag, int restartIndex,
			BasicUtilityFunctions BUF) throws SuspendableException{
		
		Random _rng = new Random();
        // fn Parameters, learn from the input value
//        double fnPriorMean = BUF.getFNPrior(obsGenotypeMatrix, fnADO);
        double fnPriorMean = fnADO;
        double fnPriorSD = fnPriorMean*0.1;
        double betaPriora = ((1 - fnPriorMean)*fnPriorMean*fnPriorMean/(fnPriorSD*fnPriorSD)) - fnPriorMean;
        double betaPriorb = betaPriora * ((1/fnPriorMean) - 1);
        
        // Parameterize prior distribution for delProb
        double delProbPriora = BUF.getBetaPriora(delProb, delProb*0.1);
        double delProbPriorb = delProbPriora * ((1/delProb) - 1);
        
        // Parameterize prior distribution for LOHProb
        double LOHProbPriora = BUF.getBetaPriora(LOHProb, LOHProb*0.1);
        double LOHProbPriorb = LOHProbPriora * ((1/LOHProb) - 1);

       
        

        
        // Obtain matrix for running parsimony
        ArrayList<Integer[]> obsCellGenomes = BUF.getCellGenomes(obsGenotypeMatrix, nMut);
        
		// Generate random tree
		STITree<Double> randTreeMaker = BUF.generateRandomTree(singleCellNames);
		STITree<Double> randTree = BUF.getTree(randTreeMaker.toNewick());
		

		// For even numbered restarts, Parsimony algorithm 
		// is run to obtain a good initial tree
		FitchParsimonyCalc FitchObj = new FitchParsimonyCalc();		
		if (restartIndex % 2 == 0){
//			int parsimonyIter = n_iter/10;
			// Now running parsimony for 10000 iterations
			int parsimonyIter = 10000;
			randTree = FitchObj.findParsimonyTree(randTreeMaker.toNewick(), singleCellNames, obsCellGenomes, parsimonyIter, nCell, BUF);			
		}
		

		// Data Structures to store tree, error rate and likelihood values
        ArrayList<STITree<Double>> treeList = new ArrayList<STITree<Double>>();
        double[] treeLoglikelihoodArr = new double[n_iter];
        double[] fnArr = new double[n_iter];
        double delProbEst = delProb;
        double LOHProbEst = LOHProb;
        double recurProbEst = recurProb;
        treeList.add(randTree);
        fnArr[0] = fnPriorMean;
        
//        TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fnPriorMean, mutRate, BUF, dataFlag, modelFlag);
        TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fnPriorMean, delProb, LOHProb, recurProb, BUF, dataFlag, modelFlag);
        treeLoglikelihoodArr[0] = jointSampleHelper.computeTreeLikelihood(obsGTObsMatrix, randTree, fnPriorMean);

        
        double bestLikelihood = Double.NEGATIVE_INFINITY;
        STITree<Double> bestTree = randTree;
        STITree<Double> previousTree = randTree;
        int bestTreeIndex = 0;
        double bestFn = fnPriorMean;
        double bestDelProb = delProbEst;
        double bestLOHProb = LOHProbEst;
        double bestRecurProb = recurProbEst;
        
        
        for (int i = 1; i < n_iter; i++){
        	STITree<Double> currTree = previousTree;
        	double loglikelihood_currTree = treeLoglikelihoodArr[i-1];
        	double currFn = fnArr[i-1];
        	double rr = _rng.nextDouble();
        	double currDelProb = delProbEst;
        	double currLOHProb = LOHProbEst;
        	double currRecurProb = recurProbEst;
        	
        	// Propose a new error value
        	if (rr <= errorlearn){
        		double newFn = jointSampleHelper.proposeFnObj.proposeNewFn(currFn, fnPriorSD);
        		if (newFn < 0.05){
        			fnArr[i] = currFn;
        			treeLoglikelihoodArr[i] = loglikelihood_currTree;
        			continue;
        		}
        		double randTheta = _rng.nextDouble();
        		// Hill-climbing Step
        		if (randTheta <= 1 - mhProb){
        			double newFnLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newFn, currTree);
        			// new error value accepted
        			if (newFnLogLikelihoodRatio > 0){
        				fnArr[i] = newFn;
    	        		treeLoglikelihoodArr[i] = jointSampleHelper.newFnSameTreeLikelihood;
            		}
            		// new error value rejected
            		else{
            			fnArr[i] = currFn;
    	        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}
        		else{
	        		double newFnLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newFn, betaPriora, betaPriorb);
	        		double newFnPriorRatio = newFnLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currFn, betaPriora, betaPriorb);
	        		double newFnProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currFn, newFn, fnPriorSD);
	        		double newFnLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newFn, currTree);
	        		double fnAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newFnLogLikelihoodRatio, newFnPriorRatio, newFnProposalRatio);
	        		boolean fnAcceptFlag = jointSampleHelper.getAcceptanceFlag(fnAcceptanceRatio);
		        			        		
	        		// new error value accepted
	        		if(fnAcceptFlag == true){    
		        		fnArr[i] = newFn;
		        		treeLoglikelihoodArr[i] = jointSampleHelper.newFnSameTreeLikelihood;
	        		}
	        		// new error value rejected
	        		else{
	        			fnArr[i] = currFn;
		        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
	        		}
        		}
        	}
        	// Propose a new model parameter value
        	else if (rr > errorlearn & rr <= (errorlearn + modelParamLearn)){
        		double rr2 = _rng.nextDouble();
        		fnArr[i] = currFn;
        		// Propose a new value for delProb
        		if (rr2 <= 0.5){
        			double newDelProb = jointSampleHelper.modelC.proposeNewParamVal(currDelProb, currDelProb*0.1);
//        			double newDelProbLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newDelProb, delProbPriora, delProbPriorb);
//        			double newDelProbPriorRatio = newDelProbLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currDelProb, delProbPriora, delProbPriorb);
//        			double newDelProbProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currDelProb, newDelProb, currDelProb*0.1);
        			double newDelProbLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newDelProb, currTree, 0);
//        			double newDelProbAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newDelProbLogLikelihoodRatio, newDelProbPriorRatio, newDelProbProposalRatio);
//            		boolean AcceptFlag = jointSampleHelper.getAcceptanceFlag(newDelProbAcceptanceRatio);
            		if (newDelProbLogLikelihoodRatio > 0){
            			delProbEst = newDelProb;
            			treeLoglikelihoodArr[i] = jointSampleHelper.newParamSameTreeLikelihood;
            			jointSampleHelper.modelC.delProb = newDelProb;
            		}
            		else{
            			treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}
        		// Propose a new value for LOHProb
        		else{
//        		else if (rr2 <= 0.66){
        			double newLOHProb = jointSampleHelper.modelC.proposeNewParamVal(currLOHProb, currLOHProb*0.1);
//        			double newLOHProbLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newLOHProb, LOHProbPriora, LOHProbPriorb);
//        			double newLOHProbPriorRatio = newLOHProbLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currLOHProb, LOHProbPriora, LOHProbPriorb);
//        			double newLOHProbProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currLOHProb, newLOHProb, currLOHProb*0.1);
        			double newLOHProbLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newLOHProb, currTree, 1);
//        			double newLOHProbAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newLOHProbLogLikelihoodRatio, newLOHProbPriorRatio, newLOHProbProposalRatio);
//            		boolean AcceptFlag = jointSampleHelper.getAcceptanceFlag(newLOHProbAcceptanceRatio);
            		if (newLOHProbLogLikelihoodRatio > 0){
            			LOHProbEst = newLOHProb;
            			treeLoglikelihoodArr[i] = jointSampleHelper.newParamSameTreeLikelihood;
            			jointSampleHelper.modelC.LOHProb = newLOHProb;
            		}
            		else{
            			treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}

        	}
        	// Propose a new tree
        	else{
        		STITree<Double> newTree = jointSampleHelper.TBP.proposeTree(currTree, nCell);
        		fnArr[i] = currFn;
        		double randTheta = _rng.nextDouble();
        		// Hill-climbing Step
        		if (randTheta <= 1 - mhProb){
        			double newTreeLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newTree);
        			if (newTreeLogLikelihoodRatio > 0){
        				previousTree = newTree;
                		treeLoglikelihoodArr[i] = jointSampleHelper.sameFnNewTreeLikelihood;
        			}
        			else
        				treeLoglikelihoodArr[i] = loglikelihood_currTree;
        		}
        		else{
	        		double newTreeLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newTree);
	        		double treeAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioTree(newTreeLogLikelihoodRatio, jointSampleHelper.TBP.getProposalRatio());
	        		boolean treeAcceptFlag = jointSampleHelper.getAcceptanceFlag(treeAcceptanceRatio);
	        		
	        		
	        		
	        		// new tree accepted
	            	if (treeAcceptFlag == true){
//	            	if (newTreeLogLikelihoodRatio > 0){
	
	            		previousTree = newTree;
	            		treeLoglikelihoodArr[i] = jointSampleHelper.sameFnNewTreeLikelihood;
	            	}
	            	// new tree rejected
	            	else{
	//            		treeList.add(currTree);
	            		
	            		treeLoglikelihoodArr[i] = loglikelihood_currTree;
	            	}
        		}
        	}
        	if (treeLoglikelihoodArr[i] > bestLikelihood){
        		bestLikelihood = treeLoglikelihoodArr[i];

        		bestTree = previousTree;
        		bestTreeIndex = i;
        		bestFn = fnArr[i];
        		bestDelProb = delProbEst;
        		bestLOHProb = LOHProbEst;

        	}
        	if (i % p_iter == 0){
        		System.out.printf("After %d iterations, likelihood of the best configuration is %f \n", i, bestLikelihood);
        		
        	}
        	
        }
        
        
        String bestTreeNewick = bestTree.toNewick();
        System.out.printf("best delProb = %f\t best LOHProb = %f\n ", bestDelProb, bestLOHProb);

        TreeLikelihoodObj bestConfiguration = new TreeLikelihoodObj(bestFn, bestTreeNewick, bestLikelihood, bestDelProb, bestLOHProb);
		return bestConfiguration;
	}
	
	public static TreeLikelihoodObj findBestTreeCRC0827(
			ArrayList<GenotypeObservation> obsGTObsMatrix,
			ArrayList<Integer[]> obsGenotypeMatrix,
			ArrayList<String> singleCellNames,
			int nCell, int nMut, double fp, double fnADO,
			double errorlearn, double modelParamLearn, double mhProb,
			double delProb, double LOHProb, double recurProb,
			int n_iter, int p_iter,
			int dataFlag, int modelFlag, int restartIndex,
			BasicUtilityFunctions BUF) throws SuspendableException{
		
		Random _rng = new Random();
        // fn Parameters, learn from the input value
//        double fnPriorMean = BUF.getFNPrior(obsGenotypeMatrix, fnADO);
        double fnPriorMean = fnADO;
        double fnPriorSD = fnPriorMean*0.1;
        double betaPriora = ((1 - fnPriorMean)*fnPriorMean*fnPriorMean/(fnPriorSD*fnPriorSD)) - fnPriorMean;
        double betaPriorb = betaPriora * ((1/fnPriorMean) - 1);
        
        // Parameterize prior distribution for delProb
        double delProbPriora = BUF.getBetaPriora(delProb, delProb*0.1);
        double delProbPriorb = delProbPriora * ((1/delProb) - 1);
        
        // Parameterize prior distribution for LOHProb
        double LOHProbPriora = BUF.getBetaPriora(LOHProb, LOHProb*0.1);
        double LOHProbPriorb = LOHProbPriora * ((1/LOHProb) - 1);

        // Parameterize prior distribution for recurProb
        double recurProbPriora = BUF.getBetaPriora(recurProb, recurProb*0.3); // More variance for recurProb
        double recurProbPriorb = recurProbPriora * ((1/recurProb) - 1);
        

        
        // Obtain matrix for running parsimony
        ArrayList<Integer[]> obsCellGenomes = BUF.getCellGenomes(obsGenotypeMatrix, nMut);
        
		// Generate random tree
		STITree<Double> randTreeMaker = BUF.generateRandomTree(singleCellNames);
		STITree<Double> randTree = BUF.getTree(randTreeMaker.toNewick());
		

		// For even numbered restarts, Parsimony algorithm 
		// is run to obtain a good initial tree
		FitchParsimonyCalc FitchObj = new FitchParsimonyCalc();		
		if (restartIndex % 2 == 0){
//			int parsimonyIter = n_iter/10;
			// Now running parsimony for 10000 iterations
			int parsimonyIter = 10000;
			randTree = FitchObj.findParsimonyTree(randTreeMaker.toNewick(), singleCellNames, obsCellGenomes, parsimonyIter, nCell, BUF);			
		}
		

		// Data Structures to store tree, error rate and likelihood values
        ArrayList<STITree<Double>> treeList = new ArrayList<STITree<Double>>();
        double[] treeLoglikelihoodArr = new double[n_iter];
        double[] fnArr = new double[n_iter];
        double delProbEst = delProb;
        double LOHProbEst = LOHProb;
        double recurProbEst = recurProb;
        treeList.add(randTree);
        fnArr[0] = fnPriorMean;
        
//        TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fnPriorMean, mutRate, BUF, dataFlag, modelFlag);
        TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fnPriorMean, delProb, LOHProb, recurProb, BUF, dataFlag, modelFlag);
        treeLoglikelihoodArr[0] = jointSampleHelper.computeTreeLikelihood(obsGTObsMatrix, randTree, fnPriorMean);

        
        double bestLikelihood = Double.NEGATIVE_INFINITY;
        STITree<Double> bestTree = randTree;
        STITree<Double> previousTree = randTree;
        int bestTreeIndex = 0;
        double bestFn = fnPriorMean;
        double bestDelProb = delProbEst;
        double bestLOHProb = LOHProbEst;
        double bestRecurProb = recurProbEst;
        
        
        for (int i = 1; i < n_iter; i++){
        	STITree<Double> currTree = previousTree;
        	double loglikelihood_currTree = treeLoglikelihoodArr[i-1];
        	double currFn = fnArr[i-1];
        	double rr = _rng.nextDouble();
        	double currDelProb = delProbEst;
        	double currLOHProb = LOHProbEst;
        	double currRecurProb = recurProbEst;
        	
        	// Propose a new error value
        	if (rr <= errorlearn){
        		double newFn = jointSampleHelper.proposeFnObj.proposeNewFn(currFn, fnPriorSD);
//        		double randTheta = _rng.nextDouble();
        		double newFnLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newFn, currTree);
    			// new error value accepted
    			if (newFnLogLikelihoodRatio > 0){
    				fnArr[i] = newFn;
	        		treeLoglikelihoodArr[i] = jointSampleHelper.newFnSameTreeLikelihood;
        		}
        		// new error value rejected
        		else{
        			fnArr[i] = currFn;
	        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
        		}
        		// Hill-climbing Step
//        		if (randTheta <= 1 - mhProb){
//        			double newFnLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newFn, currTree);
//        			// new error value accepted
//        			if (newFnLogLikelihoodRatio > 0){
//        				fnArr[i] = newFn;
//    	        		treeLoglikelihoodArr[i] = jointSampleHelper.newFnSameTreeLikelihood;
//            		}
//            		// new error value rejected
//            		else{
//            			fnArr[i] = currFn;
//    	        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
//            		}
//        		}
//        		else{
//	        		double newFnLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newFn, betaPriora, betaPriorb);
//	        		double newFnPriorRatio = newFnLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currFn, betaPriora, betaPriorb);
//	        		double newFnProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currFn, newFn, fnPriorSD);
//	        		double newFnLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newFn, currTree);
//	        		double fnAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newFnLogLikelihoodRatio, newFnPriorRatio, newFnProposalRatio);
//	        		boolean fnAcceptFlag = jointSampleHelper.getAcceptanceFlag(fnAcceptanceRatio);
//		        			        		
//	        		// new error value accepted
//	        		if(fnAcceptFlag == true){    
//		        		fnArr[i] = newFn;
//		        		treeLoglikelihoodArr[i] = jointSampleHelper.newFnSameTreeLikelihood;
//	        		}
//	        		// new error value rejected
//	        		else{
//	        			fnArr[i] = currFn;
//		        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
//	        		}
//        		}
        	}
        	// Propose a new model parameter value
        	else if (rr > errorlearn & rr <= (errorlearn + modelParamLearn)){
        		double rr2 = _rng.nextDouble();
        		fnArr[i] = currFn;
        		// Propose a new value for delProb
        		if (rr2 <= 0.5){
        			double newDelProb = jointSampleHelper.modelC.proposeNewParamVal(currDelProb, currDelProb*0.1);
//        			double newDelProbLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newDelProb, delProbPriora, delProbPriorb);
//        			double newDelProbPriorRatio = newDelProbLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currDelProb, delProbPriora, delProbPriorb);
//        			double newDelProbProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currDelProb, newDelProb, currDelProb*0.1);
        			double newDelProbLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newDelProb, currTree, 0);
//        			double newDelProbAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newDelProbLogLikelihoodRatio, newDelProbPriorRatio, newDelProbProposalRatio);
//            		boolean AcceptFlag = jointSampleHelper.getAcceptanceFlag(newDelProbAcceptanceRatio);
            		if (newDelProbLogLikelihoodRatio > 0){
            			delProbEst = newDelProb;
            			treeLoglikelihoodArr[i] = jointSampleHelper.newParamSameTreeLikelihood;
            			jointSampleHelper.modelC.delProb = newDelProb;
            		}
            		else{
            			treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}
        		// Propose a new value for LOHProb
        		else{
//        		else if (rr2 <= 0.66){
        			double newLOHProb = jointSampleHelper.modelC.proposeNewParamVal(currLOHProb, currLOHProb*0.1);
//        			double newLOHProbLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newLOHProb, LOHProbPriora, LOHProbPriorb);
//        			double newLOHProbPriorRatio = newLOHProbLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currLOHProb, LOHProbPriora, LOHProbPriorb);
//        			double newLOHProbProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currLOHProb, newLOHProb, currLOHProb*0.1);
        			double newLOHProbLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newLOHProb, currTree, 1);
//        			double newLOHProbAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newLOHProbLogLikelihoodRatio, newLOHProbPriorRatio, newLOHProbProposalRatio);
//            		boolean AcceptFlag = jointSampleHelper.getAcceptanceFlag(newLOHProbAcceptanceRatio);
            		if (newLOHProbLogLikelihoodRatio > 0){
            			LOHProbEst = newLOHProb;
            			treeLoglikelihoodArr[i] = jointSampleHelper.newParamSameTreeLikelihood;
            			jointSampleHelper.modelC.LOHProb = newLOHProb;
            		}
            		else{
            			treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}

        	}
        	// Propose a new tree
        	else{
        		STITree<Double> newTree = jointSampleHelper.TBP.proposeTree(currTree, nCell);
        		fnArr[i] = currFn;
//        		double randTheta = _rng.nextDouble();
        		double newTreeLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newTree);
    			if (newTreeLogLikelihoodRatio > 0){
    				previousTree = newTree;
            		treeLoglikelihoodArr[i] = jointSampleHelper.sameFnNewTreeLikelihood;
    			}
    			else
    				treeLoglikelihoodArr[i] = loglikelihood_currTree;
        		// Hill-climbing Step
//        		if (randTheta <= 1 - mhProb){
//        			double newTreeLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newTree);
//        			if (newTreeLogLikelihoodRatio > 0){
//        				previousTree = newTree;
//                		treeLoglikelihoodArr[i] = jointSampleHelper.sameFnNewTreeLikelihood;
//        			}
//        			else
//        				treeLoglikelihoodArr[i] = loglikelihood_currTree;
//        		}
//        		else{
//	        		double newTreeLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newTree);
//	        		double treeAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioTree(newTreeLogLikelihoodRatio, jointSampleHelper.TBP.getProposalRatio());
//	        		boolean treeAcceptFlag = jointSampleHelper.getAcceptanceFlag(treeAcceptanceRatio);
//	        		
//	        		
//	        		
//	        		// new tree accepted
//	            	if (treeAcceptFlag == true){
////	            	if (newTreeLogLikelihoodRatio > 0){
//	
//	            		previousTree = newTree;
//	            		treeLoglikelihoodArr[i] = jointSampleHelper.sameFnNewTreeLikelihood;
//	            	}
//	            	// new tree rejected
//	            	else{
//	//            		treeList.add(currTree);
//	            		
//	            		treeLoglikelihoodArr[i] = loglikelihood_currTree;
//	            	}
//        		}
        	}
        	if (treeLoglikelihoodArr[i] > bestLikelihood){
        		bestLikelihood = treeLoglikelihoodArr[i];

        		bestTree = previousTree;
        		bestTreeIndex = i;
        		bestFn = fnArr[i];
        		bestDelProb = delProbEst;
        		bestLOHProb = LOHProbEst;

        	}
        	if (i % p_iter == 0){
        		System.out.printf("After %d iterations, likelihood of the best configuration is %f \n", i, bestLikelihood);
        		
        	}
        	
        }
        
        
        String bestTreeNewick = bestTree.toNewick();
        System.out.printf("best delProb = %f\t best LOHProb = %f\n ", bestDelProb, bestLOHProb);

        TreeLikelihoodObj bestConfiguration = new TreeLikelihoodObj(bestFn, bestTreeNewick, bestLikelihood, bestDelProb, bestLOHProb);
		return bestConfiguration;
	}
	
	public static PriorityQueue<TreeLikelihoodObj> findTopBestTrees(
							ArrayList<GenotypeObservation> obsGTObsMatrix,
							ArrayList<Integer[]> obsGenotypeMatrix,
							ArrayList<String> singleCellNames,
							String goodtree,
							int nCell, int nMut, double fp, double fnADO,
							double errorlearn, double modelParamLearn, 
							double delProb, double LOHProb, double recurProb,
							int n_iter, int p_iter,
							int dataFlag, int modelFlag, int restartIndex,
							BasicUtilityFunctions BUF) throws SuspendableException{
		
		Random _rng = new Random();
        // fn Parameters, learn from the input value
//        double fnPriorMean = BUF.getFNPrior(obsGenotypeMatrix, fnADO);
        double fnPriorMean = fnADO;
        double fnPriorSD = fnPriorMean*0.1;
        double betaPriora = ((1 - fnPriorMean)*fnPriorMean*fnPriorMean/(fnPriorSD*fnPriorSD)) - fnPriorMean;
        double betaPriorb = betaPriora * ((1/fnPriorMean) - 1);
        
        // Parameterize prior distribution for delProb
        double delProbPriora = BUF.getBetaPriora(delProb, delProb*0.1);
        double delProbPriorb = delProbPriora * ((1/delProb) - 1);
        
        // Parameterize prior distribution for LOHProb
        double LOHProbPriora = BUF.getBetaPriora(LOHProb, LOHProb*0.1);
        double LOHProbPriorb = LOHProbPriora * ((1/LOHProb) - 1);
        
		// Generate random tree
//		STITree<Double> randTreeMaker = BUF.generateRandomTree(singleCellNames);
//		STITree<Double> randTree = BUF.getTree(randTreeMaker.toNewick());
		
                
		STITree<Double> randTree = BUF.getTree(goodtree);
		
		// Data Structures to store tree, error rate and likelihood values
        ArrayList<STITree<Double>> treeList = new ArrayList<STITree<Double>>();
        double[] treeLoglikelihoodArr = new double[n_iter];
        double[] fnArr = new double[n_iter];
        double delProbEst = delProb;
        double LOHProbEst = LOHProb;
        double recurProbEst = recurProb;
        treeList.add(randTree);
        fnArr[0] = fnPriorMean;
        // Comparator for using with the priority queue
        Comparator<TreeLikelihoodObj> comparator = new Comparator<TreeLikelihoodObj>() {
			
			@Override
			public int compare(TreeLikelihoodObj o1, TreeLikelihoodObj o2) {
				// TODO Auto-generated method stub
				if (o1.bestTreeLikelihood > o2.bestTreeLikelihood)
					return 1;
				else if (o1.bestTreeLikelihood == o2.bestTreeLikelihood)
					return 0;
				else 
					return -1;
			}
		};
		// This is the list of best 100 configurations
        PriorityQueue<TreeLikelihoodObj> topTreeLikelihoodConfigs = 
	            new PriorityQueue<TreeLikelihoodObj>(100, comparator);
        
        TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fnPriorMean, delProb, LOHProb, recurProb, BUF, dataFlag, modelFlag);
        treeLoglikelihoodArr[0] = jointSampleHelper.computeTreeLikelihood(obsGTObsMatrix, randTree, fnPriorMean);
        System.out.println(treeLoglikelihoodArr[0]);
        
        TreeLikelihoodObj startingTL = new TreeLikelihoodObj(fnPriorMean, randTree.toNewick(), treeLoglikelihoodArr[0], delProbEst, LOHProbEst);
        topTreeLikelihoodConfigs.add(startingTL);
        double bestLikelihood = Double.NEGATIVE_INFINITY;
        STITree<Double> bestTree = randTree;
        STITree<Double> previousTree = randTree;
        int bestTreeIndex = 0;
        double bestFn = fnPriorMean;
        double bestDelProb = delProbEst;
        double bestLOHProb = LOHProbEst;
        double bestRecurProb = recurProbEst;
        
        for (int i = 1; i < n_iter; i++){
        	STITree<Double> currTree = previousTree;
        	double loglikelihood_currTree = treeLoglikelihoodArr[i-1];
        	double currFn = fnArr[i-1];
        	double rr = _rng.nextDouble();
        	double currDelProb = delProbEst;
        	double currLOHProb = LOHProbEst;

        	// Propose a new error value
        	if (rr <= errorlearn){
        		double newFn = jointSampleHelper.proposeFnObj.proposeNewFn(currFn, fnPriorSD);
//        		if (newFn < 0.05){
//        			fnArr[i] = currFn;
//	        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
//	        		continue;
//        		}
        		double newFnLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newFn, betaPriora, betaPriorb);
        		double newFnPriorRatio = newFnLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currFn, betaPriora, betaPriorb);
        		double newFnProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currFn, newFn, fnPriorSD);
        		double newFnLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newFn, currTree);
        		double fnAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newFnLogLikelihoodRatio, newFnPriorRatio, newFnProposalRatio);
        		boolean fnAcceptFlag = jointSampleHelper.getAcceptanceFlag(fnAcceptanceRatio);
//        		treeList.add(currTree);
        		
        		
        		// new error value accepted
//        		if(fnAcceptFlag == true){    
        		if (newFnLogLikelihoodRatio > 0){
//        			System.out.printf("nu fn accepted in iteration %d\n", i);
//        			System.out.println(jointSampleHelper.newFnSameTreeLikelihood);
	        		fnArr[i] = newFn;
	        		treeLoglikelihoodArr[i] = jointSampleHelper.newFnSameTreeLikelihood;
        		}
        		// new error value rejected
        		else{
        			fnArr[i] = currFn;
	        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
        		}
        	}
        	// Propose a new model parameter value
        	else if (rr > errorlearn & rr <= (errorlearn + modelParamLearn)){
        		double rr2 = _rng.nextDouble();
        		fnArr[i] = currFn;
        		if (rr2 <= 0.5){
        			double newDelProb = jointSampleHelper.modelC.proposeNewParamVal(currDelProb, currDelProb*0.1);
        			double newDelProbLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newDelProb, delProbPriora, delProbPriorb);
        			double newDelProbPriorRatio = newDelProbLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currDelProb, delProbPriora, delProbPriorb);
        			double newDelProbProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currDelProb, newDelProb, currDelProb*0.1);
        			double newDelProbLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newDelProb, currTree, 0);
        			double newDelProbAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newDelProbLogLikelihoodRatio, newDelProbPriorRatio, newDelProbProposalRatio);
            		boolean AcceptFlag = jointSampleHelper.getAcceptanceFlag(newDelProbAcceptanceRatio);
            		if (newDelProbLogLikelihoodRatio > 0){
//            		if (AcceptFlag == true){
            			delProbEst = newDelProb;
            			treeLoglikelihoodArr[i] = jointSampleHelper.newParamSameTreeLikelihood;
            			jointSampleHelper.modelC.delProb = newDelProb;
            		}
            		else{
            			treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}
        		// Propose a new value for LOHProb
        		else{
        			double newLOHProb = jointSampleHelper.modelC.proposeNewParamVal(currLOHProb, currLOHProb*0.1);
        			double newLOHProbLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newLOHProb, LOHProbPriora, LOHProbPriorb);
        			double newLOHProbPriorRatio = newLOHProbLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currLOHProb, LOHProbPriora, LOHProbPriorb);
        			double newLOHProbProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currLOHProb, newLOHProb, currLOHProb*0.1);
        			double newLOHProbLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newLOHProb, currTree, 1);
        			double newLOHProbAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newLOHProbLogLikelihoodRatio, newLOHProbPriorRatio, newLOHProbProposalRatio);
            		boolean AcceptFlag = jointSampleHelper.getAcceptanceFlag(newLOHProbAcceptanceRatio);
            		if (newLOHProbLogLikelihoodRatio > 0){
//            		if (AcceptFlag == true){
            			LOHProbEst = newLOHProb;
            			treeLoglikelihoodArr[i] = jointSampleHelper.newParamSameTreeLikelihood;
            			jointSampleHelper.modelC.LOHProb = newLOHProb;
            		}
            		else{
            			treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}
        	}
        	// Propose a new tree
        	else{
        		STITree<Double> newTree = jointSampleHelper.TBP.proposeTree(currTree, nCell);
//        		STITree<Double> newTree = jointSampleHelper.TBP.changeBranchLength(currTree);
        		double newTreeLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newTree);
//        		double treeAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioTree(newTreeLogLikelihoodRatio, jointSampleHelper.TBP.getProposalRatio());
//        		boolean treeAcceptFlag = jointSampleHelper.getAcceptanceFlag(treeAcceptanceRatio);
        		
        		fnArr[i] = currFn;
        		
        		// new tree accepted
//            	if (treeAcceptFlag == true){
            	if (newTreeLogLikelihoodRatio > 0){
//            		System.out.printf("nu tree accepted in iteration %d\n", i);
//        			System.out.println(jointSampleHelper.sameFnNewTreeLikelihood);
//            		treeList.add(newTree);
            		previousTree = newTree;
            		treeLoglikelihoodArr[i] = jointSampleHelper.sameFnNewTreeLikelihood;
            	}
            	// new tree rejected
            	else{
//            		treeList.add(currTree);
            		
            		treeLoglikelihoodArr[i] = loglikelihood_currTree;
            	}
        	}
        	TreeLikelihoodObj currentTLBound = topTreeLikelihoodConfigs.peek();
        	if (treeLoglikelihoodArr[i] > currentTLBound.bestTreeLikelihood){
        		bestLikelihood = treeLoglikelihoodArr[i];
//        		bestTree = treeList.get(i);
        		bestTree = previousTree;
        		bestTreeIndex = i;
        		bestFn = fnArr[i];
        		bestDelProb = delProbEst;
        		bestLOHProb = LOHProbEst;
        		TreeLikelihoodObj newTLObj = new TreeLikelihoodObj(bestFn, bestTree.toNewick(), bestLikelihood, bestDelProb, bestLOHProb);
        		topTreeLikelihoodConfigs.add(newTLObj);
        		if (topTreeLikelihoodConfigs.size() > 100)
        			topTreeLikelihoodConfigs.remove();
//        		bestRecurProb = recurProbEst;
        	}
        	if (i % p_iter == 0){
        		System.out.printf("After %d iterations, likelihood of the best configuration is %f\n", i, bestLikelihood);
        	}
        }
        
        return topTreeLikelihoodConfigs;
	}
	
	public static TreeLikelihoodObj findBestTreeGivenTree(
			ArrayList<GenotypeObservation> obsGTObsMatrix,
			ArrayList<Integer[]> obsGenotypeMatrix,
			ArrayList<String> singleCellNames,
			String goodtree,
			int nCell, int nMut, double fp, double fnADO,
			double errorlearn, double modelParamLearn, 
			double delProb, double LOHProb, double recurProb,
			int n_iter, int p_iter,
			int dataFlag, int modelFlag, int restartIndex,
			BasicUtilityFunctions BUF) throws SuspendableException{
		
		Random _rng = new Random();
        // fn Parameters, learn from the input value
//        double fnPriorMean = BUF.getFNPrior(obsGenotypeMatrix, fnADO);
        double fnPriorMean = fnADO;
        double fnPriorSD = fnPriorMean*0.5;
        double betaPriora = ((1 - fnPriorMean)*fnPriorMean*fnPriorMean/(fnPriorSD*fnPriorSD)) - fnPriorMean;
        double betaPriorb = betaPriora * ((1/fnPriorMean) - 1);
        
        // Parameterize prior distribution for delProb
        double delProbPriora = BUF.getBetaPriora(delProb, delProb*0.1);
        double delProbPriorb = delProbPriora * ((1/delProb) - 1);
        
        // Parameterize prior distribution for LOHProb
        double LOHProbPriora = BUF.getBetaPriora(LOHProb, LOHProb*0.1);
        double LOHProbPriorb = LOHProbPriora * ((1/LOHProb) - 1);

        // Parameterize prior distribution for recurProb
        double recurProbPriora = BUF.getBetaPriora(recurProb, recurProb*0.3); // More variance for recurProb
        double recurProbPriorb = recurProbPriora * ((1/recurProb) - 1);
        
        // Obtain matrix for running parsimony
        ArrayList<Integer[]> obsCellGenomes = BUF.getCellGenomes(obsGenotypeMatrix, nMut);
        
		// Generate random tree
		STITree<Double> randTreeMaker = BUF.generateRandomTree(singleCellNames);
//		STITree<Double> randTree = BUF.getTree(randTreeMaker.toNewick());
		
                
		STITree<Double> randTree = BUF.getTree(goodtree);
		
//		BUF.resizeBranchLengths(randTree, 1);
		// For even numbered restarts, Parsimony algorithm 
		// is run to obtain a good initial tree
//		if (restartIndex % 2 == 0){
//		STITree<Double> randTree1 = BUF.getTree(goodtree);
////			STITree<Double> randTree1 = BUF.getTree(randTreeMaker.toNewick());
////			int parsimonyIter = n_iter/10;
//			// Now running parsimony for 10000 iterations
//			int parsimonyIter = 10000;
//			
//			// Run Parsimony to get the initial tree
//			FitchParsimonyCalc FitchObj = new FitchParsimonyCalc();
//			TopologyBranchPerturbations TBP = new TopologyBranchPerturbations();
//			String[] singleCellNameArray = new String[nCell];
//			for (int i=0; i< nCell; i++){
//				singleCellNameArray[i] = singleCellNames.get(i);
//			}
//			
//			STITree<Double> randTree = FitchObj.findParsimonyTree(randTree1, singleCellNameArray, obsCellGenomes, parsimonyIter, nCell, TBP);
//		}

		// Data Structures to store tree, error rate and likelihood values
        ArrayList<STITree<Double>> treeList = new ArrayList<STITree<Double>>();
        double[] treeLoglikelihoodArr = new double[n_iter];
        double[] fnArr = new double[n_iter];
        double delProbEst = delProb;
        double LOHProbEst = LOHProb;
        double recurProbEst = recurProb;
        treeList.add(randTree);
        fnArr[0] = fnPriorMean;
        
//        TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fnPriorMean, mutRate, BUF, dataFlag, modelFlag);
        TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fnPriorMean, delProb, LOHProb, recurProb, BUF, dataFlag, modelFlag);
        treeLoglikelihoodArr[0] = jointSampleHelper.computeTreeLikelihood(obsGTObsMatrix, randTree, fnPriorMean);
//        System.out.println(treeLoglikelihoodArr[0]);
        
        double bestLikelihood = treeLoglikelihoodArr[0]; //Double.NEGATIVE_INFINITY;
        STITree<Double> bestTree = randTree;
        STITree<Double> previousTree = randTree;
        int bestTreeIndex = 0;
        double bestFn = fnPriorMean;
        double bestDelProb = delProbEst;
        double bestLOHProb = LOHProbEst;
        double bestRecurProb = recurProbEst;
        
        
        for (int i = 1; i < n_iter; i++){
//        	STITree<Double> currTree = treeList.get(i-1);
        	STITree<Double> currTree = previousTree;
        	double loglikelihood_currTree = treeLoglikelihoodArr[i-1];
        	double currFn = fnArr[i-1];
        	double rr = _rng.nextDouble();
        	double currDelProb = delProbEst;
        	double currLOHProb = LOHProbEst;
        	double currRecurProb = recurProbEst;
        	
        	// Propose a new error value
        	if (rr <= errorlearn){
        		double newFn = jointSampleHelper.proposeFnObj.proposeNewFn(currFn, fnPriorSD);
        		if (newFn < 0.05){
        			fnArr[i] = currFn;
	        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
	        		continue;
        		}
        		double newFnLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newFn, betaPriora, betaPriorb);
        		double newFnPriorRatio = newFnLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currFn, betaPriora, betaPriorb);
        		double newFnProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currFn, newFn, fnPriorSD);
        		double newFnLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newFn, currTree);
        		double fnAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newFnLogLikelihoodRatio, newFnPriorRatio, newFnProposalRatio);
        		boolean fnAcceptFlag = jointSampleHelper.getAcceptanceFlag(fnAcceptanceRatio);
//        		treeList.add(currTree);
        		
        		
        		// new error value accepted
//        		if(fnAcceptFlag == true){
        		if (newFnLogLikelihoodRatio > 0){	
	        		fnArr[i] = newFn;
	        		treeLoglikelihoodArr[i] = jointSampleHelper.newFnSameTreeLikelihood;
        		}
        		// new error value rejected
        		else{
        			fnArr[i] = currFn;
	        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
        		}
        	}
        	// Propose a new model parameter value
        	else if (rr > errorlearn & rr <= (errorlearn + modelParamLearn)){
        		double rr2 = _rng.nextDouble();
        		fnArr[i] = currFn;
        		// Propose a new value for delProb
        		if (rr2 <= 0.5){
        			double newDelProb = jointSampleHelper.modelC.proposeNewParamVal(currDelProb, currDelProb*0.1);
        			double newDelProbLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newDelProb, delProbPriora, delProbPriorb);
        			double newDelProbPriorRatio = newDelProbLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currDelProb, delProbPriora, delProbPriorb);
        			double newDelProbProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currDelProb, newDelProb, currDelProb*0.1);
        			double newDelProbLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newDelProb, currTree, 0);
        			double newDelProbAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newDelProbLogLikelihoodRatio, newDelProbPriorRatio, newDelProbProposalRatio);
            		boolean AcceptFlag = jointSampleHelper.getAcceptanceFlag(newDelProbAcceptanceRatio);
            		if (AcceptFlag == true){
            			delProbEst = newDelProb;
            			treeLoglikelihoodArr[i] = jointSampleHelper.newParamSameTreeLikelihood;
            			jointSampleHelper.modelC.delProb = newDelProb;
            		}
            		else{
            			treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}
        		// Propose a new value for LOHProb
        		else{
//        		else if (rr2 <= 0.66){
        			double newLOHProb = jointSampleHelper.modelC.proposeNewParamVal(currLOHProb, currLOHProb*0.1);
        			double newLOHProbLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newLOHProb, LOHProbPriora, LOHProbPriorb);
        			double newLOHProbPriorRatio = newLOHProbLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currLOHProb, LOHProbPriora, LOHProbPriorb);
        			double newLOHProbProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currLOHProb, newLOHProb, currLOHProb*0.1);
        			double newLOHProbLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newLOHProb, currTree, 1);
        			double newLOHProbAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newLOHProbLogLikelihoodRatio, newLOHProbPriorRatio, newLOHProbProposalRatio);
            		boolean AcceptFlag = jointSampleHelper.getAcceptanceFlag(newLOHProbAcceptanceRatio);
            		if (AcceptFlag == true){
            			LOHProbEst = newLOHProb;
            			treeLoglikelihoodArr[i] = jointSampleHelper.newParamSameTreeLikelihood;
            			jointSampleHelper.modelC.LOHProb = newLOHProb;
            		}
            		else{
            			treeLoglikelihoodArr[i] = loglikelihood_currTree;
            		}
        		}
        		// Propose a new value for recurProb
//        		else{
//        			double newRecurProb = jointSampleHelper.modelC.proposeNewParamVal(currRecurProb, currRecurProb*0.3);
//        			double newRecurProbLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newRecurProb, recurProbPriora, recurProbPriorb);
//        			double newRecurProbPriorRatio = newRecurProbLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currRecurProb, recurProbPriora, recurProbPriorb);
//        			double newRecurProbProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currRecurProb, newRecurProb, currRecurProb*0.5);
//        			double newRecurProbLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newRecurProb, currTree, 2);
//        			double newRecurProbAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newRecurProbLogLikelihoodRatio, newRecurProbPriorRatio, newRecurProbProposalRatio);
//            		boolean AcceptFlag = jointSampleHelper.getAcceptanceFlag(newRecurProbAcceptanceRatio);
//            		if (AcceptFlag == true){
//            			recurProbEst = newRecurProb;
//            			treeLoglikelihoodArr[i] = jointSampleHelper.newParamSameTreeLikelihood;
//            			jointSampleHelper.modelC.recurProb = newRecurProb;
//            		}
//            		else{
//            			treeLoglikelihoodArr[i] = loglikelihood_currTree;
//            		}
//        		}
        	}
        	// Propose a new tree
        	else{
//        		STITree<Double> newTree = jointSampleHelper.TBP.proposeTree(currTree, nCell);
        		STITree<Double> newTree = jointSampleHelper.TBP.changeBranchLength(currTree);
        		double newTreeLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newTree);
        		double treeAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioTree(newTreeLogLikelihoodRatio, jointSampleHelper.TBP.getProposalRatio());
        		boolean treeAcceptFlag = jointSampleHelper.getAcceptanceFlag(treeAcceptanceRatio);
        		fnArr[i] = currFn;
        		
        		// new tree accepted
//            	if (treeAcceptFlag == true){
            	if (newTreeLogLikelihoodRatio > 0){
//            		System.out.printf("nu tree accepted in iteration %d\n", i);
//        			System.out.println(jointSampleHelper.sameFnNewTreeLikelihood);
//            		treeList.add(newTree);
            		previousTree = newTree;
            		treeLoglikelihoodArr[i] = jointSampleHelper.sameFnNewTreeLikelihood;
            	}
            	// new tree rejected
            	else{
//            		treeList.add(currTree);
            		
            		treeLoglikelihoodArr[i] = loglikelihood_currTree;
            	}
        	}
        	if (treeLoglikelihoodArr[i] > bestLikelihood){
        		bestLikelihood = treeLoglikelihoodArr[i];
//        		bestTree = treeList.get(i);
        		bestTree = previousTree;
        		bestTreeIndex = i;
        		bestFn = fnArr[i];
        		bestDelProb = delProbEst;
        		bestLOHProb = LOHProbEst;
//        		bestRecurProb = recurProbEst;
        	}
        	if (i % p_iter == 0){
        		System.out.printf("After %d iterations, likelihood of the best configuration is %f\n", i, bestLikelihood);
        	}
        	
        }
        
        
        String bestTreeNewick = bestTree.toNewick();
//        System.out.printf("best delProb = %f\t best LOHProb = %f\t best recurProb = %f\n ", bestDelProb, bestLOHProb, bestRecurProb);
//        TreeLikelihoodObj bestConfiguration = new TreeLikelihoodObj(bestFn, bestTreeNewick, bestLikelihood);
        TreeLikelihoodObj bestConfiguration = new TreeLikelihoodObj(bestFn, bestTreeNewick, bestLikelihood, bestDelProb, bestLOHProb);
		return bestConfiguration;
	}
	
	/**
	 * This one takes doublet information into account, not being used now.
	 * @param obsGTObsMatrix
	 * @param obsGenotypeMatrix
	 * @param singleCellNames
	 * @param doubletFlagMap
	 * @param nCell
	 * @param nMut
	 * @param fp
	 * @param fnADO
	 * @param errorlearn
	 * @param mutRate
	 * @param n_iter
	 * @param p_iter
	 * @param dataFlag
	 * @param modelFlag
	 * @param BUF
	 * @return
	 * @throws SuspendableException
	 */
	public static TreeLikelihoodObj findBestTree(
			ArrayList<GenotypeObservation> obsGTObsMatrix,
			ArrayList<Integer[]> obsGenotypeMatrix,
			ArrayList<String> singleCellNames,
			HashMap<String, Integer> doubletFlagMap,
			int nCell, int nMut, double fp, double fnADO,
			double errorlearn, double mutRate, 
			int n_iter, int p_iter,
			int dataFlag, int modelFlag,
			BasicUtilityFunctions BUF) throws SuspendableException{
		
		Random _rng = new Random();
        // fn Parameters
        double fnPriorMean = BUF.getFNPrior(obsGenotypeMatrix, fnADO);
        double fnPriorSD = fnPriorMean*0.5;
        double betaPriora = ((1 - fnPriorMean)*fnPriorMean*fnPriorMean/(fnPriorSD*fnPriorSD)) - fnPriorMean;
        double betaPriorb = betaPriora * ((1/fnPriorMean) - 1);
        
        // Obtain matrix for running parsimony
        ArrayList<Integer[]> obsCellGenomes = BUF.getCellGenomes(obsGenotypeMatrix, nMut);
        
		// Generate random tree
		STITree<Double> randTreeMaker = BUF.generateRandomTree(singleCellNames);
		STITree<Double> randTree1 = BUF.getTree(randTreeMaker.toNewick());
		BUF.resizeBranchLengths(randTree1, 1);
		
		// Run Parsimony to get the initial tree
		FitchParsimonyCalc FitchObj = new FitchParsimonyCalc();
		TopologyBranchPerturbations TBP = new TopologyBranchPerturbations();
		String[] singleCellNameArray = new String[nCell];
		for (int i=0; i< nCell; i++){
			singleCellNameArray[i] = singleCellNames.get(i);
		}
		
		STITree<Double> randTree = FitchObj.findParsimonyTree(randTree1, singleCellNameArray, obsCellGenomes, 3000, nCell, TBP);
		

		// Data Structures to store tree, error rate and likelihood values
        ArrayList<STITree<Double>> treeList = new ArrayList<STITree<Double>>();
        double[] treeLoglikelihoodArr = new double[n_iter];
        double[] fnArr = new double[n_iter];
        treeList.add(randTree);
        fnArr[0] = fnPriorMean;
        
        TreeSearchHelper jointSampleHelper = new TreeSearchHelper(fp, fnPriorMean, mutRate, BUF, dataFlag, modelFlag, doubletFlagMap);
        treeLoglikelihoodArr[0] = jointSampleHelper.computeTreeLikelihood(obsGTObsMatrix, randTree, fnPriorMean);
        
        double bestLikelihood = Double.NEGATIVE_INFINITY;
        STITree<Double> bestTree = randTree;
        STITree<Double> previousTree = randTree;
        int bestTreeIndex = 0;
        double bestFn = fnPriorMean;
        
        for (int i = 1; i < n_iter; i++){
//        	STITree<Double> currTree = treeList.get(i-1);
        	STITree<Double> currTree = previousTree;
        	double loglikelihood_currTree = treeLoglikelihoodArr[i-1];
        	double currFn = fnArr[i-1];
        	double rr = _rng.nextDouble();
        	
        	// Propose a new error value
        	if (rr <= errorlearn){
        		double newFn = jointSampleHelper.proposeFnObj.proposeNewFn(currFn, fnPriorSD);
        		double newFnLogScore = jointSampleHelper.proposeFnObj.logBetaPDF(newFn, betaPriora, betaPriorb);
        		double newFnPriorRatio = newFnLogScore - jointSampleHelper.proposeFnObj.logBetaPDF(currFn, betaPriora, betaPriorb);
        		double newFnProposalRatio = jointSampleHelper.proposeFnObj.getProposalRatio(currFn, newFn, fnPriorSD);
        		double newFnLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newFn, currTree);
        		double fnAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioFn(newFnLogLikelihoodRatio, newFnPriorRatio, newFnProposalRatio);
        		boolean fnAcceptFlag = jointSampleHelper.getAcceptanceFlag(fnAcceptanceRatio);
//        		treeList.add(currTree);
        		
        		
        		// new error value accepted
        		if(fnAcceptFlag == true){
	        		fnArr[i] = newFn;
	        		treeLoglikelihoodArr[i] = jointSampleHelper.newFnSameTreeLikelihood;
        		}
        		// new error value rejected
        		else{
        			fnArr[i] = currFn;
	        		treeLoglikelihoodArr[i] = loglikelihood_currTree;
        		}
        	}
        	// Propose a new tree
        	else{
        		STITree<Double> newTree = jointSampleHelper.TBP.proposeTree(currTree, nCell);
        		double newTreeLogLikelihoodRatio = jointSampleHelper.getLogLikelihoodRatio(obsGTObsMatrix, loglikelihood_currTree, newTree);
        		double treeAcceptanceRatio = jointSampleHelper.computeAcceptanceRatioTree(newTreeLogLikelihoodRatio, jointSampleHelper.TBP.getProposalRatio());
        		boolean treeAcceptFlag = jointSampleHelper.getAcceptanceFlag(treeAcceptanceRatio);
        		fnArr[i] = currFn;
        		
        		// new tree accepted
            	if (treeAcceptFlag == true){
//            		treeList.add(newTree);
            		previousTree = newTree;
            		treeLoglikelihoodArr[i] = jointSampleHelper.sameFnNewTreeLikelihood;
            	}
            	// new tree rejected
            	else{
//            		treeList.add(currTree);
            		
            		treeLoglikelihoodArr[i] = loglikelihood_currTree;
            	}
        	}
        	if (treeLoglikelihoodArr[i] > bestLikelihood){
        		bestLikelihood = treeLoglikelihoodArr[i];
//        		bestTree = treeList.get(i);
        		bestTree = previousTree;
        		bestTreeIndex = i;
        		bestFn = fnArr[i];
        	}
        	if (i % p_iter == 0){
        		System.out.printf("After %d iterations, likelihood of the best configuration is %f\n", i, bestLikelihood);
        	}
        	
        }
        
        
        String bestTreeNewick = bestTree.toNewick();
        TreeLikelihoodObj bestConfiguration = new TreeLikelihoodObj(bestFn, bestTreeNewick, bestLikelihood);
		return bestConfiguration;
	}
	
	public static void main(String[] args) throws IOException, SuspendableException{
		

	}
	

}
