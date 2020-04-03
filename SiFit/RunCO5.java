/**
 * 
 */
package SiFit;

import static edu.rice.hj.Module0.launchHabaneroApp;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import SiFit.algorithm.PhyloTreeSearch;
import SiFit.io.VariantMatrixReader;
import SiFit.objects.GenotypeObservation;
import SiFit.objects.TreeLikelihoodObj;

/**
 * @author hz22
 *
 */
public class RunCO5 {
	
	public static int nCell = 178;							// Number of cells
	public static int nMut = 16;								// Number of mutations
	public static double _fp = 0.0152;					// FP rate, alpha
	public static double _fnStart = 0.1;				// FN rate, beta
	public static double _errorLearnProb = 0.25;		// Probability for changing error rate
	public static double _modelParamLearnProb = 0.1;    // Probability for changing model parameter
	public static double _metropolisHastingsProb = 0.1; // Probability for using Metropolis-Hastings Steps
	public static double _delProb = 0.05;				// Deletion Rate
	public static double _LOHProb = 0.01;				// LOH Rate
	public static double _recurProb = 0.05;				// Recurrent mutation Rate
	public static int iterT = 50000;					// Total number of iterations
	public static int iterP = 1000;						// Number of iterations after which likelihood will be printed
	public static int restart = 1;						// Number of restarts
	public static double _mu = 0.1;						// Mu, parameter of the evolution model
	public static int dataFlag = 0;						// Flag to indicate the data type, binary or ternary
	public static int modelFlag = 1;					// Flag to indicate which model to use
	public static String varMatFilename = null;			// Filename of the input matrix
	public static String trueTreeFilename = null;		// File containing the true tree in newick form
	public static String cellNamesFilename = null;		// File containing the cell names
	public static String opTreeFileName = null;

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		readArguments(args);
		
		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		Random _rng = new Random();

        // Single cell info
        ArrayList<String> singleCellNames = new ArrayList<>();
        if (cellNamesFilename == null){
        	for (int i = 1; i <= nCell; i++){
            	singleCellNames.add("sc" + Integer.toString(i));
            }
        }
        else{
	        String singleCellNameString = BUF.readNewickString(cellNamesFilename);
	        String[] singleCellNameArray = singleCellNameString.split(" ");
	        for (String s: singleCellNameArray){
	        	singleCellNames.add(s);
	        }
        }
                
		// Obtain Genotype Observations from Input File
        VariantMatrixReader vr = new VariantMatrixReader(varMatFilename);
        vr.populateVarGTMatrix(varMatFilename, nCell);
        ArrayList<Integer[]> obsGenotypeMatrix = vr.obsVarGTMatrix;
		ArrayList<GenotypeObservation> obsGTObsMatrix = BUF.getGenotypeObsFromGTMatrix(singleCellNames, obsGenotypeMatrix);

		// List of likelihoods and Trees for restarts
		ArrayList<Double> listLikelihood = new ArrayList<>();
		ArrayList<String> listBestTreeNewick = new ArrayList<>();
		ArrayList<Double> listError = new ArrayList<>();
		ArrayList<TreeLikelihoodObj> listTreeLikelihoodObj = new ArrayList<>();
		
		
		// Launching HabaneroApp for parallel computation of likelihood
		launchHabaneroApp(() -> {
		
		// Looping over the number of restarts
		for (int i = 0; i < restart; i++){
			System.out.printf("Run %d\n", i);
			
			
			// Inference with the complex model of evolution
			TreeLikelihoodObj currRunTreeLikelihood = PhyloTreeSearch.findBestTreeCO5(
														obsGTObsMatrix, 
														obsGenotypeMatrix, 
														singleCellNames,
														nCell, nMut, _fp,
														_fnStart, _errorLearnProb, 
														_modelParamLearnProb, _metropolisHastingsProb,
														_delProb, _LOHProb, _recurProb,
														iterT, iterP, 
														dataFlag, modelFlag, i, BUF);
			
			
			// This is for original implementation
			listError.add(currRunTreeLikelihood.bestBeta);
			listBestTreeNewick.add(currRunTreeLikelihood.bestTreeNewickString);
			listLikelihood.add(currRunTreeLikelihood.bestTreeLikelihood);
			listTreeLikelihoodObj.add(currRunTreeLikelihood);
		}
		}); // Habanero App launch ends here
		
		double bestLikelihood = listLikelihood.get(0);
		int bestConfigIndex = 0;
		for (int i = 1; i < restart; i++){
			if (listLikelihood.get(i) > bestLikelihood){
				bestConfigIndex = i;
				bestLikelihood = listLikelihood.get(i);
			}			
		}
		
		String bestTreeNewick = listBestTreeNewick.get(bestConfigIndex);
		double bestError = listError.get(bestConfigIndex);
		
		System.out.println("\n");
		System.out.printf("best false negative rate = %f\n", bestError);
		System.out.printf("best deletion parameter = %f\n", listTreeLikelihoodObj.get(bestConfigIndex).bestDelProb);
		System.out.printf("best LOH parameter = %f\n", listTreeLikelihoodObj.get(bestConfigIndex).bestLOHProb);
		System.out.printf("best log-likelihood score = %f\n", bestLikelihood);
		System.out.printf("best tree = %s\n", bestTreeNewick);
		
		// Write best tree to output file
		opTreeFileName = varMatFilename.split(".txt")[0] + "_mlTree.newick";
		BUF.writeNewickTree(opTreeFileName, bestTreeNewick);

	}

	private static void readArguments(String[] args) {
		int nPar = args.length;
		for (int i = 0; i < nPar; i = i+2){
			if (args[i].equals("-mp") == true){
				_modelParamLearnProb = Double.parseDouble(args[i+1]);
				if (_modelParamLearnProb < 0 || _modelParamLearnProb > 1)
					throw new IllegalArgumentException("Invalid input for option -mp, should be between 0 and 1");
			}			
			else if (args[i].equals("-iter") == true){
				iterT = Integer.parseInt(args[i+1]);
			}
			else if (args[i].equals("-printIter") == true){
				iterP = Integer.parseInt(args[i+1]);
			}
			else if (args[i].equals("-r") == true){
				restart = Integer.parseInt(args[i+1]);
			}
			else if (args[i].equals("-ipMat") == true){
				varMatFilename = args[i+1];
			}
			else if (args[i].equals("-trueTree") == true){
				trueTreeFilename = args[i+1];
			}
			else if (args[i].equals("-cellNames") == true){
				cellNamesFilename = args[i+1];
			}
		}
		
	}
}
