/**
 * 
 */
package SiFit.algorithm;

import static edu.rice.hj.Module1.launchHabaneroApp;
import static edu.rice.hj.Module1.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import SiFit.BasicUtilityFunctions;
import SiFit.io.VariantMatrixReader;
import SiFit.objects.GenotypeObservation;
import cern.colt.Arrays;
import jeigen.DenseMatrix;
import SiFit.model.ComplexEvolutionModel;
import SiFit.model.JCModelSingleCell;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.NucleotideProbabilityAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * @author hz22
 *
 */
public class ComputeLikelihood extends NucleotideProbabilityAlgorithm {
	
	/**
     *  The parameters that remain static throughout the algorithm:
     *  fp - the false positive error rate
     *  fn - the false negative error rate
     *  rates - the error rates; the i,j entry represents the probability of observing state i as state j
     *  mutation_rate - the approximated/expected mutation rate over the entire tree (needed for Jukes-Cantor)
     *  genotypes - the possible genotypes for single cell data
     */
	public Tree tree;
	public JCModelSingleCell model;
	public ComplexEvolutionModel modelC;
	private static double fp;
    private static double fn;
    private static double mu;
    private static double[][] rates;
    private static double[][] ratesD;
    public static final String genotypes = "012";
    
    public ComputeLikelihood(double fpositive, double fnegative, double mut_rate, Tree t, JCModelSingleCell model){
    	fp = fpositive; fn = fnegative; mu = mut_rate;
        rates = new double[][]
                {{(1 - fp*fn/2.0 - fp),     fp,         fp*fn/2.0},
                 {fn/2.0,                   1 - fn,     fn/2.0},
                 {0.0,                      0.0,        1.0}};

        this.tree = t;
        this.model = model;
    }
    
    // This constructor uses the complex evolution model
    public ComputeLikelihood(double fpositive, double fnegative, double mut_rate, Tree t, ComplexEvolutionModel model){
    	fp = fpositive; fn = fnegative; mu = mut_rate;
        rates = new double[][]
                {{(1 - fp*fn/2.0 - fp),     fp,         fp*fn/2.0},
                 {fn/2.0,                   1 - fn,     fn/2.0},
                 {0.0,                      0.0,        1.0}};
        this.tree = t;
        this.modelC = model;
//        System.out.println(modelC.delProb);
    }
    
    /**
     * Sets a new false negative error rate
     * @param newFn, the new false negative error rate
     * @return true if error rate is valid and set, false otherwise
     */
    @SuppressWarnings("static-access")
	public boolean setFn(double newFn) {
        if (newFn > 1.0 || newFn < 0.0) {
            return false;
        }
        this.fn = newFn;
        this.rates[0][0] = (1 - fp*fn/2.0 - fp);
        this.rates[0][2] = fp*fn/2.0;
        this.rates[1][0] = fn/2.0;
        this.rates[1][1] = 1 - fn;
        this.rates[1][2] = fn/2.0;
        return true;
    }
    
    /**
     * Create the base case leaf matrix for the given genotype
     * @param genotype - single-character string representing the genotype of leaf 
     * @return DenseMatrix for the base cases of the leaf
     */
    private static DenseMatrix createLeafMatrix(Integer genotype){
    	/* Handle missing data. Missing data will be 
    	 * encoded by '-' or '3' and for that
    	 * each genotype is 1.
    	 */
    	if (genotype == 3){
            double [][] t = {{1},{1},{1}};
            DenseMatrix D = new DenseMatrix(t);
            return D;
    	}
    	double [][] temp = {
                {0},
                {0},
                {0}            
            };
    	if (genotypes.indexOf(Integer.toString(genotype)) == -1)
    		throw new RuntimeException(genotype +" is not a valid genotype");

    	
    	
    	// SCS Error model
    	for (int i = 0; i <= 2; i++){
    		temp[i][0] = rates[i][genotype];
    	}
    	return new DenseMatrix(temp);    		    	   	
    }
    
    private static DenseMatrix createLeafMatrixDoublet(Integer genotype){
    	/* Handle missing data. Missing data will be 
    	 * encoded by '-' or '3' and for that
    	 * each genotype is 1.
    	 */
    	if (genotype == 3){
            double [][] t = {{1},{1},{1}};
            DenseMatrix D = new DenseMatrix(t);
            return D;
    	}
    	double [][] temp = {
                {0},
                {0},
                {0}            
            };
    	if (genotypes.indexOf(Integer.toString(genotype)) == -1)
    		throw new RuntimeException(genotype +" is not a valid genotype");

    	
    	
    	// SCS Error model
    	for (int i = 0; i <= 2; i++){
    		temp[i][0] = ratesD[i][genotype];
    	}
    	return new DenseMatrix(temp);    		    	   	
    }
    
    /**
     * This returns the likelihood of a node having each genotype
     * @param node - The node to find the likelihood of.
     * @param obs_genotypes - A map of genotypes of the single cells.
     * @return A column matrix representing the likelihood. Reference {@link #nucleotides} to see which index corresponds to each base.
     */
    private DenseMatrix probabilityOfNode(TNode node, GenotypeObservation obs_genotypes, int modelFlag){
    	if (node.isLeaf()){    	
    		
    		return createLeafMatrix(obs_genotypes.cellGenotypeMap.get(node.getName()));
    	}
    	else{
    		double[][] probs = new double[][]{{1.0}, {1.0}, {1.0}};
            for (TNode child : node.getChildren()) {
                double dist = child.getParentDistance();
                if (Double.isNaN(dist) || Double.isInfinite(dist))
                    throw new IllegalArgumentException("Node : " + child.getName() + " has non-finite dist: "+dist);
                double[] edgeP = new double[3];
                DenseMatrix childProbs = probabilityOfNode(child, obs_genotypes, modelFlag);
             
                double[][] JCTransMat = new double[3][3];
                
                // Now just testing the new complex evolution model
                
                JCTransMat = modelC.getTransitionMatrix(dist);

                for (int i = 0; i < 3; i++){
                	edgeP[0] += JCTransMat[0][i] * childProbs.get(i, 0);
                	edgeP[1] += JCTransMat[1][i] * childProbs.get(i, 0);
                	edgeP[2] += JCTransMat[2][i] * childProbs.get(i, 0);
                	
                }
                probs[0][0] *= edgeP[0]; probs[1][0] *= edgeP[1]; probs[2][0] *= edgeP[2];
            }
            DenseMatrix probMat = new DenseMatrix(probs);

            return probMat;
    	}
    }
    
    /**
     * Function that returns the likelihood of genotype observation at given locus
     * @param arg0 - A map of genotypes of the single cells
     * @return the likelihood
     */
	public double getGenotypeProbability(GenotypeObservation arg0, int modelFlag) {
		
		return probabilityOfNode(tree.getRoot(), arg0, modelFlag).mul(modelC.getEquilibriumVector()).sum().s(); // returns densematrix.get(0,0);
	}
	
	public double getGenotypeProbability(GenotypeObservation arg0,
			int modelFlag, HashMap<String, Integer> _doubletFlagMap) {
		return probabilityOfNode(tree.getRoot(), arg0, modelFlag, _doubletFlagMap).mul(model.getEquilibriumVector()).sum().s();
	}

    private DenseMatrix probabilityOfNode(TNode node, GenotypeObservation obs_genotypes, int modelFlag, HashMap<String, Integer> _doubletFlagMap){
    	if (node.isLeaf()){
    		if (_doubletFlagMap.get(node.getName()) == 0)
    			return createLeafMatrix(obs_genotypes.cellGenotypeMap.get(node.getName()));
    		else
    			return createLeafMatrixDoublet(obs_genotypes.cellGenotypeMap.get(node.getName()));
    	}
    	else{
    		double[][] probs = new double[][]{{1.0}, {1.0}, {1.0}};
            for (TNode child : node.getChildren()) {
                double dist = child.getParentDistance();
                if (Double.isNaN(dist) || Double.isInfinite(dist))
                    throw new IllegalArgumentException("Node : " + child.getName() + " has non-finite dist: "+dist);
                double[] edgeP = new double[3];
                DenseMatrix childProbs = probabilityOfNode(child, obs_genotypes, modelFlag);
             
                double[][] JCTransMat = new double[3][3];
//                JCTransMat = new BasicUtilityFunctions().computeTransMat(0.025, 0.01, 0.01);
                if (modelFlag == 0){
                	JCTransMat = model.jukesCantor012Matrix(dist);
                }
                else{
                	JCTransMat = model.jukesCantor012AllMatrix(dist);
                }
                for (int i = 0; i < 3; i++){
                	edgeP[0] += JCTransMat[0][i] * childProbs.get(i, 0);
                	edgeP[1] += JCTransMat[1][i] * childProbs.get(i, 0);
                	edgeP[2] += JCTransMat[2][i] * childProbs.get(i, 0);
                	
                }
                probs[0][0] *= edgeP[0]; probs[1][0] *= edgeP[1]; probs[2][0] *= edgeP[2];
            }
            DenseMatrix probMat = new DenseMatrix(probs);
            return probMat;
    	}
    }


	/* (non-Javadoc)
	 * @see edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.NucleotideProbabilityAlgorithm#getProbability(edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation)
	 */
	@Override
	public double getProbability(NucleotideObservation arg0) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	public static void main(String[] args) throws IOException{
		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		String t = "(((((((sc44:0.09267023264241692,sc86:0.7106264668386636):0.8193475095058459,sc34:0.6045190667188953):0.14953103222906605,(sc35:0.745147167878201,sc6:0.2061639133327101):0.9203722061518866):0.46492352836969997,(((sc45:0.5395378113075877,sc37:0.5844042282429451):0.37438313211457463,sc88:0.0016406267297417454):0.032435450509692876,(sc65:0.34936439387004814,(sc73:0.9013300778200583,(sc2:0.1216193434713968,sc12:0.5053563416219522):0.16718374502158095):0.6951562782217547):0.2600348306642277):0.4064122322397178):0.36384127820832357,((((sc19:0.29080347938106976,sc11:0.2449621762602091):0.14869372015568383,sc84:0.33758840840971804):0.04970748393769653,(sc30:0.8544932182594116,sc38:0.820111624656587):0.26134283329587416):0.8295508288903154,(sc96:0.2926567007183887,sc22:0.9831647875314927):0.7711404676296971):0.1563762609615572):0.6969689113428985,((((((sc72:0.5920872415996888,sc14:0.09113051506121439):0.4970590121693742,sc29:0.5311758181322215):0.7116304294879077,(sc83:0.41247815467420124,sc17:0.6503973690532682):0.9422272471528415):0.08476015478247001,(sc81:0.292964604332248,(sc74:0.22147365620284176,sc13:0.9210047290888839):0.8307619514734793):0.7696177324618156):0.2006021320533018,((sc80:0.9604162810362551,sc55:0.4899691468725482):0.07910197027688348,((sc51:0.017501505998699707,(sc66:0.6754805984774501,sc36:0.6120272050506667):0.6052407132546068):0.9763867307773165,(sc89:0.3785516886482517,sc7:0.19116926821004465):0.39671243735715234):0.05135903835865441):0.39318257056173944):0.2422888597888898,((((sc98:0.31277465074103494,sc47:0.6905864522774174):0.3120442220463494,sc39:0.5124660872494996):0.5827497609044572,sc60:0.06196762774663822):0.030703608262264592,((sc54:0.9758697995063188,sc85:0.17777807276241608):0.6147082523584194,((sc71:0.4607374036582691,sc92:0.5744274862385146):0.4779244397419118,sc75:0.36957597682368937):0.06305279105877748):0.9128248017014047):0.7994951694072885):0.8609911840660774):0.05374871844770279,(((((((sc67:0.988210621254911,sc97:0.2384982664199664):0.7296538667908803,(sc53:0.6404175419859323,sc95:0.2886059615530461):0.3078389820763713):0.9203196456517028,sc25:0.4063214548938363):0.06333499772010864,(sc61:0.2500163080138659,(sc62:0.9674127797997226,sc94:0.962988756963935):0.16488701320015509):0.36695303366073795):0.7652658541784468,((((sc23:0.19286037405624357,sc24:0.8700254840130007):0.07222577258599794,sc21:0.9166616449187812):0.15828591391296243,sc100:0.40643920400082856):0.0893538913891525,(sc77:0.19430131630030267,(sc63:0.9973465515711889,sc70:0.4729544973761104):0.2278943333881912):0.9202876030255625):0.3079646551507301):0.025287292175400378,((((sc1:0.5649562454268737,sc69:0.014236563710952721):0.5490272985336153,sc91:0.3847044527312089):0.5526184785759859,((sc32:0.3598157752938751,sc50:0.45545265485535247):0.582309709440751,((sc16:0.5462621819892782,sc27:0.858050568545695):0.7865991435453384,sc57:0.7616061759471058):0.719572643451926):0.22113374580766776):0.21402416121210044,((((sc58:0.6738588509402227,(sc18:0.8099887838597847,sc40:0.5282100590370509):0.5904278736298433):0.1450060025433475,sc33:0.8691676807580687):0.31173086944457373,sc79:0.9447414154680907):0.7886185152173748,(((sc90:0.8379285161549896,sc64:0.07958356802535849):0.718485755434668,sc59:0.6844605670224873):0.5906099256317233,sc41:0.6532666129156959):0.672572197599557):0.07114422188054748):0.5373249912794722):0.9151734892115296,((((sc93:0.5014674484883928,(sc8:0.5227224389568337,sc42:0.7713765102027619):0.6802057457996086):0.29994194839621,(sc26:0.36827102586321847,sc9:0.06513504749243482):0.10975731700353808):0.749404121246939,((sc56:0.31665428402557705,sc46:0.004209901766304425):0.19090330407482992,(sc4:0.625504370219896,sc28:0.8751181382674112):0.9762422498725231):0.2356115790797817):0.24480675524574258,((((sc43:0.3701821857773967,sc31:0.5109198909476993):0.24681713653109083,sc49:0.2233560585255886):0.45327333544919024,((((sc15:0.21742158057284422,(sc99:0.8228521440579092,sc5:0.3724840828685849):0.43149584620634085):0.3015204422000698,sc20:0.4268227351124053):0.39243795459332675,sc10:0.5165875049183466):0.4595156436219041,sc78:0.244785986997905):0.16694597722228333):0.4038406173647606,((sc3:0.7661485994473264,(sc87:0.5928594949511611,sc82:0.42578552441466755):0.7974935064675907):0.04894536208209821,((sc76:0.48834090490733506,sc52:0.3363253970669098):0.38124267991401284,(sc48:0.2145705080692446,sc68:0.12249045865988195):0.24892768980248803):0.01644362573758895):0.16762254819397526):0.4599983296732344):0.4182456117657932):0.8419343495905224);";
		STITree<Double> randTree = BUF.getTree(t);
		String varMatFilename = "/Users/hz22/Acads/Research/Cancer_Genome_Research/Research Ideas/Nu_Phylo_SNV_g_m/Datasets/nuDatasets/finiteSites/100_cells/200_sites/noDoublet/dataset1/noisy_genotype_mat_dataset1.txt";
        VariantMatrixReader vr = new VariantMatrixReader(varMatFilename);
        vr.populateVarGTMatrix(varMatFilename, 100);
        ArrayList<Integer[]> obsGenotypeMatrix = vr.obsVarGTMatrix;
     // Single cell info
        ArrayList<String> singleCellNames = new ArrayList<>();
        
        	for (int i = 1; i <= 100; i++){
            	singleCellNames.add("sc" + Integer.toString(i));
            }
        
		ArrayList<GenotypeObservation> obsGTObsMatrix = BUF.getGenotypeObsFromGTMatrix(singleCellNames, obsGenotypeMatrix);
		
		// serial compute
		JCModelSingleCell model = new JCModelSingleCell(0.1);
		ComputeLikelihood CL = new ComputeLikelihood(0.001, 0.2, 0.1, randTree, model);
		int nMut = obsGTObsMatrix.size();
		double[] val = new double[nMut];
		for(int i=0;i<=nMut -1;i++)
			val[i] = CL.getGenotypeProbability(obsGTObsMatrix.get(i), 1);
		double serial_likelihood = 0;
		for (double d: val){
			double log_d;
			if (d == 0)
				log_d = -10000;
			else
				log_d = Math.log(d);
			serial_likelihood += log_d;
		}
		System.out.println(serial_likelihood);
		
		// parallel compute
		double[] pval = new double[nMut];
		double par_likelihood = 0;
		launchHabaneroApp(() -> {
			forallChunked(0, nMut -1, (i) -> {
				pval[i] = CL.getGenotypeProbability(obsGTObsMatrix.get(i), 1);
			});
		});
		for (double d: pval){
			double log_d;
			if (d == 0)
				log_d = -10000;
			else
				log_d = Math.log(d);
			par_likelihood += log_d;
		}
		System.out.println(par_likelihood);
		
//		String pt = "((((((sc44:0.09267023264241692,sc86:0.7106264668386636,sc34:0.6045190667188953):0.14953103222906605,(sc35:0.745147167878201,sc6:0.2061639133327101):0.9203722061518866):0.46492352836969997,(((sc45:0.5395378113075877,sc37:0.5844042282429451):0.37438313211457463,sc88:0.0016406267297417454):0.032435450509692876,(sc65:0.34936439387004814,(sc73:0.9013300778200583,(sc2:0.1216193434713968,sc12:0.5053563416219522):0.16718374502158095):0.6951562782217547):0.2600348306642277):0.4064122322397178):0.36384127820832357,((((sc19:0.29080347938106976,sc11:0.2449621762602091):0.14869372015568383,sc84:0.33758840840971804):0.04970748393769653,(sc30:0.8544932182594116,sc38:0.820111624656587):0.26134283329587416):0.8295508288903154,(sc96:0.2926567007183887,sc22:0.9831647875314927):0.7711404676296971):0.1563762609615572):0.6969689113428985,((((((sc72:0.5920872415996888,sc14:0.09113051506121439):0.4970590121693742,sc29:0.5311758181322215):0.7116304294879077,(sc83:0.41247815467420124,sc17:0.6503973690532682):0.9422272471528415):0.08476015478247001,(sc81:0.292964604332248,(sc74:0.22147365620284176,sc13:0.9210047290888839):0.8307619514734793):0.7696177324618156):0.2006021320533018,((sc80:0.9604162810362551,sc55:0.4899691468725482):0.07910197027688348,((sc51:0.017501505998699707,(sc66:0.6754805984774501,sc36:0.6120272050506667):0.6052407132546068):0.9763867307773165,(sc89:0.3785516886482517,sc7:0.19116926821004465):0.39671243735715234):0.05135903835865441):0.39318257056173944):0.2422888597888898,((((sc98:0.31277465074103494,sc47:0.6905864522774174):0.3120442220463494,sc39:0.5124660872494996):0.5827497609044572,sc60:0.06196762774663822):0.030703608262264592,((sc54:0.9758697995063188,sc85:0.17777807276241608):0.6147082523584194,((sc71:0.4607374036582691,sc92:0.5744274862385146):0.4779244397419118,sc75:0.36957597682368937):0.06305279105877748):0.9128248017014047):0.7994951694072885):0.8609911840660774):0.05374871844770279,(((((((sc67:0.988210621254911,sc97:0.2384982664199664):0.7296538667908803,(sc53:0.6404175419859323,sc95:0.2886059615530461):0.3078389820763713):0.9203196456517028,sc25:0.4063214548938363):0.06333499772010864,(sc61:0.2500163080138659,(sc62:0.9674127797997226,sc94:0.962988756963935):0.16488701320015509):0.36695303366073795):0.7652658541784468,((((sc23:0.19286037405624357,sc24:0.8700254840130007):0.07222577258599794,sc21:0.9166616449187812):0.15828591391296243,sc100:0.40643920400082856):0.0893538913891525,(sc77:0.19430131630030267,(sc63:0.9973465515711889,sc70:0.4729544973761104):0.2278943333881912):0.9202876030255625):0.3079646551507301):0.025287292175400378,((((sc1:0.5649562454268737,sc69:0.014236563710952721):0.5490272985336153,sc91:0.3847044527312089):0.5526184785759859,((sc32:0.3598157752938751,sc50:0.45545265485535247):0.582309709440751,((sc16:0.5462621819892782,sc27:0.858050568545695):0.7865991435453384,sc57:0.7616061759471058):0.719572643451926):0.22113374580766776):0.21402416121210044,((((sc58:0.6738588509402227,(sc18:0.8099887838597847,sc40:0.5282100590370509):0.5904278736298433):0.1450060025433475,sc33:0.8691676807580687):0.31173086944457373,sc79:0.9447414154680907):0.7886185152173748,(((sc90:0.8379285161549896,sc64:0.07958356802535849):0.718485755434668,sc59:0.6844605670224873):0.5906099256317233,sc41:0.6532666129156959):0.672572197599557):0.07114422188054748):0.5373249912794722):0.9151734892115296,((((sc93:0.5014674484883928,(sc8:0.5227224389568337,sc42:0.7713765102027619):0.6802057457996086):0.29994194839621,(sc26:0.36827102586321847,sc9:0.06513504749243482):0.10975731700353808):0.749404121246939,((sc56:0.31665428402557705,sc46:0.004209901766304425):0.19090330407482992,(sc4:0.625504370219896,sc28:0.8751181382674112):0.9762422498725231):0.2356115790797817):0.24480675524574258,((((sc43:0.3701821857773967,sc31:0.5109198909476993):0.24681713653109083,sc49:0.2233560585255886):0.45327333544919024,((((sc15:0.21742158057284422,(sc99:0.8228521440579092,sc5:0.3724840828685849):0.43149584620634085):0.3015204422000698,sc20:0.4268227351124053):0.39243795459332675,sc10:0.5165875049183466):0.4595156436219041,sc78:0.244785986997905):0.16694597722228333):0.4038406173647606,((sc3:0.7661485994473264,(sc87:0.5928594949511611,sc82:0.42578552441466755):0.7974935064675907):0.04894536208209821,((sc76:0.48834090490733506,sc52:0.3363253970669098):0.38124267991401284,(sc48:0.2145705080692446,sc68:0.12249045865988195):0.24892768980248803):0.01644362573758895):0.16762254819397526):0.4599983296732344):0.4182456117657932):0.8419343495905224);";
		String pt = "(((((sc44:0.09267023264241692,sc86:0.7106264668386636,sc34:0.6045190667188953,sc35:0.745147167878201,sc6:0.2061639133327101):0.46492352836969997,(((sc45:0.5395378113075877,sc37:0.5844042282429451):0.37438313211457463,sc88:0.0016406267297417454):0.032435450509692876,(sc65:0.34936439387004814,(sc73:0.9013300778200583,(sc2:0.1216193434713968,sc12:0.5053563416219522):0.16718374502158095):0.6951562782217547):0.2600348306642277):0.4064122322397178):0.36384127820832357,((((sc19:0.29080347938106976,sc11:0.2449621762602091):0.14869372015568383,sc84:0.33758840840971804):0.04970748393769653,(sc30:0.8544932182594116,sc38:0.820111624656587):0.26134283329587416):0.8295508288903154,(sc96:0.2926567007183887,sc22:0.9831647875314927):0.7711404676296971):0.1563762609615572):0.6969689113428985,((((((sc72:0.5920872415996888,sc14:0.09113051506121439):0.4970590121693742,sc29:0.5311758181322215):0.7116304294879077,(sc83:0.41247815467420124,sc17:0.6503973690532682):0.9422272471528415):0.08476015478247001,(sc81:0.292964604332248,(sc74:0.22147365620284176,sc13:0.9210047290888839):0.8307619514734793):0.7696177324618156):0.2006021320533018,((sc80:0.9604162810362551,sc55:0.4899691468725482):0.07910197027688348,((sc51:0.017501505998699707,(sc66:0.6754805984774501,sc36:0.6120272050506667):0.6052407132546068):0.9763867307773165,(sc89:0.3785516886482517,sc7:0.19116926821004465):0.39671243735715234):0.05135903835865441):0.39318257056173944):0.2422888597888898,((((sc98:0.31277465074103494,sc47:0.6905864522774174):0.3120442220463494,sc39:0.5124660872494996):0.5827497609044572,sc60:0.06196762774663822):0.030703608262264592,((sc54:0.9758697995063188,sc85:0.17777807276241608):0.6147082523584194,((sc71:0.4607374036582691,sc92:0.5744274862385146):0.4779244397419118,sc75:0.36957597682368937):0.06305279105877748):0.9128248017014047):0.7994951694072885):0.8609911840660774):0.05374871844770279,(((((((sc67:0.988210621254911,sc97:0.2384982664199664):0.7296538667908803,(sc53:0.6404175419859323,sc95:0.2886059615530461):0.3078389820763713):0.9203196456517028,sc25:0.4063214548938363):0.06333499772010864,(sc61:0.2500163080138659,(sc62:0.9674127797997226,sc94:0.962988756963935):0.16488701320015509):0.36695303366073795):0.7652658541784468,((((sc23:0.19286037405624357,sc24:0.8700254840130007):0.07222577258599794,sc21:0.9166616449187812):0.15828591391296243,sc100:0.40643920400082856):0.0893538913891525,(sc77:0.19430131630030267,(sc63:0.9973465515711889,sc70:0.4729544973761104):0.2278943333881912):0.9202876030255625):0.3079646551507301):0.025287292175400378,((((sc1:0.5649562454268737,sc69:0.014236563710952721):0.5490272985336153,sc91:0.3847044527312089):0.5526184785759859,((sc32:0.3598157752938751,sc50:0.45545265485535247):0.582309709440751,((sc16:0.5462621819892782,sc27:0.858050568545695):0.7865991435453384,sc57:0.7616061759471058):0.719572643451926):0.22113374580766776):0.21402416121210044,((((sc58:0.6738588509402227,(sc18:0.8099887838597847,sc40:0.5282100590370509):0.5904278736298433):0.1450060025433475,sc33:0.8691676807580687):0.31173086944457373,sc79:0.9447414154680907):0.7886185152173748,(((sc90:0.8379285161549896,sc64:0.07958356802535849):0.718485755434668,sc59:0.6844605670224873):0.5906099256317233,sc41:0.6532666129156959):0.672572197599557):0.07114422188054748):0.5373249912794722):0.9151734892115296,((((sc93:0.5014674484883928,(sc8:0.5227224389568337,sc42:0.7713765102027619):0.6802057457996086):0.29994194839621,(sc26:0.36827102586321847,sc9:0.06513504749243482):0.10975731700353808):0.749404121246939,((sc56:0.31665428402557705,sc46:0.004209901766304425):0.19090330407482992,(sc4:0.625504370219896,sc28:0.8751181382674112):0.9762422498725231):0.2356115790797817):0.24480675524574258,((((sc43:0.3701821857773967,sc31:0.5109198909476993):0.24681713653109083,sc49:0.2233560585255886):0.45327333544919024,((((sc15:0.21742158057284422,(sc99:0.8228521440579092,sc5:0.3724840828685849):0.43149584620634085):0.3015204422000698,sc20:0.4268227351124053):0.39243795459332675,sc10:0.5165875049183466):0.4595156436219041,sc78:0.244785986997905):0.16694597722228333):0.4038406173647606,((sc3:0.7661485994473264,(sc87:0.5928594949511611,sc82:0.42578552441466755):0.7974935064675907):0.04894536208209821,((sc76:0.48834090490733506,sc52:0.3363253970669098):0.38124267991401284,(sc48:0.2145705080692446,sc68:0.12249045865988195):0.24892768980248803):0.01644362573758895):0.16762254819397526):0.4599983296732344):0.4182456117657932):0.8419343495905224);";
		STITree<Double> randTreeP = BUF.getTree(pt);
		ComputeLikelihood CLP = new ComputeLikelihood(0.001, 0.2, 0.1, randTreeP, model);
		double[] poly_val = new double[nMut];
		for(int i=0;i<=nMut -1;i++)
			poly_val[i] = CLP.getGenotypeProbability(obsGTObsMatrix.get(i), 1);
		double poly_likelihood = 0;
		for (double d: poly_val){
			double log_d;
			if (d == 0)
				log_d = -10000;
			else
				log_d = Math.log(d);
			poly_likelihood += log_d;
		}
		System.out.println(poly_likelihood);
	}

}
