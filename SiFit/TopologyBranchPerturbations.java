/**
 * 
 */
package SiFit;

import java.util.ArrayList;
import java.util.Random;

import cern.jet.random.tdouble.Normal;
import cern.jet.random.tdouble.engine.DoubleMersenneTwister;
import cern.jet.random.tdouble.engine.DoubleRandomEngine;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

/**
 * @author hz22
 *
 */
public class TopologyBranchPerturbations extends BasicUtilityFunctions {
	private Random _rng;
	private boolean _topologyHasChanged = false;
	private double _extension = 0.5;
	private double _tuning = 2.0 * Math.log(1.6);
	private double _tuningNNI = 2.0 * Math.log(1.20);
	private double _local = 2.0 * Math.log(1.05);
	private double _sigma = 0.03;
	private double _proposalRatio = 0.0;
	public double _maxV = 20.0;
	
	
	/**
	 * Variety of constructors. Prefer the first or the fourth one.
	 */
	public TopologyBranchPerturbations(Random rng) {
		super(rng);
		this._rng = rng;
	}

	public TopologyBranchPerturbations() {
		this._rng = new Random();
	}
	
	/**
	 * Bunch of setter functions follow:
	 */
	
	/*
	 * maxBranchLength: Maximum allowed branch lengths for consideration.
	 */
	public void setMaxBranchLength(double maxBranchLength) {
		this._maxV = maxBranchLength;
	}
	
	/*
	 * LOCAL tuning: Used in LOCAL(Larget and Simon)
	 */	
	public void setLocalParameter(double local) {
		this._local = local;
	}

	/*
	 * tuningParameter: Used in Multiplier
	 */
	public void setTuningParameter(double tuningParameter) {
		this._tuning = tuningParameter;
	}
	
	/*
	 * sigmaParameter:  Used in CC.
	 */
	public void setSigmaParameter (double sigmaParameter) {
		this._sigma = sigmaParameter;
	}
	
	
	/**
	 * Returns proposal ratio for the last perturbation.
	 */
	public double getProposalRatio() {
		return _proposalRatio;
	}
	
	/**
	 * Indicator of whether last perturbation changed topology.
	 */
	public boolean hasTopologyChanged(){
		return _topologyHasChanged;
	}

	
	/**
	 *  Tree Perturbation Code Follows. 
	 *  Implementation of stNNI, SPR, SS.
	 */	

	/*
	 * stNNI: stochastic NNI.
	 */
	public STITree<Double> stNNI (STITree<Double> otree) {

		//Create copy of the original tree.
		STITree<Double> tree = new STITree<Double> (otree);
		
		STINode<Double> u = selectRandomNode(tree, false, false, _rng);
		STINode<Double> v = u.getParent();
		
		// Get children of "u".
		STINode<Double> a = getRandomChild(u);
		STINode<Double> b = getOtherChild(u,a);
		
		// Get other child of "v".
		STINode<Double> c = getOtherChild(v,u);
		
		double randDouble = _rng.nextDouble();
		//Swap "a" with "c".
		if (randDouble < 1.0/3) {
			v.adoptChild(a);
			u.adoptChild(c);	
			_proposalRatio = 1;
			_topologyHasChanged = true;
		}
		//Swap "b" with "c".
		else if (randDouble < 2.0/3){
			v.adoptChild(b);
			u.adoptChild(c);
			_proposalRatio = 1;
			_topologyHasChanged = true;
		}
		else {
			_topologyHasChanged = false;
//			ArrayList<STINode<Double>> nodes = new ArrayList<STINode<Double>>();
//			nodes.add(u);
//			nodes.add(a);
//			nodes.add(b);
//			nodes.add(c);
//			tree = multiplier (tree, nodes, _tuningNNI);
		}
//		System.out.println(tree.toNewick());
//		System.out.println(_topologyHasChanged);
		ArrayList<STINode<Double>> nodes = new ArrayList<STINode<Double>>();
		nodes.add(u);
		nodes.add(a);
		nodes.add(b);
		nodes.add(c);
		tree = multiplier (tree, nodes, _tuningNNI);

		return tree;
	}
	
	/*
	 * Changes a list of branch lengths.
	 */
	public STITree<Double> multiplier (STITree<Double> otree) {
		//Create copy of the original tree.
		STITree<Double> tree = new STITree<Double> (otree);
		
		STINode<Double> u = selectRandomNode(tree, false, false, _rng);
		STINode<Double> v = u.getParent();	
		
		// Get children of "u".
		STINode<Double> a = getRandomChild(u);
		STINode<Double> b = getOtherChild(u,a);

		// Get other child of "v".
		STINode<Double> c = getOtherChild(v,u);		
		
		ArrayList<STINode<Double>> nodes = new ArrayList<STINode<Double>>();
		nodes.add(u);
		nodes.add(a);
		nodes.add(b);
		nodes.add(c);
		tree = multiplier (tree, nodes, _tuning);
		
		return tree;
	}
	
	public STITree<Double> changeBranchLength(STITree<Double> otree){
		//Create copy of the original tree.
		STITree<Double> tree = new STITree<Double> (otree);
		STINode<Double> u = tree.selectRandomNode(true, false);
//		STINode<Double> u = selectRandomNode(tree, true, false, _rng);
		double oldB = u.getParentDistance();
		double newB = 0.0, m = 0.0;
		do {
			m = Math.exp(_tuning * (_rng.nextDouble()-0.5));
			newB = oldB * m;
		} while (newB >= _maxV);
		u.setParentDistance(newB);
		return tree;
	}
	
	/**
	 * Function that performs "multiplier"
	 */
	private STITree<Double> multiplier (STITree<Double> otree, ArrayList<STINode<Double>> nodes, double tuningMultiplier) {
		//Create copy of the original tree.
		STITree<Double> tree = new STITree<Double> (otree);
		
		_proposalRatio = 1.0;
//		_topologyHasChanged = false;
		for (int i = 0; i < nodes.size(); i++) {
			STINode<Double> curNode = tree.getNode(nodes.get(i).getID());
			double oldB = curNode.getParentDistance();
			double newB = 0.0, m = 0.0;
			do {
				m = Math.exp(tuningMultiplier * (_rng.nextDouble()-0.5));
				newB = oldB * m;
			} while (newB >= _maxV);
			
			curNode.setParentDistance(newB);
			_proposalRatio *= m;
		}

		return tree;
	}
	
	/*
	 * Continuous Change: Jow et al. 2002
	 */
	public STITree<Double> CC (STITree<Double> otree) {
		//Create copy of the original tree.
		STITree<Double> tree = new STITree<Double> (otree);
	
		STINode<Double> u = selectRandomNode(tree, false, false, _rng);
		STINode<Double> v = getRandomChild(u);
		
		//Normal Distribution: N(0, sigma)
		DoubleRandomEngine rngEngine = new DoubleMersenneTwister(_rng.nextInt());
		Normal rngN = new Normal(0.0,_sigma,rngEngine);
		
		double newV = v.getParentDistance() + rngN.nextDouble();
		
		//Branch length should not exceed maximum allowed value.
		if (Math.abs(newV) < _maxV) {
			v.setParentDistance(Math.abs(newV));
		}

		_proposalRatio = 1;
		_topologyHasChanged = false;
		
		//Change Topology
		if (newV < 0 && !v.isLeaf()) {
			_topologyHasChanged = true;
			STINode<Double> a = getRandomChild(v);
			STINode<Double> c = getOtherChild(u,v);
			
			v.adoptChild(c);
			u.adoptChild(a);
		}
			
		return tree;
	}
	
	/**
	 * Change topology using Random Subtree Pruning and Regrafting (rSPR) move.
	 * @param otree
	 * @return
	 */
	public STITree<Double> rSPR (STITree<Double> otree) {
		//Create copy of the original tree.
		STITree<Double> tree = new STITree<Double> (otree);
		
		STINode<Double> u = selectRandomNode(tree, false, false, _rng);
		//Node "u" has two children: "v" and "w".
		STINode<Double> v = getRandomChild(u);
		STINode<Double> w = getOtherChild(u,v);
		
		//We want to prune subtree rooted at "v".
		//We want to re-graft the pruned subtree between node "x" and it's child "cx".
		STINode<Double> x = null;
		STINode<Double> cx = null;
		
		//"x" cannot be any node in the pruned subtree OR the same node from which it was pruned.
		do {
			x = selectRandomNode(tree, false, true, _rng);
		} while(v.isAncestor(x) || (u.getID() == x.getID()));

		u.getParent().adoptChild(w);
		x.adoptChild(u);

		//"cx" cannot be: "w" as it gives rise to the same topology again; "u" as it is part of re-grafting.
		do {
			cx = getRandomChild(x);
		} while (cx.getID() == u.getID() || cx.getID() == w.getID());
		u.adoptChild(cx);

		//Use Multiplier:
		ArrayList<STINode<Double>> nodes = new ArrayList<STINode<Double>>();
		nodes.add(v);
		nodes.add(cx);
		tree = multiplier (tree, nodes, _tuning);	
		
		_proposalRatio *= 1.0;
		_topologyHasChanged = true;
		return tree;
	}
	
	/**
	 * Change topology using Extending Subtree Pruning and Regrafting (eSPR) move.
	 * @param otree
	 * @return
	 */
	public STITree<Double> eSPR (STITree<Double> otree) {
		STITree<Double> tree = new STITree<Double> (otree);
		
		//New change: Prevent waste of computation when "w" is leaf.
		STINode<Double> u, v, w;
		do {
			u = selectRandomNode(tree, false, false, _rng);
			//Node "u" has two children: "v" and "w".
			v = getRandomChild(u);
			w = getOtherChild(u,v);			
		} while (w.isLeaf());		
		
		boolean isConstrainedPrunedBranch = v.isLeaf();
		boolean isConstrainedRegraftingBranch = false;
		
		STINode<Double> x = null;
		STINode<Double> cx = null;
		
		_proposalRatio = 1;
		_topologyHasChanged = false;
		//When "w" is a leaf, no change possible.
		if (!w.isLeaf()) {
			x = w;
			cx = getRandomChild(w); 
			while (!cx.isLeaf()) {
				if (_rng.nextDouble() < (1 - _extension)) {
					break;
				}
				x = cx;
				cx = getRandomChild(cx);
			}
			
			u.getParent().adoptChild(w);
			x.adoptChild(u);
			u.adoptChild(cx);			
		
			isConstrainedRegraftingBranch = cx.isLeaf();	
			//Set proposal ratio.
			if (isConstrainedPrunedBranch && !isConstrainedRegraftingBranch) {
				_proposalRatio = 1 - _extension;
			}
			else if (!isConstrainedPrunedBranch && isConstrainedRegraftingBranch) {
				_proposalRatio = 1 / (1 - _extension);
			}
			
			
			//Use Multiplier:
			double topologyRatio = _proposalRatio;
			ArrayList<STINode<Double>> nodes = new ArrayList<STINode<Double>>();
			nodes.add(v);
			nodes.add(cx);
			tree = multiplier (tree, nodes, _tuning);	
			_proposalRatio *= topologyRatio;
			_topologyHasChanged = true;
			
		}
		
		return tree;
	}	
	
	/* 
	 * Change topology using Random Subtree Swapping (rSTS) move.
	 */
	public STITree<Double> rSTS (STITree<Double> otree) {
		//Create copy of the original tree.
		STITree<Double> tree = new STITree<Double> (otree);
		
		STINode<Double> u = selectRandomNode(tree, false, false, _rng);
		//Node "u" has two children: "v" and "w". As "w" is not involved, it is not defined.
		STINode<Double> v = getRandomChild(u);
		
		// One of the two subtrees taking part in swap is rooted at "v".
		STINode<Double> x = null;
		STINode<Double> cx = null;
		
		//"x" cannot be any node in the pruned subtree OR the same node from which it was pruned.
		do {
			x = selectRandomNode(tree, false, true, _rng);
		} while(v.isAncestor(x) || (u.getID() == x.getID()));

		do {
			cx = getRandomChild(x);
		} while (cx.getID() == u.getID() || cx.isAncestor(v));
		
		
		//We want to swap "v" (child of "u") with "cx" (child of "x")
		u.adoptChild(cx);
		x.adoptChild(v);

		
		//Use Multiplier:
		ArrayList<STINode<Double>> nodes = new ArrayList<STINode<Double>>();
		nodes.add(v);
		nodes.add(cx);
		tree = multiplier (tree, nodes, _tuning);	
		
		_proposalRatio *= 1.0;
		_topologyHasChanged = true;
		return tree;			
	}

	
	/* 
	 * Change topology using Extending Subtree Swapping (eSTS) move.
	 */
	public STITree<Double> eSTS (STITree<Double> otree) {	
		//Create copy of the original tree.
		STITree<Double> tree = new STITree<Double> (otree);		
		
		//New change: Prevent waste of computation when "w" is leaf.
		STINode<Double> u, v, w;
		do {
			u = selectRandomNode(tree, false, false, _rng);
			//Node "u" has two children: "v" and "w".
			v = getRandomChild(u);
			w = getOtherChild(u,v);			
		} while (w.isLeaf());
		
		boolean isConstrainedBranchOne = v.isLeaf();
		boolean isConstrainedBranchTwo = false;
		
		STINode<Double> x = null;
		STINode<Double> cx = null;
		
		_proposalRatio = 1;
		//When "w" is a leaf, no change possible.
		if (!w.isLeaf()) {
			x = w;
			cx = getRandomChild(w); 
			while (!cx.isLeaf()) {
				if (_rng.nextDouble() < (1 - _extension)) {
					break;
				}
				x = cx;
				cx = getRandomChild(cx);
			}
			
			u.adoptChild(cx);
			x.adoptChild(v);
			
			isConstrainedBranchTwo = cx.isLeaf();	
			//Set proposal ratio.
			if (isConstrainedBranchOne && !isConstrainedBranchTwo) {
				_proposalRatio = 1 - _extension;
			}
			else if (!isConstrainedBranchOne && isConstrainedBranchTwo) {
				_proposalRatio = 1 / (1 - _extension);
			}

			//Use Multiplier:
			double topologyRatio = _proposalRatio;
			ArrayList<STINode<Double>> nodes = new ArrayList<STINode<Double>>();
			nodes.add(v);
			nodes.add(cx);
			tree = multiplier (tree, nodes, _tuning);	
			_proposalRatio *= topologyRatio;
			_topologyHasChanged = true;				
		}

		return tree;
	}	
	
	public STITree<Double> addEdge (STITree<Double> intree) {
		//Create copy of the original tree.
		STITree<Double> tree = new STITree<Double> (intree);
		ArrayList<STINode<Double>> nonBinaryPars = new ArrayList<>();
		for (STINode<Double> intNode : tree.getNodes()){
			if (!intNode.isLeaf()){
				if (intNode.getChildCount() > 2)
					nonBinaryPars.add(intNode);
			}
		}
		for (STINode<Double> node: nonBinaryPars)
			System.out.println(node.getName());
		return null;
		
	}
	
	/**
	 * proposes Tree 
	 * @param curTree
	 * @param nCell
	 * @return
	 */
	public STITree<Double> proposeTree(STITree<Double> curTree, int nCell){
		double randDouble = _rng.nextDouble();
		STITree<Double> newTree = curTree;
		
		//For smaller problems rely on small subset of perturbations:
		if (nCell <= 7) {			
			if (randDouble < 1.0/2) {
				newTree = this.stNNI(curTree);
			}
			else {
				newTree = this.multiplier(curTree);
			}			
		}
		
		//For larger problems rely on mixture of all perturbations:
		else{
			if (randDouble < 1.0/4) {
				newTree = this.multiplier(curTree);
			}
			else if (randDouble < 3.0/5) {
				newTree = this.stNNI(curTree);
			}
			else if (randDouble < 0.7) {
//				System.out.println("eSPR happens");
				newTree = this.eSPR(curTree);
			}
			else if (randDouble < 0.8) {
				newTree = this.rSPR(curTree);
			}
			else if (randDouble < 0.9) {
				newTree = this.eSTS(curTree);
			}
			else {
				newTree = this.rSTS(curTree);
			}			
		}
		return newTree;
	}
	
	public STITree<Double> proposeTreeParsimony(STITree<Double> curTree, int nCell){
		double randDouble = _rng.nextDouble();
		STITree<Double> newTree = curTree;
		if (randDouble < 0.25){
//			System.out.println("stNNI happens");
			newTree = this.stNNI(curTree);
		}
		else if (randDouble < 0.5) {
//			System.out.println("eSPR happens");
			newTree = this.eSPR(curTree);
		}
		else if (randDouble < 0.75) {
			newTree = this.rSPR(curTree);
		}
		else {
			newTree = this.eSTS(curTree);
		}
		return newTree;
	}
	
	public static void main(String[] args){
		String t = "((((a:1,b:2,c:3):4,d:1.1):6,e:2.5):2.9,f:8):0.5;";
//		BasicUtilityFunctions BUF = new BasicUtilityFunctions();
		TopologyBranchPerturbations TBP = new TopologyBranchPerturbations();
		STITree<Double> polytree = TBP.getTree(t);
		STITree<Double> nuTree = TBP.changeBranchLength(polytree);
		System.out.println(nuTree);
		
//		int i = 1;
//		for (STINode<Double> node : polytree.getNodes()){
//			if (node.getName().equals("")){
//				node.setName("IN" + Integer.toString(i));
//				i++;
//			}
//		}
//		STITree<Double> otree = TBP.addEdge(polytree);
//		System.out.println(polytree.toNewick());
	}


}
