package org.panda.proteinsitemotif;

import org.panda.utility.Tuple;
import org.panda.utility.graph.DirectedGraph;

import java.util.Comparator;
import java.util.Map;

public class DeviationGraphGenerator
{
	private Sequences sequences;
	private int selectStatus;
	private Character centralAA;
	private double fdrThr;

	public DeviationGraphGenerator(Sequences sequences, int selectStatus, Character centralAA, double fdrThr)
	{
		this.sequences = sequences;
		this.selectStatus = selectStatus;
		this.centralAA = centralAA;
		this.fdrThr = fdrThr;
	}

	public void run()
	{
		Motif base = new Motif(0, centralAA, true);

		FreqMatrix bgM = new FreqMatrix(sequences.seqWidth, centralAA);
		bgM.addSequences(sequences, base, null);

		FreqMatrix selM = new FreqMatrix(sequences.seqWidth, centralAA);
		selM.addSequences(sequences, base, selectStatus);

		DeviationDetector dd = new DeviationDetector(bgM, selM, fdrThr, Deviation.BOTHWAYS);
		Map<AAInPos, Tuple> signif = dd.getSignificantDeviations();

		signif.keySet().stream().sorted(Comparator.comparingDouble(o -> signif.get(o).p)).forEach(ap-> System.out.println(ap + "  \t" + signif.get(ap).p));

		DirectedGraph graph = new DirectedGraph("Pattern", "controls-state-change-of");

		for (AAInPos aap : signif.keySet())
		{
			System.out.println("\nTesting: " + aap);
			Motif motifPos = new Motif(base);
			motifPos.add(aap);

			bgM = new FreqMatrix(sequences.seqWidth, centralAA);
			bgM.addSequences(sequences, motifPos, null);

			selM = new FreqMatrix(sequences.seqWidth, centralAA);
			selM.addSequences(sequences, motifPos, selectStatus);

			dd = new DeviationDetector(bgM, selM, fdrThr, Deviation.BOTHWAYS);
			Map<AAInPos, Tuple> sigPos = dd.getSignificantDeviations();
			System.out.println("------Pos-----");
			sigPos.keySet().stream().sorted(Comparator.comparingDouble(o -> sigPos.get(o).p)).forEach(ap-> System.out.println(ap + "  \t" + sigPos.get(ap).p));
			sigPos.keySet().forEach(ap -> graph.putRelation(aap.toString(), ap.toString()));

			Motif motifNeg = new Motif(base);
			motifNeg.add(aap.pos, aap.aa, !aap.prefer);

			bgM = new FreqMatrix(sequences.seqWidth, centralAA);
			bgM.addSequences(sequences, motifNeg, null);

			selM = new FreqMatrix(sequences.seqWidth, centralAA);
			selM.addSequences(sequences, motifNeg, selectStatus);

			dd = new DeviationDetector(bgM, selM, fdrThr, Deviation.BOTHWAYS);
			Map<AAInPos, Tuple> sigNeg = dd.getSignificantDeviations();
			System.out.println("------Neg-----");
			sigNeg.keySet().stream().sorted(Comparator.comparingDouble(o -> sigNeg.get(o).p)).forEach(ap-> System.out.println(ap + "  \t" + sigNeg.get(ap).p));
		}

		graph.write("/Users/ozgun/Documents/Analyses/motif/binary-graph.sif");
	}
}
