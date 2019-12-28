package org.panda.proteinsitemotif;

import org.panda.utility.Tuple;

import java.io.IOException;
import java.util.*;

public class ExhaustiveMotifGenerator
{
	private Sequences origSeq;
	private Sequences sequences;
	private int selectStatus;
	private Character centralAA;
	private double fdrThr;

	public ExhaustiveMotifGenerator(Sequences sequences, int selectStatus, Character centralAA, double fdrThr)
	{
		this.sequences = sequences.copy();
		this.origSeq = sequences;
		this.selectStatus = selectStatus;
		this.centralAA = centralAA;
		this.fdrThr = fdrThr;
	}

	public Set<Motif> run()
	{
		Set<Motif> motifs = runRecursive(new Motif(0, centralAA, true), Deviation.POSITIVE);

		Set<Motif> extended = new HashSet<>();
		motifs.stream().filter(m -> motifs.stream().noneMatch(m::isGeneralOf)).forEach(m ->
			extended.addAll(runRecursive(m, Deviation.NEGATIVE)));
		motifs.addAll(extended);

		return motifs;
	}

	private Set<Motif> runRecursive(Motif motif, Deviation dev)
	{
		System.out.println("motif = " + motif);
		FreqMatrix bgM = new FreqMatrix(sequences.seqWidth, centralAA);
		bgM.addSequences(sequences, motif, null);

		FreqMatrix selM = new FreqMatrix(sequences.seqWidth, centralAA);
		selM.addSequences(sequences, motif, selectStatus);

		if (bgM.seqCnt == selM.seqCnt) return Collections.emptySet();

		DeviationDetector dd = new DeviationDetector(bgM, selM, fdrThr, dev);
		Map<AAInPos, Tuple> signif = dd.getSignificantDeviations();
		System.out.println("signif.size() = " + signif.size());

		Set<Motif> results = new HashSet<>();

		for (AAInPos aap : signif.keySet())
		{
			Motif sub = new Motif(motif);
			sub.add(aap);
			results.add(sub);

			results.addAll(runRecursive(sub, dev));
		}

		return results;
	}
}
