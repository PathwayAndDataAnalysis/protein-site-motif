package org.panda.proteinsitemotif;

import org.panda.utility.Tuple;

import java.io.IOException;
import java.util.*;

public class GreedyIterativeMotifGenerator
{
	Sequences origSeq;
	Sequences sequences;
	int selectStatus;
	Character centralAA;
	double fdrThr;

	public GreedyIterativeMotifGenerator(Sequences sequences, int selectStatus, Character centralAA, double fdrThr)
	{
		this.sequences = sequences.copy();
		this.origSeq = sequences;
		this.selectStatus = selectStatus;
		this.centralAA = centralAA;
		this.fdrThr = fdrThr;
	}

	public Set<Motif> run()
	{
		Set<Motif> motifs = new HashSet<>();
		Motif x;
		do
		{
			x = generateOne();

			if (x != null)
			{
				motifs.add(x);
				sequences.removeMatching(x);
			}
		}
		while (x != null);

		return motifs;
	}

	private Motif generateOne()
	{
		Motif motif = new Motif(0, centralAA, true);

		searchIteratively(motif, Deviation.POSITIVE);
		if (motif.size() > 1) searchIteratively(motif, Deviation.NEGATIVE);

		if (motif.size() > 1) return motif;

		return null;
	}

	private void searchIteratively(Motif motif, Deviation devDir)
	{
		int prevSize = 0;
		int currentSize = motif.size();

		while (prevSize < currentSize)
		{
			prevSize = currentSize;

			FreqMatrix bgM = new FreqMatrix(sequences.seqWidth, centralAA);
			bgM.addSequences(sequences, motif, null);

			FreqMatrix selM = new FreqMatrix(sequences.seqWidth, centralAA);
			selM.addSequences(sequences, motif, selectStatus);

			DeviationDetector dd = new DeviationDetector(bgM, selM, fdrThr, devDir);
			Map<AAInPos, Tuple> signif = dd.getSignificantDeviations();

			Optional<AAInPos> first = signif.keySet().stream().min(Comparator.comparing(o -> signif.get(o).p));
			if (first.isPresent())
			{
				AAInPos aap = first.get();
				motif.add(aap.pos, aap.aa, signif.get(aap).v > 0);
			}

			currentSize = motif.size();
		}
	}
}
