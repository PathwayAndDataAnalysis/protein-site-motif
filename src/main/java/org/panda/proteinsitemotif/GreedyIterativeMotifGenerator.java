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

	public void run() throws IOException
	{
		while (generateOne());
	}

	private boolean generateOne() throws IOException
	{
		Motif motif = new Motif(0, centralAA, true);

		int prevSize = 0;
		int currentSize = motif.size();

		while (prevSize < currentSize)
		{
			prevSize = currentSize;

			FreqMatrix bgM = new FreqMatrix(sequences.seqWidth, centralAA);
			bgM.addSequences(sequences, motif, null);

			FreqMatrix selM = new FreqMatrix(sequences.seqWidth, centralAA);
			selM.addSequences(sequences, motif, selectStatus);

			DeviationDetector dd = new DeviationDetector(bgM, selM, fdrThr, Deviation.POSITIVE);
			Map<AAInPos, Tuple> signif = dd.getSignificantDeviations();

			Optional<AAInPos> first = signif.keySet().stream()
//				.filter(o -> signif.get(o).v > 0)
				.min(Comparator.comparing(o -> signif.get(o).p));
			if (first.isPresent())
			{
				AAInPos aap = first.get();
				motif.add(aap.pos, aap.aa, signif.get(aap).v > 0);
			}

			currentSize = motif.size();
		}

		if (motif.size() > 1)
		{
			System.out.print(motif.toString(sequences.seqWidth, null));//*/Constants.REV_REDUC_MAP));
			EnrichmentTester et = new EnrichmentTester();
			Map<String, Double> enrichedTFs = et.findEnrichedTFs(origSeq, motif, new Motif(0, centralAA, true), 0.1, 10000, 5);
			enrichedTFs.keySet().stream().sorted(Comparator.comparing(tf -> Math.abs(enrichedTFs.get(tf)))).forEach(tf -> System.out.print(" " + (enrichedTFs.get(tf) < 0 ? "-" : "") + tf));// + ("|" + enrichedTFs.get(tf))));
			System.out.println();
			sequences.removeMatching(motif);
			return true;
		}

		return false;
	}
}
