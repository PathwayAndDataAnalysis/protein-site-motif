package org.panda.proteinsitemotif;

import org.panda.utility.Tuple;
import org.panda.utility.statistics.ChiSquare;
import org.panda.utility.statistics.FDR;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DeviationDetector
{
	FreqMatrix allFM;
	FreqMatrix selFM;

	private double fdrThr;
	private Deviation deviationType;

	public DeviationDetector(FreqMatrix allFM, FreqMatrix selFM, double fdrThr, Deviation deviationType)
	{
		this.allFM = allFM;
		this.selFM = selFM;
		this.fdrThr = fdrThr;
		this.deviationType = deviationType;
	}

	public Map<AAInPos, Tuple> getSignificantDeviations()
	{
		int size = allFM.seqCnt;
		int select = selFM.seqCnt;
		Map<AAInPos, Double> pvals = new HashMap<>();
		Map<AAInPos, Double> limits = new HashMap<>();

		for (Integer pos : allFM.matrix.keySet())
		{
			Map<Character, Integer> aaMap = allFM.matrix.get(pos);

			for (Character aa : aaMap.keySet())
			{
				AAInPos ap = new AAInPos(pos, aa);

				int feat = allFM.matrix.get(pos).get(aa);
				int featSel = selFM.matrix.get(pos).getOrDefault(aa, 0);

				double pval = deviationType == Deviation.BOTHWAYS ? ChiSquare.testDependence(size, feat, select, featSel) :
					deviationType == Deviation.POSITIVE ? ChiSquare.testEnrichment(size, feat, select, featSel) :
						ChiSquare.testExclusivity(size, feat, select, featSel);

				if (Double.isNaN(pval)) continue;

				double expected = feat * (select / (double) size);
				double limit = 0;

				if ((deviationType == Deviation.BOTHWAYS && expected < featSel) || deviationType == Deviation.POSITIVE)
				{
					limit = deviationType == Deviation.BOTHWAYS ? ChiSquare.testDependence(size, feat, select, Math.min(select, feat)) :
						ChiSquare.testEnrichment(size, feat, select, Math.min(select, feat));

					ap.prefer = true;
				}
				else if ((deviationType == Deviation.BOTHWAYS && expected >= featSel) || deviationType == Deviation.NEGATIVE)
				{
					limit = deviationType == Deviation.BOTHWAYS ? ChiSquare.testDependence(size, feat, select, Math.max(0, select - (size - feat))) :
						ChiSquare.testExclusivity(size, feat, select, Math.max(0, select - (size - feat)));
					ap.prefer = false;
				}

				pvals.put(ap, pval);
				limits.put(ap, limit);
			}
		}

		Map<AAInPos, Tuple> results = new HashMap<>();

		List<AAInPos> signif = FDR.select(pvals, limits, fdrThr);
		signif.forEach(ap -> results.put(ap, new Tuple(ap.prefer ? 1: -1, pvals.get(ap))));

		return results;
	}

}
