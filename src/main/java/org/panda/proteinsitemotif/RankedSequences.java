package org.panda.proteinsitemotif;

import org.panda.resource.HGNC;
import org.panda.resource.UniProtSequence;
import org.panda.utility.ArrayUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.TermCounter;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class RankedSequences
{
	List<Sequence> seqList;
	int seqWidth;
	int halfW;

	public RankedSequences(int seqWidth)
	{
		this.seqWidth = seqWidth;
		this.halfW = seqWidth / 2;
		seqList = new ArrayList<>();
	}

	public void loadFromCPFiles(String file, String symColName, String siteColName, String signedPColName) throws IOException
	{
		Map<Sequence, Double> seqToP = new HashMap<>();

		String[] header = FileUtil.readHeader(file);
		int symInd = ArrayUtil.indexOf(header, symColName);
		int siteInd = ArrayUtil.indexOf(header, siteColName);
		int pInd = ArrayUtil.indexOf(header, signedPColName);

		int maxColInd = Math.max(Math.max(symInd, siteInd), pInd);

		Files.lines(Paths.get(file)).skip(1).map(l -> l.split("\t")).filter(t -> t.length > maxColInd && !t[siteInd].isEmpty()).forEach(t ->
		{
			String sym = t[symInd].split(" ")[0];
			String upID = HGNC.get().getUniProt(sym);
			if (upID == null) return;

			String sites = t[siteInd].split(" ")[0];

			double val = Double.valueOf(t[pInd]);

			if (val == 0)
			{
				val = 1e-20;
				if (t[pInd].startsWith("-")) val = -val;
			}

			for (String siteStr : sites.split("\\|"))
			{
				String aa = siteStr.substring(0, 1);
				int site = Integer.valueOf(siteStr.substring(1));

				if (site > halfW)
				{
					String aaFound = UniProtSequence.get().getAminoacidAt(upID, site);
					if (aaFound != null && aaFound.equals(aa))
					{
						String seq = UniProtSequence.get().getSeqAround(upID, site, seqWidth);

						if (seq != null)
						{
							seqToP.put(new Sequence(seq, sym, siteStr), val);
						}
					}
				}
			}
		});

		// Sort from most significantly upregulated to least sigificantly upregulated
		seqList = seqToP.keySet().stream().sorted((s1, s2) ->
		{
			double p1 = seqToP.get(s1);
			double p2 = seqToP.get(s2);

			if (p1 == p2)
			{
				return 0;
			}
			else if (Math.signum(p1) != Math.signum(p2))
			{
				return Math.signum(p1) > 0 ? -1 : 1;
			}
			else if (Math.signum(p1) > 0)
			{
				return p1 < p2 ? -1 : 1;
			}
			else
			{
				return p1 > p2 ? -1 : 1;
			}
		}).collect(Collectors.toList());
	}

	public RankedSequences reduce(AAReducer reducer)
	{
		RankedSequences copy = new RankedSequences(seqWidth);
		seqList.forEach(cs -> copy.seqList.add(new Sequence(reducer.reduce(cs.seq), cs.gene, cs.site)));
		return copy;
	}

	public RankedSequences copy()
	{
		RankedSequences copy = new RankedSequences(seqWidth);
		seqList.forEach(cs -> copy.seqList.add(new Sequence(cs.seq, cs.gene, cs.site)));
		return copy;
	}

	public List<Sequence> removeMatching(Motif motif)
	{
		List<Sequence> match = seqList.stream().filter(cs -> cs.matches(motif)).collect(Collectors.toList());
		seqList.removeAll(match);
		return match;
	}

	public void printAAFreqs()
	{
		TermCounter tc = new TermCounter();
		seqList.forEach(cs ->
		{
			for (int i = 0; i < cs.seq.length(); i++)
			{
				tc.addTerm(cs.seq.substring(i, i+1));
			}
		});
		tc.print();
	}
}
