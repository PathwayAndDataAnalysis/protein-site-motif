package org.panda.proteinsitemotif;

import org.panda.resource.HGNC;
import org.panda.resource.UniProtSequence;
import org.panda.utility.TermCounter;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Sequences
{
	List<ClassifiedSequence> seqList;
	int seqWidth;
	int halfW;

	public Sequences(int seqWidth)
	{
		this.seqWidth = seqWidth;
		this.halfW = seqWidth / 2;
		seqList = new ArrayList<>();
	}

	public void loadFromAslanCPFiles(String file) throws IOException
	{
		Files.lines(Paths.get(file)).skip(1).map(l -> l.split("\t")).filter(t -> t.length > 4 && !t[2].isEmpty()).forEach(t ->
		{
			String sym = t[1].split(" ")[0];
			String upID = HGNC.get().getUniProt(sym);
			if (upID == null) return;

			String sites = t[2].split(" ")[0];

			double val = Double.valueOf(t[4]);

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
							int status = (int) Math.signum(val);
							seqList.add(new ClassifiedSequence(seq, status, sym, siteStr));
						}
					}
				}
			}
		});
	}

	public Stream<String> getSeqStream(Motif motif, Integer selectionStatus)
	{
		return getClassifiedSeqStream(motif, selectionStatus).map(cs -> cs.seq);
	}

	public Stream<ClassifiedSequence> getClassifiedSeqStream(Motif motif, Integer selectionStatus)
	{
		return seqList.stream().filter(cs -> selectionStatus == null || cs.status == selectionStatus)
			.filter(cs -> motif == null || cs.matches(motif));
	}

	public Stream<ClassifiedSequence> getRandomizedClassifiedSeqStream(Motif motif, Motif base, Integer selectionStatus)
	{
		long count = seqList.stream().filter(cs -> selectionStatus == null || cs.status == selectionStatus)
			.filter(cs -> cs.matches(motif)).count();

		List<ClassifiedSequence> list = seqList.stream().filter(cs -> selectionStatus == null || cs.status == selectionStatus)
			.filter(cs -> base == null || cs.matches(base)).collect(Collectors.toList());

		Collections.shuffle(list);
		list = list.subList(0, (int) count);
		return list.stream();
	}

	public Sequences reduce(AAReducer reducer)
	{
		Sequences copy = new Sequences(seqWidth);
		seqList.forEach(cs -> copy.seqList.add(new ClassifiedSequence(reducer.reduce(cs.seq), cs.status, cs.gene, cs.site)));
		return copy;
	}

	public Sequences copy()
	{
		Sequences copy = new Sequences(seqWidth);
		seqList.forEach(cs -> copy.seqList.add(new ClassifiedSequence(cs.seq, cs.status, cs.gene, cs.site)));
		return copy;
	}

	public List<ClassifiedSequence> removeMatching(Motif motif)
	{
		List<ClassifiedSequence> match = seqList.stream().filter(cs -> cs.matches(motif)).collect(Collectors.toList());
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

	public Set<String> getGenesAndSites(Motif motif, int status)
	{
		return getClassifiedSeqStream(motif, status).map(cs -> cs.gene + "-" + cs.site).collect(Collectors.toSet());
	}

	public void randomizeSelection()
	{
		Set<Integer> set = seqList.stream().map(s -> s.status).collect(Collectors.toSet());
		System.out.println("set = " + set);

		long upCnt = seqList.stream().filter(s -> s.status == 1).count();
		long dwCnt = seqList.stream().filter(s -> s.status == -1).count();

		long chCnt = upCnt + dwCnt;

		List<ClassifiedSequence> list = new ArrayList<>(seqList);
		Collections.shuffle(list);
		List<ClassifiedSequence> chList = list.subList(0, (int) chCnt);
		List<ClassifiedSequence> upList = chList.subList(0, (int) upCnt);

		upList.forEach(s -> s.status = 1);
		chList.stream().filter(s -> !upList.contains(s)).forEach(s -> s.status = -1);
		list.stream().filter(s -> !chList.contains(s)).forEach(s -> s.status = 0);
	}
}
