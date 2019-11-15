package org.panda.proteinsitemotif;

import org.panda.utility.Progress;
import org.panda.utility.statistics.ChiSquare;
import org.panda.utility.statistics.FDR;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class EnrichmentTester
{
	Map<String, Set<String>> posCtrl;
	Map<String, Set<String>> negCtrl;

	public EnrichmentTester() throws IOException
	{
		load();
	}

	public Map<String, Double> findEnrichedTFs(Sequences seq, Motif motif, Motif base, double fdrThr, int iterations, int minTotalTrg)
	{
		Map<String, Set<String>>[] baseTotals = readEffectorCounts(seq, base, null);
		Set<String>[] consider = new Set[2];
		for (int i = 0; i < 2; i++)
		{
			int index = i;
			consider[i] = baseTotals[i].keySet().stream().filter(eff -> baseTotals[index].get(eff).size() >= minTotalTrg).collect(Collectors.toSet());
		}

		Map<String, List<Integer>>[] baseMaps = readBackgroundDistribution(seq, motif, base, iterations, consider);
		Map<String, Set<String>>[] selectMaps = readEffectorCounts(seq, motif, consider);

		Map<String, Double> results = new HashMap<>();

		for (int i = 0; i < 2; i++)
		{
			Map<String, Double> pvals = new HashMap<>();

			for (String eff : selectMaps[i].keySet())
			{
				int seen = selectMaps[i].get(eff).size();
				long rndCnt = baseMaps[i].get(eff).stream().filter(num -> num >= seen).count();
				if (rndCnt == 0) rndCnt++;

				double p = rndCnt / (double) iterations;
				pvals.put(eff, p);
			}

			List<String> signif = FDR.select(pvals, null, fdrThr);
			int index = i;
			signif.forEach(eff -> results.put(eff, index == 0 ? pvals.get(eff) : -pvals.get(eff)));

			if (i == 1 && !signif.isEmpty())
			{
				System.out.println("Negative effectors found!");
			}
		}

		return results;
	}

	private Map<String, Set<String>>[] readEffectorCounts(Sequences seq, Motif motif, Set<String>[] consider)
	{
		Map<String, Set<String>> posMap = new HashMap<>();
		Map<String, Set<String>> negMap = new HashMap<>();

		seq.getClassifiedSeqStream(motif, null).forEach(cs ->
		{
			String key = cs.gene + " " + cs.site;
			if (posCtrl.containsKey(key)) updateCounts(posCtrl.get(key), cs.gene, posMap, consider == null ? null : consider[0]);
			if (negCtrl.containsKey(key)) updateCounts(negCtrl.get(key), cs.gene, negMap, consider == null ? null : consider[1]);
		});

		return new Map[]{posMap, negMap};
	}

	private Map<String, List<Integer>>[] readBackgroundDistribution(Sequences seq, Motif motif, Motif base, int iter, Set<String>[] consider)
	{
		Map<String, List<Integer>> posMap = new HashMap<>();
		Map<String, List<Integer>> negMap = new HashMap<>();

		List<ClassifiedSequence> bgList = seq.getClassifiedSeqStream(base, null).collect(Collectors.toList());

		bgList.forEach(cs ->
		{
			String key = cs.gene + " " + cs.site;
			if (posCtrl.containsKey(key)) posCtrl.get(key).stream().filter(consider[0]::contains).forEach(eff -> posMap.put(eff, new ArrayList<>()));
			if (negCtrl.containsKey(key)) negCtrl.get(key).stream().filter(consider[1]::contains).forEach(eff -> negMap.put(eff, new ArrayList<>()));
		});

		int selectSize = (int) bgList.stream().filter(cs -> cs.matches(motif)).count();

//		Progress prg = new Progress(iter, "Generating background");
		for (int i = 0; i < iter; i++)
		{
			Collections.shuffle(bgList);
			Map<String, Set<String>> posTrgGenes = new HashMap<>();
			Map<String, Set<String>> negTrgGenes = new HashMap<>();

			bgList.stream().limit(selectSize).forEach(cs ->
			{
				String key = cs.gene + " " + cs.site;
				if (posCtrl.containsKey(key)) updateCounts(posCtrl.get(key), cs.gene, posTrgGenes, consider[0]);
				if (negCtrl.containsKey(key)) updateCounts(negCtrl.get(key), cs.gene, negTrgGenes, consider[1]);
			});

			addToDistribution(posMap, posTrgGenes);
			addToDistribution(negMap, negTrgGenes);
//			prg.tick();
		}

		return new Map[]{posMap, negMap};
	}

	private void addToDistribution(Map<String, List<Integer>> distMap, Map<String, Set<String>> currentTargs)
	{
		for (String eff : currentTargs.keySet())
		{
			List<Integer> dist = distMap.get(eff);
			if (dist == null)
			{
				dist = new ArrayList<>();
				dist.add(currentTargs.get(eff).size());
				distMap.put(eff, dist);
			}
			else
			{
				distMap.get(eff).add(currentTargs.get(eff).size());
			}
		}
	}

	private void updateCounts(Set<String> effs, String gene, Map<String, Set<String>> map, Set<String> consider)
	{
		for (String eff : effs)
		{
			if (consider == null || consider.contains(eff))
			{
				if (!map.containsKey(eff)) map.put(eff, new HashSet<>());
				map.get(eff).add(gene);
			}
		}
	}

	private void load() throws IOException
	{
		posCtrl = new HashMap<>();
		negCtrl = new HashMap<>();

		Files.lines(Paths.get("/Users/ozgun/Documents/Data/causal-priors.txt")).map(l -> l.split("\t"))
			.filter(t -> t[1].startsWith("p") || t[1].startsWith("de")).filter(t -> t.length >= 5 && !t[4].isEmpty())
			.forEach(t ->
		{
			String eff = t[0];
			String trg = t[2];

			String sites = t[4];

			Map<String, Set<String>> map = t[1].startsWith("p") ? posCtrl : negCtrl;

			for (String site : sites.split(";"))
			{
				String key = trg + " " + site;
				if (!map.containsKey(key)) map.put(key, new HashSet<>());
				map.get(key).add(eff);
			}
		});
	}
}
