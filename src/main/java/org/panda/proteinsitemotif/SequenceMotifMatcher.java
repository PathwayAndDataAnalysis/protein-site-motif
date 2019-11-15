package org.panda.proteinsitemotif;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class SequenceMotifMatcher
{
	private Sequences sequences;
	private int selection;
	private String cpResultFile;
	private String motifFile;
	Set<Motif> motifs;
	Map<Motif, String> origLineMap;

	public SequenceMotifMatcher(Sequences sequences, int selection, String cpResultFile, String motifFile)
	{
		this.sequences = sequences;
		this.selection = selection;
		this.cpResultFile = cpResultFile;
		this.motifFile = motifFile;
	}

	public void run(String outFile) throws IOException
	{
		if (motifs == null) readMotifsFromFile();

		Map<String, Set<String>> explained = readExplained();

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		motifs.stream().sorted(Comparator.comparing(Motif::size)).forEach(m ->
		{
			Set<String> expSet = new HashSet<>();
			Set<String> unexpSet = new HashSet<>();
			sequences.getClassifiedSeqStream(m, selection).forEach(cs ->
			{
				String key = cs.gene + "-" + cs.site;
				if (explained.containsKey(key)) expSet.add(key);
				else unexpSet.add(key);
			});
			FileUtil.lnwrite(origLineMap.get(m), writer);
			FileUtil.lnwrite("Explained:", writer);
			expSet.stream().sorted().forEach(trg ->
			{
				FileUtil.write(" " + trg + "(", writer);
				List<String> effList = explained.get(trg).stream().sorted().collect(Collectors.toList());
				FileUtil.write(CollectionUtil.merge(effList, ",") + ")", writer);
			});
			FileUtil.lnwrite("Unexplained:", writer);
			unexpSet.stream().sorted().forEach(trg -> FileUtil.write(" " + trg, writer));
			FileUtil.write("\n", writer);
		});
		writer.close();
	}

	private void readMotifsFromFile() throws IOException
	{
		motifs = new HashSet<>();
		origLineMap = new HashMap<>();

		Files.lines(Paths.get(motifFile)).filter(l -> !l.isEmpty()).forEach(l ->
		{
			Motif motif = new Motif(l.split("\t")[0]);
			motifs.add(motif);
			origLineMap.put(motif, l);
		});
	}

	private Map<String, Set<String>> readExplained() throws IOException
	{
		Set<String> genesAndSites = sequences.getGenesAndSites(null, selection);
		Map<String, Set<String>> map = new HashMap<>();
		Files.lines(Paths.get(cpResultFile)).map(l -> l.split("\t"))
			.filter(t -> t[1].startsWith("p") || t[1].startsWith("de"))
			.filter(t -> ((int) Math.signum(Double.valueOf(t[8]))) == selection).forEach(t ->
		{
			String trg = t[2];
			Set<String> sites = new HashSet<>(Arrays.asList(t[7].split("_")));
			sites.retainAll(Arrays.asList(t[3].split(";")));
			for (String site : sites)
			{
				String key = trg + "-" + site;
				if (genesAndSites.contains(key))
				{
					String eff = t[0];
					int sign = t[1].startsWith("p") ? 1 : -1;
					sign *= selection == 1 ? 1 : -1;
					if (sign < 0) eff = "-" + eff;
					if (!map.containsKey(key)) map.put(key, new HashSet<>());
					map.get(key).add(eff);
				}
			}
		});
		return map;
	}
}
