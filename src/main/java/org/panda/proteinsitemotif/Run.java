package org.panda.proteinsitemotif;

import org.panda.utility.FileUtil;
import org.panda.utility.Tuple;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class Run
{
	public static void main(String[] args) throws IOException
	{
//		load();
//		findMotifs();
//		generateBinaryGraph();
		runExhaustive();
//		findEnrichmentsOfMotifs();
	}

	public static void load() throws IOException
	{
		int seqWidth = 17;
		int selectStatus = 1;
		double fdrThr = 0.1;
		Character aa = 'S';
		String file = "/Users/ozgun/Documents/Analyses/Aslan-platelet/with-feedback/data-fdr0.1.txt";

		Sequences sq = new Sequences(seqWidth);
		sq.loadFromAslanCPFiles(file);
//		sq.printAAFreqs();
//		System.out.println();

		sq = sq.reduce(new AAReducer());

		Motif motif = new Motif(0, 'S', true);

		FreqMatrix bgM = new FreqMatrix(seqWidth, aa);
		bgM.addSequences(sq, motif, null);
//		System.out.println(bgM);
//		System.out.println();

		FreqMatrix selM = new FreqMatrix(seqWidth, aa);
		selM.addSequences(sq, motif, selectStatus);
//		System.out.println(selM);
//		System.out.println();

		DeviationDetector dd = new DeviationDetector(bgM, selM, fdrThr, Deviation.POSITIVE);
		Map<AAInPos, Tuple> signif = dd.getSignificantDeviations();
		System.out.println("signif.size() = " + signif.size());

		signif.keySet().stream().sorted(Comparator.comparingDouble(o -> signif.get(o).p)).forEach(ap->
		{
			System.out.println(ap + "\t" + signif.get(ap));
		});
	}

	public static void findMotifs() throws IOException
	{
		int seqWidth = 13;
		int selectStatus = 1;
		double fdrThr = 0.1;
		Character aa = 'S';
		String file = "/Users/ozgun/Documents/Analyses/Aslan-platelet/with-feedback/data-fdr0.1.txt";

		Sequences sq = new Sequences(seqWidth);
		sq.loadFromAslanCPFiles(file);

		AAReducer reducer = new AAReducer();
		sq = sq.reduce(reducer);

		new GreedyIterativeMotifGenerator(sq, selectStatus, aa, fdrThr).run();
	}

	public static void generateBinaryGraph() throws IOException
	{
		int seqWidth = 13;
		int selectStatus = 1;
		double fdrThr = 0.1;
		Character aa = 'S';
		String file = "/Users/ozgun/Documents/Analyses/Aslan-platelet/with-feedback/data-fdr0.1.txt";

		Sequences sq = new Sequences(seqWidth);
		sq.loadFromAslanCPFiles(file);

		AAReducer reducer = new AAReducer();
		sq = sq.reduce(reducer);

		new DeviationGraphGenerator(sq, selectStatus, aa, fdrThr).run();
	}

	public static void runExhaustive() throws IOException
	{
		int seqWidth = 13;
		int selectStatus = 1;
		double fdrThr = 0.1;
		Character aa = 'S';
		String caseName = "with-feedback";
//		String inDir = "/Users/ozgun/Documents/Analyses/Aslan-platelet/" + caseName + "/";
		String inDir = "/home/ozgun/Analyses/Aslan-platelet/" + caseName + "/";
		String file = inDir + "data-fdr0.1.txt";
//		String outDir = "/Users/ozgun/Documents/Analyses/motif/";
		String outDir = "/home/ozgun/Analyses/ProteinMotif/exhaustive/";
		int iterations = 10000;
		int minTotalTrg = 3;

		caseName += "-" + aa + "-w" + seqWidth + (selectStatus == 1 ? "-upreg" : "-dwreg");

		Sequences sq = new Sequences(seqWidth);
		sq.loadFromAslanCPFiles(file);

		sq = sq.reduce(new AAReducer());

		ExhaustiveMotifGenerator mg = new ExhaustiveMotifGenerator(sq, selectStatus, aa, fdrThr);
//		GreedyIterativeMotifGenerator mg = new GreedyIterativeMotifGenerator(sq, selectStatus, aa, fdrThr);
		Set<Motif> motifs = mg.run();

		Motif.write(motifs, outDir + caseName + "-motif.txt");
		MotifDAG dag = new MotifDAG(motifs);
		dag.writeGraph(outDir + caseName + "-motif-DAG.sif");

		String motEnrichFile = outDir + caseName + "-motifs-with-enrichment.txt";
		BufferedWriter writer1 = Files.newBufferedWriter(Paths.get(motEnrichFile));
		EnrichmentTester et = new EnrichmentTester();
		Map<String, String> tooltipMap = new HashMap<>();

		for (Motif motif : motifs)
		{
			Map<String, Double> signif = et.findEnrichedTFs(sq, motif, new Motif(0, aa, true), fdrThr, iterations, minTotalTrg);
			FileUtil.lnwrite(motif.toShortString() + "\t", writer1);
			StringBuilder sb = new StringBuilder();
			signif.keySet().stream().sorted(Comparator.comparing(tf -> Math.abs(signif.get(tf)))).forEach(tf ->
			{
				int sign = (signif.get(tf) < 0 ? -1 : 1) * selectStatus;
				sb.append(sign < 0 ? "-" : "").append(tf).append(", ").append(signif.get(tf)).append("\\n");
				int tfSign = (int) Math.signum(signif.get(tf));
				FileUtil.write(" " + (tfSign * selectStatus < 0 ? "-" : "") + tf, writer1);
			});
			if (sb.length() > 0)
			{
				String tooltip = sb.toString();
				tooltip = tooltip.substring(0, tooltip.length() - 2);
				tooltipMap.put(motif.toShortString(), tooltip);
			}
		}

		writer1.close();

		BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(outDir + caseName + "-motif-DAG.format"));

		for (String motifID : tooltipMap.keySet())
		{
			FileUtil.writeln("node\t" + motifID + "\tcolor\t255 255 200", writer2);
			FileUtil.writeln("node\t" + motifID + "\ttooltip\t" + tooltipMap.get(motifID), writer2);
		}

		writer2.close();

		new SequenceMotifMatcher(sq, selectStatus, inDir + "results.txt", motEnrichFile).run(outDir + caseName + "-mot-seq-assoc.txt");
	}

	public static void findEnrichmentsOfMotifs() throws IOException
	{
		int seqWidth = 13;
		int selectStatus = 1;
		double fdrThr = 0.1;
		Character aa = 'S';
		String file = "/Users/ozgun/Documents/Analyses/Aslan-platelet/with-feedback/data-fdr0.1.txt";
		String outDir = "/Users/ozgun/Documents/Analyses/motif/";

		Set<Motif> motifs = Motif.loadSet(outDir + "exhaustive.txt");

		Sequences sq = new Sequences(seqWidth);
		sq.loadFromAslanCPFiles(file);

		AAReducer reducer = new AAReducer();
		sq = sq.reduce(reducer);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outDir + "with-enrichments.txt"));
		EnrichmentTester et = new EnrichmentTester();

		for (Motif motif : motifs)
		{
			Map<String, Double> signif = et.findEnrichedTFs(sq, motif, new Motif(0, 'S', true), fdrThr, 10000, 5);
			System.out.print("\n" + motif.toShortString() + "\t");
			FileUtil.lnwrite(motif.toShortString() + "\t", writer);
			signif.keySet().stream().sorted(Comparator.comparing(tf -> Math.abs(signif.get(tf)))).forEach(tf -> FileUtil.write(" " + (signif.get(tf) < 0 ? "-" : "") + tf, writer));
			signif.keySet().stream().sorted(Comparator.comparing(tf -> Math.abs(signif.get(tf)))).forEach(tf -> System.out.print(" " + (signif.get(tf) < 0 ? "-" : "") + tf));
		}

		System.out.println();

		for (Motif base : motifs)
		{
			for (Motif motif : motifs)
			{
				if (base.isGeneralOf(motif) && motif.size() == base.size() + 1)
				{
					Map<String, Double> signif = et.findEnrichedTFs(sq, motif, base, fdrThr, 10000, 5);
					System.out.print("\n" + base.toShortString() + " --> " + motif.toShortString() + "\t");
					FileUtil.lnwrite(base.toShortString() + " --> " + motif.toShortString() + "\t", writer);
					signif.keySet().stream().sorted(Comparator.comparing(tf -> Math.abs(signif.get(tf)))).forEach(tf -> FileUtil.write(" " + (signif.get(tf) < 0 ? "-" : "") + tf, writer));
					signif.keySet().stream().sorted(Comparator.comparing(tf -> Math.abs(signif.get(tf)))).forEach(tf -> System.out.print(" " + (signif.get(tf) < 0 ? "-" : "") + tf));
				}
			}
		}

		writer.close();
	}
}
