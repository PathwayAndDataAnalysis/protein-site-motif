package org.panda.proteinsitemotif;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

public class SimulatedDataGenerator
{
	public static void main(String[] args) throws IOException
	{
		generate();
	}

	public static void generate() throws IOException
	{
		String inFile = "../Finkle-PHYS-479/test_data/sample-phosphoproteomic-data.txt";
		String outFile = "../Finkle-PHYS-479/test_data/simulated-phosphoproteomic-data.txt";

		RankedSequences rs = new RankedSequences(11);
		rs.loadFromCPFiles(inFile, "Symbols", "Sites", "SignedP", "Feature");
		rs = rs.reduce(new AAReducer());

		Motif motif = new Motif();
		motif.add(new AAInPos(0, 'S', true));
		motif.add(new AAInPos(1, 'P', true));
		motif.add(new AAInPos(-4, 'S', false));
		motif.add(new AAInPos(5, 'K', true));

		BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
		writer.write("ID\tSymbols\tSites\tFeature\tEffect\tSignedP");

		Set<String> memory = new HashSet<>();

		Random rand = new Random();

		for (Sequence seq : rs.seqList)
		{
			String id = seq.gene + "-" + seq.site;
			if (!memory.contains(id)) memory.add(id);
			else continue;

			writer.write("\n" + id + "\t" + seq.gene + "\t" + seq.site + "\tP\t");

			double p = rand.nextDouble();

			if (seq.matches(motif))
			{
				p *= rand.nextDouble();
			}
			else if (rand.nextDouble() < 0.5)
			{
				p = -p;
			}
			writer.write("\t" + p);
		}

		writer.close();
	}
}
