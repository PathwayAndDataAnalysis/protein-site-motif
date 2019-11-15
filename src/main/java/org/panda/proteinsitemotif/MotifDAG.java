package org.panda.proteinsitemotif;

import org.panda.utility.FileUtil;
import org.panda.utility.graph.DirectedGraph;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class MotifDAG
{
	Map<String, Motif> nameToMotif;
	Set<Motif> motifs;
	DirectedGraph graph;

	public MotifDAG(Set<Motif> motifs)
	{
		this.motifs = motifs;
		this.nameToMotif = motifs.stream().collect(Collectors.toMap(Motif::toShortString, m -> m));
		buildGraph();
	}

	private void buildGraph()
	{
		this.graph = new DirectedGraph("Motif DAG", "controls-state-change-of");

		for (Motif m1 : motifs)
		{
			for (Motif m2 : motifs)
			{
				if (m1 != m2)
				{
					if (m1.isGeneralOf(m2) && m2.size() == m1.size() + 1)
					{
						graph.putRelation(m1.toShortString(), m2.toShortString());
					}
				}
			}
		}
	}

	public void writeGraph(String file) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(file));
		this.graph.write(writer);
		for (Motif motif : motifs)
		{
			if (!graph.hasNode(motif.toShortString()))
			{
				FileUtil.writeln(motif.toShortString(), writer);
			}
		}
		writer.close();
	}

}
