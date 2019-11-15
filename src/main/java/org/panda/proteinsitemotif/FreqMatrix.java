package org.panda.proteinsitemotif;

import java.util.HashMap;
import java.util.Map;


public class FreqMatrix
{
	int seqWidth;
	int halfW;
	Map<Integer, Map<Character, Integer>> matrix;
	int seqCnt;
	Character centralAA;

	public FreqMatrix(int seqWidth, Character centralAA)
	{
		this.seqWidth = seqWidth;

		if (seqWidth % 2 != 1)
			throw new IllegalArgumentException("Sequence length has to be an odd number. Given: " + seqWidth);

		this.halfW = seqWidth / 2;
		matrix = new HashMap<>();

		for (int i = 1; i <= halfW; i++)
		{
			matrix.put(i, new HashMap<>());
			matrix.put(-i, new HashMap<>());
		}

		this.seqCnt = 0;
		this.centralAA = centralAA;
	}

	public void addSequences(Sequences sq, Motif motif, Integer status)
	{
		sq.getSeqStream(motif, status).forEach(this::addSequence);
	}

	public void addSequence(String seq)
	{
		if (seq.length() != seqWidth) throw new IllegalArgumentException("Each sequence must have the same width: " +
			seqWidth + ". Found = " + seq.length());

		if (seq.charAt(halfW) != centralAA) throw new IllegalArgumentException(
			"Central amino acid does not match: " + seq);

		for (int i = 0; i < seq.length(); i++)
		{
			Character aa = seq.charAt(i);
			int index = -halfW + i;

			if (index != 0)
			{
				Map<Character, Integer> cnts = matrix.get(index);
				cnts.put(aa, cnts.getOrDefault(aa, 0) + 1);
			}
		}

		seqCnt++;
	}

	@Override
	public String toString()
	{
		StringBuilder sb = new StringBuilder("Sequences = " + seqCnt);

		matrix.keySet().stream().sorted().forEach(p ->
		{
			sb.append("\n").append(p).append(" --");
			matrix.get(p).keySet().stream().sorted().forEach(a -> sb.append(" ").append(a).append(":").append(matrix.get(p).get(a)));
		});

		return sb.toString();
	}
}
