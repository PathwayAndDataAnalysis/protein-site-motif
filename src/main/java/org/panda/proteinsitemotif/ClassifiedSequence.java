package org.panda.proteinsitemotif;

public class ClassifiedSequence extends Sequence
{
	public int status;

	public ClassifiedSequence(String seq, int status, String gene, String site)
	{
		super(seq, gene, site);
		this.status = status;
	}

	public ClassifiedSequence(String seq, int status)
	{
		super(seq);
		this.status = status;
	}
}
