package org.panda.proteinsitemotif;

public class Sequence
{
	public String seq;
	public String gene;
	public String site;

	public Sequence(String seq, String gene, String site)
	{
		this.seq = seq;
		this.gene = gene;
		this.site = site;
	}

	public Sequence(String seq)
	{
		this.seq = seq;
	}

	public boolean matches(Motif motif)
	{
		int halfW = seq.length() / 2;

		return motif.stream().allMatch(aap -> aap.prefer ?
				aap.aa == seq.charAt(aap.pos+halfW) :
				aap.aa != seq.charAt(aap.pos+halfW));
	}
}
