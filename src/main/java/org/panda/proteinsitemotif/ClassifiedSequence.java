package org.panda.proteinsitemotif;

public class ClassifiedSequence
{
	public String seq;
	public int status;
	public String gene;
	public String site;

	public ClassifiedSequence(String seq, int status, String gene, String site)
	{
		this.seq = seq;
		this.status = status;
		this.gene = gene;
		this.site = site;
	}

	public ClassifiedSequence(String seq, int status)
	{
		this.seq = seq;
		this.status = status;
	}

	public boolean matches(Motif motif)
	{
		int halfW = seq.length() / 2;

		return motif.stream().allMatch(aap -> aap.prefer ?
				aap.aa == seq.charAt(aap.pos+halfW) :
				aap.aa != seq.charAt(aap.pos+halfW));
	}
}
