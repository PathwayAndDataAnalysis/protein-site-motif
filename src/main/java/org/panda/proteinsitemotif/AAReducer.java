package org.panda.proteinsitemotif;

public class AAReducer
{
	public String reduce(String seq)
	{
		StringBuilder sb = new StringBuilder();

		for (int i = 0; i < seq.length(); i++)
		{
			sb.append(Constants.REDUC_MAP.get(seq.charAt(i)));
		}

		return sb.toString();
	}
}
