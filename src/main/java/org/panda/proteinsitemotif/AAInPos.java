package org.panda.proteinsitemotif;

public class AAInPos
{
	public int pos;
	public Character aa;
	public boolean prefer;

	public AAInPos(int pos, Character aa, boolean prefer)
	{
		this.pos = pos;
		this.aa = aa;
		this.prefer = prefer;
	}

	public AAInPos(int pos, Character aa)
	{
		this(pos, aa, true);
	}

	@Override
	public String toString()
	{
		return pos + (prefer ? "" : "!") + aa;
	}

	@Override
	public int hashCode()
	{
		return aa.hashCode() + Integer.hashCode(pos) - Boolean.hashCode(prefer);
	}

	@Override
	public boolean equals(Object obj)
	{
		if (obj instanceof AAInPos)
		{
			AAInPos aap = (AAInPos) obj;

			return aap.pos == pos && aap.aa.equals(aa) && aap.prefer == prefer;
		}
		return false;
	}
}
