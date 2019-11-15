package org.panda.proteinsitemotif;

import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Motif
{
	private List<AAInPos> members;

	public Motif(int pos, Character aa, boolean prefer)
	{
		this();
		members.add(new AAInPos(pos, aa, prefer));
	}

	public Motif(Motif base)
	{
		this();
		members.addAll(base.members);
	}

	public Motif(String shortStr)
	{
		this();
		for (String s : shortStr.split(" "))
		{
			boolean prefer = true;
			if (s.contains("!"))
			{
				prefer = false;
				s = s.replace("!", "");
			}
			Character aa = s.charAt(s.length() - 1);
			s = s.substring(0, s.length() - 1);
			int pos = Integer.valueOf(s);

			members.add(new AAInPos(pos, aa, prefer));
		}
	}

	public Motif()
	{
		this.members = new ArrayList<>();
	}

	public void add(AAInPos aap)
	{
		this.members.add(aap);
	}

	public void add(int pos, Character aa, boolean prefer)
	{
		add(new AAInPos(pos, aa, prefer));
	}

	public boolean contains(AAInPos aap)
	{
		return this.members.contains(aap);
	}

	public Stream<AAInPos> stream()
	{
		return this.members.stream();
	}

	public int size()
	{
		return members.size();
	}

	public boolean isGeneralOf(Motif m)
	{
		return m.members.containsAll(this.members) && m.size() > size();
	}

	@Override
	public String toString()
	{
		return members.toString();
	}

	public String toString(int width, Map<Character, List<Character>> convMap)
	{
		int halfW = width / 2;
		StringBuilder sb = new StringBuilder();

		for (int i = -halfW; i <= halfW; i++)
		{
			int ind = i;
			if (members.stream().anyMatch(ap -> ap.pos == ind))
			{
				members.stream().filter(ap -> ap.pos == ind).forEach(ap ->
				{
					if (!ap.prefer) sb.append("!");
					List<Character> subs = convMap == null ? Collections.singletonList(ap.aa) : convMap.get(ap.aa);
					if (subs.size() == 1) sb.append(ap.aa);
					else sb.append("[").append(CollectionUtil.merge(subs, "|")).append("]");
				});
			}
			else sb.append("-");

			sb.append(" ");
		}

		return sb.toString();
	}

	public String toShortString()
	{
		return CollectionUtil.merge(members.stream().sorted(Comparator.comparing(a -> a.pos)).collect(Collectors.toList()), " ");
	}

	@Override
	public int hashCode()
	{
		int h = 0;
		for (AAInPos aap : members)
		{
			h += aap.hashCode();
		}
		return h;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (obj instanceof Motif)
		{
			Motif m = (Motif) obj;
			return m.size() == size() && m.members.containsAll(members);
		}
		return false;
	}

	public static void write(Collection<Motif> motifs, String file) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(file));
		motifs.forEach(m -> FileUtil.writeln(m.toShortString(), writer));
		writer.close();
	}

	public static Set<Motif> loadSet(String file) throws IOException
	{
		return Files.lines(Paths.get(file)).map(l -> new Motif(l.split("\t")[0])).collect(Collectors.toSet());
	}

	public static List<Motif> loadList(String file) throws IOException
	{
		return Files.lines(Paths.get(file)).map(l -> new Motif(l.split("\t")[0])).collect(Collectors.toList());
	}
}
