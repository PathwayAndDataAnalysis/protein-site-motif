package org.panda.proteinsitemotif;

import java.util.*;
import java.util.stream.Collectors;

public class Constants
{
	public static final List<Character> AA_LIST = Arrays.asList(
		'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V');

	public static final Set<Character> AA_SET = new HashSet<>(AA_LIST);

	public static final Map<Character, Character> REDUC_MAP = new HashMap<>();
	public static final Map<Character, List<Character>> REV_REDUC_MAP = new HashMap<>();

	static
	{
		REDUC_MAP.put('A', 'A');
		REDUC_MAP.put('R', 'K');
		REDUC_MAP.put('N', 'Q');
		REDUC_MAP.put('D', 'D');
		REDUC_MAP.put('C', 'C');
		REDUC_MAP.put('E', 'D');
		REDUC_MAP.put('Q', 'Q');
		REDUC_MAP.put('G', 'A');
		REDUC_MAP.put('H', 'H');
		REDUC_MAP.put('I', 'I');
		REDUC_MAP.put('L', 'I');
		REDUC_MAP.put('K', 'K');
		REDUC_MAP.put('M', 'I');
		REDUC_MAP.put('F', 'F');
		REDUC_MAP.put('P', 'P');
		REDUC_MAP.put('S', 'S');
		REDUC_MAP.put('T', 'S');
		REDUC_MAP.put('W', 'W');
		REDUC_MAP.put('Y', 'F');
		REDUC_MAP.put('V', 'I');

		REV_REDUC_MAP.put('A', Arrays.asList('A', 'G'));
		REV_REDUC_MAP.put('K', Arrays.asList('R', 'K'));
		REV_REDUC_MAP.put('Q', Arrays.asList('N', 'Q'));
		REV_REDUC_MAP.put('D', Arrays.asList('D', 'E'));
		REV_REDUC_MAP.put('C', Collections.singletonList('C'));
		REV_REDUC_MAP.put('H', Collections.singletonList('H'));
		REV_REDUC_MAP.put('I', Arrays.asList('I', 'L', 'M', 'V'));
		REV_REDUC_MAP.put('F', Arrays.asList('F', 'Y'));
		REV_REDUC_MAP.put('P', Collections.singletonList('P'));
		REV_REDUC_MAP.put('S', Arrays.asList('S', 'T'));
		REV_REDUC_MAP.put('W', Collections.singletonList('W'));
	}

	public static final List<Character> AA_LIST_REDUC = REDUC_MAP.values().stream().distinct().sorted()
		.collect(Collectors.toList());


}
