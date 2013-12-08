package utils;

public class Utils {

	/**
	 * Return GC bias of A on positions 1 and 3 of each codon. Assumes the string is only codons and already ordered.
	 * @param A
	 * @return
	 */
	public static double GCPercent13(String A)
	{
		String B = "";
		int count = 1;
		for (char a : A.toCharArray())
		{
			if (count%3!=2)
				B += a;
			count++;
		}
		return GCPercent(B);
	}
	
	/**
	 * Returns a value between 0 and 1000. Max 1000. Min 0. Low value means less GC, high value means a lot of GC. 0 and 1000
	 * are special if no GC and no AT exist and should be ignored from mean value calculations.
	 * @param A
	 * @return
	 */
	public static double GCPercent(String A)
	{
		int ATcount = 0;
		int GCcount = 0;
		for (char a : A.toUpperCase().toCharArray())
		{
			if (a=='T' || a=='A')
			{
				ATcount++;
			}
			else if (a=='G' || a=='C')
			{
				GCcount++;
			}
		}
		if (GCcount==0)
			return 0;
		else if (ATcount==0)
			return 100;
		else
		{
			double bias = ((double) GCcount) / ((double) (GCcount+ATcount));
			if (bias>100)
				return 100;
			else
				return bias*100;
		}
	}

}
