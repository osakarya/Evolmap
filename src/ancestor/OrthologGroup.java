package ancestor;

import java.util.ArrayList;
import java.util.HashMap;

public class OrthologGroup {

	public ArrayList<String> genes;
	public String type;
	public int score;
	public int no;
	
	public OrthologGroup (String line, int no)
	{
		genes = new ArrayList<String>();
		String[] items = line.split("\t");
		this.type = items[0];
		this.no = no;
		if (this.type.equals("SYM-BET"))
		{
			this.score = Integer.parseInt(items[1]);
			
			for (int i=2; i<items.length; i++)
			{
				genes.add(items[i]);
			}
		}
		else
		{
			if (this.type.equals("DIVERGED"))
			{
				for (int i=2; i<items.length; i++)
				{
					genes.add(items[i]);
				}
			}
			else
			{
				for (int i=1; i<items.length; i++)
				{
					genes.add(items[i]);
				}
			}
		}
	}
	
	public String toString()
	{
		String r = score + "\t";
		for (String g : genes)
		{
			r+=g+"\t";
		}
		return (r.substring(0, r.length()-1));
	}
}
