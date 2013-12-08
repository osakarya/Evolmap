package ancestor;

import java.util.HashMap;
import java.util.Map.Entry;

public class GeneDB {

	HashMap<String, Gene> genes;
	
	public GeneDB()
	{
		this.genes = new HashMap<String, Gene>();
	}
	
	private void add(Gene gene)
	{
		Gene existingEntry = null;
		if ((existingEntry = genes.get(gene.name))!=null)
		{
			if (existingEntry.originAncestor.equals("DIVERGED")
					|| existingEntry.originAncestor.equals("SINGULAR"))
			{
				genes.put(gene.name, gene);
				//existingEntry = gene;
			}
		}
		else
		{
			genes.put(gene.name, gene);
		}
	}
	
	private void add(OrthologGroup og, String ancestorName, HashMap<String, String> indexDB)
	{
		if (og.type.equals("DIVERGED") || og.type.equals("SINGULAR"))
		{
			ancestorName = og.type;
		}
		for (String gene : og.genes)
		{
			Gene g = new Gene(gene, og, indexDB.get(gene).split("_")[0], ancestorName, og.score);
			this.add(g);
		}
	}
	
	public void add(Ancestor a, HashMap<String, String> indexDB)
	{
		for (Entry<String, OrthologGroup> og : a.groups.entrySet())
		{
			OrthologGroup g = og.getValue();
			this.add(g, a.name, indexDB);
		}
	}
	
	public void printStatus()
	{
		int presentCount = 0;
		int divergedCount = 0;
		int singularCount = 0;
		
		for (Entry<String, Gene> g : genes.entrySet())
		{
			String t = g.getValue().originAncestor;
			if (t.equals("SINGULAR"))
			{
				singularCount++;
			}
			else if (t.equals("DIVERGED"))
			{
				divergedCount++;
			}
			else {
				presentCount++;
			}
		}
		int total = presentCount + divergedCount + singularCount;
		System.out.println("Present: " + presentCount);
		System.out.println("Diverged: " + divergedCount);
		System.out.println("Singular: " + singularCount);
		System.out.println("Total: " + total);
	}
}
