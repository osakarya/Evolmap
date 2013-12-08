package ancestor;

public class Gene {

	public String name;
	public String species;
	public String originAncestor;
	public OrthologGroup og;
	public int score;
	
	
	public Gene (String name)
	{
		this.name = name;
	}
	
	public Gene (String name, String type)
	{
		this.name = name;
		this.originAncestor = originAncestor;
	}
	
	public Gene (String name, String originAncestor, int score)
	{
		this.name = name;
		this.originAncestor = originAncestor;
		this.score = score;
	}
	
	public Gene (String name, String species, String originAncestor, int score)
	{
		this.name = name;
		this.originAncestor = originAncestor;
		this.score = score;
		this.species = species;
	}
	
	public Gene (String name, OrthologGroup og, String species, String originAncestor, int score)
	{
		this.name = name;
		this.originAncestor = originAncestor;
		this.score = score;
		this.species = species;
		this.og = og;
	}
	
	public boolean notAssigned()
	{
		if (this.originAncestor.equals("DIVERGED") || this.originAncestor.equals("SINGULAR"))
			return true;
		else
			return false;
	}
	
	public String printHeader(String speciesHeader)
	{
		return "GeneName\tSpecies\tOrigin\tOrtholog-no\tScore" + "\t" + speciesHeader;		
	}
	
	public String toString(String speciesCount)
	{
		return this.name + "\t" + this.species + "\t" + this.originAncestor 
				+ "\t" + this.og.no +"\t" + this.score + "\t" + speciesCount;
	}
}
