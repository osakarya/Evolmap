package ancestor;

public class HgtGene {

	public String name;
	public String origins;
	public int score;
	
	public HgtGene(String name, String origins, int score)
	{
		this.name = name;
		this.origins = origins;
		this.score = score;
	}
	
	public String toString()
	{
		return this.name + "\t" + this.origins + "\t" + this.score;
	}
}
