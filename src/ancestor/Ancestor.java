package ancestor;

import java.util.HashMap;
import java.util.Map.Entry;

public class Ancestor 
{
	String name;
	HashMap<String, OrthologGroup> groups;
	
	public Ancestor(String line)
	{
		String[] items = line.split("\t");
		this.name = items[1];
		groups = new HashMap<String, OrthologGroup>();
		changeNames();
	}
	
	public void add(OrthologGroup og)
	{
		this.groups.put(og.genes.get(0), og);
	}
	
	public int totalScore()
	{
		int sum = 0;
		for (Entry<String, OrthologGroup> og : groups.entrySet())
		{
			sum += og.getValue().score;
		}
		return sum;
	}
	
	public void changeNames()
	{
		if (this.name.equalsIgnoreCase("Homo_Mus_Ciona_Branchiostoma_Drosophila_Nematostella_Hydra_Trichoplax_Amphimedon"))
		{
			name = "Metazoa";
		}
		else if (this.name.equalsIgnoreCase("Homo_Mus_Ciona_Branchiostoma_Drosophila_Nematostella_Hydra_Trichoplax_Amphimedon_Monosiga"))
		{
			name = "Choanozoa";
		}
		else if (this.name.equalsIgnoreCase("Homo_Mus_Ciona_Branchiostoma_Drosophila_Nematostella_Hydra_Trichoplax_Amphimedon_Monosiga_Saccharomyces_Neurospora_Schizosaccharomyces"))
		{
			name = "Opisthokont";
		}
		else if (this.name.equalsIgnoreCase("Homo_Mus_Ciona_Branchiostoma_Drosophila_Nematostella_Hydra_Trichoplax_Amphimedon_Monosiga_Saccharomyces_Neurospora_Schizosaccharomyces_Dictyostelium"))
		{
			name = "Eukaryote-1";
		}
		else if (this.name.equalsIgnoreCase("Homo_Mus_Ciona_Branchiostoma_Drosophila_Nematostella_Hydra_Trichoplax_Amphimedon_Monosiga_Saccharomyces_Neurospora_Schizosaccharomyces_Dictyostelium_Arabidopsis_Volvox_Chlamydomonas"))
		{
			name = "Eukaryote-2";
		}
		else if (this.name.equalsIgnoreCase("Homo_Mus_Ciona_Branchiostoma_Drosophila_Nematostella_Hydra_Trichoplax_Amphimedon_Monosiga_Saccharomyces_Neurospora_Schizosaccharomyces_Dictyostelium_Arabidopsis_Volvox_Chlamydomonas_Thalassiosira_Phaeodactylum"))
		{
			name = "Eukaryote-3";
		}
		else if (this.name.equalsIgnoreCase("Homo_Mus_Ciona_Branchiostoma_Drosophila_Nematostella_Hydra_Trichoplax_Amphimedon_Monosiga_Saccharomyces_Neurospora_Schizosaccharomyces_Dictyostelium_Arabidopsis_Volvox_Chlamydomonas_Thalassiosira_Phaeodactylum_Archae"))
		{
			name = "Archae";
		}
		else if (this.name.equalsIgnoreCase("Homo_Mus_Ciona_Branchiostoma_Drosophila_Nematostella_Hydra_Trichoplax_Amphimedon_Monosiga_Saccharomyces_Neurospora_Schizosaccharomyces_Dictyostelium_Arabidopsis_Volvox_Chlamydomonas_Thalassiosira_Phaeodactylum_Archae_Bacteria"))
		{
			name = "Bacteria";
		}
	}
}
