package ancestor;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map.Entry;

public class CalculateOrigins {

	public static void main(String args[]) throws IOException
	{
		String filename1 = "Run_with_21_species_results/all_genomes.ancestors_pass2.rn";
		String filename2 = "Run_with_21_species_results/all_genomes.ancestors_pass1";
		String filename3 = "Run_with_21_species_results/all_genomes.gd.fa";
		String outname = "Run_with_21_species_results/gene.origins.txt";
		String outname2 = "Run_with_21_species_results/amphimedon.origins.txt";
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(outname));
		BufferedWriter writer2 = new BufferedWriter(new FileWriter(outname2));
		Pass2File ad1 = new Pass2File(filename1);
		Pass1File ad2 = new Pass1File(filename2);
		IndexDB id = new IndexDB(filename3);
		
		//System.out.println(ad1.ancestors.get(1).groups.get("ENSP00000313272").score);
		//System.out.println(ad2.ancestors.get(ad2.ancestors.size()-2).groups.get(id.index2.get("ENSP00000313272")).score);
		//System.out.println(ad1.ancestors.get(1).totalScore());
		//System.out.println(ad2.ancestors.get(ad2.ancestors.size()-2).totalScore());
		for (Ancestor anc1 : ad1.ancestors)
		{
			Ancestor anc2 = ad2.findAncestor(anc1.name);
			if (anc2==null)
			{
				System.out.println("This is not supposed to happen " + anc1.name);
			}
			for (Entry<String, OrthologGroup> og : anc1.groups.entrySet())
			{
				String key = og.getKey();
				String index_key = id.index2.get(key);
				OrthologGroup og1 = og.getValue();
				//System.out.println(og1);
				if (index_key!=null)
				{
					OrthologGroup og2 = anc2.groups.get(index_key);
					if (og2!=null)
					{
						og1.score = og2.score;
					}
					else
					{
					}
				}
				else
				{
				}
			}
		}
		
		GeneDB gd = new GeneDB();
		for (Ancestor anc1 : ad1.ancestors)
		{
			gd.add(anc1, id.index2);
			//System.out.println(anc1.name);
			//gd.printStatus();
			//System.out.println();
		}
		//System.out.println(ad1.ancestors.get(1).groups.get("ENSP00000313272").score);
		//System.out.println(ad1.ancestors.get(1).totalScore());
		
		
		boolean headerWrote = false;
		int count = 0;
		for (Entry<String, Gene> g : gd.genes.entrySet())
		{
			Gene gene = g.getValue();
			if (!headerWrote)
			{
				writer2.write(gene.printHeader(id.speciesHeader()) + "\n");
				headerWrote=true;
			}
			if (gene.name.contains("Aqu1"))
			{
				writer2.write(gene.toString(id.speciesCount(gene.og.genes)) + "\n");
			}
			if (count++ % 1000 == 0)
			{
				System.out.println(count);
			}
		}
		
		headerWrote = false;
		count = 0;
		for (Entry<String, Gene> g : gd.genes.entrySet())
		{
			Gene gene = g.getValue();
			if (!headerWrote)
			{
				writer.write(gene.printHeader(id.speciesHeader()) + "\n");
				headerWrote=true;
			}
			if (count++ % 1000 == 0)
			{
				System.out.println(count);
			}
			writer.write(gene.toString(id.speciesCount(gene.og.genes)) + "\n");
		}
		
		writer.close();
		writer2.close();
	}
}
