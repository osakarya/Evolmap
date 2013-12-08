package ancestor;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;


import utils.Utils;
import fasta.FastaRecord;


public class AmphimedonOriginsDB {

	public HashMap<String, Gene> genes = new HashMap<String, Gene>();
	
	public AmphimedonOriginsDB(String filename) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		
		String line = "";
		while (reader.ready())
		{
			line = reader.readLine();
			String[] items = line.split("\t");
			String[] items2 = items[0].split("PACid");
			String geneName = items2[0].substring(0, items2[0].length()-1);
			Gene gene = new Gene(items[0],items[1], items[2], Integer.parseInt(items[3]));
			genes.put(geneName, gene);
		}
		
		reader.close();
	}
	
	public void printContents()
	{
		for (Entry<String, Gene> g : genes.entrySet())
		{
			System.out.println(g.getKey() + "\t" + g.getValue());
		}
	}
}
