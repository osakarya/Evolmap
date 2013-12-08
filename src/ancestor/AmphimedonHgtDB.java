package ancestor;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map.Entry;



public class AmphimedonHgtDB {
	
	public HashMap<String, HgtGene> genes = new HashMap<String, HgtGene>();
	
	public AmphimedonHgtDB(String filename) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		
		String line = "";
		while (reader.ready())
		{
			line = reader.readLine();
			String[] items = line.split("\t");
			String[] items2 = items[2].split("PACid");
			String geneName = items2[0].substring(0, items2[0].length()-1);
			HgtGene gene = new HgtGene(items[2],items[3], Integer.parseInt(items[1]));
			genes.put(geneName, gene);
		}
		
		reader.close();
	}
	
	public void printContents()
	{
		for (Entry<String, HgtGene> g : genes.entrySet())
		{
			System.out.println(g.getKey() + "\t" + g.getValue());
		}
	}
	
	
}
