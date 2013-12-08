package ancestor;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

public class IndexDB 
{
	public HashMap<String, String> index1 = new HashMap<String, String>(); // Index, RN
	public HashMap<String, String> index2 = new HashMap<String, String>(); // RN, Index
	public HashMap<String, Integer> species = new HashMap<String, Integer>();
	int size;
	
	public IndexDB(String filename) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		
		String line = reader.readLine();
		size = Integer.parseInt(line);
		String index = "";
		String rn = "";
		while (reader.ready())
		{
			line = reader.readLine();
			if (!line.equals("") && line.startsWith(">"))
			{
				String[] items = line.substring(1).split("\t");
				index = items[0] + "_" + items[1];
				String specie = items[0];
				if (species.get(specie)==null)
				{
					species.put(specie, 0);
					System.out.println(specie);
				}
				line = reader.readLine();
				String[] items2 = line.substring(1).split(" ");
				rn = items2[0];
				line = reader.readLine(); // Move one ahead, sequence must be at least 1 line
				index1.put(index, rn);
				index2.put(rn,index);
			}
		}
	}
	
	public String speciesHeader()
	{
		String r ="";
		for (Entry<String, Integer> e : this.species.entrySet())
		{
			r += e.getKey() + "\t";
		}
		return r.substring(0, r.length()-1);
	}
	
	public String speciesCount(ArrayList<String> genes)
	{
		for (String g: genes)
		{
			String otherID = index2.get(g); 
			String specie = otherID.substring(0, otherID.indexOf("_"));
			Integer count = 0;
			if ((count = species.get(specie))!=null)
			{
				species.remove(specie);
				species.put(specie, count+1);
			}
			else
				System.out.println("This is not supposed to happen. Species: " + specie);
		}
		
		String r ="";
		
		for (Entry<String, Integer> e : this.species.entrySet())
		{
			r += e.getValue() + "\t";
		}
		
		for (Entry<String, Integer> e : this.species.entrySet())
		{
			e.setValue(0); // set all values back to zero
		}
		
		return r.substring(0, r.length()-1);
	}

}
