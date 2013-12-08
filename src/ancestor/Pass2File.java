package ancestor;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;


public class Pass2File {
	
	ArrayList<Ancestor> ancestors = new ArrayList<Ancestor>();
	
	public Pass2File(String filename) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		
		String line = "";
		Ancestor currentAnc = null;
		int no = 0; // Ortholog group no
		while (reader.ready())
		{
			line = reader.readLine();
			if (line.startsWith("ANCESTOR"))
			{
				if (currentAnc!=null)
				{
					ancestors.add(currentAnc);
				}
				currentAnc = new Ancestor(line);
				System.out.println(currentAnc.name);
				no = 0;
			}
			else
			{
				no++;
				currentAnc.add(new OrthologGroup(line, no));
			}
		}
		if (currentAnc!=null)
		{
			ancestors.add(currentAnc);
		}
	}
	
	public Ancestor findAncestor(String name) {
		for (Ancestor a : ancestors)
		{
			if (a.name.equals(name))
				return a;
		}
		return null;
	}
}
