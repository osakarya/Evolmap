package ancestor;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Does not record singulars
 * Records scores in OrthologGroup
 * @author osakarya
 *
 */
public class Pass1File {

	ArrayList<Ancestor> ancestors = new ArrayList<Ancestor>();
	
	public Pass1File(String filename) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		
		String line = "";
		Ancestor currentAnc = null;
		int no = 0;
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
				no=0;
			}
			else
			{
				no++;
				// Might think to change below, currently only reading SYM-BETs
				if (!line.startsWith("SINGULAR"))
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
