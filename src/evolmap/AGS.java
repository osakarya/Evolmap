
package EvolMAP;
import java.util.*;
import java.io.*;
import javax.naming.Name;
/**
 * AGS = Ancestral gene set
 *
 * This object holds a vector of ortholog genes.
 * @author onur
 */
public class AGS 
{
    Vector<OrthologGenes> AGS;
    String name;
    int avg_symbetScore;
    int std_symbetScore;
    
    /** Creates a new instance of AGS */
    public AGS() 
    {
        this.name = "NO_NAME";
        this.AGS = new Vector<OrthologGenes>();
    }
    public AGS(Vector<OrthologGenes> AGS)
    {
        this.name = "NO_NAME";
        this.AGS = AGS;
    }
    public AGS(GeneDatabase db)
    {
        this.name = "NO_NAME";
        this.AGS = db.getOrthologGenesVector();
    }    
    public AGS(String name, GeneDatabase db)
    {
        this.name = name;
        this.AGS = db.getOrthologGenesVector();
    }
    public AGS(String name, Vector<OrthologGenes> AGS)
    {
        this.name = name;
        this.AGS = AGS;
    }
      
    public OrthologGenes findGene(String geneName)
    {
        if (geneName==null || geneName.equals("Unknown"))
            return null;
        for (int i=0; i<this.AGS.size(); i++)
        {
            OrthologGenes og = this.AGS.get(i);
            if (og.contains(geneName))
                return og;
        }
        return null;
    } 
    
    public void setAvgandSTDsymbetScore()
    {
        this.avg_symbetScore = 0;
        for (int i=0; i<this.AGS.size(); i++)
        {
            if (this.AGS.get(i).symbet)
                this.avg_symbetScore += this.AGS.get(i).symbet_score;
        }
        if (this.symbetCount()>0)
            this.avg_symbetScore /= this.symbetCount();
        
        this.std_symbetScore = 0;
        
        for (int i=0; i<this.AGS.size(); i++)
        {
            if (this.AGS.get(i).symbet)
            {
                this.std_symbetScore += Math.pow((this.AGS.get(i).symbet_score - this.avg_symbetScore), 2);
            }
        }
        if (this.symbetCount()>0)
        {
            this.std_symbetScore /= this.symbetCount();
            this.std_symbetScore = (int)Math.sqrt(this.std_symbetScore);
        }
    }
    public int symbetCount()
    {
        int count = 0;
        for (int i=0; i<this.AGS.size(); i++)
        {
            if (this.AGS.get(i).symbet)
                count++;
        }
        return count;
    }
    public int presentCount()
    {
        int count = 0;
        for (int i=0; i<this.AGS.size(); i++)
        {
            if (this.AGS.get(i).present)
                count++;
        }
        return count;
    }
    public void resetProcessed()
    {
        for (int i=0; i<this.AGS.size(); i++)
        {
            this.AGS.get(i).resetProcessed();
        }
    }
    
    public int size()
    {
        return this.AGS.size();
    }
    
    public void addOrthologGene(OrthologGenes orthologGenes)
    {
        this.AGS.add(orthologGenes);
    }
    
    public void addOrthologGenesSet(Vector<OrthologGenes> orthologGeneSet)
    {
        for (int i=0; i<orthologGeneSet.size(); i++)
        {
            this.AGS.add(orthologGeneSet.get(i));
        }
    }
    public Vector<OrthologGenes> getAGS()
    {
        return this.AGS;
    }
    public String toString()
    {
        String output_string = "";
        for (int i=0; i<this.AGS.size(); i++)
        {
            OrthologGenes nextOrthologGenes = this.AGS.get(i);
            output_string += nextOrthologGenes.toString();     
            
            // Write final empty space
            output_string += System.getProperty("line.separator");
        }
        return output_string;
    }
    
    // PASS 1 printout
    public String toStringNames()
    {
        String output_string = "";
        for (int i=0; i<this.AGS.size(); i++)
        {
            OrthologGenes nextOrthologGenes = this.AGS.get(i);
            if (nextOrthologGenes.possible_source!=null) // If it is a singular with no source it will be "Unknown" as source
                output_string += "SINGULAR\t" + nextOrthologGenes.possible_source + "\t" + nextOrthologGenes.singular_score + "\t";
            else
                output_string += "SYM-BET\t" + nextOrthologGenes.symbet_score + "\t";
            for (int j=0; j<nextOrthologGenes.genes.size(); j++)
            {
                String seq = nextOrthologGenes.genes.get(j);
                output_string += seq + "\t";
            }
            // Write final empty space
            output_string += System.getProperty("line.separator");
        }
        return output_string;        
    }
    
    private static BufferedReader open_file(File filename)
    {
            try {
                    return new BufferedReader(new FileReader(filename));
            }
            catch (IOException e) {return null;}
    }

    private static FileWriter create_file(File filename)
    {
            try {
                    return new FileWriter(filename);
            }
            catch (IOException e) {return null;}
    }

    private static boolean more_records(Reader r)
    {
            try {
                    return r.ready();
            }
            catch (IOException e) {return false;}
    }

    private static String read_a_line(BufferedReader r) 
    {
            try {
                    return r.readLine();
            }
            catch (IOException e) {return null;}
    }      
}
