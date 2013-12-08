
package EvolMAP;
import java.util.*;

/**
 *
 * @author onur
 */
public class OrthologGenes 
{
    Vector<String> genes; // Stores the indexed name of gene
    Vector<String> divergedInparalogs; // Stores the indexed name of diverged inparalogs of this gene
    
    String possible_source;
    boolean symbet;
    int symbet_score;
    int singular_score;
    boolean singularFirstLineage;
    boolean singularSecondLineage;
    boolean present; // Every gene is originally present. An outside method must change this if gene is singular and is not found in ancestral lineages
    boolean processed = false;
    boolean divergedParalog;
    int paralogSource; // Paralog source for ancestor i
    
    /** Creates a new instance of OrthologGenes */
    public OrthologGenes() 
    {
        this.genes = new Vector<String>();
        this.symbet = true; // Set each gene as symbets initially, then check to change them with outside methods
        this.symbet_score = 0; // This is the ortholog symbet score (best match-- set by AncestorDatabase.java)
        this.present = true;  // Set each gene as present initially, then check to change them with outside methods, if it is singular it is initially set as not present
        this.singularFirstLineage = false; // Check this if need be with outside methods
        this.singularSecondLineage = false; // Check this if need be with outside methods
        this.divergedParalog = false;
        this.paralogSource = -1;
    }
    public OrthologGenes(Vector<String> sequences)
    {
        this.genes = new Vector<String>();
        this.symbet = true; // Set each gene as symbets initially, then check to change them with outside methods
        this.present = true;  // Set each gene as present initially, then check to change them with outside methods, if it is singular it is initially set as not present
        this.singularFirstLineage = false; // Check this if need be with outside methods
        this.singularSecondLineage = false; // Check this if need be with outside methods
        this.divergedParalog = false;
        this.paralogSource = -1;
        
        for (int i=0; i<sequences.size(); i++)
        {
            this.addOrtholog(sequences.get(i));
        }
    }
    // Create and add one gene. Used when reading GeneDatabase files
    public OrthologGenes(String gene)
    {
        this.genes = new Vector<String>();
        
        this.addOrtholog(gene);
        this.symbet = true;
        this.present = true;
        this.singularFirstLineage = false;
        this.singularSecondLineage = false;
        this.divergedParalog = false;
        this.paralogSource = -1;
    }
    public void setsymbet()
    {
        this.symbet = true;
    }
    public void setsymbetScore(int score)
    {
        this.symbet_score = score;
    }
    public void setSingularFirstLineage()
    {
        this.symbet = false;
        this.singularFirstLineage = true;
        this.present = false;
    }
    public void setSingularSecondLineage()
    {
        this.symbet = false;
        this.singularSecondLineage = true;
        this.present = false;
    }
    // Used for gain-loss count methods in AGS_analyzer.java
    public void setPresent()
    {
        this.present = true;
    }
    // Used for gain-loss count methods in AGS_analyzer.java
    public void setProcessed()
    {
        this.processed = true;
    }
    public void setDivergedParalog()
    {
        this.divergedParalog = true;
    }
    public void setParalogSource(int ancestor)
    {
        this.paralogSource = ancestor;
    }
    public void resetProcessed()
    {
        this.processed = false;
    }
    public void setPossibleSource(String source, int score)
    {
        this.possible_source = source;
        this.singular_score = score;
    }

    public void addDivergedInParalogs(OrthologGenes orthologGenes)
    {
        if (this.divergedInparalogs==null)
            divergedInparalogs = new Vector<String>();
        for (int i=0; i<orthologGenes.getGenes().size(); i++)
            divergedInparalogs.add(orthologGenes.getGenes().get(i));
    }
    
    public void addOrtholog(String gene)
    {
        this.genes.add(gene);
    }
    
    public void addOrthologGenes(OrthologGenes orthologGenes)
    {
        for (int i=0; i<orthologGenes.numberOfSequences(); i++)
        {
            String nextSequence = orthologGenes.getGenes().get(i);
            this.addOrtholog(nextSequence);
        }
    }
    public boolean contains(String checkSequence)
    {
        for (int i=0; i<this.genes.size(); i++)
        {
            String otherSequence = this.genes.get(i);
            if (otherSequence.equals(checkSequence))
                return true;
        }
        return false;
    }    
    
    public int numberOfSequences()
    {
        return genes.size();
    }
    
    public Vector<String> getGenes()
    {
        return this.genes;
    } 
}
