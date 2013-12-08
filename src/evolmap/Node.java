
package EvolMAP;

/**
 *
 * @author onur
 */
public class Node 
{
    String name;
    Node ancestor;
    Node descendant1;
    Node descendant2;
    boolean species_node;
    AGS ags = null;
    
    int total_genes; // Total present genes at this ancestor
    int total_gains; // Total genes gained in the branch leading towards this ancestor
    int total_losses; // Total genes lost in the branch leading towards this ancestor
          
    /** Creates a new instance of Node */
    public Node(String name, Node descendant1, Node descendant2) 
    {
        this.name = this.convertName(name);
        this.ancestor = null;
        this.descendant1 = descendant1;
        this.descendant2 = descendant2;
        
        this.total_genes = 0;
        
        if (descendant1==null || descendant2 == null)
            species_node = true;
        else
            species_node = false;
    }
    
    public AGS getAGS()
    {
        return this.ags;
    }
    public void setAGS(AGS ags)
    {
        this.ags = ags;
    }
    
    public void setAncestor(Node ancestor)
    {
        this.ancestor = ancestor;
    }
    public void setTotalGeneCount(int geneCount)
    {
        this.total_genes = geneCount;
    }
    public String convertName(String name)
    {
        name = name.replaceAll(",", "_");
        String noParenthesisName = "";
        for (int i=0; i<name.length(); i++)
        {
            if (name.charAt(i)!=')' && name.charAt(i)!='(')
            {
                noParenthesisName += name.charAt(i);
            }
        }
        return noParenthesisName;
    }
    public String getTreePrintWithsymbetDistances()
    {
        if (this.species_node)
            return this.name;
        else
        {
            return ("(" + this.descendant1.getTreePrintWithsymbetDistances() + ":" + this.getDescendant1BranchLength() + "," +
                          this.descendant2.getTreePrintWithsymbetDistances() + ":" + this.getDescendant2BranchLength() + ")");
        }
    }   
    public String getTreePrintWithGeneExpansionDistances()
    {
        if (this.species_node)
            return this.name;
        else
        {
            return ("(" + this.descendant1.getTreePrintWithGeneExpansionDistances() + ":" + this.getDescendant1GeneExpansions() + "," +
                          this.descendant2.getTreePrintWithGeneExpansionDistances() + ":" + this.getDescendant2GeneExpansions() + ")");
        }        
    }
    public String getTreePrintWithGeneGainDistances()
    {
        if (this.species_node)
            return this.name;
        else
        {
            return ("(" + this.descendant1.getTreePrintWithGeneGainDistances() + ":" + this.getDescendant1GeneGains() + "," +
                          this.descendant2.getTreePrintWithGeneGainDistances() + ":" + this.getDescendant2GeneGains() + ")");
        }        
    }
    public String getTreePrintWithGeneLossDistances()
    {
        if (this.species_node)
            return this.name;
        else
        {
            return ("(" + this.descendant1.getTreePrintWithGeneLossDistances() + ":" + this.getDescendant1GeneLosses() + "," +
                          this.descendant2.getTreePrintWithGeneLossDistances() + ":" + this.getDescendant2GeneLosses() + ")");
        }        
    }      
    public double getDescendant1BranchLength()
    {
        double total_distance = (double)((double)(1000 - this.ags.avg_symbetScore)/2);
        double descendant_distance = 0;
        if (!this.descendant1.species_node)
            descendant_distance = (double)((double)(1000 - this.descendant1.ags.avg_symbetScore)/2);
        double returnDistance = total_distance - descendant_distance;
        if (returnDistance>0)
            return returnDistance;
        else
            return 1;
    }
    
    public double getDescendant2BranchLength()
    {
        double total_distance = (double)((double)(1000 - this.ags.avg_symbetScore)/2);
        double descendant_distance = 0;
        if (!this.descendant2.species_node)
            descendant_distance = (double)((double)(1000 - this.descendant2.ags.avg_symbetScore)/2);
        double returnDistance = total_distance - descendant_distance;
        if (returnDistance>0)
            return returnDistance;
        else
            return 1;
    }
    public int getDescendant1GeneExpansions()
    {
        int total_expansion = this.descendant1.total_genes - this.total_genes;
        if (total_expansion < 0)
            total_expansion = 0;
        return total_expansion;
    }
    public int getDescendant2GeneExpansions()
    {
        int total_expansion = this.descendant2.total_genes - this.total_genes;
        if (total_expansion < 0)
            total_expansion = 0;
        return total_expansion;        
    }
    public int getDescendant1GeneGains()
    {
        int total_gains = this.descendant1.total_gains;
        if (total_gains < 0)
            total_gains = 0;
        return total_gains;
    }
    public int getDescendant2GeneGains()
    {
        int total_gains = this.descendant2.total_gains;
        if (total_gains < 0)
            total_gains = 0;
        return total_gains;        
    }
    public int getDescendant1GeneLosses()
    {
        int total_losses = this.descendant1.total_losses;
        if (total_losses < 0)
            total_losses = 0;
        return total_losses;
    }
    public int getDescendant2GeneLosses()
    {
        int total_losses = this.descendant2.total_losses;
        if (total_losses < 0)
            total_losses = 0;
        return total_losses;        
    }     
}
