
package EvolMAP;

import java.util.*;
import java.io.*;


public class Sequence
{
        String species;
        String header;
	String name;
        String gene_name;
        Vector<ProteinDomain> domains;
	String sequence;
        int scoreMatrixPosition;
        boolean displayFullName = false;
        boolean noScoresFound = false; // If doing sparse-scoring and no scores except itself (and perhaps lacking even itself) found, then tag it true somewhere -- currently done in CalculateBlastFirstScoreMatrices.java
	boolean divergedInparalog = false;
        
        // Copy constructor
        public Sequence(Sequence copySequence)
        {
            this.species = copySequence.species;
            this.header = copySequence.header;
            this.name = copySequence.name;
            this.gene_name = copySequence.gene_name;
            this.domains = copySequence.domains;
            this.sequence = copySequence.sequence;
            this.scoreMatrixPosition = copySequence.scoreMatrixPosition;
            this.displayFullName = copySequence.displayFullName;
            this.noScoresFound = copySequence.noScoresFound;
            this.divergedInparalog = false;
            
        }
        
	public Sequence(String header, String sequence)
	{
            // A sequence should either have a single original header line
            // Or it can have a gene database format (which gives species and order in first line followed by original line-- first line is in the form human 9 etc.)
            // Or it can have a domain format which is like human 9 3 PDZ GENENAME
            // Or it can have a gene format plus each domain in each line where each domain has 6 attributes-- check example files 

            this.species = null; // If null no species information was entered
            this.header = header;    
            StringTokenizer st = new StringTokenizer(header, System.getProperty("line.separator"));
            this.name = st.nextToken();
            if (this.name.indexOf("\t")!= -1)
            {
                this.gene_name = this.name.substring(1);
                this.species = (new StringTokenizer(gene_name).nextToken());
                if (st.hasMoreTokens())
                    this.name = st.nextToken();
            }
            if (gene_name == null)
                gene_name = "";

            StringTokenizer domainSt = new StringTokenizer(gene_name);

            if (domainSt.countTokens() == 4 || domainSt.countTokens() == 5)
            {
                this.domains = new Vector<ProteinDomain>();
                String speciesName = domainSt.nextToken();
                String geneNo = domainSt.nextToken();
                String domainNo = domainSt.nextToken();
                String domainType = domainSt.nextToken();
                String domainSequence = sequence;
                ProteinDomain nextDomain = new ProteinDomain (speciesName, geneNo, domainNo, domainType, domainSequence);
                domains.add(nextDomain);
            }
            else
            {
                // Domain partitioning
                int domain_count = 0;
                if (st.hasMoreTokens())
                {
                    this.domains = new Vector<ProteinDomain>();
                }
                while(st.hasMoreTokens())
                {
                    String nextToken = st.nextToken();
                    StringTokenizer st2 = new StringTokenizer (nextToken);
                    // This is the format we want, otherwise dont bother doing this.
                    if (st2.countTokens()==6) // Reading domains in this
                    {
                        ProteinDomain nextDomain = new ProteinDomain(this.name, nextToken, sequence, domain_count);
                        domains.add(nextDomain);
                        domain_count++;
                    }
                    else
                    {
                        break;
                    }
                }
            }
            // Remove any *s from the sequence
            String newSequence = "";
            for (int i=0; i<sequence.length(); i++)
            {
                if (!(sequence.charAt(i)=='*'))
                    newSequence += sequence.charAt(i);
            }
            this.sequence = newSequence;
	}
        // Used by ShowAncestor.java
        public void displayFullName()
        {
            this.displayFullName = true;
        }
		
	public String getName()
	{
		return name;
	}
        
	public String getGeneName()
	{
		return gene_name;
	}
        
	public String getSequence()
	{
		return sequence;
	}
        
        public String getDomainSequence()
        {
            return this.concatDomainSequences();
        }   
        
        public String concatDomainSequences()
        {
            String newSequence = "";
            for (int i=0; i<domains.size(); i++)
            {
                ProteinDomain nextDomain = (ProteinDomain) domains.get(i);
                newSequence += nextDomain.domainSequence;
            }
            
            return newSequence;
        }
        
        public String getPrintSequence()
        {
            String printSequence = "";
            
            for (int i=0; i<sequence.length(); i++)
            {
                if ((i>0) && (i%80 == 0) && i<sequence.length()-1)
                    printSequence += System.getProperty("line.separator");
                printSequence += sequence.charAt(i);

            }
            return printSequence;
        }
        
        public void setSpecies(String species)
        {
            this.species = species;
        }
        public String getSpecies()
        {
            return this.species;
        }
        
        public String getDomainPrintSequence()
        {
            String newSequence = " ";
            if (domains==null)
                return newSequence;
            for (int i=0; i<domains.size(); i++)
            {
                ProteinDomain nextDomain = domains.get(i);
                ProteinDomain followingDomain = null;
                if (i<domains.size()-1)
                    followingDomain = domains.get(i+1);
                
                if (nextDomain.domainName.length()>3 && nextDomain.domainName.charAt(nextDomain.domainName.length()-2) == '_')
                    nextDomain.domainName = nextDomain.domainName.substring(0, nextDomain.domainName.length()-2);
                if (followingDomain!=null && followingDomain.domainName.length()>3 && followingDomain.domainName.charAt(followingDomain.domainName.length()-2) == '_')
                    followingDomain.domainName = followingDomain.domainName.substring(0, followingDomain.domainName.length()-2);
                
                // If domains overlap with next domain, skip this one
                if (followingDomain!=null && (followingDomain.domainStart < (nextDomain.domainEnd - ((nextDomain.domainLength)/2))))
                {
                    continue;
                }
                // Non-overlap small repeat domains are to be represented as single domain as well 
                else if ((nextDomain.domainLength<50) && (followingDomain!=null) && (nextDomain.domainName.equals(followingDomain.domainName)))
                {
                }
                else if (followingDomain==null)
                {
                    newSequence += nextDomain.domainName; 
                }
                else
                {
                    newSequence += nextDomain.domainName + "-"; 
                }
            }
            return newSequence;            
        }   
	
        public int getLength()
	{
		return sequence.length();
	}
        
        public String getGeneName2()
        {
            StringTokenizer st = new StringTokenizer(gene_name);
            String returnString = "";
            while (st.hasMoreTokens())
            {
                returnString += st.nextToken() + "_";
            }
            returnString = returnString.substring(0, returnString.length()-1);
            return returnString;
        }
        
	public String getGeneName3()
        {
            if (this.name.indexOf(" ") < 0)
                return this.getGeneName2();
            else
                return (this.name.substring(1, this.name.indexOf(" ")));
        }
        
	public boolean equals(Object obj)
	{
		if ((((Sequence)obj).getName() == this.name)
			&& (((Sequence)obj).getSequence() == this.sequence))
			{
				return true;
			}
		else
			return false;
		
	}
        public String toString()
        {
            if (this.displayFullName)
                return this.name + " (" + this.getDomainPrintSequence() + ")";
            else
            {
                if (species == null)
                    return this.getDomainPrintSequence();
                else
                {
                    if (divergedInparalog)
                        return (this.getGeneName2() + " (DIVERGED) " + this.getDomainPrintSequence());
                    else
                        return (this.getGeneName2() + this.getDomainPrintSequence());
                }
            }
                
        }
	public String toStringWhole()
	{
		return header + System.getProperty("line.separator") + this.getPrintSequence()
		+ System.getProperty("line.separator");
	}
        public String toStringCompact()
        {
            return ">" + this.getGeneName2()+ System.getProperty("line.separator") + this.getPrintSequence()
		+ System.getProperty("line.separator"); 
        }
        public String toStringOriginal()
        {
            StringTokenizer st = new StringTokenizer(header, "\n");
            st.nextToken();
            return (st.nextToken() + System.getProperty("line.separator") + this.getPrintSequence()
            + System.getProperty("line.separator"));            
        }
        public String toStringGeneID()
        {
            return header.substring(header.indexOf("gene:")+5, header.indexOf("transcript:")-1);  
        }
}