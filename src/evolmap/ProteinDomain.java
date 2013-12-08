
package EvolMAP;

import java.util.*;
import java.io.*;

/**
 * @author onur
 */
public class ProteinDomain
{
    String geneName;
    String domainName;
    int domainPosition; // The ordered position of this domain on its scaffold
    int domainStart;
    int domainEnd;
    int domainLength;
    String domainSequence;
    String domainScore;    // Pfam score
    String domainEvalue;   // Pfam evalue -- this value changes with dataset size
    
    /** Creates a new instance of ProteinDomain */
    public ProteinDomain(String geneName, String domainHeaderLine, String sequence, int realDomainPosition) 
    {
        this.geneName = geneName;
        StringTokenizer st = new StringTokenizer(domainHeaderLine.substring(1));
        
        //System.out.println(geneName);
        
        this.domainName = st.nextToken();
        String domainPos = st.nextToken();   // skip this because it comes as extra input now
        //this.domainPosition = Integer.parseInt((String)(domainPos.substring(0, domainPos.indexOf("/")))); // Then cut the extra part
        this.domainStart = Integer.parseInt(st.nextToken());
        this.domainEnd = Integer.parseInt(st.nextToken());
        this.domainScore = st.nextToken();
        this.domainEvalue = st.nextToken();
        this.domainSequence = sequence.substring(domainStart, domainEnd);
        this.domainPosition = realDomainPosition;
        this.domainLength = domainEnd - domainStart + 1;
            
    }
    public ProteinDomain(String speciesName, String geneNo, String domainNo, String domainType, String domainSequence)
    {
        this.geneName = speciesName + "\t" + geneNo + "\t" + domainNo + "\t" + domainType;
        this.domainPosition = Integer.parseInt(domainNo);
        this.domainName = domainType;
        this.domainLength = domainSequence.length();
        this.domainSequence = domainSequence;
        this.domainEnd = domainSequence.length();
        this.domainStart = 0;
        this.domainEvalue = "";
    }
    
    public String getPrintSequence()
    {
        String printSequence = "";

        for (int i=0; i<domainSequence.length(); i++)
        {
            if ((i>0) && (i%80 == 0) && i<domainSequence.length()-1)
                printSequence += System.getProperty("line.separator");
            printSequence += domainSequence.charAt(i);

        }
        return printSequence;
    }    
}
