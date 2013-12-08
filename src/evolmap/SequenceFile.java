
package EvolMAP;

import java.util.*;
import java.io.*;

/**
 *
 * @author onur
 */
public class SequenceFile
{
	String fileName;
	SequenceList list;
	
	public SequenceFile(String fileName, SequenceList list)
	{
		this.fileName = fileName;
		this.list = list;
	}

	public void setfileName(String fileName)
	{
		this.fileName = fileName;
	}
	
	public void setSequenceList(SequenceList list)
	{
		this.list = list;
	}
			
	public String getfileName()
	{
		return fileName;
	}
	
	public SequenceList getSequenceList()
	{
		return list;
	}
        
        public Sequence getSequence(String header)
        {
            return this.list.getSequence(header);
        }
        
        public Vector getAllSequences()
        {
            return this.list.getAllSequences();
        }
	public String toString()
	{
		return list.toString();
	}      
}