
package EvolMAP;

import java.util.*;
import java.io.*;

/**
 *
 * @author onur
 */
public class SequenceList
{
	LinkedList list;
        int alignment_length;


	public SequenceList()
	{
		this.list = new LinkedList();
                this.alignment_length = 0;
	}
	
	public SequenceList(LinkedList list)
	{
		this.list = list;
                this.alignment_length = 0; // Set while reading sequences
	}
	
	// The typical constructor, takes a filename, and fetches sequences from it
	public SequenceList(String fileName)
		throws IOException
	{
		this.list = new LinkedList();
		convertFileToSequenceVector(fileName);		
	}
	
	public void convertFileToSequenceVector(String fileName)
		throws IOException
	{
		BufferedReader input;
		File inputFile = new File(fileName);
		
		// Try to open the file
		input = open_file(inputFile);
		
		if (input == null)
			throw new IOException("Could not open input file " + inputFile );			
				
		int number_of_sequences = 0;
		
		String line = read_a_line(input);
			
		// Read line by line and fetch every protein into a list	 
		while (more_records(input))
		{
			
			String current_sequence_name = "";
			String current_sequence = "";
			
			if (line == null)
				fatal("problem reading file");

			// Skip all blank lines
			while (!line.startsWith(">"))
			{
				line = read_a_line(input);
			}
			
			if (line.startsWith(">"))
			{
				current_sequence_name = line;
				line  = read_a_line(input);
				while (line!=null && (line.length()!= 0) && line.startsWith(">"))
				{
					current_sequence_name += 
							System.getProperty("line.separator") + line;					
					line  = read_a_line(input);
					
				}

				// Skip all blank lines
				while ((line.length()== 0))
				{
					line = read_a_line(input);
				}
			
				
				while (line!=null && !line.startsWith(">"))
				{                                        
                                        current_sequence += line;
					
					line  = read_a_line(input);
                                        					
				}	
				
				if (current_sequence_name.length()!=0 && current_sequence.length()!=0)
				{
					this.addElement(new Sequence(current_sequence_name, current_sequence));
                                        this.alignment_length = current_sequence.length();
					number_of_sequences++;
				}
							
			}
			else
			{
				line = read_a_line(input);
			}
		}		
		
		input.close();
		
		System.out.println ("File " + fileName + " is opened and " + number_of_sequences + " read!");
	
	}
	
	public void addElement(Sequence element)
	{
		this.list.add((Sequence)element);
	}

	public void removeElement(Sequence element)
	{
		this.list.remove((Sequence)element);
	}
        
        // Return sequence that includes given header
	public Sequence getSequence(String header)
	{
		ListIterator iterator = this.list.listIterator();
		SequenceList newList = new SequenceList();
		
		while (iterator.hasNext())
		{
			Sequence nextElement = (Sequence)iterator.next();
			String name = nextElement.getName();
			// Write them to file as well
			if (name.contains(header))
			{
				return nextElement;
			}			
		}
		
		return null;		
	}
        public Vector getAllSequences()
        {
                Vector sequences = new Vector<Sequence>();
		ListIterator iterator = this.list.listIterator();
		SequenceList newList = new SequenceList();
		
		while (iterator.hasNext())
		{
			Sequence nextElement = (Sequence)iterator.next();
			sequences.add(nextElement);			
		}
		
		return sequences;            
        }

	public String toString()
	{
		String wholeList = "";
		ListIterator iterator = this.list.listIterator();
		
		while (iterator.hasNext())
		{
			Sequence nextElement = (Sequence)iterator.next();
			wholeList += nextElement.toString();
		}
		
		return wholeList;
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
	
	private static void fatal(String message)
	{
		System.out.println(message + "\nProgram terminating");
		System.exit(1);
	}
	
}