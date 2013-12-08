
package EvolMAP;
import java.io.*;
import java.util.*;


public class CreateSequenceDatabaseFile 
{    
    public CreateSequenceDatabaseFile(String[] INPUT_NAMES, String SEQUENCE_TAG, String OUTPUT_NAME)
                                    throws IOException
    {
        FileWriter output = create_file(new File(OUTPUT_NAME));
        if (output == null) throw new IOException("Could not open output file."); 
        SequenceFile[] allSequences = new SequenceFile[INPUT_NAMES.length];
        for (int i=0; i<allSequences.length; i++)
        {
            allSequences[i] = prepareSequenceFile(INPUT_NAMES[i] + SEQUENCE_TAG);
        }
        // Calculate total size
        int total_size = 0;
        for (int i=0; i<allSequences.length; i++)
        {
            total_size += allSequences[i].getAllSequences().size();
        }
        output.write("" + total_size);
        output.write(System.getProperty("line.separator"));
        for (int i=0; i<allSequences.length; i++)
        {
            int gene_count = 0;
            Vector<Sequence> thisSequences = allSequences[i].getAllSequences();
            for(int j=0; j<thisSequences.size(); j++)
            {
                output.write(">" + INPUT_NAMES[i] + "\t" + gene_count);
                output.write(System.getProperty("line.separator"));
                output.write(thisSequences.get(j).toStringWhole());
                output.write(System.getProperty("line.separator"));
                gene_count++;
            }
        }        

        output.close();        
    }
    
    public SequenceFile prepareSequenceFile (String fileName)
    throws IOException
    {
        SequenceList sequenceList = new SequenceList(fileName);
        SequenceFile sequenceFile = new SequenceFile(fileName, sequenceList);
        return sequenceFile;
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
