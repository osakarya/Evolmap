
package EvolMAP;
import java.util.*;
import java.io.*;

/**
 *
 * @author onur
 */
public class PfamRetriever 
{
    
    /** Creates a new instance of PfamRetriever */
    public PfamRetriever(String speciesName, String fastaTag, String domainName)
                            throws IOException
    {
        GeneDatabase GD = null;
        String hmmerOutputFile = speciesName + "." + domainName + ".hmmerOut";
        String outputFastaFile = speciesName + "." + domainName + fastaTag;;
        
        FileWriter output = create_file(new File(hmmerOutputFile));
        Process theProcess = null;

        try
        {
            theProcess = Runtime.getRuntime().exec("./hmmpfam --cut_ga " + domainName+".HMM" + " " + speciesName+fastaTag);
            InputStream stderr = theProcess.getErrorStream();
            OutputStream stdout = theProcess.getOutputStream();
            BufferedReader input = new BufferedReader(new InputStreamReader(theProcess.getInputStream()), 100000);
            GD = new GeneDatabase(speciesName+fastaTag, true);
            
            String lastProteinName = "lastProteinName";
            String line = "";
            boolean more_records_found = true;
            int count = 0;
            while (more_records_found && GD.getGene(GD.allGenes.length-1).header.indexOf(lastProteinName)==-1)
            {
                more_records_found = false;
                Thread.sleep(2000);
                while (more_records(input))
                {
                    more_records_found = true;
                    line = read_a_line(input);
                    if (line.startsWith("Query sequence: "))
                    {
                        lastProteinName = line.substring(16);
                        count++;
                    }
                    /*
                    if ( count%10==0)
                        System.out.println("Gene " + count + " read.");
                     */
                    output.write(line);
                    output.write(System.getProperty("line.separator"));
                }                
            }
            // theProcess.waitFor(); // This must finish
        }
        catch(Exception e)
        {
            System.err.println("Error on exec() method");
            e.printStackTrace();  
        }
        output.close();
        BufferedReader input = open_file(new File(hmmerOutputFile));
        output = create_file(new File(outputFastaFile));
        String line = "";
        String protein_name = "";
        String processed_protein_name = "";
        int count = 0;
        while (more_records(input))
        {
            line = read_a_line(input);                    
            if (line.startsWith("Query sequence: "))
            {
                protein_name = line.substring(16);
            }
            if (line.startsWith(domainName) && !processed_protein_name.equals(protein_name))
            {
                processed_protein_name = protein_name;
                count++;
                Sequence thisSeq = GD.getGeneByAnyName(protein_name);
                output.write(thisSeq.toStringWhole());
                output.write(System.getProperty("line.separator"));             
            }           
        }
        System.out.println("HMMER finished. " + count + " genes written to file " + outputFastaFile + ".");
        // Delete hmmer output file
        input.close();
        output.close();
        new File(hmmerOutputFile).delete();
    }
    
    
    private static BufferedReader open_file(File filename)
    {
            try {
                    return new BufferedReader(new FileReader(filename));
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
    private static FileWriter create_file(File filename)
    {
            try {
                    return new FileWriter(filename);
            }catch (Exception exception) {System.out.println(exception.getStackTrace());
            return null;}
    }    
    private static String removeEmptySpace(String s)
    {
        return s.replaceAll(" ", "");
    }        
}

