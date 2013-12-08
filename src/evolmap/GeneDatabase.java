
package EvolMAP;
import java.io.*;
import java.util.*;

/**
 * @author onur
 */
public class GeneDatabase 
{
    Sequence[] allGenes;
    String species;
    String[] species_list;
    int[] species_start_indexes;
    private int LENGTH_MIN = 0; // Sequences shorter than this are removed
    private int LENGTH_MAX = 10000; // Sequences longer than this are cut to this length to prevent out of heap space alignments
    private boolean REMOVE_X_CHARACTERS = true; // Sequences with X characters are deleted.
    
    /** Creates a new instance of GeneDatabase */
    public GeneDatabase(String sequenceDatabaseFile)
                            throws IOException
    {
        this.species = "";
        readSequences(sequenceDatabaseFile, false);
        setSpeciesIndexes();
    }
    /** Creates a new instance of GeneDatabase */
    public GeneDatabase(String sequenceDatabaseFile, boolean countFirst)
                            throws IOException
    {
        this.species = "";
        readSequences(sequenceDatabaseFile, countFirst);
        setSpeciesIndexes();
    }
    
    public GeneDatabase(String species, Sequence[] genes)
    {
        this.species = species;
        this.allGenes = genes;
    }
    
    public Sequence[] getAllGenes() 
    {
        return allGenes;
    }    
    
    public void setSpeciesIndexes()
    {
        Vector<String> speciesNames = new Vector<String>();
        Vector<Integer> speciesStartIndexes = new Vector<Integer>();
        for (int i=0; i<allGenes.length; i++)
        {
            Sequence seq = allGenes[i];
            if (!speciesNames.contains(seq.species))
            {
                speciesNames.add(seq.species);
                speciesStartIndexes.add(i);
            }
        }
        this.species_list = new String[speciesNames.size()];
        this.species_start_indexes = new int[speciesStartIndexes.size()];
        
        for (int i=0; i<this.species_list.length; i++)
        {
            species_list[i] = speciesNames.get(i);
            species_start_indexes[i] = speciesStartIndexes.get(i);
        }
    }
    public int getSpeciesStartIndex(String species)
    {
        for (int i=0; i<this.species_list.length; i++)
        {
            if (species.equals(species_list[i]))
            {
                return this.species_start_indexes[i];
            }
        }
        return -1;
    }
    public Sequence getGeneByPreComputedIndex(String gene_name)
    {
        return this.getGene(this.getGeneLocationByPreComputedIndex(gene_name));
    }
            
    public int getGeneLocationByPreComputedIndex(String gene_name)
    {
        StringTokenizer st = new StringTokenizer(gene_name, "_");
        String species = st.nextToken();
        String nextToken = st.nextToken();
        while(st.hasMoreTokens())
            nextToken = st.nextToken();
        int speciesOrder = Integer.parseInt(nextToken);
        int speciesStartIndex = this.getSpeciesStartIndex(species);
        if (speciesStartIndex == -1)
            System.out.println("Something wrong going on with indexes. Check method getGeneLocationByPreComputedIndex in GeneDatabase.java");
        return (speciesStartIndex + speciesOrder);
    }
    public void setScoreMatrixPositions()
    {
         // Since everything is set in order in score matrices before, above is not necessary really
          for (int i=0; i<this.allGenes.length; i++)
          {
              this.allGenes[i].scoreMatrixPosition = i;
          }
    }
    public GeneDatabase getSpeciesGenes(String species)
    {
        Vector<Sequence> speciesGenes = new Vector<Sequence>();
        
        for (int i=0; i<this.allGenes.length; i++)
        {
            if (this.allGenes[i].species!=null && this.allGenes[i].species.equals(species))
            {
                speciesGenes.add(this.allGenes[i]);
            }
        }
        Sequence[] speciesGenesArray = new Sequence[speciesGenes.size()];
        for (int i=0; i<speciesGenes.size(); i++)
        {
            speciesGenesArray[i] = speciesGenes.get(i);
        }
        GeneDatabase speciesBase = new GeneDatabase(species, speciesGenesArray);
        return speciesBase;
    }
    
    public int size()
    {
        return allGenes.length;
    }
    
    public Sequence getGene(int i)
    {
        return allGenes[i];
    }
    
    public Sequence getGeneByAnyName(String g)
    {
        for (int i=0; i<this.allGenes.length; i++)
        {
            if ((this.allGenes[i]).header.contains(g))
                return this.allGenes[i];
        }
        return null;
    }
    public Sequence getGeneByName2(String g) // Used only in extra methods -- can be removed
    {
        for (int i=0; i<this.allGenes.length; i++)
        {
            if ((this.allGenes[i]).getGeneName2().equals(g))
                return this.allGenes[i];
        }
        return null;        
    }  
    public Vector<String> getAllGenesByName2(String[] geneNames)
    {
        Vector<String> foundGenes = new Vector<String>();
        for (int i=0; i<geneNames.length; i++)
        {
            Sequence currentSequence = this.allGenes[getGeneLocationByPreComputedIndex(geneNames[i])];
            if (currentSequence==null)
                System.out.println("Read error in getAllGenesByName2 file-- gene " +  geneNames[i] + " not found in database");
            else
                foundGenes.add(currentSequence.getGeneName2());
        }
        return foundGenes;
    }
    
    public void readSequences(String fileName, boolean countFirst)
		throws IOException
    {
        BufferedReader input = open_file(new File(fileName));
        if (input == null) 
        {
            throw new IOException("Could not open gene database file " + fileName); 
        }
        int number_of_sequences = 0;
        String line = "";
        
        if (countFirst)
        {
            while (more_records(input))
            {
                if (line.startsWith(">"))
                {
                   number_of_sequences++;
                   line = read_a_line(input);
                   while (line.startsWith(">"))
                   {
                       line = read_a_line(input);
                   }
                }
                else
                {
                    line = read_a_line(input);
                }
            }
            input.close();
            input = open_file(new File(fileName));
            line = read_a_line(input);
            this.allGenes = new Sequence[number_of_sequences];
            number_of_sequences = 0;
        }
        else
        {
            line = read_a_line(input);

            // First line is the database size
            this.allGenes = new Sequence[Integer.parseInt(line)];
            line = read_a_line(input);
        }
        while (more_records(input))
        {
            String current_sequence_name = "";
            String current_sequence = "";
            if (line == null) 
                throw new IOException("Problem reading file.");
            // Skip all blank lines
            while (line!=null && !line.startsWith(">")) 
                line = read_a_line(input);

            if (line.startsWith(">"))
            {
                current_sequence_name = line;
                line  = read_a_line(input);
                while (line!=null && (line.length()!= 0) && line.startsWith(">"))
                {
                        current_sequence_name += System.getProperty("line.separator") + line;					
                        line  = read_a_line(input);
                }
                // Skip all blank lines
                while ((line.length()== 0)) 
                    line = read_a_line(input);

                while (line!=null && !line.startsWith(">"))
                {  
                    current_sequence += line;
                    line  = read_a_line(input);
                }	

                if (current_sequence_name.length()!=0 && current_sequence.length()!=0)
                {
                    if (current_sequence.length()>LENGTH_MAX)
                        current_sequence = current_sequence.substring(0, LENGTH_MAX);
                    
                    if (REMOVE_X_CHARACTERS)
                    {
                        current_sequence = current_sequence.replaceAll("X", "");
                        if (current_sequence.length()==0)
                            current_sequence = "A";
                    }
                    
                    if (current_sequence.length()>LENGTH_MIN)
                    {
                        this.allGenes[number_of_sequences] = new Sequence(current_sequence_name, current_sequence);
                        number_of_sequences++;
                    }
                }
            }
            else
            {
                    line = read_a_line(input);
            }
        }		

        System.out.println ("Fasta file " + fileName + " is opened and " + number_of_sequences + " sequences read.");
        
        
	input.close();
    }
     
    // This method is only used to convert a single species fasta file to OrthologGenes vector
    // Therefore it includes species names.
    public Vector<OrthologGenes> getOrthologGenesVector()
    {
        Vector<OrthologGenes> og = new Vector<OrthologGenes>();        
        for (int i=0; i<this.size(); i++)
        {
            og.add(new OrthologGenes(this.getGene(i).getGeneName2()));
        }
        return og;
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

    public String toString()
    {
        String returnString = "";
        for (int i=0; i<this.allGenes.length; i++)
        {
            returnString += allGenes[i].toStringWhole() + System.getProperty("line.separator");
        }
        return returnString;
    }
    public String toStringCompact()
    {
        String returnString = "";
        for (int i=0; i<this.allGenes.length; i++)
        {
            returnString += allGenes[i].toStringCompact() + System.getProperty("line.separator");
        }
        return returnString;        
    }
}
