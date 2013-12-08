
package EvolMAP;
import java.io.*;
import java.util.*;
import org.omg.SendingContext.RunTime;
/**
 * @author onur
 */
public class BLAST_parser implements Runnable 
{
    int COMMON_BLAST_SIZE;
    int MIN_LENGTH = 50; // Alignment min length from BLAST for consideration as a score
    int MIN_BIT_SCORE = 50;
    int alignmentCount = 0;  // Count total number of alignments to do
    
    String[] currentGenes;
    //ArrayList<String>[] currentGeneMatchesNames;
    ArrayList<Integer>[] currentGeneMatchesIndexes;
    ArrayList<Integer>[] currentGeneMatchesScores;
    
    int id;
    int processors;
    String fileName;
    int alignments;
    boolean protein;
    
    public BLAST_parser (int id, int processors, String fileName, int alignments, boolean protein)
    {
        this.id = id;
        this.processors = processors;
        this.fileName = fileName;
        this.alignments = alignments;
        this.protein = protein;
    }
    public void run()
    {
        int lowerlimit = (int)(((double)((double)id/(double)processors))*(double)AGS_console.GD.getAllGenes().length);
        int upperlimit = (int)((((double)((double)(id+1)/(double)processors))*(double)AGS_console.GD.getAllGenes().length ));
        if (upperlimit > AGS_console.GD.getAllGenes().length )
            upperlimit = AGS_console.GD.getAllGenes().length;
        
        System.out.println("Thread: " + id + " LL: " + lowerlimit + " UL: " + upperlimit);
        
        String thisFileName = fileName + "_" + id;
        String outputName = thisFileName + ".out";
        // Prepeare new gene database
        Sequence[] thisSet = new Sequence[upperlimit-lowerlimit];
        for (int i=lowerlimit; i<upperlimit; i++)
        {
            thisSet[i-lowerlimit] = AGS_console.GD.getGene(i); 
        }
        GeneDatabase GD2 = new GeneDatabase("", thisSet);
        
        // Write to file
        try{
            FileWriter output = create_file(new File(thisFileName));
            for (int i=0; i<GD2.size(); i++)
            {
                output.write (GD2.allGenes[i].toStringCompact() + System.getProperty("line.separator"));
            }
            output.close();
        }
        catch (IOException e){};
        
        // Run blast on that file
        this.calculateBlastTableMatrix(fileName, thisFileName, alignments, outputName, protein);
        new File(thisFileName).delete();         
    }
    public static void merge(int processors, String fileName)
    {
        String mergedFile = fileName + ".BLAST";
        FileWriter output = null;
        try{
            output = create_file(new File(mergedFile));
            // Open files for each results, read them and write to output
            for (int i=0; i<processors; i++)
            {
                String inputName = fileName + "_" + i + ".out";
                BufferedReader input = open_file(new File(inputName));
                String line = "";
                while (more_records(input))
                {
                    line = read_a_line(input);  
                    output.write(line);
                    output.write(System.getProperty("line.separator"));;
                }
                input.close();
                new File(inputName).delete(); 
            }  
            output.close();
            
        }
        catch (IOException e){};  
        
    }    
    
    /** Creates a new instance of BLAST_parser */
    public BLAST_parser(int processors, String fileName, int ALIGNMENT_LIMIT, boolean RECALCULATE_BLAST_DATABASE, boolean DELETE_BLAST_RESULTS, boolean PROTEIN) 
                        throws IOException
    {
        String outputName = fileName + ".BLAST";
        COMMON_BLAST_SIZE = ALIGNMENT_LIMIT;
        
        if (RECALCULATE_BLAST_DATABASE)
        {
            this.formatDatabase(fileName, ALIGNMENT_LIMIT, outputName, PROTEIN);
            
            long time = -System.currentTimeMillis();
            System.out.println("Started running parallel");            
            
            Runnable runners[] = new Runnable[processors];
            Thread th[] = new Thread[processors];                
            for (int j=0; j<processors; j++)
            {
                runners[j] = new BLAST_parser(j, processors, fileName, ALIGNMENT_LIMIT, PROTEIN);
                th[j] = new Thread(runners[j]);
                th[j].start();
            }
            for(int j=0;j<processors;j++){
                try{
                th[j].join();
                }
                catch(InterruptedException e){}
            }
            merge(processors, fileName);
            
            time += System.currentTimeMillis();
            System.out.println (time + " msec"); 
        }
        
        
        currentGenes = new String[AGS_console.GD.size()];
        //currentGeneMatchesNames = new ArrayList[GD.size()];
        currentGeneMatchesIndexes = new ArrayList[AGS_console.GD.size()];
        //currentGeneMatchesIdentitites = new ArrayList[GD.size()];
        //currentGeneMatchesLengths = new ArrayList[GD.size()];
        currentGeneMatchesScores = new ArrayList[AGS_console.GD.size()]; // Bit score
        
        // Initialize array lists
        for (int i=0; i<AGS_console.GD.size(); i++)
        {
            //currentGeneMatchesNames[i] = new ArrayList<String>(COMMON_BLAST_SIZE);
            currentGeneMatchesIndexes[i] = new ArrayList<Integer>(COMMON_BLAST_SIZE/3);
            //currentGeneMatchesIdentitites[i] = new ArrayList<Double>(COMMON_BLAST_SIZE);
            //currentGeneMatchesLengths[i] = new ArrayList<Integer>(COMMON_BLAST_SIZE);
            //currentGeneMatchesScores[i] = new ArrayList<Double>(COMMON_BLAST_SIZE);
            currentGeneMatchesScores[i] = new ArrayList<Integer>(COMMON_BLAST_SIZE/3);
        }
        
        BufferedReader input = open_file(new File(outputName));
        if (input == null) 
            throw new IOException("Could not open Blast-all score files."); 
        
        int count = -1;  // First gene is counted immediately to 0
        int match_count = 0;
        String lastGene = "";
        String currentGene = "";
        String currentMatchName = "";
        int currentMatch = 0;
        String line = "";
        StringTokenizer st;
        int length = 0;
        int score = 0;
        int olderMatchIndex = 0;
        
        
        while (more_records(input))
        {
            line = read_a_line(input);

            // Get contents of line
            st = new StringTokenizer(line, "\t");
            currentGene = st.nextToken(); // query id
            currentMatchName = st.nextToken(); // subject id
            currentMatch = AGS_console.GD.getGeneLocationByPreComputedIndex(currentMatchName);
            st.nextToken(); //double identity = Double.parseDouble(st.nextToken()); // % identity
            length = Integer.parseInt(st.nextToken()); // al. length
            
            st.nextToken(); // mismatches
            st.nextToken(); //int gap_opens = Integer.parseInt(st.nextToken()); // gap openings
            st.nextToken(); //int qstart = Integer.parseInt(st.nextToken()); // query start
            st.nextToken(); //int qend = Integer.parseInt(st.nextToken()); // query end
            st.nextToken(); //int mstart = Integer.parseInt(st.nextToken()); // match - subject start
            st.nextToken(); //int mend = Integer.parseInt(st.nextToken()); // match - subject end
            st.nextToken(); // e-value
            score = (int)Double.parseDouble(st.nextToken());//double bit_score = Double.parseDouble(st.nextToken()); // bit score
            
            if (!PROTEIN || length < MIN_LENGTH)
                continue;   // Do not consider small alignments at all
            if (score < MIN_BIT_SCORE) // Do not consider scores less than 50(?)
                continue;
            
            // If the last gene is not same with to current gene, do not increment current count
            if (currentGene.equals(lastGene))
            {
                if(currentGeneMatchesIndexes[count].indexOf(currentMatch) != -1) // We got this score against this gene already so pick the better of them!
                {
                    // If new match has better identity than the first one, then replace this one
                    olderMatchIndex = currentGeneMatchesIndexes[count].indexOf(currentMatch);
                    if (currentGeneMatchesScores[count].get(olderMatchIndex) < score) //Takes the final score as current score
                    {
                        //currentGeneMatchesIdentitites[count].set(olderMatchIndex, identity);
                        //currentGeneMatchesLengths[count].set(olderMatchIndex, length);
                        //currentGeneMatchesScores[count].set(olderMatchIndex, bit_score);
                        currentGeneMatchesScores[count].set(olderMatchIndex, score);
                    }
                    
                }
                else // This gene's hit is being added for first time
                {
                    alignmentCount++;
                    currentGeneMatchesIndexes[count].add(currentMatch);
                    //currentGeneMatchesIdentitites[count].add(identity);
                    //currentGeneMatchesLengths[count].add(length);
                    //currentGeneMatchesScores[count].add(bit_score);
                    currentGeneMatchesScores[count].add(score);
                }
            }
            else
            {
                count++;
                
                // Not the next gene in the database is our match!!!! -- So skip this and all its arrays will be empty
                if (!AGS_console.GD.allGenes[count].getGeneName2().equals(currentGene)) 
                {
                    System.out.println("Warning: A gene was skipped by blast: " + AGS_console.GD.allGenes[count].getGeneName2());
                    // Increment count till you find the right gene
                    while(!AGS_console.GD.allGenes[count].getGeneName2().equals(currentGene))
                    {
                        currentGenes[count] = AGS_console.GD.allGenes[count].getGeneName2();
                        count++;
                    }
                }
                
                currentGenes[count] = currentGene;
                
                alignmentCount++;
                currentGeneMatchesIndexes[count].add(currentMatch);
                //currentGeneMatchesIdentitites[count].add(identity);
                //currentGeneMatchesLengths[count].add(length);
                //currentGeneMatchesScores[count].add(bit_score);
                currentGeneMatchesScores[count].add(score); 
                
                lastGene = currentGene;
                
            }
         }
        count++;
        // Add genes missing from the file back to the list finally--
        while(count<AGS_console.GD.size())
        {
            currentGenes[count] = AGS_console.GD.allGenes[count].getGeneName2();
            count++;
        }
        
        input.close();
        
        // Delete blast scores, temp database file and blast database files
        if (DELETE_BLAST_RESULTS)
            this.deleteExtraFiles(fileName, outputName);
        
    }
    public static void formatDatabase (String fileName, int ALIGNMENT_LIMIT, String outputName, boolean PROTEIN)
    {
      System.out.print("Creating BLAST database...");
      Process theProcess = null;
      
      try
      {   
          if (PROTEIN)
          {
              if (System.getProperty("os.name").contains("Windows"))
              {
                  theProcess = Runtime.getRuntime().exec("formatdb -i " + fileName);  
              }             
              else
              {
                  theProcess = Runtime.getRuntime().exec("./formatdb -i " + fileName);  
              }
          }
          else
          {
              if (System.getProperty("os.name").contains("Windows"))
              {
                  theProcess = Runtime.getRuntime().exec("formatdb -p F -i " + fileName);  
              }             
              else
              {
                  theProcess = Runtime.getRuntime().exec("./formatdb -p F -i " + fileName);  
              }
          }
          theProcess.waitFor();
      }
      catch(Exception e)
      {
         System.err.println("Error on exec() method");
         e.printStackTrace();  
      }  
      System.out.println("complete.");
    }
    
    public static void calculateBlastTableMatrix(String databaseName, String fileName, int ALIGNMENT_LIMIT, String outputName, boolean PROTEIN)
    {
      Process theProcess = null;
 
      System.out.print("Running all-to-all BLAST...");
      
      String runString;
      if (PROTEIN)
      {
          if (System.getProperty("os.name").contains("Windows"))
          {
            runString = "blastall -p blastp -d " + databaseName + " -i " + fileName + " -o " + outputName + " -m 8" + " -b " + ALIGNMENT_LIMIT;
          }
          else
          {
            runString = "./blastall -p blastp -d " + databaseName + " -i " + fileName + " -o " + outputName + " -m 8" + " -b " + ALIGNMENT_LIMIT; 
          }
      
      }
      else
      {
          if (System.getProperty("os.name").contains("Windows"))
          {
            runString = "blastall -p blastn -d " + databaseName + " -i " + fileName + " -o " + outputName + " -m 8" + " -b " + ALIGNMENT_LIMIT;
          }
          else
          {
            runString = "./blastall -p blastn -d " + databaseName + " -i " + fileName + " -o " + outputName + " -m 8" + " -b " + ALIGNMENT_LIMIT;  
          }
      }
      
      try
      {
          theProcess = Runtime.getRuntime().exec(runString);  
          theProcess.waitFor(); 
      }
      catch(Exception e)
      {
         System.err.println("Error on exec() method");
         e.printStackTrace();  
      }   
      System.out.println("complete.");
      // read from the called program's standard output stream        
    }
   
    public static void deleteExtraFiles(String fileName, String outputName)
    {
        new File(fileName).delete(); // Delete temp database
        new File(outputName).delete(); // Delete blast hit file
        new File(fileName + ".phr").delete(); // Delete blast database file
        new File(fileName + ".pin").delete(); // Delete blast database file
        new File(fileName + ".psq").delete(); // Delete blast database file
        new File(fileName + ".psd").delete(); // Delete blast database file
        new File(fileName + ".psi").delete(); // Delete blast database file
    }
    
    public ArrayList<Integer> getAllHitsForGene(int i)
    {
        return currentGeneMatchesIndexes[i];
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

