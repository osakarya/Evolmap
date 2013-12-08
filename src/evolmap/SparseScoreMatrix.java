
package EvolMAP;

import java.util.*;
import java.io.*;

/**
 *
 * @author onur
 */
public class SparseScoreMatrix extends ScoreMatrix
{
    //ArrayList<String>[] matchNames;
    ArrayList<Integer>[] matchIndexes;
    ArrayList<Integer>[] scores;
    
    /** Creates a new instance of ScoreMatrix */
    public SparseScoreMatrix(String fileName, int alignmentLimit)
                        throws IOException
    {
            System.out.println("Reading sparse matrix.");
            
            BufferedReader input = open_file(new File(fileName));
            if (input == null) 
                throw new IOException("Could not open sparse matrix file " + fileName); 
            
            String line = read_a_line(input);
            StringTokenizer st = new StringTokenizer(line.substring(1), "\t");
            
            this.names = new String[AGS_console.GD.size()];
            //this.matchNames = new ArrayList[GD.size()];
            this.matchIndexes = new ArrayList[AGS_console.GD.size()];
            this.scores = new ArrayList[AGS_console.GD.size()];
            
            for (int i=0; i<AGS_console.GD.size(); i++)
            {
                //this.matchNames[i] = new ArrayList<String>(alignmentLimit);
                this.matchIndexes[i] = new ArrayList<Integer>(alignmentLimit/3);
                this.scores[i] = new ArrayList<Integer>(alignmentLimit/3);
            }
            
            int count = 0;
            while (st.hasMoreTokens())
            {
                this.names[count] = st.nextToken();
                count++;
            }
            
           
            // Read line by line and fetch every score to list	 
            int row_count = 0;
            String nextHit = "";
            while (more_records(input))
            {
                line = read_a_line(input);
                if ((line.length()==0))
                {
                    break;
                }
                st = new StringTokenizer(line, "\t");
                
                int column_count = 0;
                while (st.hasMoreTokens())
                {
                    nextHit = st.nextToken();
                    //this.matchNames[row_count].add(nextHit);
                    this.matchIndexes[row_count].add(AGS_console.GD.getGeneLocationByPreComputedIndex(nextHit));
                    int nextNumber = Integer.parseInt(st.nextToken());
                    this.scores[row_count].add(nextNumber);
                    column_count++;
                }
                row_count++;
            }
            input.close();
            
            // Set ONLY self-scoring genes
            for (int i=0; i<this.matchIndexes.length; i++)
            {
                if (this.matchIndexes[i].size()<=1) // Scores only itself or (nothing--which should not happen)
                {
                    AGS_console.GD.allGenes[i].noScoresFound = true;
                }
            } 
    }
    
    public int getElement(int i, int j)
    {
        String secondGene = this.names[j];
        //int k = this.matchNames[i].indexOf(secondGene);
        int k = this.matchIndexes[i].indexOf(j);
        if (k>=0)
            return this.scores[i].get(k);
        else
            return 0;
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

