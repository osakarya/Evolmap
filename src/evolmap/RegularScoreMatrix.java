
package EvolMAP;

import java.util.*;
import java.io.*;

/**
 *
 * @author onur
 */
public class RegularScoreMatrix extends ScoreMatrix
{
    int[][] scores;
    
    /** Creates a new instance of ScoreMatrix */
    public RegularScoreMatrix(String fileName)
                        throws IOException
    {
            BufferedReader input = open_file(new File(fileName));
            if (input == null) 
                throw new IOException("Could not open score matrix file " + fileName); 
                   
            String line = read_a_line(input);
            StringTokenizer st = new StringTokenizer(line.substring(1), "\t");
            
            this.names = new String[st.countTokens()];
            this.scores = new int[st.countTokens()][st.countTokens()];
            
            int count = 0;
            while (st.hasMoreTokens())
            {
                this.names[count] = st.nextToken();
                count++;
            }
           
            // Read line by line and fetch every score to list	 
            int row_count = 0;
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
                    int nextNumber = Integer.parseInt(st.nextToken());
                    this.scores[row_count][column_count] = nextNumber;
                    column_count++;
                }
                row_count++;
            }
            input.close();
    }
    
    public int getElement(int i, int j)
    {
        return this.scores[i][j];
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
