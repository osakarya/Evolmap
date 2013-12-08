
package EvolMAP;
import java.io.*;
import java.util.*;

/**
 * @author onur
 */
public class CalculateBlastFirstScoreMatrices implements Runnable
{
    static private ScoringMatrix scoring; // The scoring matrix
    static private boolean PROTEIN_ALIGNMENT; // If protein alignment, it is made faster by first preparing int vectors
    static private boolean LENGTH_NORMALIZED_SCORES = true;
    static private int[][] ALL_INT_ARRAYS; // Used for fast alignments
    static int[] MAX_SCORES;
    static String outputName;
    static int ALIGNMENT_LIMIT;
    
    int id;
    int processors;
    boolean protein;   
    boolean mode;    
    static int totalAlignments;
    
    public CalculateBlastFirstScoreMatrices (int id, int processors, boolean protein, boolean mode)
    {
        this.id = id;
        this.processors = processors;
        this.protein = protein;
        this.mode = mode;
        totalAlignments = 0;
    }
    
    public void run()
    {
        int lowerlimit = (int)(((double)((double)id/(double)processors))*(double)AGS_console.GD.getAllGenes().length);
        int upperlimit = (int)((((double)((double)(id+1)/(double)processors))*(double)AGS_console.GD.getAllGenes().length ));
        if (upperlimit > AGS_console.GD.getAllGenes().length )
            upperlimit = AGS_console.GD.getAllGenes().length;
        
        if (mode==true)
        {
            try {
                // Calculate max scores
                for (int i=lowerlimit; i<upperlimit; i++)
                {
                    Sequence currentSequence = AGS_console.GD.getGene(i);
                    MAX_SCORES[i] = getScore(currentSequence, currentSequence, i, i);
                }
            }
            catch (Exception e) {};

            System.out.println("Max scores calculated.");
        }

        if (mode==false)
        {
            ArrayList<Integer>[] ALL_SCORES = new ArrayList[upperlimit-lowerlimit];
            for (int i=0; i<ALL_SCORES.length; i++)
            {
                ALL_SCORES[i] = new ArrayList<Integer>(ALIGNMENT_LIMIT/3);
            }            
            FileWriter output = create_file(new File(outputName + "_" + id));     
            
            for (int i=lowerlimit; i<upperlimit; i++)
            {
                Sequence currentSequence = AGS_console.GD.getGene(i);
                int maxScore1 = MAX_SCORES[i];
                for (int j=0; j<AGS_console.blast_parser.currentGeneMatchesIndexes[i].size(); j++)
                {
                    totalAlignments++;
                    Sequence otherSequence = AGS_console.GD.getGene(AGS_console.blast_parser.currentGeneMatchesIndexes[i].get(j));
                    int maxScore2 = MAX_SCORES[AGS_console.blast_parser.currentGeneMatchesIndexes[i].get(j)];
                    int score = 0;
                    try {
                        score = getScore(currentSequence, otherSequence, i, AGS_console.blast_parser.currentGeneMatchesIndexes[i].get(j));
                    }
                    catch (Exception e) {};
                    if (score > 0)
                    {
                        if (maxScore1 < maxScore2)
                            score = score*1000 / maxScore1;
                        else
                            score = score*1000 / maxScore2;

                        if (score>1000)
                            score = 1000;

                        if (score<0)
                            score = 0;

                        ALL_SCORES[i-lowerlimit].add(score);
                    }
                }
                
                if (i%1000==0)
                    System.out.println("Gene " + i + " alignments are calculated.");
   
            }
            // Write to file
            try
            {
                for (int i=0; i<ALL_SCORES.length; i++)
                {
                    if (ALL_SCORES[i].size()==0) // No hits found, add a self hit
                    {
                        String currentGene = AGS_console.blast_parser.currentGenes[lowerlimit+i];
                        int self_score = 1000;
                        output.write(currentGene + "\t" + self_score + "\t");                    
                    }
                    for (int j=0; j<ALL_SCORES[i].size(); j++)
                    {
                        output.write((AGS_console.GD.getGene(AGS_console.blast_parser.currentGeneMatchesIndexes[lowerlimit+i].get(j))).getGeneName2() + "\t" + ALL_SCORES[i].get(j) + "\t");
                    }
                    output.write(System.getProperty("line.separator"));
                }   
                output.close(); 
            } catch (Exception e){}
        }
        
    }
    /** Creates a new instance of CalculateAllScoreMatrices */
    public CalculateBlastFirstScoreMatrices(int processors, String outputName, boolean BLAST_ONLY,
                           boolean PROTEIN_ALIGNMENT, int ALIGNMENT_LIMIT, boolean LENGTH_NORMALIZED_SCORES)
                                    throws IOException, InvalidScoringMatrixException, IncompatibleScoringSchemeException
    {
        this.PROTEIN_ALIGNMENT = PROTEIN_ALIGNMENT;
        this.LENGTH_NORMALIZED_SCORES = LENGTH_NORMALIZED_SCORES;
        totalAlignments = 0;
        this.outputName = outputName;
        this.ALIGNMENT_LIMIT = ALIGNMENT_LIMIT;
        
        // Matrix file for protein inputs
        File blosum62Matrix = new File("blosum62.txt");
        FileReader matrix_file = new FileReader(blosum62Matrix);
        scoring = new ScoringMatrix (matrix_file);            

        FileWriter output = create_file(new File(outputName));
        if (output == null) throw new IOException("Could not open output file."); 

        if (BLAST_ONLY)
        {
            // Write header of file
            output.write(">");
            for (int i=0; i<AGS_console.GD.size(); i++)
            {
                output.write(AGS_console.GD.getGene(i).getGeneName2() +"\t");
            }
            output.write(System.getProperty("line.separator"));            
            
            for (int i=0; i<AGS_console.blast_parser.currentGenes.length; i++)
            {
                String currentGene = AGS_console.blast_parser.currentGenes[i];
                boolean selfHitFound = false;
                // No hits -- probably insignificant sequence (doesnt even score itself), add itself as a hit, give error message, deal with it later (i.e. remove from dataset)
                // For mammalian genomes sequences less than length 50 and all X chars were eliminated using tags in GeneDatabase.java so this error message should never happen
                if (AGS_console.blast_parser.currentGeneMatchesIndexes[i].size()==0)
                {
                    System.out.println("Warning: No matches found for gene " + currentGene);
                    int self_score = 1000;
                    output.write(currentGene + "\t" + self_score +"\t");
                    output.write(System.getProperty("line.separator"));
                }
                else // There is at least one hit (probably at least itself) for this gene
                {
                    // First find max scores
                    int maxScore = 0;
                    
                    for (int j=0; j<AGS_console.blast_parser.currentGeneMatchesScores[i].size(); j++)
                    {
                        int thisScore = AGS_console.blast_parser.currentGeneMatchesScores[i].get(j);
                        
                        if (maxScore < thisScore)
                        {
                            maxScore = thisScore;
                        }                    
                    }
                    
                    System.out.println("Number of matches for gene " + currentGene + " are " + AGS_console.blast_parser.currentGeneMatchesIndexes[i].size());
                    
                    for (int j=0; j<AGS_console.blast_parser.currentGeneMatchesIndexes[i].size(); j++)
                    {
                        
                        if (AGS_console.blast_parser.currentGeneMatchesIndexes[i].get(j).equals(i))
                            selfHitFound = true;
                        
                        output.write(AGS_console.GD.getGene(AGS_console.blast_parser.currentGeneMatchesIndexes[i].get(j)).getGeneName2() + "\t" +
                                    (int)(AGS_console.blast_parser.currentGeneMatchesScores[i].get(j)*1000/maxScore)+ "\t");
    
                    }
                    if (!selfHitFound)
                    {
                        System.out.println("Warning: No self hits found for gene (probably length is less than max allowed): " + currentGene);
                        int self_score = 1000;
                        output.write(currentGene + "\t" + self_score +"\t");
                    }
                    output.write(System.getProperty("line.separator"));
                }
            }            
        }
        else // Blast first then calculate alignments
        {
            
            // For faster alignments
            if (PROTEIN_ALIGNMENT)
                ALL_INT_ARRAYS = prepareSequenceIndex(AGS_console.GD);
            System.out.println("Sequence indexes calculated.");

            MAX_SCORES = new int[AGS_console.GD.size()];
            
            
            long time = -System.currentTimeMillis();
            System.out.println("Started running parallel");            
                            // Calculate max scores
            
            // Calculate max scores
            Runnable runners[] = new Runnable[processors];
            Thread th[] = new Thread[processors];                
            for (int j=0; j<processors; j++)
            {
                runners[j] = new CalculateBlastFirstScoreMatrices(j, processors, PROTEIN_ALIGNMENT, true);
                th[j] = new Thread(runners[j]);
                th[j].start();
            }
            for(int j=0;j<processors;j++){
                try{
                th[j].join();
                }
                catch(InterruptedException e){}
            }
    
            // Calculate regular scores
            for (int j=0; j<processors; j++)
            {
                runners[j] = new CalculateBlastFirstScoreMatrices(j, processors, PROTEIN_ALIGNMENT, false);
                th[j] = new Thread(runners[j]);
                th[j].start();
            }
            for(int j=0;j<processors;j++){
                try{
                th[j].join();
                }
                catch(InterruptedException e){}
            }            
            //merge(processors, fileName);
            
            time += System.currentTimeMillis();
            System.out.println (time + " msec"); 
            
            
            // Write header of file
            output.write(">");
            for (int i=0; i<AGS_console.GD.size(); i++)
            {
                output.write(AGS_console.GD.getGene(i).getGeneName2() +"\t");
            }
            output.write(System.getProperty("line.separator"));
            
            // Merge all scores here
            for (int i=0; i<processors; i++)
            {
                String inputName = outputName + "_" + i;
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

        }
        output.close();              
    }

    /*
     * Gets two sequences, aligns them and returns the resulting alignment string.
     * Set whether protein or nucleotide alignment.
     * 
     */
    public AlignmentString alignTwoSequences (Sequence seq1, Sequence seq2, int i, int j)
     						throws InvalidScoringMatrixException,
							IncompatibleScoringSchemeException
    {
        String S1, S2;
        S1 = seq1.getSequence();
        S2 = seq2.getSequence();
        
        NeedlemanWunsch nw1;
        if (PROTEIN_ALIGNMENT)
        {
            int[] index1 = ALL_INT_ARRAYS[i]; // Faster alingnments
            int[] index2 = ALL_INT_ARRAYS[j]; // Faster alignments
            nw1 = new NeedlemanWunsch(index1, index2, S1, S2, scoring);
            return nw1.computeFastPairwiseAlignment();
        }
        else
        {
            nw1 = new NeedlemanWunsch(S1, S2, PROTEIN_ALIGNMENT, scoring); // This is only for nucleotide alignments -- without array
            return nw1.computePairwiseAlignment();
        }
            
    }
    
    public int getScore(Sequence seq1, Sequence seq2, int i, int j)
     						throws InvalidScoringMatrixException,
							IncompatibleScoringSchemeException    
    {    
        AlignmentString as = alignTwoSequences(seq1, seq2, i, j);
        if (LENGTH_NORMALIZED_SCORES)
        {
            return (int)(as.normalized_score*1000);
        }
        else
        {
            return (int)(as.score*1000);
        }
    }
    private ArrayList<Integer> getAllGeneHitLocations(GeneDatabase ALL_GENES, ArrayList<String> hits, int limit)
    {
        ArrayList<Integer> intHits = new ArrayList<Integer>(limit);
        for (int i=0; i<hits.size();i++)
        {
            if (i>limit)
                break;
            String currentHit = hits.get(i);
            //intHits.add(ALL_GENES.getGeneLocationByName2(currentHit));
            intHits.add(ALL_GENES.getGeneLocationByPreComputedIndex(currentHit));
        }
        return intHits;
    }   
    
    
    private int[][] prepareSequenceIndex(GeneDatabase ALL_GENES)
    {
        int[][] ALL_INT_ARRAYS = new int[ALL_GENES.size()][];
        
        for (int i=0; i<ALL_INT_ARRAYS.length; i++)
        {
            Sequence seq = ALL_GENES.getGene(i);
            ALL_INT_ARRAYS[i] = scoring.getIndexesOfSequences(seq.sequence);
        }
        return ALL_INT_ARRAYS;
    }
    
    private static FileWriter create_file(File filename)
    {
            try {
                    return new FileWriter(filename);
            }
            catch (IOException e) {return null;}
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
}
