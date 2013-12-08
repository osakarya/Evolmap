
package EvolMAP;
import java.io.*;
import java.util.*;

/**
 * @author onur
 */
public class CalculateAllScoreMatrices 
{
    
    private ScoringMatrix scoring; // The scoring matrix
    private boolean PROTEIN_ALIGNMENT; // If protein alignment, it is made faster by first preparing int vectors
    private boolean LENGTH_NORMALIZED_SCORES;
    private int[][] ALL_INT_ARRAYS; // Used for fast alignments
    
    /** Creates a new instance of CalculateAllScoreMatrices */
    public CalculateAllScoreMatrices(GeneDatabase ALL_GENES, String outputName, boolean PROTEIN_ALIGNMENT, boolean LENGTH_NORMALIZED_SCORES)
                                    throws IOException, InvalidScoringMatrixException, IncompatibleScoringSchemeException
    {
        this.PROTEIN_ALIGNMENT = PROTEIN_ALIGNMENT;
        this.LENGTH_NORMALIZED_SCORES = LENGTH_NORMALIZED_SCORES;
        
        // Matrix file for protein inputs
        File blosum62Matrix = new File("blosum62.txt");
        FileReader matrix_file = new FileReader(blosum62Matrix);
        scoring = new ScoringMatrix (matrix_file);            

        FileWriter output = create_file(new File(outputName));
        if (output == null) throw new IOException("Could not open output file.");
        
        // For faster alignments
        if (PROTEIN_ALIGNMENT)
            ALL_INT_ARRAYS = prepareSequenceIndex(ALL_GENES);
        
        int[][] ALL_SCORES = new int[ALL_GENES.size()][ALL_GENES.size()]; // This might go out of memory
        int[] MAX_SCORES = new int[ALL_GENES.size()];

        // Calculate max scores
        for (int i=0; i<MAX_SCORES.length; i++)
        {
            Sequence currentSequence = ALL_GENES.getGene(i);
            MAX_SCORES[i] = getScore(currentSequence, currentSequence, i, i);
        }

        // Calculate lower triangular scores matrix
        for (int i=0; i<ALL_SCORES.length; i++)
        {
            Sequence currentSequence = ALL_GENES.getGene(i);
            int maxScore1 = MAX_SCORES[i];
            for (int j=i; j<ALL_SCORES[i].length; j++)
            {
                Sequence otherSequence = ALL_GENES.getGene(j);
                int maxScore2 = MAX_SCORES[j];
                int score = getScore(currentSequence, otherSequence, i, j);

                if (maxScore1 < maxScore2)
                    ALL_SCORES[i][j] = score*1000 / maxScore1;
                else
                    ALL_SCORES[i][j] = score*1000 / maxScore2;

                if (ALL_SCORES[i][j]>1000)
                    ALL_SCORES[i][j] = 1000;
                
                if (ALL_SCORES[i][j]<0)
                    ALL_SCORES[i][j] = 0;
            }
            System.out.println("Gene " + i + " alignments are calculated.");
        }

        // Fill full matrix by transposing
        for (int i=ALL_SCORES.length-1; i>0; i--)
        {
            for (int j=i-1; j>=0; j--)
            {
                ALL_SCORES[i][j] = ALL_SCORES[j][i];
            }
        }

        // Write header of file
        output.write(">");
        for (int i=0; i<ALL_GENES.size(); i++)
        {
            Sequence currentSequence = ALL_GENES.getGene(i);
            output.write(ALL_GENES.getGene(i).getGeneName2() +"\t");
        }
        output.write(System.getProperty("line.separator"));

        for (int i=0; i<ALL_SCORES.length; i++)
        {
            for (int j=0; j<ALL_SCORES[i].length; j++)
            {
                output.write(ALL_SCORES[i][j] + "\t");
            }
            output.write(System.getProperty("line.separator"));
        }

        output.close();              
    }

    /*
     * Gets two sequences, aligns them and returns the resulting alignment string.
     * Set whether protein or nucleotide alignment.
     * 
     * BE CAREFUL CURRENTLY ALIGNS DOMAIN SEQUENCES ONLY, CHANGE THAT TO SEQUENCE FOR REGULAR ALIGNMENT
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
            int[] indexed1 = ALL_INT_ARRAYS[i]; // Faster alingnments
            int[] indexed2 = ALL_INT_ARRAYS[j]; // Faster alignments
            nw1 = new NeedlemanWunsch(indexed1, indexed2, S1, S2, scoring);
        }
        else
            nw1 = new NeedlemanWunsch(S1, S2, PROTEIN_ALIGNMENT, scoring); // This is only for nucleotide alignments -- without array
        
        if (PROTEIN_ALIGNMENT)
            return nw1.computeFastPairwiseAlignment();
        else    
            return nw1.computePairwiseAlignment();
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
}
