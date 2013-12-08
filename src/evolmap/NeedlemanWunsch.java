
package EvolMAP;

public class NeedlemanWunsch
{
	
	final static int MATCH_SCORE = 5;
	final static int MISMATCH_SCORE = -4;
	final static int GAP_START = -1;
	final static int GAP = 0;
        final static int GAP_NUCLEOTIDE = -1;
	
	final String GAP_CHARACTER = "-";
	final String GAP_TAG = "*";
	final String MATCH_TAG = "|";
	final String MISMATCH_TAG = " ";
	
	protected String seq1;
	protected String seq2;
	
        protected int[] seq_int_1;
        protected int[] seq_int_2;
        
        boolean FAST_ALIGNMENT = false; // Change this to true if called with the right constructor
	boolean protein_alignment;
	ScoringScheme scoring;

	/**
	 * The dynamic programming matrix. Each position (i, j) represents the best score
	 * between the firsts i characters of seq1 and j characters of seq2.
	 */
	protected int[][] matrix;
	
	// Shows if a gap was already started at the given matrix location
	protected boolean[][][] gapMatrix;
	
	/**
	 * Constructor 1
	 */
	public NeedlemanWunsch()
	{
		this.seq1 = null;
		this.seq2 = null;
		protein_alignment = false;
	}
	
	/**
	 * Constructor 2
	 */	
	public NeedlemanWunsch(String seq1, String seq2)
	{
		this.seq1 = seq1;
		this.seq2 = seq2;
		protein_alignment = false;
	}	
	
	/**
	 * Constructor 3
	 */	
	public NeedlemanWunsch(String seq1, String seq2, boolean protein_alignment, ScoringScheme scoring)
	{
		this.seq1 = seq1;
		this.seq2 = seq2;
		this.protein_alignment = protein_alignment;
                this.setScoringScheme(scoring);
	}	
        
        /**
         * Constructor 4 -- only for protein alignments so don't worry about nucleotide scoring if FAST_ALIGNMENT is true
         */
        public NeedlemanWunsch(int[] seq_int_1, int[] seq_int_2, String seq1, String seq2, ScoringScheme scoring)
        {
                this.FAST_ALIGNMENT = true;
		this.seq_int_1 = seq_int_1;
                this.seq1 = seq1;
		this.seq_int_2 = seq_int_2;
                this.seq2 = seq2;
		this.protein_alignment = true;
                this.setScoringScheme(scoring);            
        }
	
	/**
	 * Load the sequences if they are created with the first Constructor
	 */
	public void loadSequences (String input1, String input2)
	{
		this.seq1 = input1;
		this.seq2 = input2;
	}
	
	/**
	 * Call this method to directly get the alignment
	 * Not to be used in this project.
         * This is not used by fast alignments...
	 */
	public AlignmentString computePairwiseAlignment ()
		throws IncompatibleScoringSchemeException
	{
		// compute the matrix
		computeMatrix();

		// build and return an optimal global alignment
		AlignmentString alignment = buildOptimalAlignment(seq2.length(), seq1.length());

		return alignment;
	}
	
        public AlignmentString computeFastPairwiseAlignment()
                throws IncompatibleScoringSchemeException
        {
            computeMatrix();
            return (new AlignmentString(this.getMaxScore(seq2.length(), seq1.length()), this.getSingleGapAlignmentLength(seq2.length(), seq1.length())));
            
        }
	public void unloadSequences()
	{
		this.seq1 = null;
		this.seq2 = null;
		this.matrix = null;
		this.gapMatrix = null;
	}

	/**
	 * Computes the dynamic programming matrix.
	 */
	protected void computeMatrix ()
			throws IncompatibleScoringSchemeException
	{
		int r, c, rows, cols, gap1, gap2, sub;
                
                if (FAST_ALIGNMENT)
                {
                    rows = seq_int_2.length+1;
                    cols = seq_int_1.length+1;
                }
                else
                {
                    rows = seq2.length()+1;
                    cols = seq1.length()+1;
                }
                
		matrix = new int [rows][cols];
		
		// gapMatrix[0] holds if there is an ongoing gap for seq1
		// gapMatrix[1] holds if there is an ongoing gap for seq2
		gapMatrix = new boolean [2][rows][cols];

		// initiate first row
		matrix[0][0] = 0;
		gapMatrix[1][0][0] = false;
		for (c = 1; c < cols; c++)
		{
			matrix[0][c] = matrix[0][c-1] + gapPenalty(1, 0, c-1);
			gapMatrix[1][0][c] = true;
		}
		// calculates the similarity matrix (row-wise)
		for (r = 1; r < rows; r++)
		{
			// initiate first column
			matrix[r][0] = matrix[r-1][0] + gapPenalty(0, r-1, 0);
			gapMatrix[0][r][0] = true;
			for (c = 1; c < cols; c++)
			{
                                if (FAST_ALIGNMENT)
                                {
                                    sub = matrix[r-1][c-1] + scoreSubstitution(seq_int_1[c-1], seq_int_2[r-1]);
                                }
                                else
                                {
                                    sub = matrix[r-1][c-1] + scoreSubstitution(seq1.charAt(c-1), seq2.charAt(r-1));
                                }
                                gap1 = matrix[r][c-1] + gapPenalty(1, r, c-1);
				gap2 = matrix[r-1][c] + gapPenalty(0, r-1, c);

				// choose the greatest
				matrix[r][c] = max (sub, gap1, gap2);
			
			
				if (sub == max(sub, gap1, gap2))
				{
					gapMatrix[1][r][c] = false;
					gapMatrix[0][r][c] = false;
				}		
				
				else if (gap1 == max(sub,gap1, gap2))
				{
					gapMatrix[1][r][c] = true;
					gapMatrix[0][r][c] = false;
				}
				
				else if (gap2 == max(sub, gap1, gap2))
				{
					gapMatrix[1][r][c] = false;
					gapMatrix[0][r][c] = true;
				}
			}
		}
	}
	
	/**
	 * Helper method to compute the the greater of three values.
	 */
	protected int max (int v1, int v2, int v3)
	{
		return (v1 >= v2) ? ((v1 >= v3)? v1 : v3) : ((v2 >= v3)? v2 : v3);
	}
	
	/**
	 * Helper method to invoke the scoreSubstitution of the scoring
	 * scheme set to this algorithm.
	 *
	 * a and b must be in the same case.
	 */
	protected int scoreSubstitution (char a, char b)
		throws IncompatibleScoringSchemeException
	{
		if (!protein_alignment)
		{
			if (a==b)
				return this.MATCH_SCORE;
			else
				return this.MISMATCH_SCORE;
		}
		else
		{
			return scoring.scoreSubstitution (a, b);                     
		}
	}
        protected int scoreSubstitution(int r, int c)
            throws IncompatibleScoringSchemeException
        {
            return scoring.scoreSubstitution (r, c);  
        }
	
	/* Set the scoring scheme */
	public void setScoringScheme (ScoringScheme scoring)
	{

		this.scoring = scoring;
		
	}
	/**
	 * Helper method to get the gapPenalty
	 */
	protected int gapPenalty (int k, int i, int j)
	{
		if (gapMatrix[k][i][j] == false)
			return this.GAP_START;
		else if (this.protein_alignment)
			return this.GAP;
                else
                        return this.GAP_NUCLEOTIDE;
	}
	
	protected int getMaxScore(int r, int c)
        {
            return matrix[r][c];
        }
        protected int getSingleGapAlignmentLength(int r, int c)
                    throws IncompatibleScoringSchemeException
        {
            int sub = 0;
            int total_length = 0;
            boolean in_gap = false;
            
            while ((r > 0) || (c > 0))
            {
                if (c > 0)
                {
                        if (matrix[r][c] == (matrix[r][c-1] + gapPenalty(1, r, c-1)))
                        {
                                c = c - 1;
                                if (!in_gap)
                                    total_length++;
                                in_gap = true;
                                // skip to the next iteration
                                continue;
                        }
                }

                if ((r > 0) && (c > 0))
                {
                    if (protein_alignment)
                        sub = scoreSubstitution(seq_int_1[c-1], seq_int_2[r-1]);
                    else
                        sub = scoreSubstitution(seq1.charAt(c-1), seq2.charAt(r-1));

                    if (matrix[r][c] == matrix[r-1][c-1] + sub)
                    {
                        r = r - 1; c = c - 1;
                        total_length++;
                        in_gap = false;
                        // skip to the next iteration
                        continue;
                    }
                }

                // must be a deletion
                r = r - 1;
                if (!in_gap)
                    total_length++;
                in_gap = true;
            }    
            return total_length;
        }
	/**
	 * Builds an optimal global alignment between the loaded sequences. Before it is
	 * executed, the dynamic programming matrix must already have been computed.
	 *
	 */
	protected AlignmentString buildOptimalAlignment (int r, int c)
		throws IncompatibleScoringSchemeException
	{
		
		StringBuffer	gapped_seq1, score_tag_line, gapped_seq2;
		int				sub, max_score;

		gapped_seq1 	= new StringBuffer();
		score_tag_line	= new StringBuffer();
		gapped_seq2 	= new StringBuffer();
		
		/*
		 *  Find the max score anywhere in the matrix. In a copy event, it is
		 *  expected not to happen in the beginning.
		 *
		 *  IF PROTEIN ALIGNMENT Find the max_score according to global alignment.
		 */ 
		max_score = matrix[r][c];
		
		/* SWITCHED OFF
		if (!protein_alignment)
			for (int i=1; i<x; i++)
				for (int j=1; j<y; j++)
				{
					if (matrix[i][j] > max_score)
					{
						max_score = matrix[i][j];
						r=i;
						c=j;
					}
				}
		
		*/

		while ((r > 0) || (c > 0))
		{
			if (c > 0)
			{
				if (matrix[r][c] == (matrix[r][c-1] + gapPenalty(1, r, c-1)))
				{
					// gap at second sequence
					gapped_seq2.insert (0, GAP_CHARACTER);
					score_tag_line.insert (0, GAP_TAG);
					gapped_seq1.insert (0, seq1.charAt(c-1));
					c = c - 1;

					// skip to the next iteration
					continue;
				}
			}
			
			if ((r > 0) && (c > 0))
			{
				sub = scoreSubstitution(seq1.charAt(c-1), seq2.charAt(r-1));

				if (matrix[r][c] == matrix[r-1][c-1] + sub)
				{
					// substitution was used
					gapped_seq1.insert (0, seq1.charAt(c-1));
					if (seq1.charAt(c-1) == seq2.charAt(r-1))
							score_tag_line.insert (0, MATCH_TAG);
					else
						score_tag_line.insert (0, MISMATCH_TAG);
					gapped_seq2.insert (0, seq2.charAt(r-1));
					r = r - 1; c = c - 1;

					// skip to the next iteration
					continue;
				}
			}
			
			// must be a deletion
			gapped_seq1.insert (0, GAP_CHARACTER);
			score_tag_line.insert (0, GAP_TAG);
			gapped_seq2.insert (0, seq2.charAt(r-1));
			r = r - 1;
		}
		
		return new AlignmentString (gapped_seq1.toString(), score_tag_line.toString(),
										gapped_seq2.toString(), max_score);
	}

    public void setScoring(ScoringScheme scoring) {
        this.scoring = scoring;
    }
		
}