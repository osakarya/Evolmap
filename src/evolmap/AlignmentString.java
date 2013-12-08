
package EvolMAP;

public class AlignmentString

{
	String line1;
	String line2;
	String tag_line;
	int alignmentLength;
        int gapless_alignment_length; // This counts gaps as single characters and returns the length of the alignment without any gaps
        int score;
        double normalized_score;

	// These vectors holds the real positions of the sequences when you are
	// on a particular position of the alignment (Simply skips gaps)
	int[] line1RealPos;
	int[] line2RealPos;
	
	public AlignmentString(String line1, String tag_line, String line2, int score)
	{
		this.line1 = line1;
		this.line2 = line2;
		this.tag_line = tag_line;
		this.alignmentLength = line1.length();
		if (line1.length()<line2.length())
			this.alignmentLength = line2.length();
                
                this.score = score;
                this.setGaplessAlignmentLength();
                this.setNormalizedScore();
		      
	}
        public AlignmentString(int score, int gapless_alignment_length)
        {
            this.score = score;
            this.gapless_alignment_length = gapless_alignment_length;
            this.setNormalizedScore();
        }
	
	public int length()
	{
		return this.alignmentLength;
	}
        
        public int getGaplessLength()
        {
            return gapless_alignment_length;
        }
        
        public void setGaplessAlignmentLength()
        {
            int start_length = alignmentLength;
            boolean currently_in_gap = false;
            for (int i=0; i<tag_line.length(); i++)
            {
                if (!currently_in_gap)
                {
                    if (tag_line.charAt(i)=='*')
                    {
                        currently_in_gap = true;
                    }
                    else
                    {
                        currently_in_gap = false;
                    }
                }
                else
                {
                    if (tag_line.charAt(i)=='*')
                    {
                        start_length--;
                    }
                    else
                    {
                        currently_in_gap = false;
                    }
                }
            }
            this.gapless_alignment_length = start_length;
        }
        
        public void setNormalizedScore()
        {
            this.normalized_score = (double)((double)score/(double)gapless_alignment_length);
        }		

}