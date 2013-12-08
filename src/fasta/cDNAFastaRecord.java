package fasta;

public class cDNAFastaRecord extends FastaRecord {

	public String contig;
	public String contigStart;
	public String contigEnd;
	public String contigFrame;
	public String gene;
	public String transcript;
	public String species;
	public String origin;
	public int symmetryScore;
	public String hgt;
	public int hgtScore; 
	
	public cDNAFastaRecord(String header, String sequence) {
		super(header, sequence);
		
		String[] items = header.split(" ");
		String contigInfo = items[2];
		String geneInfo = items[3];
		String transcriptInfo = items[4];
		
		String[] contigItems = contigInfo.split(":");
		this.contig = contigItems[2];
		this.contigStart = contigItems[3];
		this.contigEnd = contigItems[4];
		this.contigFrame = contigItems[5];
		
		String[] geneItems = geneInfo.split(":");
		this.gene = geneItems[1];
		
		String[] transcriptItems = transcriptInfo.split(":");
		this.transcript = transcriptItems[1] + ":" + transcriptItems[2];
	}
	
	public String toString()
	{
		double normGCBias = ((double)((int)(this.percentGC * 100))) / ((double) 100.0);
		double normGCBias13 = ((double)((int)(this.percentGC13 * 100))) / ((double) 100.0);
		return this.gene + "\t" + this.transcript + "\t" + this.contig + "\t" 
				+ this.contigStart + ":" + this.contigEnd + ":" + this.contigFrame + "\t"
				+ normGCBias + "\t"  + normGCBias13 + "\t"
				+ (this.hgt!=null) + "\t" + ((this.hgt==null) ? "NA" : this.hgt);
	}
	
	public static String getHeader() {
		return("#GeneName\tTranscript\tContigName\tContigLocation\t" +
				"GC\tGC13\t" +
				"HGT\tHGTgenes"
				);
	}

}
