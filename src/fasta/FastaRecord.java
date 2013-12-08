package fasta;

public class FastaRecord {

	String header;
	String sequence;
	public double percentGC;
	public double percentGC13;
	public boolean GCbiasSignificance;
	public boolean GC13biasSignificance;
	
	public FastaRecord(String header, String sequence)
	{
		this.header = header;
		this.sequence = sequence;
		percentGC = 50;
		percentGC13 = 50;
	}
}
