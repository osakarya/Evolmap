package fasta;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import utils.Utils;


public class FastaDB {

	public ArrayList<FastaRecord> seqs;
	public double meanGC;
	public double meanGC13;
	public double stdGC;
	public double stdGC13;
	
	public FastaDB()
	{
		this.seqs = new ArrayList<FastaRecord>();
	}
	
	public FastaDB(String filename) throws IOException
	{
		this.seqs = new ArrayList<FastaRecord>();
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		
		String line = "";
		String header = "";
		String sequence = "";
		while (reader.ready())
		{
			line = reader.readLine();
			if (line.startsWith(">"))
			{
				if (header.equals(""))
				{
					header = line;
				}
				else
				{
					FastaRecord fr = new FastaRecord(header, sequence);
					fr.percentGC = Utils.GCPercent(fr.sequence);
					seqs.add(fr);
					header = line;
					sequence = "";
				}
			}
			else
			{
				sequence += line;
			}
		}
	}
	public void add(FastaRecord fr)
	{
		this.seqs.add(fr);
	}
	
	public void calculateGCStats()
	{
		double sum = 0;
		double sum13 = 0;
		double n = 0;
		for (FastaRecord fr : seqs)
		{
			sum += fr.percentGC;
			sum13 += fr.percentGC13;
			n+=1;
		}
		this.meanGC = sum / n;
		this.meanGC13 = sum13 / n;
		
		double diffsqsum = 0;
		double diffsqsum13 = 0;
		for (FastaRecord fr: seqs)
		{
			diffsqsum += Math.pow((fr.percentGC - meanGC),2); 
			diffsqsum13 += Math.pow((fr.percentGC13 - meanGC13),2); 
		}
		
		this.stdGC = Math.sqrt(diffsqsum/n);
		this.stdGC13 = Math.sqrt(diffsqsum13/n);
	}

	public int size() {
		return seqs.size();
	}
}
