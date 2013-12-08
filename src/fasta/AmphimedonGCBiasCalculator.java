package fasta;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import ancestor.AmphimedonHgtDB;
import ancestor.AmphimedonOriginsDB;



public class AmphimedonGCBiasCalculator 
{
	
	public static void main(String args[]) throws IOException
	{
		String filename = "Amphimedon_queenslandica.Aqu1.14.cdna.all.fa";
		String filename3 = "all_genomes.ancestors_pass1.HGS.out";
		String outFilename = "amphimedon.gcstats.txt";
		BufferedWriter outFile = new BufferedWriter(new FileWriter(outFilename));
		
		AmphimedonHgtDB hgtDB = new AmphimedonHgtDB(filename3);
		//hgtDB.printContents();
		
		cDNAFastaDB fd = new cDNAFastaDB(filename);
		fd.annotateWithHGT(hgtDB);
		fd.calculateGCStats();
	
		// Write header info
		double normMeanGC = ((double)((int)(fd.meanGC * 100))) / ((double) 100.0);
		double normStdGC = ((double)((int)(fd.stdGC * 100))) / ((double) 100.0);
		double normMeanGC12 = ((double)((int)(fd.meanGC13 * 100))) / ((double) 100.0);
		double normStdGC13 = ((double)((int)(fd.stdGC13 * 100))) / ((double) 100.0);
		outFile.write("#size: " + fd.size() +"\t" 
				+ "GC-mean-%: " + normMeanGC 
				+ "\tGC-std-%: " + normStdGC
				+ "\tGC13-mean-%: " + normMeanGC12 
				+ "\tGC13-std-%: " + normStdGC13 + "\n");
		
		// Write header
		outFile.write(cDNAFastaRecord.getHeader() +"\n");
		
		// Write file
		int different = 0;
		int different13 = 0;
		for (FastaRecord fr : fd.seqs)
		{
			if ((fr.percentGC > (fd.meanGC+(2.0)*fd.stdGC))
					|| (fr.percentGC < (fd.meanGC-(2.0)*fd.stdGC)))
			{
				fr.GCbiasSignificance = true;
				different++;
			}
			else
			{
				fr.GCbiasSignificance = false;
			}
			
			if ((fr.percentGC13 > (fd.meanGC13+(2.0)*fd.stdGC13))
					|| (fr.percentGC13 < (fd.meanGC13-(2.0)*fd.stdGC13)))
			{
				fr.GC13biasSignificance = true;
				different13++;
			}
			else
			{
				fr.GC13biasSignificance = false;
			}
			
			outFile.write(((cDNAFastaRecord) fr).toString() + "\n");
		}
		System.out.println("Number of records with significant GC bias (2std): " + different);
		System.out.println("Number of records with significant 1st and 3rd codon GC bias (2std): " + different13);
		
		outFile.close();
	}

}
