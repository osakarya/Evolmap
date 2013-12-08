package fasta;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import ancestor.AmphimedonHgtDB;
import ancestor.Gene;
import ancestor.AmphimedonOriginsDB;
import ancestor.HgtGene;


import utils.Utils;


public class cDNAFastaDB extends FastaDB {

	public cDNAFastaDB()
	{
		super();
	}
	
	public cDNAFastaDB(String filename) throws IOException
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
					cDNAFastaRecord fr = new cDNAFastaRecord(header, sequence);
					fr.percentGC = Utils.GCPercent(fr.sequence);
					fr.percentGC13 = Utils.GCPercent13(fr.sequence);
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

	public void annotateWithAncestors(AmphimedonOriginsDB odb) {
		
		for (FastaRecord fr : this.seqs)
		{
			Gene g = odb.genes.get(((cDNAFastaRecord)fr).gene);
			if (g!=null)
			{
				((cDNAFastaRecord)fr).species = g.species;
				((cDNAFastaRecord)fr).origin = g.originAncestor;
				((cDNAFastaRecord)fr).symmetryScore = g.score;
			}
		}
	}

	public void annotateWithHGT(AmphimedonHgtDB hgtDB) {
		for (FastaRecord fr : this.seqs)
		{
			HgtGene g = hgtDB.genes.get(((cDNAFastaRecord)fr).gene);
			if (g!=null)
			{
				((cDNAFastaRecord)fr).hgt = g.origins;
				((cDNAFastaRecord)fr).hgtScore = g.score;
			}
			else
			{
				((cDNAFastaRecord)fr).hgt = null;
				((cDNAFastaRecord)fr).hgtScore = 0;
			}
		}
	}
}
