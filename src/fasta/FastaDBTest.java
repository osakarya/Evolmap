package fasta;
import static org.junit.Assert.*;

import org.junit.Test;


public class FastaDBTest {

	@Test
	public void test() {
		FastaRecord fr1 = new FastaRecord("a", "b");
		fr1.percentGC = 1;
		
		FastaRecord fr2 = new FastaRecord("a", "b");
		fr2.percentGC = 2;
		
		FastaDB fd = new FastaDB();
		fd.add(fr1);
		fd.add(fr2);
		fd.calculateGCStats();
		
		assertEquals((Double)fd.meanGC, new Double(1.5));
		assertEquals((Double)fd.stdGC, new Double(0.5));
		
		
		FastaRecord fr3 = new FastaRecord("a", "b");
		fr3.percentGC = 3;
		fd.add(fr3);
		fd.calculateGCStats();
		
		assertEquals((Double)fd.meanGC, new Double(2));
		assertEquals((Double)fd.stdGC, new Double(Math.sqrt(2.0/3.0)));
		
	}

}
