package utils;
import static org.junit.Assert.*;

import org.junit.Test;


public class UtilsTests {

	@Test
	public void test() {
		assertEquals((Double)Utils.GCPercent("AG"), new Double(50));
		assertEquals((Double)Utils.GCPercent("AA"), new Double(0));
		assertEquals((Double)Utils.GCPercent("TT"), new Double(0));
		assertEquals((Double)Utils.GCPercent("CG"), new Double(100));
		assertEquals((Double)Utils.GCPercent("GG"), new Double(100));
		assertEquals((Double)Utils.GCPercent("cc"), new Double(100));
		assertEquals((Double)Utils.GCPercent("AGGG"), new Double(75));
		assertEquals((Double)Utils.GCPercent("AGCG"), new Double(75));
		assertEquals((Double)Utils.GCPercent("AAAAG"), new Double(20));
		assertEquals((Double)Utils.GCPercent("AAAAAAATAC"), new Double(10));
		
		assertEquals((Double)Utils.GCPercent13("AG"), new Double(0));
		assertEquals((Double)Utils.GCPercent13("AA"), new Double(0));
		assertEquals((Double)Utils.GCPercent13("TT"), new Double(0));
		assertEquals((Double)Utils.GCPercent13("CG"), new Double(100));
		assertEquals((Double)Utils.GCPercent13("GG"), new Double(100));
		assertEquals((Double)Utils.GCPercent13("cc"), new Double(100));
		assertEquals((Double)Utils.GCPercent13("AGG"), new Double(50));
		assertEquals((Double)Utils.GCPercent13("AGCGAG"), new Double(75));
		assertEquals((Double)Utils.GCPercent13("AAAAG"), new Double(0));
		assertEquals((Double)Utils.GCPercent13("AAAAAAATACGGAAA"), new Double(20));
		
	}

}
