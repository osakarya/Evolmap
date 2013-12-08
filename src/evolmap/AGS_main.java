
package EvolMAP;
import java.util.*;
import java.io.*;
import javax.swing.JFrame;
import javax.swing.UIManager;
/**
 *
 * @author onur
 */
public class AGS_main 
{
    /**
     * Create the GUI and show it.  For thread safety,
     * this method should be invoked from the
     * event-dispatching thread.
     */
    public static void createAndShowGUI(NewickTree tree, GeneDatabase GD, boolean displayFullName) 
                throws IOException
    {
        try {
                UIManager.setLookAndFeel(
                    UIManager.getSystemLookAndFeelClassName());
            } catch (Exception e){}
        //Create and set up the window.
        JFrame frame = new JFrame("Gene Expansion Viewer");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        //Create and set up the content pane.
        ShowAncestor ancestor = new ShowAncestor(tree, GD, displayFullName);
        ancestor.setOpaque(true); //content panes must be opaque
        frame.setContentPane(ancestor);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }
    
    private static void outputDomains(GeneDatabase GD, AncestorDatabase ancestorDatabase, String DOMAIN_NAME, String ANCESTOR, int MAX_DOMAIN_LENGTH, int symbet_TOLERANCE)
                                    throws IOException
    {
        FileWriter output = create_file(new File(ANCESTOR + "." + DOMAIN_NAME + ".fa"));
        if (output == null) throw new IOException();
        
        FileWriter output2 = create_file(new File(ANCESTOR + "." + DOMAIN_NAME + ".groups"));
        if (output2 == null) throw new IOException();
        
        AGS ancestor_AGS = ancestorDatabase.tree.findNode(ANCESTOR).ags;
        int domainCount=0;
        int symbetCount=0;
        for (int i=0; i<ancestor_AGS.AGS.size(); i++)
        {
            OrthologGenes og = ancestor_AGS.AGS.get(i);
            if (GD.getGeneByPreComputedIndex(og.genes.get(0)).name.contains(DOMAIN_NAME))
            {
                if (!og.symbet)
                    continue;
                if (og.symbet_score < symbet_TOLERANCE)
                    continue;
                
                for (int j=0; j<og.genes.size(); j++)
                {
                   if (GD.getGeneByPreComputedIndex(og.genes.get(j)).sequence.length() > MAX_DOMAIN_LENGTH)
                       continue;
                   output.write(GD.getGeneByPreComputedIndex(og.genes.get(j)).toStringCompact());
                   output.write(System.getProperty("line.separator"));
                   
                   output2.write(og.genes.get(j) + "\t");
                }
                output2.write(System.getProperty("line.separator"));
                domainCount += og.genes.size();
                symbetCount++;
            }
        }
        System.out.println("Total " + domainCount + " many " + DOMAIN_NAME + " domains from " + symbetCount + " sym-bets from "+ ANCESTOR + " are written to file.");
        output.close();
        output2.close();
        
    }

    private static FileWriter create_file(File filename)
    {
            try {
                    return new FileWriter(filename);
            }
            catch (IOException e) {return null;}
    }    
}
