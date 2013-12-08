
package EvolMAP;

import java.util.*;
import java.io.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.tree.*;
import javax.swing.event.*;

/**
 *
 * @author onur
 */
public class ShowAncestor extends JPanel 
                             implements ActionListener, TreeSelectionListener
{
    private DynamicTree treePanel;
    private JEditorPane genePane;
    private JTextField searchBox;
    private int newNodeSuffix = 1;
    private static String REMOVE_COMMAND = "remove";
    private static String EXPORT_COMMAND = "export";
    private static String SEARCH_COMMAND = "search";
    
    NewickTree tree;
    GeneDatabase GD;
    boolean DISPLAY_FULL_NAME;
    
    /** Creates a new instance of SpeciesTree */
    public ShowAncestor(NewickTree tree, GeneDatabase GD, boolean displayFullName)
            throws IOException
    {
        this.tree = tree;
        this.GD = GD;
        this.DISPLAY_FULL_NAME = displayFullName;
        
        //Create the components.
        treePanel = new DynamicTree();
        treePanel.tree.addTreeSelectionListener(this);   
        
        for (int i=0; i<tree.allNodes.size(); i++)
        {
            Node node = tree.getNode(i);
            AGS ancestor = node.getAGS();
            
            if (!node.species_node)
            {
                DefaultMutableTreeNode ANC = treePanel.addObject(null, node.name + " (SYM-BETS: " + node.getAGS().symbetCount() + ") " + 
                                                                        "(PRESENT-LOCI: " + node.getAGS().presentCount() + ")" );
                int presentCount = 0;
                for (int j=0; j<ancestor.size(); j++)
                {
                    OrthologGenes og = ancestor.getAGS().get(j);
                    if (og.present) 
                        presentCount++;
                    if (!og.divergedParalog)
                        addNewGene(presentCount, treePanel, ANC, og); 
                }
            }
        }
        
        //Lay everything out.
        genePane = new JEditorPane();
        genePane.setEditable(false);
        treePanel.setPreferredSize(new Dimension(1000, 700));
        JScrollPane geneView = new JScrollPane(genePane);
        geneView.setPreferredSize(new Dimension(1000, 300));
        JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
        splitPane.setTopComponent(treePanel);
        splitPane.setBottomComponent(geneView);
        add(splitPane, BorderLayout.CENTER);
        JPanel panel = new JPanel(new GridLayout(0,1));  
        
        searchBox = new JTextField();
        panel.add(searchBox);
        
        JButton searchButton = new JButton("Search in Ancestor");
        searchButton.setActionCommand(SEARCH_COMMAND);
        searchButton.addActionListener(this);  
        panel.add(searchButton);
        
        JButton removeButton = new JButton("Remove");
        removeButton.setActionCommand(REMOVE_COMMAND);
        removeButton.addActionListener(this);  
        panel.add(removeButton);
        
        JButton exportButton = new JButton("Export selected genes");
        exportButton.setActionCommand(EXPORT_COMMAND);
        exportButton.addActionListener(this);    
        panel.add(exportButton);
        
        add(panel, BorderLayout.LINE_END);
    }
   
    public void addNewGene(int geneCount, DynamicTree treePanel, DefaultMutableTreeNode Split_Ancestor, OrthologGenes og) 
    {
        String geneNo = "";
        String presentText = "";
        String name = "";
        DefaultMutableTreeNode Gene_Ancestor;
        
        if (og.symbet) 
        {
            presentText += " (SYM-BET with score " + og.symbet_score + ")";
            geneNo += geneCount;
        }
        else if (og.present) 
        {
            presentText += " (PRESENT)";
            geneNo += geneCount;
            
        }
        else 
        {
            presentText += " (SINGULAR with unknown source)";
            geneNo += "-";
        }
        
        if (DISPLAY_FULL_NAME)
            name = GD.getGeneByPreComputedIndex(og.getGenes().get(0)).name;
        else
            name = og.getGenes().get(0);
        
        int total_sequences = og.numberOfSequences();
        if (og.divergedInparalogs!=null)
        {
            total_sequences += og.divergedInparalogs.size();
        }
        
        Gene_Ancestor = treePanel.addObject(Split_Ancestor, "( " + geneNo + " ) " + "( " + name + " ) " +
                                                     "    (" + total_sequences + " genes)" + presentText);     
        for (int i=0; i<og.numberOfSequences(); i++)
        {
            Sequence sequence = GD.getGeneByPreComputedIndex(og.getGenes().get(i));
            String speciesName = sequence.species;
            sequence.setSpecies(speciesName);
            if (DISPLAY_FULL_NAME)
                sequence.displayFullName();            
            treePanel.addObject(Gene_Ancestor, sequence);
        }
        if (og.divergedInparalogs!=null)
        {
            for (int i=0; i<og.divergedInparalogs.size(); i++)
            {
                String nextDiverged = og.divergedInparalogs.get(i);
                Sequence sequence = new Sequence(GD.getGene(GD.getGeneLocationByPreComputedIndex(nextDiverged)));
                String speciesName = sequence.species;
                sequence.divergedInparalog = true;
                if (DISPLAY_FULL_NAME)
                    sequence.displayFullName();            
                treePanel.addObject(Gene_Ancestor, sequence);   
            }
        }       
    }

    
    public void actionPerformed(ActionEvent e) {
        String command = e.getActionCommand();
        
        if (SEARCH_COMMAND.equals(command))
        {
            String searchText = this.searchBox.getText();
            if (searchText.length()==0)
                return;
            TreePath currentSelection = treePanel.tree.getSelectionPath();
            if (currentSelection != null) 
            {
                DefaultMutableTreeNode currentNode = (DefaultMutableTreeNode)
                             (currentSelection.getLastPathComponent());
                MutableTreeNode parent = (MutableTreeNode)(currentNode.getParent());
                if (parent!= treePanel.rootNode) 
                {
                    JOptionPane.showMessageDialog(null, "Please select an ancestral node to search for", "Selection error", JOptionPane.ERROR_MESSAGE);
                }
                else
                {
                    boolean found = false;

                    for (int i=0; i<currentNode.getChildCount(); i++)
                    {
                        if (!found && !currentNode.getChildAt(i).isLeaf())
                        {
                            DefaultMutableTreeNode currentChild = (DefaultMutableTreeNode) currentNode.getChildAt(i);
                            for (int j=0; j<currentChild.getChildCount(); j++)
                            {
                                if (currentChild.getChildAt(j).isLeaf())
                                {
                                    Sequence retrievedSequence = (Sequence)((DefaultMutableTreeNode)(currentChild.getChildAt(j))).getUserObject();
                                    if (retrievedSequence.header.contains(searchText))
                                    {
                                        treePanel.tree.setSelectionPath(new TreePath(((DefaultMutableTreeNode)currentChild.getChildAt(j)).getPath()));
                                        found = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if (!found)
                    {
                        JOptionPane.showMessageDialog(null, "Search sequence " + searchText + " not found in headers.", "Selection error", JOptionPane.ERROR_MESSAGE);
                    } 
                }
            }
        }
        else if (REMOVE_COMMAND.equals(command)) {
            //Remove button clicked
            treePanel.removeCurrentNode();
        }
        else if (EXPORT_COMMAND.equals(command))
        {
            TreePath currentSelection = treePanel.tree.getSelectionPath();
            if (currentSelection != null) {
                DefaultMutableTreeNode currentNode = (DefaultMutableTreeNode)
                             (currentSelection.getLastPathComponent());
                MutableTreeNode parent = (MutableTreeNode)(currentNode.getParent());
                if (parent != null && parent!= treePanel.rootNode && !currentNode.isLeaf()) 
                {
                    String fileName = "anc_" + parent.getParent().getIndex(parent) + "_" + (parent.getIndex(currentNode)+1) + ".fa";
                    FileWriter output = create_file(new File(fileName));
                    if (output == null) fatal("could not open output file");  
                   
                    for (int i=0; i<currentNode.getChildCount(); i++)
                    {
                        if (currentNode.getChildAt(i).isLeaf())
                        {
                            Sequence retrievedSequence = (Sequence)((DefaultMutableTreeNode)(currentNode.getChildAt(i))).getUserObject();
                            try {
                                output.write(retrievedSequence.toStringWhole() + System.getProperty("line.separator"));
                            } catch(IOException ex) {}
                        }
                    }
                    try {
                        output.close();
                    } catch (IOException ex){}
                    JOptionPane.showMessageDialog(null, "" + currentNode.getLeafCount() + " descendant genes written to file " + fileName , "Fasta exported!", JOptionPane.INFORMATION_MESSAGE);
                }
                else
                {
                    JOptionPane.showMessageDialog(null, "Please select a (ancestral) gene", "Selection error", JOptionPane.ERROR_MESSAGE);
                }
            }
        }
    }
    /** Required by TreeSelectionListener interface. */
    public void valueChanged(TreeSelectionEvent e) 
    {
        DefaultMutableTreeNode node = (DefaultMutableTreeNode)
                           (treePanel.tree).getLastSelectedPathComponent();

        if (node == null) return;

        Object nodeInfo = node.getUserObject();
        if (node.isLeaf()) {
            Sequence gene = (Sequence)nodeInfo;
            try{displayGene(gene.toStringWhole());}
            catch (IOException et) {}
            }
    }    
    
    private void displayGene(String gene) 
                    throws IOException
    {
        if (gene != null) 
        {
            genePane.setText(gene);
        } 
        else 
        {
            genePane.setText("File Not Found");
        }
    }    
    
    private static FileWriter create_file(File filename)
    {
            try {
                    return new FileWriter(filename);
            }
            catch (IOException e) {return null;}
    }

    private static void fatal(String message)
    {
            System.out.println(message + "\nProgram terminating");
            System.exit(1);
    }       
}
