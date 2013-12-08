/*
 * AGS_console.java
 *
 * Created on August 13, 2007, 1:14 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package EvolMAP;
import java.util.*;
import java.io.*;
/**
 *
 * @author onur
 */
public class AGS_console {
    
    static GeneDatabase GD = null;
    static BLAST_parser blast_parser = null;
    static ScoreMatrix scoreMatrix = null; 
    
    public AGS_console(int processors, NewickTree tree, String tag, boolean protein, boolean Blastall, int alignments, boolean bit_scores, boolean stopAfterPass1,
                    String analysis_name, boolean read_database, boolean read_scores, boolean read_blast_scores, boolean readPass1, boolean view_ancestors,
                    int ortholog_threshold, int diverged_threshold, int diverged_std, boolean avg_of_paralogs, boolean sfa, ArrayList<String> subtree_names,
                    ArrayList<String> subtrees)
                    throws Exception
    {
        String database_name = analysis_name + ".gd" + tag;
        String sfa_name = analysis_name + tag;
        String SCORE_MATRIX_NAME = analysis_name + ".scores";
        String SPARSE_SCORE_MATRIX_NAME = analysis_name + ".sparse_scores";
        String ANCESTOR_FILE_NAME = analysis_name + ".ancestors_pass1";
        String ANCESTOR_FILE_NAME2 = analysis_name + ".ancestors_pass2";
        String GAIN_LOSS_FILE_NAME = analysis_name + ".ancestor_results";
        int intOrtholog_threshold = ortholog_threshold;
        int intDiverged_threshold = diverged_threshold;
       
        
        if (!view_ancestors && !read_database)
        {
            new CreateSequenceDatabaseFile(tree.getAllSpeciesNames(), tag, database_name);
            System.out.println("Gene database written to file " + database_name);
        }   
        this.GD = new GeneDatabase(database_name);
        
        
        if (!view_ancestors && !readPass1 && !read_scores && Blastall)   
        {
            if (!read_blast_scores)
            {
                FileWriter output2 = create_file(new File("temp." + SPARSE_SCORE_MATRIX_NAME));
                if (output2 == null) throw new IOException();         

                for (int i=0; i<GD.allGenes.length; i++)
                {
                    output2.write(GD.allGenes[i].toStringCompact() + System.getProperty("line.separator"));
                }
                output2.close();
            }
            this.blast_parser = new BLAST_parser(processors, ("temp." + SPARSE_SCORE_MATRIX_NAME), alignments, !read_blast_scores, false, protein);
            System.out.println("Blast parser finished reading hits for " + blast_parser.currentGenes.length + " many genes. There are " + blast_parser.alignmentCount + " many alignments.");

            new CalculateBlastFirstScoreMatrices(processors, SPARSE_SCORE_MATRIX_NAME, bit_scores, protein, alignments, true);
        }
        else if (!view_ancestors && !readPass1 && !read_scores)
        {
            new CalculateAllScoreMatrices(GD, SCORE_MATRIX_NAME, protein, true);
        }
        
        if (!view_ancestors && !readPass1 && Blastall)
            scoreMatrix = new SparseScoreMatrix(SPARSE_SCORE_MATRIX_NAME, alignments);
        else if (!view_ancestors && !readPass1)
            scoreMatrix = new RegularScoreMatrix(SCORE_MATRIX_NAME);
        else if (!view_ancestors && readPass1 && Blastall) // Read and close scoreMatrix to detect no scoring genes
        {
            scoreMatrix = new SparseScoreMatrix(SPARSE_SCORE_MATRIX_NAME, alignments);
            scoreMatrix = null;
        }
        else if (!view_ancestors && readPass1 && !Blastall) // Read and close scoreMatrix to detect no scoring genes
        {
            scoreMatrix = new RegularScoreMatrix(SCORE_MATRIX_NAME);
            scoreMatrix = null;
        }        
        
        if (subtree_names!=null)
        {
            String mainTreeName = analysis_name.substring(0, analysis_name.indexOf("_final"));
            AncestorDatabase ancestorDatabase = new AncestorDatabase(tree, GD, ANCESTOR_FILE_NAME, avg_of_paralogs, mainTreeName, subtree_names, subtrees);
        }
        else if (!view_ancestors && !readPass1)
        {
            AncestorDatabase ancestorDatabase = new AncestorDatabase(tree, GD, ANCESTOR_FILE_NAME, intOrtholog_threshold, Blastall, alignments, avg_of_paralogs, protein, sfa_name, sfa);
        }
        else
        {
            AncestorDatabase ancestorDatabase = new AncestorDatabase(tree, GD, ANCESTOR_FILE_NAME, avg_of_paralogs);
        }
             
        if (!view_ancestors && !stopAfterPass1)
            new AGS_analyzer(tree, GD, GAIN_LOSS_FILE_NAME, ANCESTOR_FILE_NAME2, true, true, intDiverged_threshold, diverged_std);
        else if (view_ancestors && !stopAfterPass1)
        {
             new AGS_analyzer(tree, GD, GAIN_LOSS_FILE_NAME, ANCESTOR_FILE_NAME2, false, false, intDiverged_threshold, diverged_std);
            System.out.println("Starting display...");
            AGS_main.createAndShowGUI(tree, GD, true);
        }    
    }
    public static void main(String args[])
            throws Exception
    {
        ArrayList<String> TRUE = new ArrayList<String>();
        ArrayList<String> FALSE = new ArrayList<String>();
        TRUE.add("true"); TRUE.add("yes"); TRUE.add("t"); TRUE.add("y");
        TRUE.add("True"); TRUE.add("Yes");
        TRUE.add("TRUE"); TRUE.add("YES"); TRUE.add("T"); TRUE.add("Y");
        
        FALSE.add("false"); FALSE.add("no"); FALSE.add("f"); FALSE.add("n");
        FALSE.add("False"); FALSE.add("No");
        FALSE.add("FALSE"); FALSE.add("NO"); FALSE.add("F"); FALSE.add("N");
        
        if (args.length==1)
        {
            String optionsFile = args[0];
            
            BufferedReader input = open_file(new File(optionsFile));
            if (input == null) 
            {
                System.out.println("Could not open input file " + optionsFile);
                System.exit(1);
            }
            String treeString = "";
            String firstTag = "";
            boolean protein = true;
            boolean Blastall = true;
            int alignments = 250;
            boolean bit_scores = true;
            String database_name = "";
            boolean read_database = false;
            boolean read_scores = false;
            boolean read_blast_scores = false;
            boolean read_ancestors = false;
            boolean view_ancestors = false;
            int ortholog_threshold = 250;
            int diverged_threshold = 250;
            int diverged_std = 3;
            boolean avg_of_paralogs = true;
            boolean sfa = false;
            String score_file_name;
            String blast_database_name;
            boolean domainAnalysis = false;
            ArrayList<String> domains = new ArrayList<String>();
            boolean recalculate_sfa_subtrees = false;
            boolean recalculate_final_tree = false;
            ArrayList<String> subtree_names = new ArrayList<String>();
            ArrayList<String> subtrees = new ArrayList<String>();
            int processors = 1;
            
            
            while (more_records(input))
            {
                String line = read_a_line(input);
                if (line.length()==0)
                    continue;
                StringTokenizer st = new StringTokenizer(line, "=");
                
                String command = removeEmptySpace(st.nextToken());  // Command
                String commandText = removeEmptySpace(st.nextToken()); // Command's text
                String moreText = "";
                if (st.hasMoreTokens())
                    moreText = removeEmptySpace(st.nextToken()); // Command's text
                
                if (commandText.indexOf("//")!=-1)
                    commandText = commandText.substring(0, commandText.indexOf("//"));
                
                if (command.equalsIgnoreCase("processors"))
                {
                    processors = Integer.parseInt(commandText);
                    if (processors < 1)
                    {
                        System.out.println("Number of processors error: " + processors);
                        System.exit(1);
                    }
                }
                else if (command.equalsIgnoreCase("subtree"))
                {
                    subtree_names.add(commandText);
                    subtrees.add(moreText);
                }
                else if (command.equalsIgnoreCase("recalculate_sfa_subtree"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        recalculate_sfa_subtrees = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        recalculate_sfa_subtrees = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }
                }
                else if (command.equalsIgnoreCase("recalculate_final_tree"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        recalculate_final_tree = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        recalculate_final_tree = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }                    
                }                
                else if (command.equalsIgnoreCase("retrieve_domains"))
                {
                    StringTokenizer st2 = new StringTokenizer(commandText, ",");
                    while (st2.hasMoreTokens())
                    {
                        domains.add(st2.nextToken());
                    }
                    if (domains.size()>0)
                    {
                        domainAnalysis = true;
                    }
                }
                else if (command.equalsIgnoreCase("tree")) // Required option
                {
                    treeString = commandText;
                }
                else if (command.equalsIgnoreCase("tag"))
                {
                    firstTag = commandText;
                }
                else if (command.equalsIgnoreCase("protein"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        protein = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        protein = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }                       
                }
                else if (command.equalsIgnoreCase("Blastall"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        Blastall = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        Blastall = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }                       
                }
                else if (command.equalsIgnoreCase("alignments"))
                {
                    alignments = Integer.parseInt(commandText);
                    
                }
                else if (command.equalsIgnoreCase("bit_scores"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        bit_scores = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        bit_scores = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }                       
                }
                else if (command.equalsIgnoreCase("database_name"))
                {
                    database_name = commandText;
                }
                else if (command.equalsIgnoreCase("read_database"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        read_database = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        read_database = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }                       
                }
                else if (command.equalsIgnoreCase("read_scores"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        read_scores = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        read_scores = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }                       
                }
                else if (command.equalsIgnoreCase("read_blast_scores"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        read_blast_scores = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        read_blast_scores = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }                       
                    
                }
                else if (command.equalsIgnoreCase("read_ancestors"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        read_ancestors = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        read_ancestors = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }                     
                }
                else if (command.equalsIgnoreCase("view_ancestors"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        view_ancestors = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        view_ancestors = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }                       
                    
                }                
                else if (command.equalsIgnoreCase("ortholog_threshold"))
                {
                    ortholog_threshold = Integer.parseInt(commandText); 
                }
                else if (command.equalsIgnoreCase("diverged_threshold"))
                {
                    diverged_threshold = Integer.parseInt(commandText);
                    
                }
                else if (command.equalsIgnoreCase("diverged_std"))
                {
                    diverged_std = Integer.parseInt(commandText);
                    
                }
                else if (command.equalsIgnoreCase("avg_of_paralogs"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        avg_of_paralogs = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        avg_of_paralogs = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }                       
                    
                }
                else if (command.equalsIgnoreCase("sfa"))
                {
                    if (TRUE.indexOf(commandText)!=-1)
                    {
                        sfa = true;
                    }
                    else if (FALSE.indexOf(commandText)!=-1)
                    {
                        sfa = false;
                    }
                    else
                    {
                        System.out.println("Invalid input " + line);
                        System.exit(1);
                    }                       
                }
                else
                {
                    System.out.println("Unknown command: " + command);
                }
            }
            
            // Control if every command is in order
            if (treeString.length()==0 || firstTag.length()==0 || database_name.length()==0)
            {
                System.out.println("One of tree, tag or database_name commands are missing.");
                System.exit(1);
            }
            
            if (domainAnalysis)
            {
                // We do multiple EvolMAP runs for each domain
                for (int i=0; i<domains.size(); i++)
                {
                    String currentDomain = domains.get(i);
                    NewickTree tree = new NewickTree(treeString);
                    // First we retrieve genes containing the domain for each species
                    String[] speciesNames = tree.getAllSpeciesNames();
                    if (!read_blast_scores && !read_scores)
                    {
                        for (int j=0; j<speciesNames.length; j++)
                        {
                            // This creates a fasta file genes containing selected domain for the selected species
                            System.out.print("Retrieving genes containing domain " + currentDomain + " from " + speciesNames[j] + firstTag + "...");
                            new PfamRetriever(speciesNames[j], firstTag, currentDomain);
                        }
                    }
                    String tag = "." + currentDomain + firstTag;
                    String analysis_name = database_name + "." + currentDomain;
                    System.out.println("Options parsed successfully. Running domain-based EvolMAP analysis for domain " + currentDomain + "...");
                    new AGS_console(processors, tree, tag, protein, Blastall, alignments, bit_scores, false,
                        analysis_name, read_database, read_scores, read_blast_scores, false, view_ancestors,
                         ortholog_threshold, diverged_threshold, diverged_std, avg_of_paralogs, sfa, null, null);                    
                }
            }
            else if (subtrees.size()>0) // SFA subtree based analysis
            {
                String realTreeString = treeString;
                
                // Do regular analysis for all subtrees first if this is true -- no pass2
                if (recalculate_sfa_subtrees==true)
                {
                    for (int i=0; i<subtrees.size(); i++)
                    {
                        NewickTree tree = new NewickTree(subtrees.get(i));
                        new AGS_console(processors, tree, firstTag, protein, Blastall, alignments, bit_scores, true,
                            subtree_names.get(i), read_database, read_scores, read_blast_scores, false, view_ancestors,
                             ortholog_threshold, diverged_threshold, diverged_std, avg_of_paralogs, sfa, null, null);  
                    }
                    
                }
                
                // Do regular anaylsis for final tree with SFA files used for analysis -- no pass2
                if (recalculate_final_tree==true)
                {
                    NewickTree tree = new NewickTree(treeString);
                    new AGS_console(processors, tree, firstTag, protein, Blastall, alignments, bit_scores, true,
                            database_name, read_database, read_scores, read_blast_scores, false, view_ancestors,
                             ortholog_threshold, diverged_threshold, diverged_std, avg_of_paralogs, sfa, null, null);                          
                }
                
                // Get complete tree string
                for (int i=0; i<subtree_names.size(); i++)
                {
                    if (realTreeString.indexOf(subtree_names.get(i))!=-1)
                    {
                        realTreeString = 
                                realTreeString.substring(0, realTreeString.indexOf(subtree_names.get(i))) +
                                subtrees.get(i) +
                                realTreeString.substring(realTreeString.indexOf(subtree_names.get(i))+subtree_names.get(i).length(), realTreeString.length());
                    }
                    else
                    {
                        System.out.println("Invalid input: " + subtree_names.get(i) + " not found in main tree.");
                        System.exit(1);
                    }
                }
                System.out.println("Final tree: " + realTreeString);
                
                // Run AGS_analyzer with the gene database and pass1 file that will calculate pass2 file and finish
                NewickTree tree = new NewickTree(realTreeString);
                
                // Replace subtrees with their nodal representation
                for (int i=0; i<subtrees.size(); i++)
                {
                    String thisTree = subtrees.get(i);
                    thisTree = thisTree.replaceAll(",", "_");
                    String noParenthesisName = "";
                    for (int j=0; j<thisTree.length(); j++)
                    {
                        if (thisTree.charAt(j)!=')' && thisTree.charAt(j)!='(')
                        {
                            noParenthesisName += thisTree.charAt(j);
                        }
                    }                     
                    subtrees.set(i, noParenthesisName);
                }
                
                new AGS_console(processors, tree, firstTag, protein, Blastall, alignments, bit_scores, false,
                    database_name + "_final", read_database, read_scores, read_blast_scores, true, view_ancestors,
                     ortholog_threshold, diverged_threshold, diverged_std, avg_of_paralogs, false, subtree_names, subtrees);
                
                
            }
            else // Regular single analysis
            {
                System.out.println("Options parsed successfully. Running regular EvolMAP analysis...");
                NewickTree tree = new NewickTree(treeString);
                new AGS_console(processors, tree, firstTag, protein, Blastall, alignments, bit_scores, false,
                    database_name, read_database, read_scores, read_blast_scores, read_ancestors, view_ancestors,
                     ortholog_threshold, diverged_threshold, diverged_std, avg_of_paralogs, sfa, null, null);
            }
        }
        else
        {
            System.out.println("Too many arguments. Only argument should be the options file.");
        }
                              
    }
    
    private static BufferedReader open_file(File filename)
    {
            try {
                    return new BufferedReader(new FileReader(filename));
            }
            catch (IOException e) {return null;}
    }

    private static boolean more_records(Reader r)
    {
            try {
                    return r.ready();
            }
            catch (IOException e) {return false;}
    }

    private static String read_a_line(BufferedReader r) 
    {
            try {
                    return r.readLine();
            }
            catch (IOException e) {return null;}
    }
    private static FileWriter create_file(File filename)
    {
            try {
                    return new FileWriter(filename);
            }catch (Exception exception) {System.out.println(exception.getStackTrace());
            return null;}
    }    
    private static String removeEmptySpace(String s)
    {
        return s.replaceAll(" ", "");
    }    
}
