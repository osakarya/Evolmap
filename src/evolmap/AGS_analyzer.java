/*
 * AGS_analyzer.java
 *
 * Created on May 24, 2007, 4:54 PM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package EvolMAP;
import java.util.*;
import java.io.*;


public class AGS_analyzer 
{
    public NewickTree tree;
    private int DIVERGED_PARALOG_MIN_THRESHOLD;
    private int STD_TOLERANCE_FOR_DIVERGED_PARALOGS;
    
    /*
     * Takes as input an ancestor tree with all ancestors calculated by AncestorDatabase.java
     * as well as a gene database with all species genes in it, and an output file for gain/loss counts.
     *
     */
    public AGS_analyzer(NewickTree ancestorTree, GeneDatabase GD, String OutputFile, String presentOutputFile,
            boolean reCalculatePresentGenes, boolean reCalculateAGSfile, int DIVERGED_PARALOG_THRESHOLD, int STD_TOLERANCE_FOR_DIVERGED_PARALOGS)
                                    throws IOException
    {
        System.out.println("AGS_analyzer started.");
        this.tree = ancestorTree;
        this.DIVERGED_PARALOG_MIN_THRESHOLD = DIVERGED_PARALOG_THRESHOLD;
        this.STD_TOLERANCE_FOR_DIVERGED_PARALOGS = STD_TOLERANCE_FOR_DIVERGED_PARALOGS;
      
        Vector<Node> allNodes = reOrderNodes(tree.root); // Order nodes for analysis
        
        // For each ancestor, check every gene and set whether it is a sym-bet, or singular or not -- this finds which clade a singular belongs to
        for (int i=0; i<allNodes.size(); i++)
        {
            Node currentNode = allNodes.get(i);
            AGS currentAncestor = currentNode.getAGS();
            if (!currentNode.species_node)
            {
                Vector<String> firstDescent = getAllEndSpecies(currentNode.descendant1);
                Vector<String> secondDescent = getAllEndSpecies(currentNode.descendant2);            
                for (int m=0; m<currentAncestor.size(); m++)
                {
                    OrthologGenes currentGenes = currentAncestor.AGS.get(m);            
                    boolean thereAreFirstLineageMembers = geneIncludesMembersFromAnyOfTheSpecies(currentGenes, firstDescent);
                    boolean thereAreSecondLineageMembers = geneIncludesMembersFromAnyOfTheSpecies(currentGenes, secondDescent); 
                    if (thereAreFirstLineageMembers && !thereAreSecondLineageMembers)
                    {
                        currentGenes.setSingularFirstLineage();
                    }
                    else if (!thereAreFirstLineageMembers && thereAreSecondLineageMembers)
                    {
                        currentGenes.setSingularSecondLineage();
                    }
                    else
                    {
                        currentGenes.setsymbet();
                    }
                }
            }
        }
        System.out.println("Sym-bets and singulars are found");

        // For each non-species node of the tree, calculate avg and std of symbet scores
        for (int i=0; i<allNodes.size(); i++)
        {        
            Node currentNode = allNodes.get(i);
            AGS currentAncestor = currentNode.getAGS();

            // Find avg. and std. deviation of sym-bet scores
            if (!currentNode.species_node)
            {
                currentAncestor.setAvgandSTDsymbetScore();
            }
        }
        System.out.println("Average and STD of sym-bet scores set.");
        
        
        if (reCalculatePresentGenes)
        {
            // For each ancestor, find whether gene was present or not (i.e. lost in the other branch or not), goes from top to bottom
            for (int i=0; i<allNodes.size(); i++)
            {
                Node currentNode = allNodes.get(i);
                AGS currentAncestor = currentNode.getAGS();
                if (!currentNode.species_node)
                {
                    Vector<Node> earlierNodes = getEarlierNodes(currentNode);
                    
                    for (int l=0;l<earlierNodes.size(); l++)
                    {
                        earlierNodes.get(l).getAGS().resetProcessed();
                    }
                    
                    Vector<AGS> earlierAncestors = getAllAGS(earlierNodes);

                    Vector<String> firstDescent = getAllEndSpecies(currentNode.descendant1);
                    Vector<String> secondDescent = getAllEndSpecies(currentNode.descendant2);   

                    for (int m=0; m<currentAncestor.size(); m++)
                    {
                        OrthologGenes currentGenes = currentAncestor.AGS.get(m);
                        if (currentGenes.symbet)
                        {
                            currentGenes.setPresent();
                        }
                        else if (currentGenes.singularFirstLineage 
                            && geneIncludesMembersFromOlderSpeciesButNotFromOtherLineage(currentGenes, earlierNodes, secondDescent))
                        {
                            currentGenes.setPresent();
                        }
                        else if (currentGenes.singularSecondLineage
                            && geneIncludesMembersFromOlderSpeciesButNotFromOtherLineage(currentGenes, earlierNodes, firstDescent))
                        {
                            currentGenes.setPresent();
                        }  

                        if (m%1000==0)
                            System.out.println("Gene " + m + " from ancestor " + i + " is finished.");
                    }
                }
            }
            System.out.println("Present genes found with dollo parsimony.");
        }
        else // Read present genes from file
        {
            BufferedReader input;
            File inputFile = new File(presentOutputFile);
            input = open_file(inputFile);        
            int currentAncestor = -1;
            int currentGene = -1;
            Node currentNode = null;
            
            while (more_records(input))
            {
                String line = read_a_line(input);
                if (line.equals(""))
                {
                    continue;
                }
                else if (line.startsWith("ANCESTOR"))
                {
                    currentGene = -1;
                    currentAncestor++;
                    while (allNodes.get(currentAncestor).species_node)
                        currentAncestor++;
                    currentNode = allNodes.get(currentAncestor);
                    System.out.println("Ancestor " + currentAncestor + " being read from file...");
                }
                else if (line.startsWith("PRESENT"))
                {
                    currentGene++;
                    allNodes.get(currentAncestor).ags.AGS.get(currentGene).setPresent();
                }
                else if (line.startsWith("DIVERGED"))
                {
                    currentGene++;
                    StringTokenizer st = new StringTokenizer(line, "\t");
                    st.nextToken();
                    String singularSourceName = st.nextToken();
                    String divergedGene = st.nextToken();
                    OrthologGenes currentOrthologGenes = currentNode.ags.AGS.get(currentGene);
                    if (!currentOrthologGenes.genes.get(0).equals(divergedGene))
                        System.out.println("Index error happened...");
                    OrthologGenes singularSource = currentNode.ags.findGene(singularSourceName);                    
                    singularSource.addDivergedInParalogs(currentOrthologGenes);
                    currentOrthologGenes.setDivergedParalog();
                    
                }
                else if (line.startsWith("SINGULAR"))
                {
                    currentGene++;
                }
                
            }
            input.close();
            System.out.println("Present genes read from file.");
        }
        
        // For each species and ancestor count the number genes present or not and loss/inparalog/ambiguous_gain counts for each incoming branch
        if (reCalculateAGSfile)
        {
            FileWriter output = create_file(new File(OutputFile));
            if (output == null) throw new IOException();   

            

            // Print header
            output.write("Ancestor name\tSym-bets\tPresent loci\t" +
                    "Loss\tParalogs\tDiverged Paralogs\tAmbiguous-gains\tTotal Gains\tNo scoring genes\t");
            output.write("AVG_Sym-bet\tSTD_Sym-bet");
            output.write(System.getProperty("line.separator"));

            for (int i=0; i<allNodes.size(); i++)
            {
                allNodes.get(i).getAGS().resetProcessed();
            }

            for (int i=0; i<allNodes.size(); i++)
            {
                Node currentNode = allNodes.get(i);
                AGS currentAncestor = currentNode.getAGS();     

                // Root print
                if (currentNode.ancestor == null)
                {
                    int symbetCount = 0;
                    for (int n=0; n<currentAncestor.size(); n++)
                    {

                        OrthologGenes currentGenes = currentAncestor.AGS.get(n);
                        if (currentGenes.present)
                        {
                            symbetCount++;
                        }
                    }    
                    // Print output for this ancestor
                    output.write(currentNode.name + "\t");
                    output.write(symbetCount + "\t\t\t\t\t\t\t\t\t");

                    currentNode.setTotalGeneCount(symbetCount);
  
                    output.write(System.getProperty("line.separator"));                
                    continue;
                }

                // Other prints
                AGS earlierAncestor = currentNode.ancestor.getAGS();

                // Gene family counts
                int symbetCount = 0;
                int presentCount = 0;
                int lossCount = 0;
                int paralogCount = 0;
                int divergedParalogCount = 0;
                int ambiguousGainCount = 0;

                // Loss and paralog counts
                for (int m=0; m<earlierAncestor.size(); m++)
                {
                    OrthologGenes earlierGenes = earlierAncestor.AGS.get(m);
                    if (earlierGenes.present) // Count paralogs
                    {
                        int thisParalogCount = 0;
                        for (int n=0; n<currentAncestor.size(); n++)
                        {
                            OrthologGenes currentGenes = currentAncestor.AGS.get(n);
                            if (currentGenes.present && earlierGenes.contains(currentGenes.genes.get(0)))
                            {
                                currentGenes.setProcessed();
                                thisParalogCount++;
                            }
                        }
                        if (thisParalogCount==0) lossCount++;
                        if (thisParalogCount>1) 
                        {
                            paralogCount += (thisParalogCount-1);
                            earlierGenes.setParalogSource(i);
                        }
                    }
                    else if (!earlierGenes.present) // Count diverged paralogs
                    {
                        for (int n=0; n<currentAncestor.size(); n++)
                        {
                            OrthologGenes currentGenes = currentAncestor.AGS.get(n);
                            if (currentGenes.present && !currentGenes.processed && earlierGenes.contains(currentGenes.genes.get(0)))
                            {
                                OrthologGenes singularSource = earlierAncestor.findGene(earlierGenes.possible_source);
                                int tolerance = STD_TOLERANCE_FOR_DIVERGED_PARALOGS * currentAncestor.std_symbetScore;
                                if (tolerance < DIVERGED_PARALOG_MIN_THRESHOLD)
                                    tolerance = DIVERGED_PARALOG_MIN_THRESHOLD;
                                if (singularSource!=null && singularSource.symbet_score-tolerance < earlierGenes.singular_score)
                                {
                                    singularSource.addDivergedInParalogs(currentGenes);                                     
                                    singularSource.setParalogSource(i);
                                    earlierGenes.setDivergedParalog();
                                    divergedParalogCount++;
                                    currentGenes.setProcessed();
                                }
                            }
                        }                
                    }
                }

                // Sym-bet, presence and ambiguous gain counts
                for (int n=0; n<currentAncestor.size(); n++)
                {
                    OrthologGenes currentGenes = currentAncestor.AGS.get(n);            

                    if (currentGenes.symbet) symbetCount++;
                    if (currentGenes.present)  presentCount++;
                    if (currentGenes.present && (!currentGenes.processed)) ambiguousGainCount++;
                }

                int noScoringGenes = 0;
                if (currentNode.species_node)
                {
                    for (int t=0; t<currentNode.ags.AGS.size(); t++)
                    {
                        if (GD.getGene(GD.getGeneLocationByPreComputedIndex(currentNode.ags.AGS.get(t).genes.get(0))).noScoresFound)
                        {
                            noScoringGenes++;
                        }
                    }
                }                
                
                presentCount = presentCount - noScoringGenes;
                ambiguousGainCount = ambiguousGainCount - noScoringGenes;
                currentNode.setTotalGeneCount(presentCount);
                currentNode.total_gains = paralogCount+divergedParalogCount+ambiguousGainCount;
                currentNode.total_losses = lossCount;
                
                // Print output for this ancestor
                output.write(currentNode.name + "\t");
                output.write(symbetCount + "\t" + presentCount + "\t");
                output.write(lossCount + "\t" + paralogCount + "\t" + divergedParalogCount + "\t" + ambiguousGainCount + "\t" + (paralogCount+divergedParalogCount+ambiguousGainCount + "\t"));
                output.write(noScoringGenes + "\t");


                if (!currentNode.species_node)
                {
                    output.write(currentAncestor.avg_symbetScore + "\t" + currentAncestor.std_symbetScore);
                }
                output.write(System.getProperty("line.separator"));

                System.out.println("Ancestor " + i + " finished writing to AGS file."); 
            }

             // Print tree in newick format with branch lengths calculated from avg. symbet distances
            output.write(System.getProperty("line.separator"));
            output.write(System.getProperty("line.separator"));
            output.write("Ortholog avg. divergence tree: " + System.getProperty("line.separator"));
            output.write(this.tree.root.getTreePrintWithsymbetDistances() +";");  // Average ortholog divergence tree
            output.write(System.getProperty("line.separator"));
            output.write("Gene expansion tree: " + System.getProperty("line.separator"));
            output.write(this.tree.root.getTreePrintWithGeneExpansionDistances() +":"+ this.tree.root.ags.symbetCount() + ";"); // Gene expansion tree
            output.write(System.getProperty("line.separator"));
            output.write("Gene gain tree: " + System.getProperty("line.separator"));
            output.write(this.tree.root.getTreePrintWithGeneGainDistances() + ";"); // Gene gain tree
            output.write(System.getProperty("line.separator"));
            output.write("Gene loss tree: " + System.getProperty("line.separator"));
            output.write(this.tree.root.getTreePrintWithGeneLossDistances() + ";"); // Gene loss tree
            output.write(System.getProperty("line.separator"));            

            output.close(); 
            
            // Write present genes to file
            FileWriter output2 = create_file(new File(presentOutputFile)); // File with evolmap gene ID
            FileWriter output3 = create_file(new File(presentOutputFile + ".rn")); // File with non-unique(?) fasta header identifier lines [i.e. first word of fasta till empty space]
            FileWriter output4 = create_file(new File(presentOutputFile + ".sum")); // File with summary of descent sizes of each node
            if (output2 == null) throw new IOException("Could not open output file.");  
            if (output3 == null) throw new IOException("Could not open output file.");  
            if (output4 == null) throw new IOException("Could not open output file.");
            
            int DESCENT_MAX_SIZE = 200;
            int[] descentSizeCountArray = new int[DESCENT_MAX_SIZE];
            output4.write("Ancestor name\t");
            for (int i=1; i<descentSizeCountArray.length; i++)
            {
                descentSizeCountArray[i] = 0;
                output4.write("" + i + "\t");
            }
            output4.write(System.getProperty("line.separator")); 
            
            for (int t=0; t<allNodes.size(); t++)
            {
                Node nextNode = allNodes.get(t);
                if (!nextNode.species_node)
                {
                    output2.write("ANCESTOR\t" + nextNode.name);
                    output2.write(System.getProperty("line.separator"));
                    output3.write("ANCESTOR\t" + nextNode.name);
                    output3.write(System.getProperty("line.separator"));    
                    
                    for (int i=0; i<nextNode.getAGS().size(); i++)
                    {
                        String output_string = "";
                        String output_string2 = "";
                        OrthologGenes nextOrthologGenes = nextNode.getAGS().AGS.get(i);
                        
                        if (nextOrthologGenes.present)
                        {
                            int size = nextOrthologGenes.genes.size();
                            if (nextOrthologGenes.divergedInparalogs!=null)
                                size += nextOrthologGenes.divergedInparalogs.size();
                            if (size<DESCENT_MAX_SIZE)
                            {
                                descentSizeCountArray[size]++;
                            }
                            else
                            {
                                descentSizeCountArray[DESCENT_MAX_SIZE-1]++; 
                            }
                        }
                        
                        if (nextOrthologGenes.present)
                            output_string += "PRESENT\t";
                        else if (nextOrthologGenes.divergedParalog)
                            output_string += "DIVERGED\t";
                        else
                            output_string += "SINGULAR\t";
                        
                        output_string2 = output_string;
                        
                        if (nextOrthologGenes.divergedParalog)
                        {
                            String source = nextOrthologGenes.possible_source;
                            output_string += source + "\t";
                            
                            Sequence source_gene = GD.getGeneByName2(source);
                            String source_header = source_gene.header;
                            source_header = source_header.substring(source_header.indexOf("\n")+1);
                            int first_space = source_header.indexOf(" ");
                            if (first_space==-1) first_space = source_header.length();
                            String source_fasta_id = source_header.substring(1, first_space);
                            output_string2 += source_fasta_id + "\t";
                            
                            for (int p=0; p<nextOrthologGenes.genes.size(); p++)
                            {
                                String seq = nextOrthologGenes.genes.get(p);
                                output_string += seq + "\t";
                                
                                Sequence seq_real = GD.getGeneByName2(seq);
                                String seq_header = seq_real.header;
                                seq_header = seq_header.substring(seq_header.indexOf("\n")+1);                                
                                first_space = seq_header.indexOf(" ");
                                if (first_space==-1) first_space = seq_header.length();
                                String seq_fasta_id = seq_header.substring(1, first_space);
                                output_string2 += seq_fasta_id + "\t";    
                            }
                        }
                        else
                        {
                            for (int j=0; j<nextOrthologGenes.genes.size(); j++)
                            {
                                String seq = nextOrthologGenes.genes.get(j);
                                output_string += seq + "\t";
                                
                                Sequence seq_real = GD.getGeneByName2(seq);
                                String seq_header = seq_real.header;
                                seq_header = seq_header.substring(seq_header.indexOf("\n")+1);                                    
                                int first_space = seq_header.indexOf(" ");
                                if (first_space==-1) first_space = seq_header.length();
                                String seq_fasta_id = seq_header.substring(1, first_space);
                                output_string2 += seq_fasta_id + "\t";                                   
                            }
                        }
                        
                        output_string += System.getProperty("line.separator");
                        output_string2 += System.getProperty("line.separator");
                        
                        output2.write(output_string);
                        output3.write(output_string2);                      
                    }
                    output4.write(nextNode.name + "\t");
                    for (int i=1; i<descentSizeCountArray.length; i++)
                    {
                        output4.write("" + descentSizeCountArray[i] + "\t");
                    }
                    output4.write(System.getProperty("line.separator"));
                    
                    // Re-initialize descentSizeCountArray
                    for (int i=0; i<descentSizeCountArray.length; i++)
                    {
                        descentSizeCountArray[i] = 0;
                    }                    
                }
            }
            output2.close(); 
            output3.close();
            output4.close();

            System.out.println("Present genes written to file.");            
        }     
    }
    public static Vector<Node> reOrderNodes(Node root)
    {
        if (root.species_node)
        {
            Vector<Node> thisNode = new Vector<Node>();
            thisNode.add(root);
            return thisNode;
        }
        else // Not a species node
        {
            Vector<Node> allNodes = new Vector<Node>();
            allNodes.add(root);
            allNodes.addAll(reOrderNodes(root.descendant1));
            allNodes.addAll(reOrderNodes(root.descendant2));
            return allNodes;
        }        
    }
       
    public static boolean DA1isCoveredByDA2(String currentDA, String otherDA)
    {
        StringTokenizer current_st = new StringTokenizer(currentDA, "-");
        StringTokenizer other_st = new StringTokenizer(otherDA, "-");

        boolean DA1_covers_DA2 = true;

        String currentDomain = current_st.nextToken();
        String otherDomain = other_st.nextToken();

        while (true)
        {
            if (currentDomain.equals(otherDomain))
            {
                if (current_st.hasMoreTokens())
                {
                    currentDomain = current_st.nextToken();
                }
                else // Finished and found all domains
                {
                    break;
                }
                if (other_st.hasMoreTokens())
                {
                    otherDomain = other_st.nextToken();
                }
                else // We got a new first domain but no domains left in second sequence
                {
                    DA1_covers_DA2 = false;
                    break;
                }
            }
            else // Elements not equal
            {
                if (other_st.hasMoreTokens())
                {
                    otherDomain = other_st.nextToken();
                }
                else // We could not find a match to all domains
                {
                    DA1_covers_DA2 = false;
                    break;
                }
            }
        }
        return DA1_covers_DA2;
    }
    public static Vector<Node> getEarlierNodes(Node node)
    {
        Vector<Node> earlierNodes = new Vector<Node>();
        Node earlierNode = node.ancestor;
        while (earlierNode!=null)
        {
            earlierNodes.add(earlierNode);
            earlierNode = earlierNode.ancestor;;
        }
        return earlierNodes;        
    }
    public static Vector<AGS> getAllAGS(Vector<Node> nodes)
    {
        Vector<AGS> allAGS = new Vector<AGS>();
        for (int i=0; i<nodes.size(); i++)
        {
            allAGS.add(nodes.get(i).getAGS());
        }
        return allAGS;
    }
    
    public static Vector<String> getAllEndSpecies(Node node)
    {
        Vector<String> speciesList = new Vector<String>();
        String name = node.name;
        StringTokenizer st = new StringTokenizer(name, "_");
        while(st.hasMoreTokens())
        {
            speciesList.add(st.nextToken());
        }
        return speciesList;
    }
    
    // Method to count symbets
    public static boolean geneIncludesMembersFromAnyOfTheSpecies(OrthologGenes og, Vector<String> speciesList)
    {
        for (int i=0; i<og.genes.size(); i++)
        {
            String currentSpecies = (new StringTokenizer(og.genes.get(i), "_")).nextToken();
            for (int j=0; j<speciesList.size(); j++)
            {
                String otherSpecies = speciesList.get(j);
                if (currentSpecies.equals(otherSpecies))
                    return true;
            }
        }
        return false;
    }

    public boolean geneIncludesMembersFromOlderSpeciesButNotFromOtherLineage(OrthologGenes currentGenes, Vector<Node> earlierNodes, Vector<String> otherDescent)
    {
        boolean lossEvent = false;
        boolean alreadyProcessed = false;
        
        for (int i=0; i<earlierNodes.size(); i++)
        {
            AGS earlierAncestor = earlierNodes.get(i).getAGS();
            for (int j=0; j<earlierAncestor.size(); j++)
            {
                OrthologGenes otherGenes = earlierAncestor.getAGS().get(j);
                
                if (otherGenes.contains(currentGenes.genes.get(0)))
                {
                    Vector<String> firstDescent = getAllEndSpecies(earlierNodes.get(i).descendant1);
                    Vector<String> secondDescent = getAllEndSpecies(earlierNodes.get(i).descendant2);                   
                    
                    if (geneIncludesMembersFromAnyOfTheSpecies(otherGenes, firstDescent) && 
                            geneIncludesMembersFromAnyOfTheSpecies(otherGenes, secondDescent) &&
                            !geneIncludesMembersFromAnyOfTheSpecies(otherGenes, otherDescent))
                    {
                        if (otherGenes.processed)
                        {
                            alreadyProcessed = true;
                            continue;
                        }
                        else
                        {
                            otherGenes.setProcessed();
                            lossEvent = true;
                        }
                    }
                }
            }
        }
        if (alreadyProcessed)
            return false;
        else
            return lossEvent;
    }  
    
    private static BufferedReader open_file(File filename)
    {
            try {
                    return new BufferedReader(new FileReader(filename));
            }
            catch (IOException e) {return null;}
    }

    private static FileWriter create_file(File filename)
    {
            try {
                    return new FileWriter(filename);
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
        
}
