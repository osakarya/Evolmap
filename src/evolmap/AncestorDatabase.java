
package EvolMAP;
import java.util.*;
import java.io.*;

/**
 * This file takes as input a tree, gene database from end-species and a score matrix
 * It calculates ancestors and writes them to file, also sets each node in tree with their ancestors
 * 
 * @author onur
 */
public class AncestorDatabase 
{
    public boolean USE_AVG_SCORES; // If true, uses avg score of members of ortholog families, otherwise uses best score of members of families
    private int[][] ALL_INT_ARRAYS; // Used for faster alignments -- used only for calculating SFA files
    private boolean PROTEIN;
    private ScoringMatrix scoring; // The scoring matrix
    
    public NewickTree tree;
    public GeneDatabase GD;
    public int threshold; // Score threshold for sym-bets, any scores below this level are not considered sym-bets

   
    // This reads a pre-calculated ancestor and sets up tree ancestor nodes
    public AncestorDatabase(NewickTree tree, GeneDatabase gd, String ancestorFile, boolean USE_AVG_SCORES)
                                    throws IOException
    {
        this.USE_AVG_SCORES = USE_AVG_SCORES;
        this.tree = tree;
        this.GD = gd;
        this.setSpeciesGenes(); // For all species nodes on tree, find genes in the gene database and set them
        System.out.println("Reading pre-calculated ancestor pass1 data from file " + ancestorFile);
        this.readAncestors(ancestorFile);
    }
    
    // This reads a pre-calculated ancestor and sets up tree ancestor nodes
    public AncestorDatabase(NewickTree tree, GeneDatabase gd, String outputName, boolean USE_AVG_SCORES, String mainTreeName, 
            ArrayList<String> subtree_names, ArrayList<String> subtrees)
                                    throws IOException, InvalidScoringMatrixException, IncompatibleScoringSchemeException
    {
        this.USE_AVG_SCORES = USE_AVG_SCORES;
        this.tree = tree;
        this.GD = gd;
        this.setSpeciesGenes(); // For all species nodes on tree, find genes in the gene database and set them
        this.readAncestors(mainTreeName, subtree_names, subtrees);
        this.writeAncestors(outputName, false, "");// If we like to calculate SFA file for this too, then need to write ancestors here
    }
        
    // MODE 0 -- score matrix file is regular (either full sequence or domain full sequence) but it is not partitioned or is domain-wise
    /** Creates a new instance of AncestorDatabase */
    public AncestorDatabase(NewickTree tree, GeneDatabase gd, String outputName, int threshold, boolean SPARSE_MATRIX, 
            int ALIGNMENT_LIMIT, boolean USE_AVG_SCORES, boolean PROTEIN, String sfaFileName, boolean CALCULATE_SFA_FILE) 
                                            throws IOException, InvalidScoringMatrixException, IncompatibleScoringSchemeException
    {
        this.USE_AVG_SCORES = USE_AVG_SCORES;
        this.tree = tree;
        this.GD = gd;
        this.threshold = threshold;
        this.PROTEIN = PROTEIN;
        
        // Matrix file for protein inputs
        File blosum62Matrix = new File("blosum62.txt");
        FileReader matrix_file = new FileReader(blosum62Matrix);
        scoring = new ScoringMatrix (matrix_file); 
        
        this.GD.setScoreMatrixPositions(); // Set positions of each sequence on the score matrix
        System.out.println("Score matrix positions are set.");
        this.setSpeciesGenes(); // For all species nodes on tree, find genes in the gene database and set them
        System.out.println("Species ancestors are set.");
        this.setAGS(tree.root, SPARSE_MATRIX, ALIGNMENT_LIMIT);
        this.writeAncestors(outputName, CALCULATE_SFA_FILE, sfaFileName);
    }
    
    private void setSpeciesGenes()
    {
        for (int i=0; i<this.tree.size(); i++)
        {
            Node node = this.tree.getNode(i);
            if (node.species_node)
            {
                GeneDatabase db_new = this.GD.getSpeciesGenes(node.name);
                AGS new_ags = new AGS(db_new); 
                node.setAGS(new_ags);
            }
        }
    }
    
    private void setAGS(Node node, boolean SPARSE_MATRIX, int ALIGNMENT_LIMIT)
                        throws IOException
    {
        if (node.descendant1.getAGS()==null)
            this.setAGS(node.descendant1, SPARSE_MATRIX, ALIGNMENT_LIMIT);
        if (node.descendant2.getAGS()==null)
            this.setAGS(node.descendant2, SPARSE_MATRIX, ALIGNMENT_LIMIT);       
        
        if (!SPARSE_MATRIX)
            node.setAGS(calculateAncestor(node.descendant1, node.descendant2));
        else
            node.setAGS(calculateSparseAncestor(node.descendant1, node.descendant2, ALIGNMENT_LIMIT));
        
        System.out.println("AGS for ancestor " + node.name + " is calculated");
    }
    
    private AGS calculateAncestor(Node node1, Node node2)
    {
        Vector<OrthologGenes> seqSet1 = node1.getAGS().getAGS();
        Vector<OrthologGenes> seqSet2 = node2.getAGS().getAGS();

        Vector<OrthologGenes> ancestralGeneSet = new Vector<OrthologGenes>();
        Vector<OrthologGenes> orthologSet = new Vector<OrthologGenes>();
        Vector<OrthologGenes> singularSet1 = new Vector<OrthologGenes>();
        Vector<OrthologGenes> singularSet2 = new Vector<OrthologGenes>();        
        
        int[][] set1_scores_towards_set1 = new int[seqSet1.size()][seqSet1.size()];
        int[][] set1_scores_towards_set2 = new int[seqSet1.size()][seqSet2.size()];
        int[][] set2_scores_towards_set1 = new int[seqSet2.size()][seqSet1.size()];
        int[][] set2_scores_towards_set2 = new int[seqSet2.size()][seqSet2.size()];
        
        // Finds scores of each gene in set1 against set1
        for (int i=0; i<seqSet1.size(); i++)
        {    
            OrthologGenes set1currentGenes = seqSet1.get(i);

            for (int j=0; j<seqSet1.size(); j++)
            {
                if (i==j) 
                {
                    set1_scores_towards_set1[i][j] = 0;
                    continue;
                }
                OrthologGenes set1currentGenes_2 = seqSet1.get(j);
                int total_score = 0;
                int best_score = 0;
                
                for (int m=0; m<set1currentGenes.numberOfSequences(); m++)
                {
                    Sequence set1currentSequence = GD.getGeneByPreComputedIndex(set1currentGenes.getGenes().get(m));
                    for (int n=0; n<set1currentGenes_2.numberOfSequences(); n++)
                    {
                        Sequence set1currentSequence_2 = GD.getGeneByPreComputedIndex(set1currentGenes_2.getGenes().get(n));
                        int current_score = getScore(set1currentSequence, set1currentSequence_2);
                        if (current_score > best_score)
                            best_score = current_score;
                        total_score += current_score;
                    }
                }
                int normalized_score = total_score / (set1currentGenes.numberOfSequences() * set1currentGenes_2.numberOfSequences());

                if (USE_AVG_SCORES)
                    set1_scores_towards_set1[i][j]  = normalized_score;
                else
                    set1_scores_towards_set1[i][j]  = best_score;   
            }
        }        
        
        // Finds scores of each gene in set1 against set2
        for (int i=0; i<seqSet1.size(); i++)
        {    
            OrthologGenes set1currentGenes = seqSet1.get(i);
            for (int j=0; j<seqSet2.size(); j++)
            {
                OrthologGenes set2currentGenes = seqSet2.get(j);
                int total_score = 0;
                int best_score = 0;
                
                for (int m=0; m<set1currentGenes.numberOfSequences(); m++)
                {
                    Sequence set1currentSequence = GD.getGeneByPreComputedIndex(set1currentGenes.getGenes().get(m));
                    for (int n=0; n<set2currentGenes.numberOfSequences(); n++)
                    {
                        Sequence set2currentSequence = GD.getGeneByPreComputedIndex(set2currentGenes.getGenes().get(n));
                        int current_score = getScore(set1currentSequence, set2currentSequence);
                        if (current_score > best_score)
                            best_score = current_score;
                        total_score += current_score;
                    }
                }
                int normalized_score = total_score / (set1currentGenes.numberOfSequences() * set2currentGenes.numberOfSequences());
                
                //System.out.print(normalized_score + "\t");
                
                if (USE_AVG_SCORES)
                    set1_scores_towards_set2[i][j] = normalized_score;
                else
                    set1_scores_towards_set2[i][j] = best_score;                
            }
            //System.out.println();
        }
        
        // Finds scores of each gene in set2 against set2
        for (int i=0; i<seqSet2.size(); i++)
        {    
            OrthologGenes set2currentGenes = seqSet2.get(i);

            for (int j=0; j<seqSet2.size(); j++)
            {
                if (i==j)
                {
                    set2_scores_towards_set2[i][j] = 0;
                    continue;
                }
                OrthologGenes set2currentGenes_2 = seqSet2.get(j);
                int total_score = 0;
                int best_score = 0;
                
                for (int m=0; m<set2currentGenes.numberOfSequences(); m++)
                {
                    Sequence set2currentSequence = GD.getGeneByPreComputedIndex(set2currentGenes.getGenes().get(m));
                    for (int n=0; n<set2currentGenes_2.numberOfSequences(); n++)
                    {
                        Sequence set2currentSequence_2 = GD.getGeneByPreComputedIndex(set2currentGenes_2.getGenes().get(n));
                        int current_score = getScore(set2currentSequence, set2currentSequence_2);
                        if (current_score > best_score)
                            best_score = current_score;                        
                        total_score += current_score;
                    }
                }
                int normalized_score = total_score / (set2currentGenes.numberOfSequences() * set2currentGenes_2.numberOfSequences());
                
                if (USE_AVG_SCORES)
                    set2_scores_towards_set2[i][j]  = normalized_score;
                else
                    set2_scores_towards_set2[i][j]  = best_score;
            }
        }  
        
        // Finds scores of each gene in set2 against set1
        for (int i=0; i<seqSet2.size(); i++)
        {    
            OrthologGenes set2currentGenes = seqSet2.get(i);

            for (int j=0; j<seqSet1.size(); j++)
            {
                OrthologGenes set1currentGenes = seqSet1.get(j);
                int total_score = 0;
                int best_score = 0;
                for (int m=0; m<set2currentGenes.numberOfSequences(); m++)
                {
                    Sequence set2currentSequence = GD.getGeneByPreComputedIndex(set2currentGenes.getGenes().get(m));
                    for (int n=0; n<set1currentGenes.numberOfSequences(); n++)
                    {
                        Sequence set1currentSequence = GD.getGeneByPreComputedIndex(set1currentGenes.getGenes().get(n));
                        int current_score = getScore(set2currentSequence, set1currentSequence);
                        
                        if (current_score > best_score)
                            best_score = current_score;
                        total_score += current_score;
                    }
                }
                int normalized_score = total_score / (set2currentGenes.numberOfSequences() * set1currentGenes.numberOfSequences());
                
                if (USE_AVG_SCORES)
                    set2_scores_towards_set1[i][j] = normalized_score;
                else
                    set2_scores_towards_set1[i][j] = best_score;
            }
        }
        
        int[] set1_towards_set2_best_scores = new int[seqSet1.size()];
        int[] set2_towards_set1_best_scores = new int[seqSet2.size()];
        int[] set1_towards_set1_best_scores = new int[seqSet1.size()];
        int[] set2_towards_set2_best_scores = new int[seqSet2.size()];
        
        int[] set1_towards_set2_best_hit_locations = new int[seqSet1.size()]; // This is the position of the gene in set2 that has a best hit towards set1
        int[] set2_towards_set1_best_hit_locations = new int[seqSet2.size()]; // This is the position of the gene in set1 that has a best hit towards set2
        int[] set1_towards_set1_best_hit_locations = new int[seqSet1.size()]; // This is used to predict where a singular gene might be from
        int[] set2_towards_set2_best_hit_locations = new int[seqSet2.size()];
        
        // Find all best-hits of set 1 to set 2
        for (int i=0; i<set1_scores_towards_set2.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set1_scores_towards_set2[i].length; j++)
            {
                if (set1_scores_towards_set2[i][j] > best_score)
                {
                    best_hit_location = j;
                    best_score = set1_scores_towards_set2[i][j];
                }
            }
            set1_towards_set2_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set1_towards_set2_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set1_towards_set2_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }
        
        // Find all best-hits of set 2 to set 1
        for (int i=0; i<set2_scores_towards_set1.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set2_scores_towards_set1[i].length; j++)
            {
                if (set2_scores_towards_set1[i][j] > best_score)
                {
                    best_hit_location = j;
                    best_score = set2_scores_towards_set1[i][j];
                }
            }
            set2_towards_set1_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set2_towards_set1_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set2_towards_set1_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }
        
        // Find all best-hits of set 1 to set 1
        for (int i=0; i<set1_scores_towards_set1.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set1_scores_towards_set1[i].length; j++)
            {
                if (set1_scores_towards_set1[i][j] > best_score)
                {
                    best_hit_location = j;
                    best_score = set1_scores_towards_set1[i][j];
                }
            }
            set1_towards_set1_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set1_towards_set1_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set1_towards_set1_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }        
        
        
        // Find all best-hits of set 2 to set 2
        for (int i=0; i<set2_scores_towards_set2.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set2_scores_towards_set2[i].length; j++)
            {
                if (set2_scores_towards_set2[i][j] > best_score)
                {
                    best_hit_location = j;
                    best_score = set2_scores_towards_set2[i][j];
                }
            }
            set2_towards_set2_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set2_towards_set2_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set2_towards_set2_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }
        
        // Find all sym-bets
        boolean[] set1_symbet_genes = new boolean[seqSet1.size()];
        boolean[] set2_symbet_genes = new boolean[seqSet2.size()];
        
        for(int i=0; i<set1_symbet_genes.length; i++)
        {
            set1_symbet_genes[i] = false;
        }
        
        for(int i=0; i<set2_symbet_genes.length; i++)
        {
            set2_symbet_genes[i] = false;
        }        
        
        for (int i=0; i<set1_towards_set2_best_hit_locations.length; i++)
        {
            int thisHitLocation = set1_towards_set2_best_hit_locations[i];
            if (thisHitLocation == -1)
                continue;
            int otherHitLocation = set2_towards_set1_best_hit_locations[thisHitLocation];
            if (i == otherHitLocation)
            {
                set1_symbet_genes[i] = true;
            }
        }
        
        for (int i=0; i<set2_towards_set1_best_hit_locations.length; i++)
        {
            int thisHitLocation = set2_towards_set1_best_hit_locations[i];
            if (thisHitLocation == -1)
                continue;           
            int otherHitLocation = set1_towards_set2_best_hit_locations[thisHitLocation];
            if (i == otherHitLocation)
            {
                set2_symbet_genes[i] = true;
            }
        }
        
        // Find all inparalogs of sym-bets and add them to output
        boolean[] set1_gene_used = new boolean[seqSet1.size()];
        boolean[] set2_gene_used = new boolean[seqSet2.size()];
        
        for(int i=0; i<set1_gene_used.length; i++)
        {
            set1_gene_used[i] = false;
        }
        
        for(int i=0; i<set2_gene_used.length; i++)
        {
            set2_gene_used[i] = false;
        }
        
        for (int i=0; i<set1_symbet_genes.length; i++)
        {
            if (set1_symbet_genes[i] == true)
            {
                OrthologGenes orthologGene1 = new OrthologGenes();
                orthologGene1.addOrthologGenes(seqSet1.get(i));
                int best_hit_location = set1_towards_set2_best_hit_locations[i];
                OrthologGenes orthologGene2 = new OrthologGenes();
                orthologGene2.addOrthologGenes(seqSet2.get(best_hit_location));
                set1_gene_used[i] = true;
                set2_gene_used[best_hit_location] = true;
                
                Vector<Integer> ortholog1_locations = new Vector<Integer>();
                ortholog1_locations.add(i);
                
                // Find inparalogs of set 1 -- use some sort of clustering algorithm here to find inparalogs of inparalogs
                for (int j=0; j<set1_scores_towards_set1.length; j++)
                {
                    if (set1_gene_used[j] == false &&
                            set1_symbet_genes[j] == false)
                    {
                        for (int t=0; t<ortholog1_locations.size(); t++)
                        {
                            if ((set1_towards_set2_best_scores[i]) < set1_scores_towards_set1[ortholog1_locations.get(t)][j])
                            {
                                orthologGene1.addOrthologGenes(seqSet1.get(j));
                                set1_gene_used[j] = true;
                                ortholog1_locations.add(j);
                                j=0;
                                break;
                            }
                        }
                    }
                }
                
                Vector<Integer> ortholog2_locations = new Vector<Integer>();
                ortholog2_locations.add(best_hit_location);
                
                // Find inparalogs of set 2
                for (int j=0; j<set2_scores_towards_set2.length; j++)
                {
                    if (set2_gene_used[j] == false &&
                            set2_symbet_genes[j] == false)
                    {
                        for (int t=0; t<ortholog2_locations.size(); t++)
                        {                        
                            if ((set2_towards_set1_best_scores[best_hit_location]) < set2_scores_towards_set2[ortholog2_locations.get(t)][j])
                            {
                                orthologGene2.addOrthologGenes(seqSet2.get(j));
                                set2_gene_used[j] = true;
                                ortholog2_locations.add(j);
                                j=0;
                                break;
                            }
                        }
                    }
                }           
                
                OrthologGenes ancestralGene = new OrthologGenes();
                ancestralGene.addOrthologGenes(orthologGene1);
                ancestralGene.addOrthologGenes(orthologGene2);
                ancestralGene.setsymbetScore(set1_scores_towards_set2[i][best_hit_location]);
                
                orthologSet.add(ancestralGene);
            }
        }
        
        // At this point all orthologs are added to the set, and their symbet scores are set.
        // Sources of singulars should be found at this time including later-present genes (decided by AGSanalyzer.java)
        // Each singular can be sourced only by genes with a symbet score.
        
        // Reset set1 towards set1 and set2 towards set2 according to symbet scores
        // Find all best-hits of set 1 to set 1
        for (int i=0; i<set1_scores_towards_set1.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set1_scores_towards_set1[i].length; j++)
            {
                if (set1_gene_used[j]==true && set1_scores_towards_set1[i][j] > best_score)
                {
                    best_hit_location = j;
                    best_score = set1_scores_towards_set1[i][j];
                }
            }
            set1_towards_set1_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set1_towards_set1_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set1_towards_set1_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }        
        
        
        // Find all best-hits of set 2 to set 2
        for (int i=0; i<set2_scores_towards_set2.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set2_scores_towards_set2[i].length; j++)
            {
                if (set2_gene_used[j]==true && set2_scores_towards_set2[i][j] > best_score)
                {
                    best_hit_location = j;
                    best_score = set2_scores_towards_set2[i][j];
                }
            }
            set2_towards_set2_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set2_towards_set2_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set2_towards_set2_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }
        
        // Find singulars of set1 and set their possible source
        for (int i=0; i<set1_gene_used.length; i++)
        {
            if (set1_gene_used[i] == false)
            {
                OrthologGenes og = new OrthologGenes();
                og.setPossibleSource("Unknown", 0);
                if (set1_towards_set1_best_hit_locations[i] != -1)
                {
                    og.setPossibleSource(seqSet1.get(set1_towards_set1_best_hit_locations[i]).genes.get(0), set1_scores_towards_set1[i][set1_towards_set1_best_hit_locations[i]]);
                }
                og.addOrthologGenes(seqSet1.get(i));
                singularSet1.add(og);
            }
        }
        
        // Find singulars for set2 and set their possible source
        for (int i=0; i<set2_gene_used.length; i++)
        {
            if (set2_gene_used[i] == false)
            {
                OrthologGenes og = new OrthologGenes();
                og.setPossibleSource("Unknown", 0);                
                if (set2_towards_set2_best_hit_locations[i] != -1)
                {            
                    og.setPossibleSource(seqSet2.get(set2_towards_set2_best_hit_locations[i]).genes.get(0), set2_scores_towards_set2[i][set2_towards_set2_best_hit_locations[i]]);
                }
                og.addOrthologGenes(seqSet2.get(i));
                singularSet2.add(og);
            }
        }
        
        // Add all genes to one vector
        ancestralGeneSet.addAll(orthologSet);
        ancestralGeneSet.addAll(singularSet1);
        ancestralGeneSet.addAll(singularSet2);
        
        return new AGS(ancestralGeneSet);
    }
    private AGS calculateSparseAncestor(Node node1, Node node2, int ALIGNMENT_LIMIT)
    {
        Vector<OrthologGenes> seqSet1 = node1.getAGS().getAGS();
        Vector<OrthologGenes> seqSet2 = node2.getAGS().getAGS();
        
        Vector<OrthologGenes> ancestralGeneSet = new Vector<OrthologGenes>();
        Vector<OrthologGenes> orthologSet = new Vector<OrthologGenes>();
        Vector<OrthologGenes> singularSet1 = new Vector<OrthologGenes>();
        Vector<OrthologGenes> singularSet2 = new Vector<OrthologGenes>();        
        
        
        ArrayList<Integer>[] set1_scores_towards_set1 = new ArrayList[seqSet1.size()];
        ArrayList<Integer>[] set1_scores_towards_set1_locations = new ArrayList[seqSet1.size()];
        for (int i=0; i<set1_scores_towards_set1.length; i++)
        {
            set1_scores_towards_set1[i] = new ArrayList<Integer>(ALIGNMENT_LIMIT);
            set1_scores_towards_set1_locations[i] = new ArrayList<Integer>(ALIGNMENT_LIMIT);
        }
        
        ArrayList<Integer>[] set1_scores_towards_set2 = new ArrayList[seqSet1.size()];
        ArrayList<Integer>[] set1_scores_towards_set2_locations = new ArrayList[seqSet1.size()];
        for (int i=0; i<set1_scores_towards_set2.length; i++)
        {
            set1_scores_towards_set2[i] = new ArrayList<Integer>(ALIGNMENT_LIMIT);
            set1_scores_towards_set2_locations[i] = new ArrayList<Integer>(ALIGNMENT_LIMIT);
        }
        
        ArrayList<Integer>[] set2_scores_towards_set1 = new ArrayList[seqSet2.size()];
        ArrayList<Integer>[] set2_scores_towards_set1_locations = new ArrayList[seqSet2.size()];
        for (int i=0; i<set2_scores_towards_set1.length; i++)
        {
            set2_scores_towards_set1[i] = new ArrayList<Integer>(ALIGNMENT_LIMIT);
            set2_scores_towards_set1_locations[i] = new ArrayList<Integer>(ALIGNMENT_LIMIT);
        }
        
        ArrayList<Integer>[] set2_scores_towards_set2 = new ArrayList[seqSet2.size()];
        ArrayList<Integer>[] set2_scores_towards_set2_locations = new ArrayList[seqSet2.size()];
        for (int i=0; i<set2_scores_towards_set2.length; i++)
        {
            set2_scores_towards_set2[i] = new ArrayList<Integer>(ALIGNMENT_LIMIT);
            set2_scores_towards_set2_locations[i] = new ArrayList<Integer>(ALIGNMENT_LIMIT);
        }
        
        // Finds scores of each gene in set1 against set1
        for (int i=0; i<seqSet1.size(); i++)
        {    
            OrthologGenes set1currentGenes = seqSet1.get(i);

            for (int j=0; j<seqSet1.size(); j++)
            {
                if (i==j) 
                {
                    //set1_scores_towards_set1[i][j] = 0;
                    continue;
                }
                OrthologGenes set1currentGenes_2 = seqSet1.get(j);
                int total_score = 0;
                int best_score = 0;
                int total_comparisons = 0;
                
                for (int m=0; m<set1currentGenes.numberOfSequences(); m++)
                {
                    Sequence set1currentSequence = GD.getGeneByPreComputedIndex(set1currentGenes.getGenes().get(m));
                    for (int n=0; n<set1currentGenes_2.numberOfSequences(); n++)
                    {
                        Sequence set1currentSequence_2 = GD.getGeneByPreComputedIndex(set1currentGenes_2.getGenes().get(n));
                        int current_score = getScore(set1currentSequence, set1currentSequence_2);
                        if (current_score > 0)
                            total_comparisons++;
                        if (current_score > best_score)
                            best_score = current_score;
                        total_score += current_score;
                    }
                }
                if (total_comparisons==0)
                    total_comparisons = 1;
                int normalized_score = total_score / total_comparisons;
                
                if (normalized_score > 0)
                {
                    if (USE_AVG_SCORES)
                    {
                        set1_scores_towards_set1[i].add(normalized_score);
                        set1_scores_towards_set1_locations[i].add(j);
                    }
                    else
                    {
                        set1_scores_towards_set1[i].add(best_score);
                        set1_scores_towards_set1_locations[i].add(j);
                    }
                }
            }
        }        
        
        // Finds scores of each gene in set1 against set2
        for (int i=0; i<seqSet1.size(); i++)
        {    
            OrthologGenes set1currentGenes = seqSet1.get(i);
            for (int j=0; j<seqSet2.size(); j++)
            {
                OrthologGenes set2currentGenes = seqSet2.get(j);
                int total_score = 0;
                int best_score = 0;
                int total_comparisons = 0;
                
                for (int m=0; m<set1currentGenes.numberOfSequences(); m++)
                {
                    Sequence set1currentSequence = GD.getGeneByPreComputedIndex(set1currentGenes.getGenes().get(m));
                    for (int n=0; n<set2currentGenes.numberOfSequences(); n++)
                    {
                        Sequence set2currentSequence = GD.getGeneByPreComputedIndex(set2currentGenes.getGenes().get(n));
                        int current_score = getScore(set1currentSequence, set2currentSequence);
                        if (current_score > 0)
                            total_comparisons++;                        
                        if (current_score > best_score)
                            best_score = current_score;
                        total_score += current_score;
                    }
                }
                
                if (total_comparisons==0)
                    total_comparisons = 1;                
                int normalized_score = total_score / total_comparisons;
                
                //System.out.print(normalized_score + "\t");
                if (normalized_score > 0)
                {
                    if (USE_AVG_SCORES)
                    {
                        set1_scores_towards_set2[i].add(normalized_score);
                        set1_scores_towards_set2_locations[i].add(j);
                    }
                    else
                    {
                        set1_scores_towards_set2[i].add(best_score);
                        set1_scores_towards_set2_locations[i].add(j);
                    }
                }
            }
            //System.out.println();
        }
        
        // Finds scores of each gene in set2 against set2
        for (int i=0; i<seqSet2.size(); i++)
        {    
            OrthologGenes set2currentGenes = seqSet2.get(i);

            for (int j=0; j<seqSet2.size(); j++)
            {
                if (i==j)
                {
                    //set2_scores_towards_set2[i][j] = 0;
                    continue;
                }
                OrthologGenes set2currentGenes_2 = seqSet2.get(j);
                int total_score = 0;
                int best_score = 0;
                int total_comparisons = 0;
                
                for (int m=0; m<set2currentGenes.numberOfSequences(); m++)
                {
                    Sequence set2currentSequence = GD.getGeneByPreComputedIndex(set2currentGenes.getGenes().get(m));
                    for (int n=0; n<set2currentGenes_2.numberOfSequences(); n++)
                    {
                        Sequence set2currentSequence_2 = GD.getGeneByPreComputedIndex(set2currentGenes_2.getGenes().get(n));
                        int current_score = getScore(set2currentSequence, set2currentSequence_2);
                        if (current_score > 0)
                            total_comparisons++;
                        if (current_score > best_score)
                            best_score = current_score;                        
                        total_score += current_score;
                    }
                }
                if (total_comparisons==0)
                    total_comparisons = 1;                
                int normalized_score = total_score / total_comparisons;
                
                if (normalized_score > 0)
                {
                    if (USE_AVG_SCORES)
                    {
                        set2_scores_towards_set2[i].add(normalized_score);
                        set2_scores_towards_set2_locations[i].add(j);
                    }
                    else
                    {
                        set2_scores_towards_set2[i].add(best_score);
                        set2_scores_towards_set2_locations[i].add(j);
                    }
                }
            }
        }  
        
        // Finds scores of each gene in set2 against set1
        for (int i=0; i<seqSet2.size(); i++)
        {    
            OrthologGenes set2currentGenes = seqSet2.get(i);

            for (int j=0; j<seqSet1.size(); j++)
            {
                OrthologGenes set1currentGenes = seqSet1.get(j);
                int total_score = 0;
                int best_score = 0;
                int total_comparisons = 0;
                
                for (int m=0; m<set2currentGenes.numberOfSequences(); m++)
                {
                    Sequence set2currentSequence = GD.getGeneByPreComputedIndex(set2currentGenes.getGenes().get(m));
                    for (int n=0; n<set1currentGenes.numberOfSequences(); n++)
                    {
                        Sequence set1currentSequence = GD.getGeneByPreComputedIndex(set1currentGenes.getGenes().get(n));
                        int current_score = getScore(set2currentSequence, set1currentSequence);
                        if (current_score > 0)
                            total_comparisons++;
                        if (current_score > best_score)
                            best_score = current_score;
                        total_score += current_score;
                    }
                }
                
                if (total_comparisons==0)
                    total_comparisons = 1;                
                int normalized_score = total_score / total_comparisons;
                
                if (normalized_score > 0)
                {
                    if (USE_AVG_SCORES)
                    {
                        set2_scores_towards_set1[i].add(normalized_score);
                        set2_scores_towards_set1_locations[i].add(j);
                    }
                    else
                    {
                        set2_scores_towards_set1[i].add(best_score);
                        set2_scores_towards_set1_locations[i].add(j);
                    }
                }
            }
        }
        
        int[] set1_towards_set2_best_scores = new int[seqSet1.size()];
        int[] set2_towards_set1_best_scores = new int[seqSet2.size()];
        int[] set1_towards_set1_best_scores = new int[seqSet1.size()];
        int[] set2_towards_set2_best_scores = new int[seqSet2.size()];
        
        int[] set1_towards_set2_best_hit_locations = new int[seqSet1.size()]; // This is the position of the gene in set2 that has a best hit towards set1
        int[] set2_towards_set1_best_hit_locations = new int[seqSet2.size()]; // This is the position of the gene in set1 that has a best hit towards set2
        int[] set1_towards_set1_best_hit_locations = new int[seqSet1.size()]; // This is used to predict where a singular gene might be from
        int[] set2_towards_set2_best_hit_locations = new int[seqSet2.size()];
        
        // Find all best-hits of set 1 to set 2
        for (int i=0; i<set1_scores_towards_set2.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set1_scores_towards_set2[i].size(); j++)
            {
                if (set1_scores_towards_set2[i].get(j) > best_score)
                {
                    best_hit_location = set1_scores_towards_set2_locations[i].get(j);
                    best_score = set1_scores_towards_set2[i].get(j);
                }
            }
            set1_towards_set2_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set1_towards_set2_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set1_towards_set2_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }
        
        // Find all best-hits of set 2 to set 1
        for (int i=0; i<set2_scores_towards_set1.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set2_scores_towards_set1[i].size(); j++)
            {
                if (set2_scores_towards_set1[i].get(j) > best_score)
                {
                    best_hit_location = set2_scores_towards_set1_locations[i].get(j);
                    best_score = set2_scores_towards_set1[i].get(j);
                }
            }
            set2_towards_set1_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set2_towards_set1_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set2_towards_set1_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }
        
        // Find all best-hits of set 1 to set 1
        for (int i=0; i<set1_scores_towards_set1.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set1_scores_towards_set1[i].size(); j++)
            {
                if (set1_scores_towards_set1[i].get(j) > best_score)
                {
                    best_hit_location = set1_scores_towards_set1_locations[i].get(j);
                    best_score = set1_scores_towards_set1[i].get(j);
                }
            }
            set1_towards_set1_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set1_towards_set1_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set1_towards_set1_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }        
        
        
        // Find all best-hits of set 2 to set 2
        for (int i=0; i<set2_scores_towards_set2.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set2_scores_towards_set2[i].size(); j++)
            {
                if (set2_scores_towards_set2[i].get(j) > best_score)
                {
                    best_hit_location = set2_scores_towards_set2_locations[i].get(j);
                    best_score = set2_scores_towards_set2[i].get(j);
                }
            }
            set2_towards_set2_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set2_towards_set2_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set2_towards_set2_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }
        
        // Find all sym-bets
        boolean[] set1_symbet_genes = new boolean[seqSet1.size()];
        boolean[] set2_symbet_genes = new boolean[seqSet2.size()];
        
        for(int i=0; i<set1_symbet_genes.length; i++)
        {
            set1_symbet_genes[i] = false;
        }
        
        for(int i=0; i<set2_symbet_genes.length; i++)
        {
            set2_symbet_genes[i] = false;
        }        
        
        for (int i=0; i<set1_towards_set2_best_hit_locations.length; i++)
        {
            int thisHitLocation = set1_towards_set2_best_hit_locations[i];
            if (thisHitLocation == -1)
                continue;
            int otherHitLocation = set2_towards_set1_best_hit_locations[thisHitLocation];
            if (i == otherHitLocation)
            {
                set1_symbet_genes[i] = true;
            }
        }
        
        for (int i=0; i<set2_towards_set1_best_hit_locations.length; i++)
        {
            int thisHitLocation = set2_towards_set1_best_hit_locations[i];
            if (thisHitLocation == -1)
                continue;           
            int otherHitLocation = set1_towards_set2_best_hit_locations[thisHitLocation];
            if (i == otherHitLocation)
            {
                set2_symbet_genes[i] = true;
            }
        }
        
        // Find all inparalogs of sym-bets and add them to output
        boolean[] set1_gene_used = new boolean[seqSet1.size()];
        boolean[] set2_gene_used = new boolean[seqSet2.size()];
        
        for(int i=0; i<set1_gene_used.length; i++)
        {
            set1_gene_used[i] = false;
        }
        
        for(int i=0; i<set2_gene_used.length; i++)
        {
            set2_gene_used[i] = false;
        }
        
        for (int i=0; i<set1_symbet_genes.length; i++)
        {
            if (set1_symbet_genes[i] == true)
            {
                OrthologGenes orthologGene1 = new OrthologGenes();
                orthologGene1.addOrthologGenes(seqSet1.get(i));
                int best_hit_location = set1_towards_set2_best_hit_locations[i];
                OrthologGenes orthologGene2 = new OrthologGenes();
                orthologGene2.addOrthologGenes(seqSet2.get(best_hit_location));
                set1_gene_used[i] = true;
                set2_gene_used[best_hit_location] = true;
                
                Vector<Integer> ortholog1_locations = new Vector<Integer>();
                ortholog1_locations.add(i);
                
                // Find inparalogs of set 1 -- use some sort of clustering algorithm here to find inparalogs of inparalogs
                for (int j=0; j<set1_scores_towards_set1[i].size(); j++)
                {
                    int realLocation = set1_scores_towards_set1_locations[i].get(j);
                    if (set1_gene_used[realLocation] == false &&
                            set1_symbet_genes[realLocation] == false)
                    {
                        for (int t=0; t<ortholog1_locations.size(); t++)
                        {
                            int paralogLocation = set1_scores_towards_set1_locations[ortholog1_locations.get(t)].indexOf(realLocation);
                            
                            if ((paralogLocation!=-1) && 
                                    (set1_towards_set2_best_scores[i] < set1_scores_towards_set1[ortholog1_locations.get(t)].get(paralogLocation)))
                            {
                                orthologGene1.addOrthologGenes(seqSet1.get(realLocation));
                                set1_gene_used[realLocation] = true;
                                ortholog1_locations.add(realLocation);
                                j=0;
                                break;
                            }
                        }
                    }
                }
                
                Vector<Integer> ortholog2_locations = new Vector<Integer>();
                ortholog2_locations.add(best_hit_location);
                
                // Find inparalogs of set 2
                for (int j=0; j<set2_scores_towards_set2[best_hit_location].size(); j++)
                {
                    int realLocation = set2_scores_towards_set2_locations[best_hit_location].get(j);
                    if (set2_gene_used[realLocation] == false &&
                            set2_symbet_genes[realLocation] == false)
                    {
                        for (int t=0; t<ortholog2_locations.size(); t++)
                        {                        
                            int paralogLocation = set2_scores_towards_set2_locations[ortholog2_locations.get(t)].indexOf(realLocation);
                            
                            if ((paralogLocation!=-1) && 
                                    (set2_towards_set1_best_scores[best_hit_location] < set2_scores_towards_set2[ortholog2_locations.get(t)].get(paralogLocation)))
                            {
                                orthologGene2.addOrthologGenes(seqSet2.get(realLocation));
                                set2_gene_used[realLocation] = true;
                                ortholog2_locations.add(realLocation);
                                j=0;
                                break;
                            }
                        }
                    }
                }           
                
                OrthologGenes ancestralGene = new OrthologGenes();
                ancestralGene.addOrthologGenes(orthologGene1);
                ancestralGene.addOrthologGenes(orthologGene2);
    
                int symbetLocation = set1_scores_towards_set2_locations[i].indexOf(best_hit_location);
                if (symbetLocation != -1)
                    ancestralGene.setsymbetScore(set1_scores_towards_set2[i].get(symbetLocation));
                else
                {
                    // Is this always right??
                    ancestralGene.setsymbetScore(0);
                }
                
                orthologSet.add(ancestralGene);
            }
        }
        
        // At this point all orthologs are added to the set, and their symbet scores are set.
        // Sources of singulars should be found at this time including later-present genes (decided by AGSanalyzer.java)
        // Each singular can be sourced only by genes with a symbet score.
        
        // Reset set1 towards set1 and set2 towards set2 according to symbet scores
        // Find all best-hits of set 1 to set 1
        for (int i=0; i<set1_scores_towards_set1.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set1_scores_towards_set1[i].size(); j++)
            {
                int hit_location = set1_scores_towards_set1_locations[i].get(j);
                if (set1_gene_used[hit_location]==true && set1_scores_towards_set1[i].get(j) > best_score)
                {
                    best_hit_location = hit_location;
                    best_score = set1_scores_towards_set1[i].get(j);
                }
            }
            set1_towards_set1_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set1_towards_set1_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set1_towards_set1_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }        
        
        // Find all best-hits of set 2 to set 2
        for (int i=0; i<set2_scores_towards_set2.length; i++)
        {
            int best_score = this.threshold;
            int best_hit_location = -1;
            
            for (int j=0; j<set2_scores_towards_set2[i].size(); j++)
            {
                int hit_location = set2_scores_towards_set2_locations[i].get(j);
                if (set2_gene_used[hit_location]==true && set2_scores_towards_set2[i].get(j) > best_score)
                {
                    best_hit_location = hit_location;
                    best_score = set2_scores_towards_set2[i].get(j);
                }
            }
            set2_towards_set2_best_scores[i] = best_score;
            if (best_hit_location !=-1)
            {
                set2_towards_set2_best_hit_locations[i] = best_hit_location;
            }
            else
            {
                set2_towards_set2_best_hit_locations[i] = -1; // This represents all hits having 0 score, won't happen since all our genes have PDZ domains
            }
        }        
        
        
        // Find singulars of set1 and set their possible source
        for (int i=0; i<set1_gene_used.length; i++)
        {
            if (set1_gene_used[i] == false)
            {
                OrthologGenes og = new OrthologGenes();
                og.setPossibleSource("Unknown", 0);
                if (set1_towards_set1_best_hit_locations[i] != -1)
                {
                    int realLocation = set1_scores_towards_set1_locations[i].indexOf(set1_towards_set1_best_hit_locations[i]);
                    if (realLocation==-1)
                    {
                        System.out.println("Problem in ANCESTOR DATABASE with real location indexes");
                    }
                    else
                    {
                        og.setPossibleSource(seqSet1.get(set1_towards_set1_best_hit_locations[i]).genes.get(0), set1_scores_towards_set1[i].get(realLocation));
                    }
                }
                og.addOrthologGenes(seqSet1.get(i));
                singularSet1.add(og);
            }
        }
        
        // Find singulars for set2 and set their possible source
        for (int i=0; i<set2_gene_used.length; i++)
        {
            if (set2_gene_used[i] == false)
            {
                OrthologGenes og = new OrthologGenes();
                og.setPossibleSource("Unknown", 0);                
                if (set2_towards_set2_best_hit_locations[i] != -1)
                { 
                    int realLocation = set2_scores_towards_set2_locations[i].indexOf(set2_towards_set2_best_hit_locations[i]);
                    if (realLocation==-1)
                    {
                        System.out.println("Problem in ANCESTOR DATABASE with real location indexes");
                    }
                    else
                    {
                        og.setPossibleSource(seqSet2.get(set2_towards_set2_best_hit_locations[i]).genes.get(0), set2_scores_towards_set2[i].get(realLocation));
                    }
                }
                og.addOrthologGenes(seqSet2.get(i));
                singularSet2.add(og);
            }
        }
        
        // Add all genes to one vector
        ancestralGeneSet.addAll(orthologSet);
        ancestralGeneSet.addAll(singularSet1);
        ancestralGeneSet.addAll(singularSet2);
        
        return new AGS(ancestralGeneSet);
    }    
    private int getScore(Sequence seq1, Sequence seq2)
    {
        int bestScore = 0;
        
        int position1 = seq1.scoreMatrixPosition;
        int position2 = seq2.scoreMatrixPosition;
        bestScore = AGS_console.scoreMatrix.getElement(position1, position2);
      
        return bestScore;
    }
    
    public void writeAncestors(String outputName, boolean CALCULATE_SFA_FILE, String SFA_FILE_NAME)
                                throws IOException, InvalidScoringMatrixException, IncompatibleScoringSchemeException
    {
        FileWriter output = create_file(new File(outputName));
        if (output == null) throw new IOException("Could not write to file!");   

        for (int i=0; i<this.tree.size(); i++)
        {
            Node nextNode = this.tree.getNode(i);
            if (!nextNode.species_node)
            {
                output.write("ANCESTOR\t" + nextNode.name);
                output.write(System.getProperty("line.separator"));
                output.write(nextNode.getAGS().toStringNames());
            }
        }
        output.close();
        
        if (CALCULATE_SFA_FILE)
        {
            System.out.println("Preparing SFA file.");
            // Create output
            FileWriter output2 = create_file(new File(SFA_FILE_NAME));
            if (output2 == null) throw new IOException(); 
           
            if (PROTEIN)
                ALL_INT_ARRAYS = prepareSequenceIndex(GD);               
            System.out.println("SFA database indexes calculated.");
            // Get final ancestor
            Node currentNode = this.tree.root;
            AGS currentAncestor = currentNode.getAGS();   
            
            // For each gene of the final ancestor
            for (int m=0; m<currentAncestor.size(); m++)
            {
                OrthologGenes currentGenes = currentAncestor.AGS.get(m);  
                int bestScore = 0;
                int bestScoreIndex = -1;
                
                // Find representative gene
                for (int n=0; n<currentGenes.genes.size(); n++)
                {
                    String thisGene = currentGenes.genes.get(n);
                    int totalScore = 0;
                    for (int t=0; t<currentGenes.genes.size(); t++)
                    {
                        String otherGene = currentGenes.genes.get(t);
                        totalScore += this.getScore(GD.getGeneByPreComputedIndex(thisGene), GD.getGeneByPreComputedIndex(otherGene), 
                                GD.getGeneLocationByPreComputedIndex(thisGene), GD.getGeneLocationByPreComputedIndex(otherGene));
                    }
                    if (bestScore < totalScore)
                    {
                        bestScore = totalScore;
                        bestScoreIndex = n;
                    }
                }
                output2.write(GD.getGeneByPreComputedIndex(currentGenes.genes.get(bestScoreIndex)).toStringOriginal());
                output2.write(System.getProperty("line.separator"));
            }
            output2.close();
            System.out.println("SFA file written to file.");
        }   
        
    }
    
    private int[][] prepareSequenceIndex(GeneDatabase ALL_GENES)
    {
        int[][] ALL_INT_ARRAYS = new int[ALL_GENES.size()][];
        
        for (int i=0; i<ALL_INT_ARRAYS.length; i++)
        {
            Sequence seq = ALL_GENES.getGene(i);
            ALL_INT_ARRAYS[i] = scoring.getIndexesOfSequences(seq.sequence);
        }
        return ALL_INT_ARRAYS;
    }
    /*
     * Gets two sequences, aligns them and returns the resulting alignment string.
     * Set whether protein or nucleotide alignment.
     * 
     */
    public AlignmentString alignTwoSequences (Sequence seq1, Sequence seq2, int i, int j)
     						throws InvalidScoringMatrixException,
							IncompatibleScoringSchemeException
    {
        String S1, S2;
        S1 = seq1.getSequence();
        S2 = seq2.getSequence();
        
        NeedlemanWunsch nw1;
        if (PROTEIN)
        {
            int[] index1 = ALL_INT_ARRAYS[i]; // Faster alingnments
            int[] index2 = ALL_INT_ARRAYS[j]; // Faster alignments
            nw1 = new NeedlemanWunsch(index1, index2, S1, S2, scoring);
            return nw1.computeFastPairwiseAlignment();
        }
        else
        {
            nw1 = new NeedlemanWunsch(S1, S2, PROTEIN, scoring); // This is only for nucleotide alignments -- without array
            return nw1.computePairwiseAlignment();
        }
            
    }
    
    public int getScore(Sequence seq1, Sequence seq2, int i, int j)
     						throws InvalidScoringMatrixException,
							IncompatibleScoringSchemeException    
    {    
        AlignmentString as = alignTwoSequences(seq1, seq2, i, j);
        return (int)(as.normalized_score*1000);
    }      
    
    public void readAncestors(String mainTreeName, ArrayList<String> subtree_names, ArrayList<String> subtrees)
                        throws IOException
    {
        // First read all subtrees to the tree
        for (int i=0; i<subtree_names.size(); i++)
        {
            readAncestors(subtree_names.get(i)+".ancestors_pass1");
        }
        readAndConvertMainAncestor(mainTreeName, subtree_names, subtrees);
        
    }
    public void readAndConvertMainAncestor(String inputName, ArrayList<String> subtree_names, ArrayList<String> subtrees)
                                throws IOException    
    {
        Node node = null;
        BufferedReader input;
        File inputFile = new File(inputName+".ancestors_pass1");
        input = open_file(inputFile);
        
        if (input==null)
            throw new IOException("Could not read ancestor file " + inputFile);
        
        // Read line by line and fetch every protein into a list	 
        while (more_records(input))
        {
            String line = read_a_line(input);
            if (line.equals(""))
                continue;
            else if (line.startsWith("ANCESTOR"))
            {
                StringTokenizer st = new StringTokenizer(line, "\t");
                st.nextToken();
                String name = st.nextToken();
                
                // Get complete ancestor string
                for (int i=0; i<subtree_names.size(); i++)
                {
                    if (name.indexOf(subtree_names.get(i))!=-1)
                    {
                        name = 
                                name.substring(0, name.indexOf(subtree_names.get(i))) +
                                subtrees.get(i) +
                                name.substring(name.indexOf(subtree_names.get(i))+subtree_names.get(i).length(), name.length());
                    }
                }
                
                node = this.tree.findNode(name);
                if (node == null)
                {
                    System.out.println("Unknown node found: " + name);
                    break;
                }
                // Create AGS file
                node.setAGS(new AGS());
            }
            else if (line.startsWith("SINGULAR"))
            {
                StringTokenizer st = new StringTokenizer(line, "\t");
                String[] geneNames = new String[st.countTokens()-3];
                st.nextToken(); // Skip singular text
                String source = st.nextToken(); // get source
                int singular_score = Integer.parseInt(st.nextToken());
                int count = 0;
                while(st.hasMoreTokens())
                {
                    geneNames[count++] = st.nextToken();
                }
                
                // Now replace all subtree gene names with their real genes
                // Get complete ancestor string
                ArrayList<String> geneNamesFull = new ArrayList<String>();

                for (int j=0; j<geneNames.length; j++)
                {
                    boolean found_subtree = false;
                    for (int i=0; i<subtree_names.size(); i++)
                    {
                        if (geneNames[j].startsWith(subtree_names.get(i)))
                        {
                            found_subtree = true;
                            Node realNode = tree.findNode(subtrees.get(i));
                            int location = Integer.parseInt(geneNames[j].substring(geneNames[j].lastIndexOf("_")+1));
                            OrthologGenes allRelatedGenes = realNode.ags.AGS.get(location);
                            for (int t=0; t<allRelatedGenes.genes.size(); t++)
                            {
                                geneNamesFull.add(allRelatedGenes.genes.get(t));
                            }
                        }
                    }
                    if (!found_subtree)
                        geneNamesFull.add(geneNames[j]);
                }
                
                String[] realGeneNames = new String[geneNamesFull.size()];
                for (int i=0; i<geneNamesFull.size(); i++)
                {
                    realGeneNames[i] = geneNamesFull.get(i);
                }
                
                // Find source
                for (int i=0; i<subtree_names.size(); i++)
                {
                    if (source.startsWith(subtree_names.get(i)))
                    {
                        Node realNode = tree.findNode(subtrees.get(i));
                        int location = Integer.parseInt(source.substring(source.lastIndexOf("_")+1));
                        OrthologGenes allRelatedGenes = realNode.ags.AGS.get(location);
                        source = allRelatedGenes.genes.get(0);
                    }
                }
                
                OrthologGenes og = new OrthologGenes(this.GD.getAllGenesByName2(realGeneNames));
                og.setPossibleSource(source,singular_score);
                node.getAGS().addOrthologGene(og);                
            }
            else if (line.startsWith("SYM-BET"))// regular ortholog genes line
            {
                StringTokenizer st = new StringTokenizer(line, "\t");
                String[] geneNames = new String[st.countTokens()-2];
                st.nextToken(); // Skip sym-bet text
                int symbet_score = Integer.parseInt(st.nextToken());
                int count = 0;
                while(st.hasMoreTokens())
                {
                    geneNames[count++] = st.nextToken();
                }
                
                // Now replace all subtree gene names with their real genes
                // Get complete ancestor string
                ArrayList<String> geneNamesFull = new ArrayList<String>();

                for (int j=0; j<geneNames.length; j++)
                {
                    boolean found_subtree = false;
                    for (int i=0; i<subtree_names.size(); i++)
                    {
                        if (geneNames[j].startsWith(subtree_names.get(i)))
                        {
                            found_subtree = true;
                            Node realNode = tree.findNode(subtrees.get(i));
                            int location = Integer.parseInt(geneNames[j].substring(geneNames[j].lastIndexOf("_")+1));
                            OrthologGenes allRelatedGenes = realNode.ags.AGS.get(location);
                            for (int t=0; t<allRelatedGenes.genes.size(); t++)
                            {
                                geneNamesFull.add(allRelatedGenes.genes.get(t));
                            }
                        }
                    }
                    if (!found_subtree)
                        geneNamesFull.add(geneNames[j]);
                }
                
                String[] realGeneNames = new String[geneNamesFull.size()];
                for (int i=0; i<geneNamesFull.size(); i++)
                {
                    realGeneNames[i] = geneNamesFull.get(i);
                }                
                
                OrthologGenes og = new OrthologGenes(this.GD.getAllGenesByName2(realGeneNames));
                og.setsymbetScore(symbet_score);
                node.getAGS().addOrthologGene(og);
            }
            
        }        
    }
    public void readAncestors(String inputName)
                                throws IOException
    {
        Node node = null;
        BufferedReader input;
        File inputFile = new File(inputName);
        input = open_file(inputFile);
        
        if (input==null)
            throw new IOException("Could not read ancestor file " + inputFile);
        
        // Read line by line and fetch every protein into a list	 
        while (more_records(input))
        {
            String line = read_a_line(input);
            
            if (line.equals(""))
                continue;
            else if (line.startsWith("ANCESTOR"))
            {
                StringTokenizer st = new StringTokenizer(line, "\t");
                st.nextToken();
                String name = st.nextToken();
                node = this.tree.findNode(name);
                if (node == null)
                {
                    System.out.println("Unknown node found: " + line);
                    break;
                }
                // Create AGS file
                node.setAGS(new AGS());
            }
            else if (line.startsWith("SINGULAR"))
            {
                StringTokenizer st = new StringTokenizer(line, "\t");
                String[] geneNames = new String[st.countTokens()-3];
                st.nextToken(); // Skip singular text
                String source = st.nextToken(); // get source
                int singular_score = Integer.parseInt(st.nextToken());
                int count = 0;
                while(st.hasMoreTokens())
                {
                    geneNames[count++] = st.nextToken();
                }
                OrthologGenes og = new OrthologGenes(this.GD.getAllGenesByName2(geneNames));
                og.setPossibleSource(source,singular_score);
                node.getAGS().addOrthologGene(og);                
            }
            else if (line.startsWith("SYM-BET"))// regular ortholog genes line
            {
                StringTokenizer st = new StringTokenizer(line, "\t");
                String[] geneNames = new String[st.countTokens()-2];
                st.nextToken(); // Skip sym-bet text
                int symbet_score = Integer.parseInt(st.nextToken());
                int count = 0;
                while(st.hasMoreTokens())
                {
                    geneNames[count++] = st.nextToken();
                }
                OrthologGenes og = new OrthologGenes(this.GD.getAllGenesByName2(geneNames));
                og.setsymbetScore(symbet_score);
                node.getAGS().addOrthologGene(og);
            }
            
        }
        
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
