

package EvolMAP;
import java.io.*;
import java.util.*;

/**
 *
 * @author onur
 */
public class NewickTree 
{
    Vector<Node> allNodes;
    Node root = null;
    
    public Node getRoot()
    {
        return this.root;
    }
    
    public int size()
    {
        return this.allNodes.size();
    }
    
    public Node getNode(int i)
    {
        return this.allNodes.get(i);
    }
    
    // Copy tree -- if you like to recreate outside and add in
    public NewickTree(Node root, Vector<Node> allNodes)
    {
        this.allNodes = new Vector<Node>();
        this.allNodes.addAll(allNodes);
        this.root = root;
    }
    
    /** Creates a new instance of NewickTree */
    public NewickTree(String treeString) 
    {
        treeString = treeString.replaceAll(" ","");
        
        this.allNodes = new Vector<Node>();
        
        // Remove semicolon if there is one in the end
        if (treeString.charAt(treeString.length()-1) == ';')
            treeString = treeString.substring(0, treeString.length()-1);
        
        Vector<String> parseStack = new Vector<String>();
        int pointer = 0;
        String currentRead = "";
        
        // Read string and partion into stack
        while (pointer<treeString.length())
        {
            if (treeString.charAt(pointer) == '(')
            {
                if (currentRead == "")
                {
                    parseStack.add("(");
                }
                else
                {
                    System.out.println("ERROR: characters before (");
                }
            }
            else if (treeString.charAt(pointer) == ')')
            {
                if (currentRead == "")
                {
                    parseStack.add(")");
                }
                else
                {
                    parseStack.add(currentRead);
                    this.allNodes.add(new Node(currentRead, null, null));
                    currentRead = "";
                    parseStack.add(")");                    
                }                
            }
            else if (treeString.charAt(pointer) == ',')
            {
                if (currentRead != "")
                {
                    parseStack.add(currentRead);
                    this.allNodes.add(new Node(currentRead, null, null));
                    currentRead = "";
                    parseStack.add(",");
                }
                else
                {
                    parseStack.add(",");
                }
            }
            else
            {
                currentRead += treeString.charAt(pointer);
                
            }
            pointer++;
        }
        
        // Read stack and pop out to create nodes, algorithm below
        // Once you see a "(", you need to pop out till you see a ")"
        Vector <String> popStack = new Vector<String>();
        while(parseStack.size()>0)
        {
            String nextString = parseStack.remove(parseStack.size()-1);
            popStack.add(nextString);
            if (nextString.equals("("))
            {
                String newNode = "";
                boolean stop = false;
                while (popStack.size()>0 && !stop)
                {
                    String nextPoppedString = popStack.remove(popStack.size()-1);
                    newNode += nextPoppedString;
                    if (nextPoppedString.equals(")"))
                        break;
                }
                Node node = addNode(newNode, allNodes);
                popStack.add(node.name);
            }
        }
        this.root = allNodes.get(allNodes.size()-1);
        
    }
    public Node addNode(String newNode, Vector<Node> allNodes)
    {
        // Find first and second nodes, remove first and last bracket
        Node descendant1 = this.getNode(newNode.substring(1, newNode.indexOf(',')), allNodes);
        Node descendant2 = this.getNode(newNode.substring(newNode.indexOf(',')+1, newNode.length()-1), allNodes);
        
        Node ancestor = new Node(newNode, descendant1, descendant2);
        descendant1.setAncestor(ancestor);
        descendant2.setAncestor(ancestor);
        allNodes.add(ancestor);
        return ancestor;
    }
    
    // Returns null if node not found -- used when creating
    public Node getNode(String name, Vector<Node> allNodes)
    {
        // Create name
        name = name.replaceAll(",", "_");
        String noParenthesisName = "";
        for (int i=0; i<name.length(); i++)
        {
            if (name.charAt(i)!=')' && name.charAt(i)!='(')
            {
                noParenthesisName += name.charAt(i);
            }
        }        
        for (int i=0; i<allNodes.size(); i++)
        {
            if (noParenthesisName.equals(allNodes.get(i).name))
            {
                return allNodes.get(i);
            }
        }
        return null;
    }  
    public String[] getAllSpeciesNames()
    {
        Vector<String> speciesNames = new Vector<String>();
        
        for (int i=0; i<this.allNodes.size(); i++)
        {
            Node thisNode = allNodes.get(i);
            if (thisNode.species_node)
                speciesNames.add(thisNode.name);
        }
        String[] speciesNamesArray = new String[speciesNames.size()];
        for (int i=0; i<speciesNames.size(); i++)
        {
            speciesNamesArray[i] = speciesNames.get(i);
        }
        return speciesNamesArray;
    }
    // Finds and returns a null by name
    public Node findNode(String name)
    {
        for (int i=0; i<this.allNodes.size(); i++)
        {
            if (allNodes.get(i).name.equals(name))
                return this.allNodes.get(i);
        }
        return null;
    }
}
