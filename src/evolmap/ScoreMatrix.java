
package EvolMAP;

/**
 *
 * @author onur
 */
public abstract class ScoreMatrix 
{
    String[] names; // column names
    public abstract int getElement(int i, int j);    
}