package uni.tubingen.inference.proteinlasso;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DoubleValue;
import org.knime.core.data.StringValue;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnNameSelection;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelDouble;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.ColumnFilter;


/**
 * <code>NodeDialog</code> for the "ProteinLasso" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author enrique
 */
public class ProteinLassoNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the ProteinLasso node.
     */
    protected ProteinLassoNodeDialog() {
    	      super ();
    	        
    	        //fields to match with coming table...
    	        final SettingsModelString matches_peptides = new SettingsModelString(ProteinLassoNodeModel.CFGKEY_PEPTIDES, "Peptides");
                final SettingsModelString accsn_protein    = new SettingsModelString(ProteinLassoNodeModel.CFGKEY_PROTEIN, "Protein");
                final SettingsModelString probabilities    = new SettingsModelString(ProteinLassoNodeModel.CFGKEY_PROBABILITIES, "Probabilities");
                final SettingsModelString detectability    = new SettingsModelString(ProteinLassoNodeModel.CFGKEY_DETECTABILITY, "Detectabilities");
            
            	final SettingsModelDouble lammda           = new SettingsModelDouble(ProteinLassoNodeModel.CFGKEY_LAMBDA_PARAMETER, 0.5);
            	
            	
     
                addDialogComponent(new DialogComponentColumnNameSelection(accsn_protein, "Proteins Column", 0, true, StringValue.class));
                addDialogComponent(new DialogComponentColumnNameSelection(matches_peptides, "Peptides Column", 0, true, 
                		           
              		               new ColumnFilter() {

        		                	@Override
        							public boolean includeColumn(DataColumnSpec colSpec) {
        								if (colSpec.getType().isCollectionType() && colSpec.getType().getCollectionElementType().isCompatible(StringValue.class))
        									return true;
        								
        								if (colSpec.getType().isCompatible(StringValue.class)) 
        									return true;
        								
        								return false;
        							}

        							@Override
        							public String allFilteredMsg() {
        								return "No suitable columns (string or List/Set column) to select!";
        							}
                			
                		}));
                 addDialogComponent(new DialogComponentColumnNameSelection(probabilities, "Probabilities", 0, true, DoubleValue.class));
                 addDialogComponent(new DialogComponentColumnNameSelection(detectability, "Detectability", 0, true, DoubleValue.class));
                
                 addDialogComponent(new DialogComponentNumber(new SettingsModelDoubleBounded(ProteinLassoNodeModel.CFGKEY_LAMBDA_PARAMETER, 0.5, 0.0, 1.0),  "lambda",  0.01));

    }
}

