package uni.tubingen.inference.proteinlasso;

import java.io.File;
import java.io.IOException;

import org.knime.core.data.DataColumnSpec;
import org.knime.core.data.DataColumnSpecCreator;
import org.knime.core.data.DataTableSpec;
import org.knime.core.data.def.StringCell;
import org.knime.core.node.BufferedDataContainer;
import org.knime.core.node.BufferedDataTable;
import org.knime.core.node.CanceledExecutionException;
import org.knime.core.node.ExecutionContext;
import org.knime.core.node.ExecutionMonitor;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeLogger;
import org.knime.core.node.NodeModel;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * This is the model implementation of ProteinLasso.
 * 
 *
 * @author enrique
 */
public class ProteinLassoNodeModel extends NodeModel {
	
	private static final NodeLogger logger = NodeLogger.getLogger("ProteinLasso probabilities");
    
    static String CFGKEY_PEPTIDES = "peptides";
	static String CFGKEY_PROTEIN = "protein";
	static String CFGKEY_PROBABILITIES = "probabilities";
	static String CFGKEY_DETECTABILITY = "detectability";
	static String CFGKEY_LAMMDA_PARAMETER = "lammda";
	
	//fields to link execute variable with input variable...
	private final SettingsModelString m_peptide_column = new SettingsModelString(CFGKEY_PEPTIDES, "Peptides");
	private final SettingsModelString m_protein_column   = new SettingsModelString(CFGKEY_PROTEIN, "Protein");
	private final SettingsModelString m_probability_column   = new SettingsModelString(CFGKEY_PROBABILITIES, "Probabilities");
	private final SettingsModelString m_detectability_column   = new SettingsModelString(CFGKEY_PROBABILITIES, "Detectability");
	private final SettingsModelDoubleBounded lammda_parameter = new SettingsModelDoubleBounded(CFGKEY_LAMMDA_PARAMETER, 0.9, 0.0, 1.0);
	
	//fields to manage the input table...
	static int pep_idx    = 0;
	static int accsn_idx  = 0;
	static int proba_idx  = 0;
	static int detect_idx = 0;
	
	/**
     * Constructor for the node model.
     */
    protected ProteinLassoNodeModel() {
    
        // TODO: Specify the amount of input and output ports needed.
        super(1, 1);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected BufferedDataTable[] execute(final BufferedDataTable[] inData,
            final ExecutionContext exec) throws Exception {
    	
         this.checkTableConfiguraion(inData);
         
         DataTableSpec new_spec_table = new DataTableSpec(make_output_spec());  	
         BufferedDataContainer container = exec.createDataContainer(new_spec_table);
         
         int K= 100; //We choose 100 points between the maximal lamda value and the minimal lamda value.
		 double DECAY = 0.001;
		 double startTime = System.currentTimeMillis();

		 ProteinLasso prolas = new ProteinLasso();
		
	   	 prolas.buildDetectabilityFile(inData[0]);
		 prolas.buildPeptideFile(inData[0]);
		 
		 String tag="max";
		 prolas.getcoef(tag);
	     int totalProteins = prolas.get_proteinNum();
		 double lamda_max = prolas.get_Lamda_max();
		 double lamda_min = lamda_max*DECAY; 
		 System.out.println("------------------------------------------");
		 System.out.println("The maximal value of lamda="+lamda_max);
		 System.out.println("The minimal value of lamda="+lamda_min);
		 
		 
		 double lamda = lamda_max;	 
	     double[] result = new double[totalProteins];
		 double[] coef = new double[totalProteins];
		 for(int j= 0; j<totalProteins; j++){
			 result[j]=0;
			 coef[j]=0;
		 }
			
		 int i= 0;
		 while(i<K){	
			lamda = Math.log(lamda_max)-((double)i*(Math.log(lamda_max)-Math.log(lamda_min)))/(double)K;
			lamda = Math.pow(Math.E,lamda)*0.5;
			//System.out.println("lamda="+lamda);
			result = prolas.Coordinate_Descent(result, lamda);
	    	for(int j= 0; j<totalProteins; j++){
	    		coef[j] = coef[j] + result[j];
	    	}
	    	i++;	 
	    }
			
		 for(int j= 0; j<totalProteins; j++){
				coef[j]=coef[j]/100;
		 }	 
		 
		 prolas.writeContainer(container, coef);
		 
	     double endTime = System.currentTimeMillis();
		 double running_time = (endTime-startTime)/(double)1000;	
		 System.out.println("Running Time:"+running_time);
         
    
         container.close();

        // TODO: Return a BufferedDataTable for each output port 
        return new BufferedDataTable[]{ container.getTable() };
    }
    
  
    //method for checking the table configuration coming...
    private void checkTableConfiguraion(BufferedDataTable[] inData) throws Exception{
    	
    	pep_idx  = inData[0].getDataTableSpec().findColumnIndex(m_peptide_column.getStringValue());
    	accsn_idx= inData[0].getDataTableSpec().findColumnIndex(m_protein_column.getStringValue());
    	proba_idx= inData[0].getDataTableSpec().findColumnIndex(m_probability_column.getStringValue());
    	detect_idx= inData[0].getDataTableSpec().findColumnIndex(m_detectability_column.getStringValue());
    	
    	
    	if (pep_idx < 0 || accsn_idx < 0 || proba_idx < 0 || detect_idx < 0 || pep_idx == accsn_idx ) {
    		throw new Exception("Illegal columns: "+m_peptide_column+" "+m_protein_column+" "+m_probability_column+" "+m_detectability_column+", re-configure the node!");
    	}
    }
    
 
    private DataColumnSpec[]  make_output_spec() {
    	
    	DataColumnSpec cols[] = new DataColumnSpec[2];
    	cols[0] = new DataColumnSpecCreator("Protein ID", StringCell.TYPE).createSpec();
    	cols[1] = new DataColumnSpecCreator("ProteinLasso Probability", StringCell.TYPE).createSpec();
	
      return cols;
	}
    

    /**
     * {@inheritDoc}
     */
    @Override
    protected void reset() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected DataTableSpec[] configure(final DataTableSpec[] inSpecs)
            throws InvalidSettingsException {

        // TODO: generated method stub
    	return new DataTableSpec[]{new DataTableSpec(this.make_output_spec())};
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveSettingsTo(final NodeSettingsWO settings) {
    	
    	m_peptide_column.saveSettingsTo(settings);
        m_protein_column.saveSettingsTo(settings);
        m_probability_column.saveSettingsTo(settings);
        m_detectability_column.saveSettingsTo(settings);       
        lammda_parameter.saveSettingsTo(settings);
         // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadValidatedSettingsFrom(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        // TODO: generated method stub
    	m_peptide_column.loadSettingsFrom(settings);
        m_protein_column.loadSettingsFrom(settings);
        m_probability_column.loadSettingsFrom(settings);
        m_detectability_column.loadSettingsFrom(settings);       
        lammda_parameter.loadSettingsFrom(settings);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void validateSettings(final NodeSettingsRO settings)
            throws InvalidSettingsException {
        // TODO: generated method stub
    	m_peptide_column.validateSettings(settings);
        m_protein_column.validateSettings(settings);
        m_probability_column.validateSettings(settings);
        m_detectability_column.validateSettings(settings);       
        lammda_parameter.validateSettings(settings);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void loadInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected void saveInternals(final File internDir,
            final ExecutionMonitor exec) throws IOException,
            CanceledExecutionException {
        // TODO: generated method stub
    }

}

