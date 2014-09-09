package uni.tubingen.inference.proteinlasso;

import org.knime.core.node.NodeView;


/**
 * <code>NodeView</code> for the "ProteinLasso" Node.
 * 
 *
 * @author enrique
 */
public class ProteinLassoNodeView extends NodeView<ProteinLassoNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link ProteinLassoNodeModel})
     */
    protected ProteinLassoNodeView(final ProteinLassoNodeModel nodeModel) {
        super(nodeModel);
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() {
        // TODO: generated method stub
    	ProteinLassoNodeModel nodeModel = 
                getNodeModel();
            assert nodeModel != null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {
        // TODO: generated method stub
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {
        // TODO: generated method stub
    }

}

