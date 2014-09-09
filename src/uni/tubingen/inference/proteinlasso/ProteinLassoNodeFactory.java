package uni.tubingen.inference.proteinlasso;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "ProteinLasso" Node.
 * 
 *
 * @author enrique
 */
public class ProteinLassoNodeFactory 
        extends NodeFactory<ProteinLassoNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ProteinLassoNodeModel createNodeModel() {
        return new ProteinLassoNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 0;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<ProteinLassoNodeModel> createNodeView(final int viewIndex,
            final ProteinLassoNodeModel nodeModel) {
        return new ProteinLassoNodeView(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new ProteinLassoNodeDialog();
    }

}

