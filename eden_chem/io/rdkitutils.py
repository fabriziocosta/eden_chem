'''
interface between networkx graphs and the world of chemistry via rdkit

functionality:
# gernerating networkx graphs:
sdf_to_nx(file.sdf)
smi_to_nx(file.smi)
smiles_to_nx(smilesstringlist)



#bonus: garden style transformer.
class MoleculeToGraph
'''

import logging
from rdkit import Chem

import networkx as nx
from sklearn.base import BaseEstimator, TransformerMixin

logger = logging.getLogger(__name__)


class MoleculeToGraph(BaseEstimator, TransformerMixin):

    def __init__(self, file_format='sdf'):
        """Constructor.

        valid 'file_format' strings and what the transformer will expect

        smi: path to .smi file
        sdf: pat to .sdf file
        """
        self.file_format = file_format

    def transform(self, data):
        """Transform."""
        try:
            if self.file_format == 'smi':
                graphs = smi_to_nx(data)
            elif self.file_format == 'sdf':
                graphs = sdf_to_nx(data)
            else:
                raise Exception('file_format must be smi or sdf')
            for graph in graphs:
                yield graph

        except Exception as e:
            logger.debug('Failed iteration. Reason: %s' % e)
            logger.debug('Exception', exc_info=True)


################
# import to networkx graphs
###############

def sdf_to_nx(file):
    # read sdf file
    suppl = Chem.SDMolSupplier(file)
    for mol in suppl:
        if mol:
            yield rdkmol_to_nx(mol)


def smi_to_nx(file):
    # read smi file
    suppl = Chem.SmilesMolSupplier(file)
    for mol in suppl:
        yield rdkmol_to_nx(mol)


def rdkmol_to_nx(mol):
    #  rdkit-mol object to nx.graph
    graph = nx.Graph()
    for e in mol.GetAtoms():
        graph.add_node(e.GetIdx(), label=e.GetSymbol())
    for b in mol.GetBonds():
        graph.add_edge(b.GetBeginAtomIdx(), b.GetEndAtomIdx(),
                       label=str(int(b.GetBondTypeAsDouble())))
    return graph


def smiles_strings_to_nx(smileslist):
    # smiles strings
    for smile in smileslist:
        mol = Chem.MolFromSmiles(smile)
        yield rdkmol_to_nx(mol)


################
# exporting networkx graphs
###############

def nx_to_smi(graphs, file=None):
    # writes smiles strings to a file
    chem = [nx_to_rdkit(graph) for graph in graphs]
    smis = [Chem.MolToSmiles(m) for m in chem]
    if file:
        with open(file, 'w') as f:
            f.write('\n'.join(smis))
    return smis


def nx_to_rdkit(graph):
    m = Chem.MolFromSmiles('')
    mw = Chem.RWMol(m)
    atom_index = {}
    for n, d in graph.nodes(data=True):
        atom_index[n] = mw.AddAtom(Chem.Atom(d['label']))
    for a, b, d in graph.edges(data=True):
        start = atom_index[a]
        end = atom_index[b]
        bond_type = d.get("label", '1')
        if bond_type == '1':
            mw.AddBond(start, end, Chem.BondType.SINGLE)
        elif bond_type == '2':
            mw.AddBond(start, end, Chem.BondType.DOUBLE)
        elif bond_type == '3':
            mw.AddBond(start, end, Chem.BondType.TRIPLE)
        # more options:
        # http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.BondType-class.html
        else:
            raise Exception('bond type not implemented')

    mol = mw.GetMol()
    return mol
