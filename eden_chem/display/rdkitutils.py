from rdkit import Chem
from rdkit.Chem import AllChem, Draw

import networkx as nx
from eden_chem.io.rdkitutils import nx_to_rdkit


def set_coordinates(compounds):
    for m in compounds:
        if m:
            # updateprops fixes "RuntimeError: Pre-condition Violation"
            m.UpdatePropertyCache(strict=False)
            AllChem.Compute2DCoords(m)
        else:
            raise Exception('''set coordinates failed..''')


def get_smiles_strings(graphs):
    compounds = map(nx_to_rdkit, graphs)
    return list(map(Chem.MolToSmiles, compounds))


def nx_to_image(graphs, n_graphs_per_line=5, size=250, title_key=None, titles=None):
    # we want a list of graphs
    if isinstance(graphs, nx.Graph):
        raise Exception("give me a list of graphs")
    # make molecule objects
    compounds = list(map(nx_to_rdkit, graphs))
    # print compounds

    # take care of the subtitle of each graph
    if title_key:
        legend = [g.graph.get(title_key, 'N/A') for g in graphs]
    elif titles:
        legend = titles
    else:
        legend = [str(i) for i in range(len(graphs))]
    return compounds_to_image(compounds, n_graphs_per_line=n_graphs_per_line, size=size, legend=legend)


def compounds_to_image(compounds, n_graphs_per_line=5, size=250, legend=None):
    # calculate coordinates:
    set_coordinates(compounds)
    # make the image
    return Draw.MolsToGridImage(compounds, molsPerRow=n_graphs_per_line, subImgSize=(size, size), legends=legend)
