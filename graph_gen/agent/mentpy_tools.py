import mentpy as mp


"""
# MBQC Circuits
"""
def make_mbqc_circuit(graph, input_qubits, output_qubits):
    """
    The MBQCircuit class that deals with operations and manipulations of graph states

    Parameters:
    graph: mp.GraphState
        The graph state of the MBQC circuit.
    input_nodes: list
        The input nodes of the MBQC circuit.
    output_nodes: list
        The output nodes of the MBQC circuit.

    Examples:
        g = mp.GraphState()
        g.add_edges_from([(0,1), (1,2), (2,3), (3, 4)])
        state = make_mbqc_circuit(g, input_nodes=[0], output_nodes=[4])
    """
    return mp.MBQCircuit(graph, input_qubits, output_qubits)

def calculate_lie_algebra(mbqc_circuit):
    """
    Calculate the Lie algebra of the MBQC circuit

    Parameters
    ----------
    mbqc_circuit: mp.MBQCircuit
        The MBQC circuit.
    """
    return mp.utils.lie_algebra.calculate_gens_lie_algebra(mbqc_circuit)

