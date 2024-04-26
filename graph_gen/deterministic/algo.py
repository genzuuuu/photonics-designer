import mentpy as mp

def generate_circuit(pauli_op):

    circuit = None
    for op in pauli_op:
        if circuit is None:
            circuit = mp.templates.from_pauli(op)
        else:
            circuit = mp.hstack((circuit, mp.templates.from_pauli(op)))
    
    return circuit

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    op = mp.PauliOp('XIZ;IZX;ZXX')
    circuit = generate_circuit(op)
    mp.draw(circuit, figsize=(10, 10))  
    plt.show()

    print(circuit.input_nodes, circuit.output_nodes)
