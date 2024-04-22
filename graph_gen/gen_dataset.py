import torch
from torch.utils.data import Dataset
import numpy as np
import networkx as nx
import mentpy as mp

class MBQCDataset(Dataset):
    def __init__(self, n_nodes, input_nodes, output_nodes):

        flowspace = mp.utils.FlowSpace(n_nodes, input_nodes, output_nodes)

        self.data = []

        # Iterate through the nodes in the flowspace's flow graph
        for node in flowspace.flow_graph_space.nodes:
            node_data = flowspace.flow_graph_space.nodes[node]

            # Check if the mbqc_circuit attribute is not None
            if node_data['mbqc_circuit'] is not None:
                try:
                    # Calculate the Lie algebra
                    lie_algebra = mp.utils.lie_algebra.calculate_gens_lie_algebra(node_data['mbqc_circuit'])
                    
                    # Append the adjacency matrix and the corresponding Lie algebra matrix to the dataset
                    adjacency_matrix = nx.adjacency_matrix(node_data['mbqc_circuit'].graph).todense()
                    lie_algebra_matrix = np.array(lie_algebra.matrix)
                    
                    self.data.append((lie_algebra_matrix, adjacency_matrix))
                except:
                    # If an error occurs, skip this node
                    pass

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        input_matrix, target_matrix = self.data[idx]

        # Convert the numpy matrices to torch tensors
        input_tensor = torch.tensor(input_matrix, dtype=torch.float)
        target_tensor = torch.tensor(target_matrix, dtype=torch.float) 

        return input_tensor, target_tensor

if __name__ == "__main__":
    # Create dataset
    dataset = MBQCDataset(5, [0,1], [3,4])

    # Save the dataset
    torch.save(dataset, 'data/dataset.pt')
