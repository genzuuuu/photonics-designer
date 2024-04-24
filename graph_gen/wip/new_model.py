import torch
import torch.nn as nn
import torch_geometric.nn as geom_nn
import torch_geometric.data as geom_data

max_n = 100  # Maximum size for n
max_m = 50   # Maximum size for m

class GraphNetwork(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_dim=64):
        super(GraphNetwork, self).__init__()
        self.conv1 = geom_nn.GCNConv(input_dim, hidden_dim)
        self.conv2 = geom_nn.GCNConv(hidden_dim, hidden_dim)
        self.linear = nn.Linear(hidden_dim, output_dim)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = torch.relu(self.conv1(x, edge_index))
        x = torch.relu(self.conv2(x, edge_index))
        x = self.linear(x)
        return x

def pad_matrix(matrix, max_rows, max_cols):
    padded = torch.zeros((max_rows, max_cols))
    rows, cols = matrix.shape
    padded[:rows, :cols] = matrix
    return padded

def create_graph_data(lie_alg_matrix, adj_matrix, max_n, max_m):
    lie_alg_matrix_padded = pad_matrix(lie_alg_matrix, max_n, max_m)
    x = torch.tensor(lie_alg_matrix_padded, dtype=torch.float)
    edge_index = torch.nonzero(adj_matrix, as_tuple=False).t().contiguous()
    data = geom_data.Data(x=x, edge_index=edge_index)
    return data


if __name__ == "__main__":

    from gen_dataset import MBQCDataset

    dt = torch.load('data/dataset.pt')

    # Example usage
    input_dim = max_n * max_m
    output_dim = max_n * max_n  # Assuming square adjacency matrix

    lie_alg_matrix = torch.rand((10, 20))  # Example smaller matrix
    adj_matrix = torch.randint(0, 2, (10, 10))  # Example adjacency matrix
    data = create_graph_data(lie_alg_matrix, adj_matrix, max_n, max_m)


    n = dt[0][1].shape[0] # 5: number of nodes
    o = dt[0][0].shape[1] # 4: 2*number of output qubits
    m = dt[0][0].shape[0] # 3 here but may vary, is the number of pauli operator the that the graph implements

    model = GraphNetwork(input_dim=input_dim, output_dim=num_nodes**2)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    criterion = nn.MSELoss()

    for epoch in range(100):
        model.train()
        optimizer.zero_grad()
        output = model(data)
        loss = criterion(output, target_adj_matrix.view(-1))
        loss.backward()
        optimizer.step()

        if epoch % 10 == 0:
            print(f'Epoch {epoch}: Loss {loss.item()}')
