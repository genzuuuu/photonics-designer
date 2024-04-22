# Author: Marko Huang
import torch
from torch import nn
from torch.nn import functional as F


# class PauliEmbedder(nn.Module):
#     def __init__(self, o, hs):
#         super().__init__()
#         # Basically: self.pauli_embedder = nn.Embedding(o, hs) 
#         self.o = o
#         self.hs = hs
#         dummy = torch.zeros(1, hs, requires_grad=False)  # Zeros for padding
#         pauli_emb = nn.Parameter(torch.rand(o, hs), requires_grad=True) # Last row for instructions
#         self.embedding = torch.cat((dummy, pauli_emb), 0)
    
#     def forward(self, x):
#         """
#         changes from:
#         tensor([[0, 0, 0, 1],
#                 [0, 1, 1, 0],
#                 [1, 0, 0, 1]])
#         to:
#         tensor([[0, 0, 0, 4],
#                 [0, 2, 3, 0],
#                 [1, 0, 0, 4]])
#         and outputs something like:
#         tensor([[[0.0000, 0.0000],
#                  [0.0000, 0.0000],
#                  [0.0000, 0.0000],
#                  [0.5569, 0.6218]],

#                 [[0.0000, 0.0000],
#                  [0.1003, 0.4569],
#                  [0.3822, 0.6429],
#                  [0.0000, 0.0000]],

#                 [[0.4562, 0.5808],
#                  [0.0000, 0.0000],
#                  [0.0000, 0.0000],
#                  [0.5569, 0.6218]]])
#         and sums over dim=1
#         """
#         p_ops = x.clone()
#         idx = p_ops.nonzero(as_tuple=True)
#         p_ops[*idx] = idx[1].type(torch.uint8)+1
#         out = F.embedding(p_ops.reshape(-1).int(), self.embedding)\
#             .reshape(-1, self.o, self.hs).sum(dim=1)
#         return out

NUM_ITERATIONS=5

# MantillaNetv1 (Named by Marko - March 10, 2024)
# MarkoNetv1 (Re-named by Luis - March 11, 2024)
class MarkoNetv1(nn.Module):
    def __init__(self, o, n, hs, nh):
        super().__init__()
        # self.pauli_embedder = PauliEmbedder(o, hs)
        self.o, self.n, self.hs, self.nh = o, n, hs, nh
        self.pauli_embedder = nn.Embedding(o, hs) 
        self.p_attn = nn.MultiheadAttention(hs, nh, batch_first=True)
        self.pg_attn = nn.MultiheadAttention(hs, nh, batch_first=True)
        self.g_attn = nn.MultiheadAttention(hs, nh, batch_first=True, add_bias_kv=True)
        self.init_instr = nn.Parameter(torch.randn(hs), requires_grad=True)
        self.register_buffer('init_graph', torch.ones(self.n**2, self.hs))
        self.f_out = nn.Linear(hs, 1)

    def forward(self, p_ops):
        # initialize pauli embeddings with instruction token
        p_embs = [self.init_instr] # add instructions token at the start
        for p_op in p_ops:
            p_embs.append(self.pauli_embedder(torch.nonzero(p_op, as_tuple=False).flatten()).sum(dim=0))
        p_embs = torch.stack(p_embs)[None, :]
        # initialize empty graph
        # init_graph = torch.zeros(self.n**2, self.hs)
        curr_graph = self.init_graph
        
        for _ in range(NUM_ITERATIONS):
            p_attn_emb, _ = self.p_attn(p_embs, p_embs, p_embs)
            graph_emb = self.g_attn(curr_graph, curr_graph, curr_graph)[0]
            instr_emb, _ = self.pg_attn(p_attn_emb[None, :, 0], graph_emb[None, :], graph_emb[None, :])
            
            # feed instr_emb to current iteration of graph_emb
            new_graph = torch.empty(self.n**2+1, self.hs)
            new_graph[0, :] = instr_emb.squeeze()
            new_graph[1:, :] = graph_emb
            new_graph, _ = self.g_attn(new_graph[None, :, :], new_graph[None, :, :], new_graph[None, :, :])
            new_graph = new_graph.squeeze()[1:] # keep just the graph part
            
            # update p_embs
            new_p_embs = torch.empty_like(p_embs)
            new_p_embs[0, :] = instr_emb
            new_p_embs[1:, :] = p_embs
            p_embs = new_p_embs
            curr_graph = new_graph
            
        
        # Use Set Cross-Entropy maybe?
        # https://arxiv.org/pdf/1812.01217.pdf
        return self.f_out(new_graph)


if __name__ == "__main__":

    from gen_dataset import MBQCDataset

    dt = torch.load('data/dataset.pt') # List[Tuple[input, output]]
    # for our dataset
    n = dt[0][1].shape[0] # 5: number of nodes
    o = dt[0][0].shape[1] # 4: 2*number of output qubits
    m = dt[0][0].shape[0] # 3 here but may vary, is the number of pauli operator the that the graph implements

    hs = 8 # each pauli element is linear combination of 4 embeddings
    hhs = 16 # attention hidden size
    nh = 8 # number of attention heads
    # nl = 8 # number of layers
    
    sample = torch.tensor(dt[0][0])
    model = MarkoNetv1(o, n, hs, nh)
        
    print(model(sample).reshape(n,n)) # Output needs to be reshaped to n x n

    import torch.optim as optim
    from torch.utils.data import DataLoader

    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(model.parameters())

    loader = DataLoader(dt, batch_size=1, shuffle=True)

    num_epochs = 100
    for epoch in range(num_epochs):  # num_epochs should be defined by you
        for inputs, targets in loader:
            optimizer.zero_grad()
            outputs = model(inputs)
            targets = targets.reshape(-1,1)
            loss = criterion(outputs, targets)
            loss.backward()
            optimizer.step()

        if (epoch+1) % 10 == 0:
            print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item()}')
    