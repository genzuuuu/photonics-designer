> [!NOTE]
> Module under development.

## Graph structure generation

In this section, we show to approaches to generate graph states for quantum computing given a set of desired gates to implement.

The first approach is a deterministic algorithm based on [this work](https://dx.doi.org/10.14288/1.0435698). To run the algorithm, run `graph_gen/deterministic.py`. The algorithm will generate a graph state for a given set of gates.

The second approach is a machine learning algorithm. We train a trasformer-based model to generate the adjacency matrix of a graph state given a set of desired gates to implement.

Moreover, we identify viable graph states following the approach described by [JC Adcock in his PhD thesis](https://research-information.bris.ac.uk/ws/portalfiles/portal/202643120/Jeremy_C._Adcock_Ph.D_thesis_Generating_Optical_Graph_States_v2.pdf)

### Installation

To install the required packages, run `pip install -r requirements.txt`.


