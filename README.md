# TCM
Time-Constrained Continuous Subgraph Matching Using Temporal Information for Filtering and Backtracking

## Usages
```bash
./TCM <data_graph_file> <query_graph_file> <window_size>
```

## Experiments in the paper
- Each script file in `scripts` directory reproduces the results of the experiments in the paper (Figures 6~8).
- The results will be stored in `results/<dataset_name>` directory.

## Datasets
Put graph files into `datasets` directory and `querysets` directory.
- [Download](https://drive.google.com/drive/folders/1JCfeZFjAQT8_2jwrFexCiOSJarVSI0GQ?usp=sharing)

## Input File Format
Data graph file format is a text format to store a directed temporal data graph.
- The first line of the file should be "t # ID".
- Following lines of "v vertex-ID vertex-label" indicate the vertices in the graph.
  - The vertices should be written in the file in ascending order of their IDs, and a vertex ID should be in [0, #vertices - 1].
- Edges are represented as follows:
  - Following lines of "e source-ID destination-ID source-port destination-port edge-label sec msec" indicate the edges in the graph for Netflow.
  - Following lines of "e source-ID destination-ID edge-label sec" indicate the edges in the graph for Wikitalk and Superuser.

Query graph file format is a text format to store a directed temporal query graph and constraints on timing order.
- The first line of a graph should be "t # s ID".
- Following lines of "v vertex-ID vertex-label" indicate the vertices in the graph.
- Edges are represented as follows:
  - Following lines of "e source-ID destination-ID edge-label source-port destination-port" indicate the edges in the graph for Netflow.
  - Following lines of "e source-ID destination-ID edge-label" indicate the edges in the graph for Wikitalk and Superuser.
  - The ID of the $i$-th edge in the file is $i-1$.
- Following lines of "b edge-ID1 edge-ID2" indicate that the strict partial order edge-ID1 $\prec$ edge-ID2 holds.
