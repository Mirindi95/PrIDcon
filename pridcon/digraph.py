class WDiGraph:
    """Class implementing a directed multi graph with a weighted digraph,
    weights are the number of that specific edge. """

    class Node:
        """Node class to store different attributes required to perform DFS and eulerian paths."""

        def __init__(self):
            self.in_degree = 0
            self.out_degree = 0
            self.balance = 0  # balanced = 0, semi-balanced = 1 or -1
            self.color = 'white'
            self.predecessor = None
            self.dist = -1  # Book pseudocode has infinitive, but since distance/time cannot be negative, initialization is -1

        def update_balance(self):
            self.balance = self.out_degree - self.in_degree

    def __init__(self):
        """Empty digraph initialization. The data structure is a adjacency list implemented with
        nested dictionaries."""
        self.graph = dict()  # adjacency list
        self.nodes = dict()  # Dict of Node objects (label is the dict key)
        self.path = [] # list of nodes that form the eulerian path 

    def add_node(self, node_label) -> None:
        """Adds a node to the graph in case it does not exist."""
        if node_label not in self.nodes.keys():
            self.nodes[node_label] = self.Node()
            self.graph[node_label] = dict()

    def add_edge(self, source_node, target_node) -> None:
        """Add edge between source_node and target_node, if node exists, weight is increased by 1"""
        try:
            assert source_node in self.graph.keys() and target_node in self.graph.keys()
        except AssertionError:
            print('Nodes are not in the graph.')
        else:
            if target_node in self.graph[source_node].keys():  # if edge already exists
                self.graph[source_node][target_node] += 1
            else:
                self.graph[source_node][target_node] = 1

            self.nodes[source_node].out_degree += 1
            self.nodes[target_node].in_degree += 1
            self.nodes[source_node].update_balance()
            self.nodes[target_node].update_balance()

    def number_of_nodes(self) -> int:
        """Returns the number of nodes."""
        return len(self.graph.keys())

    def number_of_edges(self) -> int:
        """Returns the number of edges."""
        total_edges = 0
        for edges in self.graph.values():
            total_edges += sum(edges.values())
        return total_edges

    def is_eulerian(self) -> int:
        """Returns a Tuple. True if the graph is eulerian, False otherwise."""  # slide 10 definition
        total_nodes = self.number_of_nodes()
        pos, neg, balanced, _, _ = self.balance_stats()
        total_semibalanced = pos + neg
        return total_semibalanced in {0, 2} and pos == neg and total_semibalanced + balanced == total_nodes

    def euler_type(self)-> str or None:
        """Returns a string. Path if it has eulerian path, circuit if it has eulerian circuit and None otherwise."""
        if self.is_eulerian():
            pos, neg, _, _, _ = self.balance_stats()
            if pos == 1 and neg == 1:
                return 'path'
            else:
                return 'circuit'
        else:
            return None

    def balance_stats(self) -> list:
        """Returns the count of positive semibalanced, negative semibalanced, balanced nodes and list of heads and tails present in the graph."""
        balanced_nodes = 0
        pos_semibalanced = 0
        neg_semibalanced = 0
        head = []  # list of pos balanced nodes (heads)
        tail = []  # list of neg balanced nodes (tails)
        for label, stats in self.nodes.items():
            if stats.balance == 0:
                balanced_nodes += 1
            elif stats.balance == 1:
                pos_semibalanced += 1
                head.append(label)
            elif stats.balance == -1:
                neg_semibalanced += 1
                tail.append(label)
        return pos_semibalanced, neg_semibalanced, balanced_nodes, head, tail

    def next_node(self, node):
        """
        Returns the next reachable node from stack.

        Arguments
        ---------
            node: current node
        Returns
        -------
            next node
        """

        for next_node in self.graph[node].keys():
            if self.nodes[next_node].in_degree != 0 and self.graph[node][next_node] > 0:  # edge exist between nodes
                self.graph[node][next_node] -= 1  # remove edge
                self.nodes[next_node].in_degree -= 1  # decrease in_degree
                yield next_node

    def dfs_visit(self, node):
        """
        Depth First Search traversing eulerian graph.
        Arguments
        ---------
            node: node to explore
        """
        if self.nodes[node].out_degree != 0:
            self.nodes[node].out_degree -= 1  # remove edge
            next_node = next(self.next_node(node))
            self.dfs_visit(next_node)

        self.path.append(node)

    def eulerian_path(self, start_node, shift=1) -> str:
        """
        Builds the eulerian path from the graph.
        """
        self.dfs_visit(start_node)
        path = self.path[::-1]
        return path[0] + ''.join(map(lambda x: x[-shift:], path[1:]))


class deBruijnGraph(WDiGraph):
    """de Bruijn graph implementation. Subclass of WDiGraph."""

    def __init__(self, reads_dict, k=30, assembly_name='genome'):
        super().__init__()
        for _, read_seq in reads_dict.items():
            self.update_from_read(read_seq, k)
        self._contigs = self.build_assembly_from_sequence(assembly_name)
        
    def get_contigs(self):
        """Returns the found contigs as a dictionary."""
        return self._contigs
        
    def update_from_read(self, read_seq, k, shift=1) -> None:
        """
        updates the graph from a sequence.

        Arguments
        ---------
            read_seq: sequence
            k: size of K-mer
            shift: shift of left and right K-mer

        Returns
        -------
            None
        """
        for i in range(0, len(read_seq) - (k - 1), shift):  # k-1 to take the last k-mer
            k_mer = read_seq[i: i + k]
            L_k1_mer = k_mer[: -shift]
            R_k1_mer = k_mer[shift:]
            self.add_node(L_k1_mer)
            self.add_node(R_k1_mer)
            self.add_edge(L_k1_mer, R_k1_mer)
            
    def build_assembly_from_sequence(self, identifier: str) -> dict:
        """
        Creates a dictionary with sequence identifier as key and assembled string as value.
        
        Arguments
        ---------
            identifier: sequence identifier
        Returns
        -------
            dict
        """
              
        assembly = {}
        _, _, _, start_list, _ = self.balance_stats()
        for index, start_node in enumerate(start_list):
            contig = self.eulerian_path(start_node)
            assembly[identifier+'_contig_'+str(index)] = contig
        return assembly






