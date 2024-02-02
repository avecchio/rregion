class DeBruijnGraph:
    def __init__(self, sequences, k):
        self.graph = {}
        self.k = k

        # Construct the De Bruijn graph
        for sequence in sequences:
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i + k]
                prefix = kmer[:-1]
                suffix = kmer[1:]

                if prefix not in self.graph:
                    self.graph[prefix] = []

                self.graph[prefix].append(suffix)

    def traverse(self, start_node):
        # Traverse the De Bruijn graph starting from the given node
        path = [start_node]
        current_node = start_node

        while current_node in self.graph:
            next_node = self.graph[current_node].pop()
            path.append(next_node)
            current_node = next_node

        return ''.join(path)

# Example usage
sequences = ["ATCGATCG", "CGATCGAT", "GATCGATC"]
k = 5

# Construct De Bruijn graph
de_bruijn_graph = DeBruijnGraph(sequences, k)

# Traverse the graph from a specific node
start_node = 'ATC'
resulting_sequence = de_bruijn_graph.traverse(start_node)

# Print the De Bruijn graph
print("De Bruijn Graph:")
for node, neighbors in de_bruijn_graph.graph.items():
    print(f"{node} -> {', '.join(neighbors)}")

# Print the resulting sequence after traversal
print("\nResulting Sequence after Traversal:")
print(resulting_sequence)