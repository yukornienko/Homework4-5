from Bio import SeqIO
import Bio
from collections import defaultdict
import pygraphviz as pgv

class Vertex:
    
    def __init__(self, seq):
        self.seq = seq
        self.coverage = 1
        self.in_edges = {}
        self.out_edges = {}
        
    def increase_coverage(self):
        self.coverage += 1

class Edge:
    
    def __init__(self,k1,k2):
        self.seq = k1 + k2[-1]
        self.n = 2
        self.coverage = 0
    
    def calc_coverage(self,c1,c2):
        self.coverage = (c1+c2)/2


class Graph:

    def __init__(self,k):
        self.vertices = {}
        self.k = k
        
    def add_read(self,read):
        read_lng = len(read)
        if read_lng < self.k:
            return
            
        kmer = read[:k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)
        
        for next_kmer_indx in range(1,read_lng-k+1,1):
            next_kmer = read[next_kmer_indx:(next_kmer_indx+k)]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)
            
            new_edge = Edge(kmer,next_kmer)
            
            self.vertices[next_kmer].in_edges[kmer]  = [new_edge]
            
            self.vertices[kmer].out_edges[next_kmer] = [new_edge]

            kmer = next_kmer
    
    def calc_init_edge_coverage(self):
        
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0].calc_coverage(self.vertices[current_vertex].coverage,self.vertices[next_vertex].coverage)
    

    def visualize(self, fullness, output_name):
        self.G = pgv.AGraph(directed=True)
        if fullness == "short":
            for node in self.vertices.keys():
                for in_edge in self.vertices[node].in_edges:
                    self.G.add_node(in_edge, label = str(self.vertices[in_edge].coverage))
                    self.G.add_edge(in_edge, node, label = "Len=" + str(len(self.vertices[in_edge].seq[-1] + self.vertices[node].seq) - k + 1) + ",Cov= " + str((self.vertices[in_edge].coverage + self.vertices[node].coverage)/2))
                for out_edge in self.vertices[node].out_edges:
                    self.G.add_node(out_edge, label = str(self.vertices[out_edge].coverage))
                    self.G.add_edge(node, out_edge, label = "Len=" + str(len(self.vertices[node].seq[-1] + self.vertices[out_edge].seq) - k + 1) + ",Cov= " + str((self.vertices[node].coverage + self.vertices[out_edge].coverage)/2))
        else:
            for node in self.vertices.keys():
                for in_edge in self.vertices[node].in_edges:
                    self.G.add_node(in_edge, label = str(self.vertices[in_edge].seq))
                    self.G.add_edge(in_edge, node, label = (str(self.vertices[in_edge].seq) + str(self.vertices[node].seq[-1])))
                for out_edge in self.vertices[node].out_edges:
                    self.G.add_node(out_edge, label = str(self.vertices[out_edge].seq))  
                    self.G.add_edge(node, out_edge, label = (str(self.vertices[node].seq) + str(self.vertices[out_edge].seq[-1])))
   
        #self.G.layout()
        #self.G.draw(output_name, format='png', prog='dot')
        self.G.write(output_name)
        
k = 15
    
my_graph = Graph(k)
    
with open("/Users/yukornienko/Downloads/BI_spring_2018/Python/hw_4_5_dataset.fasta", "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        read = str(record.seq)
        my_graph.add_read(read)
        #my_graph.add_read( str(record.reverse_complement().seq) )

    my_graph.calc_init_edge_coverage()
    
    for v in my_graph.vertices:
        print('Vertex: {}, coverage: {}'.format(v,my_graph.vertices[v].coverage))
        for e in my_graph.vertices[v].out_edges:
            print('-> Out edge: {}'.format(e))
        for e in my_graph.vertices[v].in_edges:
            print('-> In edge: {}'.format(e))    
  
my_graph.visualize("short", "/Users/yukornienko/Downloads/BI_spring_2018/Python/HW4_short_15.dot")
my_graph.visualize("full", "/Users/yukornienko/Downloads/BI_spring_2018/Python/HW4_full_15.dot")

