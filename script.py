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
        
        def compress (self):
       
        self.to_be_comressed = [] #создаем список вершин, которые будем сжимать (первые из них)
        
        for node in self.vertices.keys():
            if len(list(self.vertices[node].out_edges.keys())) == 1:
                next_node = list(self.vertices[node].out_edges)[0]
                if len(list(self.vertices[node].in_edges.keys())) == 0 and len(list(self.vertices[next_node].in_edges.keys())) == 1:
                    self.to_be_comressed += [node]
                   
                elif len(list(self.vertices[node].in_edges.keys())) == 1 and len(list(self.vertices[next_node].in_edges.keys())) == 1:
                    if list(self.vertices[node].in_edges.keys())[0] not in self.to_be_comressed: 
                        self.to_be_comressed += [node]
                        
        for node in self.to_be_comressed:  #иттерируемся по созданному списку
            self.new_edges = []
            self.new_edges_in = []
    
            next_node = list(self.vertices[node].out_edges)[0]
            new_node = str(node) + str(next_node)[-1]

            if new_node in self.vertices:
                self.vertices[new_node].increase_coverage()                
            else:
                self.vertices[new_node] = Vertex(new_node)
                self.vertices[new_node].coverage = self.vertices[node].coverage

                if len(list(self.vertices[node].in_edges)) != 0:
                    
                    previous_node = list(self.vertices[node].in_edges)[0]
                    new_in_edge = Edge(previous_node, new_node)
                    del self.vertices[previous_node].out_edges[node]
                    self.vertices[new_node].in_edges[previous_node] = [new_in_edge]
                    self.vertices[previous_node].out_edges[new_node] = [new_in_edge]
                    del self.vertices[node].in_edges[previous_node]

                if len(list(self.vertices[next_node].out_edges)) != 0:
                    for out_edge_ in self.vertices[next_node].out_edges:
                        self.new_edges += [Edge(new_node, out_edge_)]
                        del self.vertices[out_edge_].in_edges[next_node]
                        self.vertices[out_edge_].in_edges[new_node] = self.new_edges 
                        self.vertices[new_node].out_edges[out_edge_] = self.new_edges
                    
                    del self.vertices[next_node].out_edges

            del self.vertices[next_node].in_edges                    
            del self.vertices[node].out_edges
            del self.vertices[node]
            del self.vertices[next_node]


        
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
  
my_graph.visualize("short", "/HW4_short_15.dot")
my_graph.visualize("full", "/HW4_full_15.dot")
my_graph.compress()
my_graph.visualize("full", "/hw_5_compressed_15full.dot")
my_graph.visualize("short", "/hw_5_compressed_15short.dot")

