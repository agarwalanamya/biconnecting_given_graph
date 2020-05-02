#!/usr/bin/env python
# coding: utf-8

# In[15]:


import networkx as nx # Networkx Library
import matplotlib.pyplot as plt # Matplot lib


# In[16]:


# Takes vertices of the graph as input
vertices =[];
print("Enter the vertex")
while True:
    vertices.append(input(" >>>"))
    print("Want to enter another vertex? Y/N")
    counter = input(" >>>")
    if(counter != 'Y'):
        break;
vertices; # Gives a list of vertices in graph


# In[13]:


edges = [];
print("Enter the edge")
while True:
    vertex1 = input(" >>>")
    vertex2 = input(" >>>")
    edges.append((vertex1,vertex2))
    print("Want to enter another vertex? Y/N")
    counter = input(" >>>")
    if(counter != 'Y'):
        break;
edges
    


# In[17]:


#Graph Initialisation(Will be converted into an interface later)
G = nx.Graph() 
G.add_nodes_from(vertices)
G.add_edges_from(edges)
#Draws our PTP Graph
plt.subplot(111)
nx.draw_planar(G, with_labels = True)
plt.show()


# In[18]:


# Conversion to directed
H = G.to_directed()


# In[14]:


# Get all triangles
all_cycles = list(nx.simple_cycles(H))
all_triangles = []
for cycle in all_cycles:
     if(len(cycle) == 3):
             all_triangles.append(cycle)
all_triangles #Contains all triangles in PTP graph


# In[20]:


# Get edges on outer boundary
outer_boundary = []
for edge in H.edges:
	count = 0
	for triangle in all_triangles:
		if(edge[0] in triangle and edge[1] in triangle):
			count += 1
	if(count == 2):
		outer_boundary.append(edge)
outer_boundary #Contains the outer boundary of PTP Graph


# In[21]:


# Get Vertex-Set of outerboundary
outer_vertices = []
for edge in outer_boundary:
	if(edge[0] not in outer_vertices):
		outer_vertices.append(edge[0])
	if(edge[1] not in outer_vertices):
		outer_vertices.append(edge[1])
outer_vertices #Contains vertices containing outer vertices of PTP Graph


# In[22]:


# Get top,left,right and bottom boundaries of graph
corner_implying_paths = []
loop_count = 0
while len(outer_vertices)>1:
	temp=[]
	temp.append(outer_vertices[0])
	outer_vertices.pop(0)
	for vertices in temp:
		for vertex in outer_vertices:
			temp1 = temp.copy()
			temp1.pop(len(temp)-1)
			if((temp[len(temp)-1],vertex) in outer_boundary):
				temp.append(vertex)
				outer_vertices.remove(vertex)		
				if(temp1 is not None):
					for vertex1 in temp1:
						if((vertex1,vertex) in H.edges):
							temp.remove(vertex)
							outer_vertices.append(vertex)
	corner_implying_paths.append(temp)
	outer_vertices.insert(0,temp[len(temp)-1])
	if(len(outer_vertices) == 1 and loop_count == 0):
		outer_vertices.append(corner_implying_paths[0][0])
		loop_count += 1

def create_cip ( index ):
	corner_implying_paths.insert(index+1,corner_implying_paths[index])
	corner_implying_paths[index] = corner_implying_paths[index][0:2]
	del corner_implying_paths[index+1][0:1]


if(len(corner_implying_paths) == 5):
	corner_implying_paths[0] = corner_implying_paths[4] + list(set(corner_implying_paths[0]) - set(corner_implying_paths[4]))
	corner_implying_paths.pop()

if(len(corner_implying_paths) == 3):
	index = corner_implying_paths.index(max(corner_implying_paths,key =len))
	create_cip(index)

if(len(corner_implying_paths) == 2):
	index = corner_implying_paths.index(max(corner_implying_paths,key =len))
	create_cip(index)
	create_cip(index+1)
corner_implying_paths    #Gives boundaries of PTP Graph


# In[23]:


#Adding north, south, east and west vertices and connects them to boundary vertices 
	
G.add_nodes_from(["North","South","East","West"])


def news_edges( cip , vertex ):
	"This adds north,east,south,west"
	for vertices in cip:
		G.add_edge(vertex,vertices)

news_edges(corner_implying_paths[0],'North')
news_edges(corner_implying_paths[1],'East')
news_edges(corner_implying_paths[2],'South')
news_edges(corner_implying_paths[3],'West')


G.add_edges_from([('North','West'),('West','South'),('South','East'),('East','North')])


# In[24]:


#Creates new graph
plt.subplot(111)
nx.draw_planar(G, with_labels=True)
plt.show()


# In[ ]:


#For example 1
# G.add_nodes_from(["a","b","c","d","e","f","g","h","i","j"])
# G.add_edges_from([("a","b"),("a","d"),("a","f"),("a","g"),("a","j"),("b","c"),("b","d"),("b","g"),("b","h"),("c","d"),("c","h"),("c","e"),("d","e"),("d","f"),("e","f"),("e","i"),("e","j"),("f","j"),("i","j")])

#For example 2
# G.add_nodes_from(["a","b","c","d","e"])
# G.add_edges_from([("a","d"),("a","b"),("a","e"),("b","c"),("b","e"),("c","d"),("c","e"),("d","e")])

