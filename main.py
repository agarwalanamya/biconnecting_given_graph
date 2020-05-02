import tkinter as tk
import turtle
import ptpg

import warnings
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

pen = None
scale = 300
origin = {'x': 300, 'y': -150}


def draw_undirected_graph(graph):
    pen.clear()
    pen.pencolor('black')
    pen.penup()
    for from_node in range(graph.matrix.shape[0]):
        pen.setposition(graph.node_position[from_node][0] * scale + origin['x'],
                        graph.node_position[from_node][1] * scale + origin['y'])
        if from_node == graph.north:
            pen.write("N")
        elif from_node == graph.south:
            pen.write("S")
        elif from_node == graph.east:
            pen.write("E")
        elif from_node == graph.west:
            pen.write("W")
        else:
            pen.write(from_node)
        for to_node in range(from_node):
            if graph.matrix[from_node][to_node] == 1:
                pen.setposition(graph.node_position[from_node][0] * scale + origin['x'],
                                graph.node_position[from_node][1] * scale + origin['y'])
                pen.pendown()
                pen.setposition(graph.node_position[to_node][0] * scale + origin['x'],
                                graph.node_position[to_node][1] * scale + origin['y'])
                pen.penup()


def draw_directed_graph(graph):
    pen.clear()
    pen.width(1)
    pen.penup()
    for from_node in range(graph.matrix.shape[0]):
        pen.setposition(graph.node_position[from_node][0] * scale + origin['x'],
                        graph.node_position[from_node][1] * scale + origin['y'])
        if from_node == graph.north:
            pen.write("N")
        elif from_node == graph.south:
            pen.write("S")
        elif from_node == graph.east:
            pen.write("E")
        elif from_node == graph.west:
            pen.write("W")
        else:
            pen.write(from_node)
        for to_node in range(graph.matrix.shape[0]):
            if graph.matrix[from_node][to_node] == 0:
                continue
            else:
                if graph.matrix[from_node][to_node] == 2:
                    pen.color('blue')
                elif graph.matrix[from_node][to_node] == 3:
                    pen.color('red')
                pen.setposition(graph.node_position[from_node][0] * scale + origin['x'],
                                graph.node_position[from_node][1] * scale + origin['y'])
                pen.pendown()
                pen.setposition(
                    ((graph.node_position[from_node][0] + graph.node_position[to_node][0]) * scale / 2) + origin['x'],
                    ((graph.node_position[from_node][1] + graph.node_position[to_node][1]) * scale / 2) + origin['y'])
                if graph.matrix[from_node][to_node] != 1:
                    pen.width(2)
                pen.setposition(graph.node_position[to_node][0] * scale + origin['x'],
                                graph.node_position[to_node][1] * scale + origin['y'])
                pen.penup()
                pen.color('black')
                pen.width(1)


def draw_rdg(graph):
    pen.clear()
    pen.width(1.5)
    pen.color('black')
    pen.hideturtle()
    pen.penup()
    scale = 75
    origin = {'x': 50, 'y': -300}
    for i in range(graph.room_x.shape[0]):
        if graph.room_width[i] == 0:
            continue
        pen.setposition(graph.room_x[i] * scale + origin['x'], graph.room_y[i] * scale + origin['y'])
        pen.pendown()
        pen.setposition((graph.room_x[i] + graph.room_width[i]) * scale + origin['x'],
                        graph.room_y[i] * scale + origin['y'])
        pen.setposition((graph.room_x[i] + graph.room_width[i]) * scale + origin['x'],
                        (graph.room_y[i] + graph.room_height[i]) * scale + origin['y'])
        pen.setposition(graph.room_x[i] * scale + origin['x'],
                        (graph.room_y[i] + graph.room_height[i]) * scale + origin['y'])
        pen.setposition(graph.room_x[i] * scale + origin['x'], graph.room_y[i] * scale + origin['y'])
        pen.penup()
        pen.setposition(((2 * graph.room_x[i] + graph.room_width[i]) * scale / 2) + origin['x'],
                        ((2 * graph.room_y[i] + graph.room_height[i]) * scale / 2) + origin['y'])
        pen.write(i)
        pen.penup()


def plot_graph(G):  # draw the graph with just G as input
    H = nx.from_numpy_matrix(G.matrix)
    nx.draw_planar(H, with_labels=True)


def main():
    warnings.filterwarnings("ignore")  # To ignore warning given by matplotlib
    G = ptpg.PTPG()
    plt.figure(1)
    plot_graph(G)
    if (G.isBiconnected()):
        print("The given graph is biconnected")
        #print("Below are the biconnected components")
        #G.print_biconnected_components()
    else:
        print("The given graph is not biconnected")
        print("Making the graph biconnected")
        #G.print_biconnected_components()
        #print(G.no_of_articulation_points)
        #print(G.articulation_points_value)
        #print(G.bcc_sets)
        #print(G.articulation_point_sets)
        #print(G.articulation_points_value)
        #print(G.belong_in_same_block(3, 1))
        #print("Below are the biconnected components")
        #G.print_biconnected_components()
        G.initialize_bcc_sets()
        G.find_articulation_points()
        G.make_biconnected()
        print("Is the graph now Biconnected : ", G.isBiconnected(), "; As shown in figure 2")
        plt.figure(2)
        plot_graph(G)



"""
    G.add_nesw_vertices()
    G.node_position = nx.planar_layout(nx.from_numpy_matrix(G.matrix))
    draw_undirected_graph(G)
    x = input()
    G.initialize_degrees()
    G.initialize_good_vertices()
    v, u = G.contract()
    while v != -1:
        print("Contracted the edge between " + str(v) + " and " + str(u))
        node_position = nx.planar_layout(nx.from_numpy_matrix(G.matrix))
        draw_undirected_graph(G)
        x = input()
        v, u = G.contract()
        
    G.get_trivial_rel()
    draw_directed_graph(G)

    # REMOVE following 2 lines later
    construct_rdg(G)
    draw_rdg(G)

    while len(G.contractions) != 0:
        x = input()
        G.expand()
        draw_directed_graph(G)

        #REMOVE following 2 lines later
        construct_rdg(G)
        draw_rdg(G)

    # G.populate_t1_matrix()
    # print(G.t1_matrix)
    # G.populate_t2_matrix()
    # print(G.t2_matrix)
    # G.get_dimensions()
    # print(G.room_x)
    # print(G.room_y)
    # print(G.room_width)
    # print(G.room_height)
    construct_rdg(G)
    draw_rdg(G)
"""


# remove this function later... adding this so that we can take screenshots of intermediate floorplans coresponding
# to the intermediate RELs in exoansion proces
def construct_rdg(G):
    G.t1_matrix = None
    G.t2_matrix = None
    G.t1_longest_distance = [-1] * (G.west + 1)
    G.t2_longest_distance = [-1] * (G.west + 1)
    G.t1_longest_distance_value = -1
    G.t2_longest_distance_value = -1
    G.n_s_paths = []
    G.w_e_paths = []

    G.room_x = np.zeros(G.west - 3)
    G.room_y = np.zeros(G.west - 3)
    G.room_height = np.zeros(G.west - 3)
    G.room_width = np.zeros(G.west - 3)
    G.populate_t1_matrix()
    G.populate_t2_matrix()
    G.get_dimensions()


if __name__ == '__main__':
    root_window = tk.Tk()
    root_window.title('Rectangular Dual')
    root_window.geometry(str(1000) + 'x' + str(600))
    root_window.resizable(0, 0)
    root_window.grid_columnconfigure(0, weight=1, uniform=1)
    root_window.grid_rowconfigure(0, weight=1)

    border_details = {'highlightbackground': 'black', 'highlightcolor': 'black', 'highlightthickness': 1}

    canvas = tk.Canvas(root_window, **border_details)
    canvas.pack_propagate(0)
    canvas.grid(row=0, column=0, sticky='nsew', padx=4, pady=4)
    pen = turtle.RawTurtle(canvas)
    pen.speed(30)
    main()
    root_window.mainloop()
