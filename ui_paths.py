import streamlit as st
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict

st.set_page_config(layout="wide")

st.title("Genome Graph Visualizer")

sequence = st.text_input("Enter Genome Sequence", "ATGCACTG").strip().upper()
k = st.slider("Select k-mer size", 2, 6, 3)

MAX_LEN = 300

if len(sequence) > MAX_LEN:
    sequence = sequence[:MAX_LEN]
    st.warning("Sequence truncated for performance")

def construct_kmers(seq, k):
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

def hamilton_edges(kmers):
    edges = []
    n = len(kmers)
    if n > 200:
        return []
    for i in kmers:
        for j in kmers:
            if i != j and i[1:] == j[:-1]:
                edges.append((i, j))
    return edges

def build_debruijn(kmers):
    graph = defaultdict(list)
    for kmer in kmers:
        graph[kmer[:-1]].append(kmer[1:])
    return graph

def find_start(graph):
    indeg = defaultdict(int)
    outdeg = defaultdict(int)
    for u in graph:
        outdeg[u] += len(graph[u])
        for v in graph[u]:
            indeg[v] += 1
    nodes = set(indeg) | set(outdeg)
    for n in nodes:
        if outdeg[n] - indeg[n] == 1:
            return n
    return list(graph.keys())[0]

def eulerian_path(graph):
    g = {u: graph[u][:] for u in graph}
    start = find_start(g)
    stack = [start]
    path = []
    while stack:
        v = stack[-1]
        if v in g and g[v]:
            stack.append(g[v].pop())
        else:
            path.append(stack.pop())
    return path[::-1]

def reconstruct(path):
    seq = path[0]
    for node in path[1:]:
        seq += node[-1]
    return seq

kmers = construct_kmers(sequence, k)

tab1, tab2, tab3 = st.tabs(["Hamiltonian", "Eulerian", "DeBruijn"])

with tab1:
    st.subheader("Hamiltonian Graph")
    edges = hamilton_edges(kmers)
    if not edges:
        st.warning("Too large for Hamiltonian graph")
    else:
        G = nx.DiGraph()
        for u, v in edges:
            G.add_edge(u, v)
        pos = nx.circular_layout(G)
        plt.figure()
        nx.draw(G, pos, with_labels=True)
        st.pyplot(plt)

with tab2:
    st.subheader("Eulerian Path")
    graph = build_debruijn(kmers)
    path = eulerian_path(graph)
    reconstructed_seq = reconstruct(path)
    st.write("Path:", path[:20])
    st.write("Reconstructed:", reconstructed_seq[:100])
    if sequence == reconstructed_seq:
        st.success("MATCH")
    else:
        if reconstructed_seq in sequence:
            st.success("MATCH")
        else:
            st.error("NOT MATCH")

with tab3:
    st.subheader("DeBruijn Graph")
    graph = build_debruijn(kmers)
    G = nx.DiGraph()
    count = 0
    for u in graph:
        for v in graph[u]:
            if count > 300:
                break
            G.add_edge(u, v)
            count += 1
    pos = nx.circular_layout(G)
    plt.figure()
    nx.draw(G, pos, with_labels=True)
    st.pyplot(plt)
