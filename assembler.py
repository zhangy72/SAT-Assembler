#!/usr/bin/env python
# this program is used to classify a certain number of reads that are
# in top K best paths.
import os
import sys
import operator
import subprocess
import networkx as nx
from pprint import pprint
import heapq
import cPickle
import random
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein
import copy
import Queue            

# define the read class.
class Read:
    def __init__(self, name):
        self.name = name
        self.seq = ''
        self.begin_state = 0
        self.end_state = 0
        self.score = 0.0
        self.length = 0
        self.size = 0
        self.members = set()
        self.descs = set() # desc should be unique
        self.cds_set = set()
  
    def __repr__(self):
        return repr((self.name, self.length, self.begin_state, self.end_state,
                     self.members, self.descs))
   
    def set_name(self, name):
        self.name = name

    def set_seq(self, seq):
        self.seq = seq
        self.length = len(seq)

    def set_size(self, size):
        self.size = size
   
    def set_begin_state(self, begin_state):
        self.begin_state = begin_state

    def set_end_state(self, end_state):
        self.end_state = end_state
 
    def set_score(self, score):
        self.score = score

    # add names of member reads belonging to the tag.
    def add_member(self, member):
        self.members.add(member)

    def add_desc(self, desc):
        self.descs.add(desc)
 
    def set_cds_set(self, cds_set):
        self.cds_set = cds_set

class Contig:
    def __init__(self, name):
        self.name = name
        nodes = [] # consider combined tags.
        tags = [] # all are single tags.
        members = set()
        descs = []
        seq = ''
        length = 0
        coverage = 0.0
        cds_set = set()
        shared_cds_set = set()
        final_cds_set = set()
    
    def __repr__(self):
        return repr((self.name, self.nodes, self.tags, self.members, self.descs, 
                     self.seq, self.length, self.coverage, self.cds_set))

    def set_nodes(self, nodes):
        self.nodes = nodes
 
    def set_tags(self, tags):
        self.tags = tags  
 
    def set_members(self, members):
        self.members = members 

    def set_descs(self, descs):
        self.descs = descs

    def set_seq(self, seq):
        self.seq = seq
        self.length = len(seq)

    def set_coverage(self, coverage):
        self.coverage = coverage

    def set_id(self):
        self.id = '%s_%d_%d_%.2f' % (self.name, len(self.members),
                                     self.length, self.coverage)
    
    def set_cds_set(self, cds_set):
        self.cds_set = cds_set

    def set_final_cds_set(self, final_cds_set):
        self.final_cds_set = final_cds_set

    def set_shared_cds_set(self, shared_cds_set):
        self.shared_cds_set = shared_cds_set

    def set_path(self, path):
        self.path = path
 
    def set_avg_overlap(self):
        length = len(self.seq)
        tag_num = len(self.tags)
        if tag_num == 1:
            self.avg_overlap = 0
        else:
          self.avg_overlap = (self.coverage*length-length)/(tag_num-1) 
  
    def set_name(self, name):
        self.name = name


# get domains with reads aligned to them.
# use read name instead of tag names
def get_hmmer_aligned_read_dict(in_file_name, target_domain,
                                EVALUE_THRES):
    aligned_read_dict = {} # reads that belong to domains.
    with open(in_file_name, 'Ur') as f:
        for line in f:
            if not line.strip():
                continue
            items = line.strip().split()
            domain = items[1][:7]
            evalue = float(items[3])
            if domain == target_domain and evalue <= EVALUE_THRES:
                read_name = items[0]
                score = float(items[2])
                begin_state = int(items[4])
                end_state = int(items[5])
                strand = items[8]
                composite_read_name = read_name + '$' + strand
                # keep the better hit of the same read.
                if composite_read_name in aligned_read_dict and \
                    score <= aligned_read_dict[composite_read_name].score:
                    continue
                read = Read(composite_read_name)
                read.set_score(score)
                read.set_begin_state(begin_state)
                read.set_end_state(end_state)
                read.add_member(read_name)
                aligned_read_dict[composite_read_name] = read
    return aligned_read_dict

def get_compressed_read_dict(target_read_dict, prefix):
    seq_read_dict = {}
    for read_name in target_read_dict:
        seq = target_read_dict[read_name].seq
        seq_read_dict.setdefault(seq, [])
        seq_read_dict[seq].append(read_name)
    compressed_read_dict = {}
    index = 1
    for seq in seq_read_dict:
        tag_name = '%s%d' % (prefix, index)
        index += 1
        # use the first read that has the same sequence.
        first_read_name = seq_read_dict[seq][0]
        first_read = target_read_dict[first_read_name] 
        tag = Read(tag_name)
        tag.set_begin_state(first_read.begin_state)
        tag.set_end_state(first_read.end_state)
        tag.set_score(first_read.score)
        tag.set_seq(first_read.seq)
        # keep the information 
        tag_size = len(seq)
        for read_name in seq_read_dict[seq]:
            read_descs = target_read_dict[read_name].descs
            tag.add_member(read_name)
            for desc in read_descs:
                tag.add_desc(desc)
        tag.set_size(tag_size)
        compressed_read_dict[tag_name] = tag
    return compressed_read_dict

# the read names end with '+' or '-'
def get_pairs_from_two_sets(read_set1, read_set2):
    pairs = set()
    for read in read_set1:
        if read == 'root':
            continue
        name = read[:-2]
        assert name[-1] == '1' or name[-1] == '2'
        if name[-1] == '1':
            other1 = name[:-1] + '2$+'
            other2 = name[:-1] + '2$-'
            if other1 in read_set2 or other2 in read_set2:
                pairs.add(name[:-2])
        else:
            other1 = name[:-1] + '1$+'
            other2 = name[:-1] + '1$-'
            if other1 in read_set2 or other2 in read_set2:
                pairs.add(name[:-2])
    return pairs
 
# get the sequence of the reads.
def set_read_seq(fasta_file_name, read_dict):
    with open(fasta_file_name, 'Ur') as f:
        for record in SeqIO.parse(f, 'fasta'):
            read_name = record.id
            composite_read_name = read_name + '$+'
            if composite_read_name in read_dict:
                seq = str(record.seq)
                read_dict[composite_read_name].set_seq(seq)
                read_dict[composite_read_name].set_size(len(seq))
            composite_read_name = read_name + '$-'
            if composite_read_name in read_dict:
                seq = str(record.seq.reverse_complement())
                read_dict[composite_read_name].set_seq(seq)
                read_dict[composite_read_name].set_size(len(seq))

# get the overlap length of two reads. return 0 if no overlap.
def get_seq_overlap_length(seq1, seq2):
    # no string should contain the other.
    hamming_thres = 3 
    len1 = len(seq1)
    len2 = len(seq2)
    max_overlap = 0 # currently maximum overlap.
    hamming_distance = 0
    # loop over number of overlappd pos.
    for i in xrange(min(len1, len2), 0, -1):
        hamming_distance = get_hamming_distance(seq1[-i:], seq2[:i])
        if hamming_distance <= hamming_thres:
            max_overlap = i
            break
    return (max_overlap, hamming_distance)   

def get_hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2)
    hamming_distance = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            hamming_distance += 1
    return hamming_distance

# get the overlap of two positions.
# if there is no overlap, return negative value.
def get_pos_overlap_length(begin1, end1, begin2, end2):
    return min(end1, end2) - max(begin1, begin2) + 1

# add root to transfer node score to edge weight.
def add_root_to_subgraph(subgraph):
    assert 'root' not in subgraph.nodes()
    subgraph.add_node('root', seq='', length=0, coverage=0.0, 
                      members=set(('root',)), descs=['root'], tags=['root'])
    for node in subgraph.nodes():
        if node != 'root' and subgraph.in_degree(node) == 0:
            subgraph.add_edge('root', node, overlap=0, hamming=0)

# convert weights of reads into read number.
def get_read_num_from_set(read_set):
    read_num = 0
    for read in read_set:
        # root is not counted as a read.
        if read == 'root':
            continue
        read_weight = int(read.split('_')[-1])
        read_num += read_weight
    return read_num

def find_complete_edges(G, read_length):
    if G.number_of_nodes() == 0:
        print 'Empty graph!'
        sys.exit()
    complete_edges = []
    for (node1, node2) in G.edges():
        if G[node1][node2]['overlap'] == read_length:
            complete_edges.append((node1, node2))

# get the contig from a path.
def get_contig_from_path(path, id, G):
    assert path[0] == 'root'
    contig = Contig(id)
    # path[0] is always root.
    contig_seq = ''
    nodes = []
    tags = []
    members = set()
    coverage = 0.0
    cds_set = set()
    contig_path = path[:]
    for i in xrange(1, len(path)):
        begin_node = path[i-1]
        end_node = path[i]
        nodes.append(end_node)
        tags += G.node[end_node]['tags']
        members |= G.node[end_node]['members']
        # concatenate sequence.
        overlap = G[begin_node][end_node]['overlap']
        seq = G.node[end_node]['seq']
        contig_seq += seq[overlap:]
        coverage += G.node[end_node]['coverage'] * G.node[end_node]['length']
        #print overlap, len(seq), seq
    #print len(contig), contig
    contig.set_seq(contig_seq)
    contig.set_nodes(nodes)
    contig.set_tags(tags)
    contig.set_members(members)
    coverage /= contig.length 
    contig.set_coverage(coverage)
    contig.set_path(contig_path)
    contig.set_avg_overlap()
    return contig
            
def get_sink_nodes(G):
    sink_nodes = []
    for node in G.nodes_iter():
        if G.out_degree(node) == 0:
            sink_nodes.append(node)
    return sink_nodes

def get_all_source_nodes(G):
    all_source_nodes = []
    for node in G.nodes_iter():
        if G.in_degree(node) == 0:
            all_source_nodes.append(node)
    return all_source_nodes

def set_overlap_graph_coverage(G):
    for node in G.nodes_iter():
        overlap_sum = 0.0
        for pred in G.predecessors_iter(node):
            overlap_sum += G[pred][node]['overlap']
        for succ in G.successors_iter(node):
            overlap_sum += G[node][succ]['overlap']
        G.node[node]['coverage'] = overlap_sum / G.node[node]['length']   
    # set coverage difference for edges.
    for node1, node2 in G.edges_iter():
        coverage1 = G.node[node1]['coverage']
        coverage2 = G.node[node2]['coverage']
        G[node1][node2]['delta_coverage'] = abs(coverage1 - coverage2)     

def set_overlap_graph_parameters(target_domain, G):
    G.graph['domain'] = target_domain
    G.graph['combined_basename'] = 'nodes'
    G.graph['combined_index'] = 1
    avg_length = 0.0
    for node in G.nodes_iter():
        avg_length += G.node[node]['length']
    G.graph['avg_length'] = avg_length / G.number_of_nodes()

# remove edges that do not help introduce any reads.
# return a transitive edge dict.
def remove_redundant_edges(G):
    # only for single nodes.
    transitive_nodes = {}
    edges = G.edges()
    for begin_node, end_node in edges:
        current_edge_data = G.edge[begin_node][end_node]
        G.remove_edge(begin_node, end_node)
        if not nx.has_path(G, begin_node, end_node):
            G.add_edge(begin_node, end_node, current_edge_data) 
        else:
            transitive_nodes[(begin_node, end_node)] = True
    return transitive_nodes

# remove directed cycles in the graph. make it DAG.
def remove_cycles(G):
    while not nx.is_directed_acyclic_graph(G):
        subgraphs = nx.strongly_connected_component_subgraphs(G)
        for subgraph in subgraphs:
            if subgraph.number_of_nodes() > 1:
                min_edge = get_shortest_edge(subgraph)
                #print subgraph.edges()
                #print min_edge
                G.remove_edge(*min_edge)

def get_shortest_edge(G):
    min_overlap = 100000
    for node1, node2 in G.edges_iter():
        if G[node1][node2]['overlap'] < min_overlap:
            min_edge = (node1, node2)
            min_overlap = G[node1][node2]['overlap']
    return min_edge

def get_sink_num(G):
    sink_num = 0
    for node in G.nodes_iter():
        if G.out_degree(node) == 0:
            sink_num += 1
    return sink_num

# get the mapping between contigs, tags, and reads.
def get_read_mapping_dict(in_file_name):
    contig_read_dict = {}
    with open(in_file_name, 'Ur') as f:
        contig_name = ''
        for line in f:
            if line.strip():
                if line[0] == '>':
                    contig_name = line[1:].rstrip()
                else:
                    contig_read_dict.setdefault(contig_name, [])
                    contig_read_dict[contig_name].append(line.rstrip())
    return contig_read_dict


def output_tag_original_read_mapping(compressed_read_dict):
    for tag in compressed_read_dict:
        print '>' + tag
        for read in compressed_read_dict[tag].original_members:
            print read

# find the nodes that have two out-going edges.
def output_bi_nodes(subgraph):
    for node in subgraph.nodes_iter():
        if subgraph.out_degree(node) > 1:
            print 'out', node, subgraph[node]
        if subgraph.in_degree(node) > 1:
            for node1, node2 in subgraph.in_edges(node):
                print 'in', node1, node2, subgraph[node1][node2]

# count overlap between the same genes and different genes.
def count_overlaps_between_genes(subgraph):
    overlaps = [[], []] # inner gene and inter gene.
    for node in subgraph.nodes_iter():
        if subgraph.out_degree(node) > 1:
            print 'out', node, subgraph.out_edges(node, data=True) 
            for node1, node2 in subgraph.out_edges(node):
                overlap = subgraph[node1][node2]['overlap']
                gene1 = get_gene_from_node(node1)
                gene2 = get_gene_from_node(node2)
                if gene1 & gene2:
                    overlaps[0].append(overlap)
                else:
                    overlaps[1].append(overlap)
        if subgraph.in_degree(node) > 1:
            print 'in', node, subgraph.in_edges(node, data=True) 
            for node1, node2 in subgraph.in_edges(node):
                overlap = subgraph[node1][node2]['overlap']
                gene1 = get_gene_from_node(node1)
                gene2 = get_gene_from_node(node2)
                if gene1 & gene2:
                    overlaps[0].append(overlap)
                else:
                    overlaps[1].append(overlap)
    return overlaps

# count the overlaps between genes and different genes.
# all edges will be considered.
def count_overlaps(G):
    overlaps = [[], []]
    for node1, node2 in G.edges_iter():
        overlap = G[node1][node2]['overlap']
        gene1 = get_gene_from_node(node1)
        gene2 = get_gene_from_node(node2)
        if gene1 & gene2:
            #if len(gene1) > 1:
            #    print node1, gene1, node2, gene2
            overlaps[0].append(overlap)
        else:
            overlaps[1].append(overlap)
    return overlaps

def output_tags_fasta(compressed_read_dict):
    for tag_name in compressed_read_dict:
        print '>' + compressed_read_dict[tag_name].name
        print compressed_read_dict[tag_name].seq

def keep_single_out_edge(G, pair_num_dict):
    for node in G.nodes_iter():
        if len(G.out_edges(node)) > 1:
            kept_succs = set()
            all_succs = G.successors(node)
            for succ in all_succs:
                if (node, succ) in pair_num_dict or (succ, node) in pair_num_dict:
                    kept_succs.add(succ)
            if not kept_succs:
                OVERLAP_THRES = 0.5
                metrics = []
                for succ in all_succs:
                    overlap = G[node][succ]['overlap']
                    hamming = G[node][succ]['hamming']
                    metrics.append((overlap, -hamming, succ))
                opt_overlap, opt_hamming, opt_succ = max(metrics)
                for overlap, hamming, succ in metrics:
                    if succ != opt_succ:
                        G.remove_edge(node, succ)
            else:
                for succ in all_succs:
                    if succ not in kept_succs:
                        G.remove_edge(node, succ)

def keep_single_in_edge(G, pair_num_dict):
    for node in G.nodes_iter():
        if len(G.in_edges(node)) > 1:
            kept_preds = set()
            all_preds = G.predecessors(node)
            for pred in all_preds:
                if (pred, node) in pair_num_dict or (node, pred) in pair_num_dict:
                    kept_preds.add(pred)
            if not kept_preds:
                OVERLAP_THRES = 0.5
                metrics = []
                for pred in all_preds:
                    overlap = G[pred][node]['overlap']
                    hamming = G[pred][node]['hamming']
                    metrics.append((overlap, -hamming, pred))
                opt_overlap, opt_hamming, opt_pred = max(metrics)
                for overlap, hamming, pred in metrics:
                    if pred != opt_pred:
                        G.remove_edge(pred, node)
            else:
                for pred in all_preds:
                    if pred not in kept_preds:
                        G.remove_edge(pred, node)

def trim_multiple_out_edges(G, pair_num_dict):
    for node in G.nodes_iter():
        if len(G.out_edges(node)) > 1:
            #print node 
            cds_set = G.node[node]['cds_set']
            kept_succs = set()
            all_succs = G.successors(node)
            for succ in all_succs:
                if (node, succ) in pair_num_dict or (succ, node) in pair_num_dict:
                    kept_succs.add(succ)
            #print 'kept_succs', kept_succs
            if not kept_succs:
                OVERLAP_THRES = 0.5
                metrics = []
                for succ in all_succs:
                    overlap = G[node][succ]['overlap']
                    hamming = G[node][succ]['hamming']
                    metrics.append((overlap, -hamming, succ))
                opt_overlap, opt_hamming, opt_succ = max(metrics)
                for overlap, hamming, succ in metrics:
                    if opt_overlap - overlap >= OVERLAP_THRES * overlap and \
                       opt_hamming > hamming: 
                        succ_cds_set = G.node[succ]['cds_set']
                        G.remove_edge(node, succ)
                        #print 'Out', node, cds_set, overlap, opt_hamming, \
                        #      succ, overlap, hamming, succ_cds_set 
            else:
                for succ in all_succs:
                    if succ not in kept_succs:
                        succ_cds_set = G.node[succ]['cds_set']
                        #print 'Out', node, cds_set, succ, succ_cds_set
                        G.remove_edge(node, succ)

def trim_multiple_in_edges(G, pair_num_dict):
    for node in G.nodes_iter():
        if len(G.in_edges(node)) > 1:
            #print node
            cds_set = G.node[node]['cds_set']
            kept_preds = set()
            all_preds = G.predecessors(node)
            for pred in all_preds:
                if (pred, node) in pair_num_dict or (node, pred) in pair_num_dict:
                    #print pred, node
                    kept_preds.add(pred)
            if not kept_preds:
                OVERLAP_THRES = 0.5
                metrics = []
                for pred in all_preds:
                    overlap = G[pred][node]['overlap']
                    hamming = G[pred][node]['hamming']
                    metrics.append((overlap, -hamming, pred))
                opt_overlap, opt_hamming, opt_pred = max(metrics)
                for overlap, hamming, pred in metrics:
                    if opt_overlap - overlap >= OVERLAP_THRES * overlap and \
                       opt_hamming > hamming: 
                        pred_cds_set = G.node[pred]['cds_set']
                        G.remove_edge(pred, node)
                        #print 'In', node, cds_set, overlap, opt_hamming, \
                        #      pred, overlap, hamming, pred_cds_set 
            else:
                for pred in all_preds:
                    if pred not in kept_preds:
                        pred_cds_set = G.node[pred]['cds_set']
                        #print 'In', node, cds_set, pred, pred_cds_set
                        G.remove_edge(pred, node)

def show_best_out_edges(G):
    correct = 0
    incorrect = 0
    for node in G.nodes_iter():
        if len(G.out_edges(node)) >= 2:
            gene_set = get_gene_set(node)
            metrics = []
            succ_gene_set = set()
            for succ in G.successors_iter(node):
                succ_gene_set |= get_gene_set(succ)
                overlap = G[node][succ]['overlap']
                hamming = G[node][succ]['hamming']
                delta_coverage = G[node][succ]['delta_coverage']
                metrics.append((overlap, -delta_coverage, -hamming, succ))
            metrics.sort()
            metrics.reverse()
            best_out_node = metrics[0][-1]
            succ_gene_set = get_gene_set(best_out_node)
            if gene_set & succ_gene_set:
                correct += 1
            else:
                incorrect += 1
    print 'correct:%d' % (correct,), 'incorrect:%d' % (incorrect,) 

def show_best_in_edges(G):
    correct = 0
    incorrect = 0
    for node in G.nodes_iter():
        if len(G.in_edges(node)) >= 2:
            gene_set = get_gene_set(node)
            metrics = []
            pred_gene_set = set()
            for pred in G.predecessors_iter(node):
                pred_gene_set |= get_gene_set(pred)
                overlap = G[pred][node]['overlap']
                hamming = G[pred][node]['hamming']
                delta_coverage = G[pred][node]['delta_coverage']
                metrics.append((overlap, -delta_coverage, -hamming, pred))
            metrics.sort()
            metrics.reverse()
            best_in_node = metrics[0][-1]
            pred_gene_set = get_gene_set(best_in_node)
            if gene_set & pred_gene_set:
                correct += 1
            else:
                incorrect += 1
    print 'correct:%d' % (correct,), 'incorrect:%d' % (incorrect,) 

def show_multiple_in_edges(G):
    for node in G.nodes_iter():
        if len(G.in_edges(node)) >= 2:
            gene_set = get_gene_set(node)
            metrics = []
            pred_gene_set = set()
            for pred in G.predecessors_iter(node):
                pred_gene_set |= get_gene_set(pred)
                overlap = G[pred][node]['overlap']
                hamming = G[pred][node]['hamming']
                delta_coverage = G[pred][node]['delta_coverage']
                metrics.append((overlap, -delta_coverage, -hamming, pred))
            metrics.sort()
            metrics.reverse()
            print node, metrics

def trim_out_edges(G):
    DELTA_OVERLAP_THRES = 0.1
    DELTA_COVERAGE_THRES = 5.0
    for node in G.nodes_iter():
        if len(G.out_edges(node)) >= 2:
            metrics = []
            for succ in G.successors_iter(node):
                overlap = G[node][succ]['overlap']
                hamming = G[node][succ]['hamming']
                delta_coverage = G[node][succ]['delta_coverage']
                metrics.append((overlap, -delta_coverage, -hamming, succ))
            metrics.sort()
            metrics.reverse()
            kept_succs = []
            kept_succs.append(metrics[0][-1])
            out_node_num = len(metrics)
            for i in range(out_node_num-1):
                overlap1 = metrics[i][0]
                delta_coverage1 = metrics[i][1]
                overlap2 = metrics[i+1][0]
                delta_coverage2 = metrics[i+1][1]
                if float(overlap1-overlap2)/overlap2 <= DELTA_OVERLAP_THRES and \
                   (delta_coverage1-delta_coverage2)/abs(delta_coverage2) <= DELTA_COVERAGE_THRES:
                    kept_succs.append(succs[i+1][-1])
                else:
                    break  
            for succ in G.successors_iter(node):
                if succ not in kept_succs:
                    G.remove_edge(node, succ)

# keep the best edge for nodes that have multiple out-going edges.
def keep_best_out_edges(G):
    for node in G.nodes_iter():
        if len(G.out_edges(node)) >= 2:
            metrics = []
            for succ in G.successors_iter(node):
                overlap = G[node][succ]['overlap']
                hamming = G[node][succ]['hamming']
                delta_coverage = G[node][succ]['delta_coverage']
                metrics.append((overlap, -delta_coverage, -hamming, succ))
            metrics.sort()
            metrics.reverse()
            out_edge_num = len(metrics)
            for i in range(1, out_edge_num):
                G.remove_edge(node, metrics[i][-1])

def keep_best_in_edges(G):
    for node in G.nodes_iter():
        if len(G.in_edges(node)) >= 2:
            metrics = []
            for pred in G.predecessors_iter(node):
                overlap = G[pred][node]['overlap']
                hamming = G[pred][node]['hamming']
                delta_coverage = G[pred][node]['delta_coverage']
                metrics.append((overlap, -delta_coverage, -hamming, pred))
            metrics.sort()
            metrics.reverse()
            in_edge_num = len(metrics)
            for i in range(1, in_edge_num):
                G.remove_edge(metrics[i][-1], node)

def show_all_edges_states(G):
    for node1, node2 in G.edges_iter():
         print '>', node1, node2, G[node1][node2]['overlap'], G[node1][node2]['hamming']
         print G.node[node1]['descs'], G.node[node1]['begin_state']
         print G.node[node2]['descs'], G.node[node2]['begin_state']
         
def get_all_paths(source_node, G):
    all_paths = []
    sink_nodes = get_sink_nodes(G)
    for node in sink_nodes:
        for path in nx.all_simple_paths(G, source_node, node):
            all_paths.append(path)
    return all_paths

def get_contigs_from_paths(paths, G):
    contig_basename = G.graph['contig_basename']
    contigs = []
    domain = G.graph['domain']
    i = 1
    for path in paths:
        name = '%s_%s%d' % (domain, contig_basename, i)
        contig = get_contig_from_path(path, name, G)
        contig.set_id()
        contigs.append(contig)
        i += 1
    return contigs

def get_other_target_read_dict(target_read_dict):
    other_target_read_dict = {}
    # remove strand.
    read_name_list = [read[:-2] for read in target_read_dict]
    target_read_set = set(read_name_list)
    target_pair_set = get_pair_set(target_read_set)
    for read in target_read_dict:
        base_name = read[:-4]
        strand = read[-1]
        if strand == '+':
            other_strand = '-'
        else:
            other_strand = '+'
        if base_name not in target_pair_set:
            if read[-3] == '1':
                other_end = base_name + '.2' + '$' + other_strand
            else:
                other_end = base_name + '.1' + '$' + other_strand
            other_target_read_dict[other_end] = Read(other_end)
    return other_target_read_dict

def set_read_dict_desc(target_read_dict, 
                       read_desc_dict):
    for read in target_read_dict:
        desc = get_read_desc(read_desc_dict, read)
        assert not target_read_dict[read].descs
        for item in desc.split('&'):
            target_read_dict[read].add_desc(item)

def build_other_overlap_graph(other_compressed_read_dict, 
                              target_domain, 
                              overlap_thres):
    other_read_num = len(other_compressed_read_dict)
    G = nx.DiGraph() 
    other_read_list = other_compressed_read_dict.values()
    for read in other_read_list:
        G.add_node(read.name,
                   seq =read.seq,
                   length=read.length,
                   coverage=float(read.size)/read.length,
                   members=copy.deepcopy(read.members), #original read names.
                   descs=['&'.join(sorted(list(read.descs)))],
                   tags=[read.name]) # only consider single node or tags.
    for i in xrange(other_read_num-1):
        for j in xrange(i+1, other_read_num):
            read1 = other_read_list[i]
            read2 = other_read_list[j]
            read1_length = len(read1.seq)
            read2_length = len(read2.seq)
            seq_overlap_length1, h1 = get_seq_overlap_length(read1.seq, read2.seq)
            seq_overlap_length2, h2 = get_seq_overlap_length(read2.seq, read1.seq)
            if max(seq_overlap_length1, seq_overlap_length2) >= overlap_thres:
                if seq_overlap_length1 >= seq_overlap_length2:  
                    overlap = seq_overlap_length1
                    G.add_edge(read1.name, read2.name, 
                               overlap=overlap, hamming=h1)
                else:
                    overlap = seq_overlap_length2
                    G.add_edge(read2.name, read1.name, 
                               overlap=overlap, hamming=h2)
    return G
    
def get_pair_set(domain_read_set):
    pair_set = set()
    for read in domain_read_set:
        if read[-1] == '1':
            other_end = read[:-1] + '2'
            if other_end in domain_read_set:
                pair_set.add(read[:-2])
    return pair_set

def combine_single_paths(G):
    while True:
        nodes_to_combine = []
        # combine backward the nodes with in-degree = 1
        for node in G.nodes_iter():
            in_deg = G.in_degree(node)
            if in_deg == 1:
                pred = G.predecessors(node)[0]
                out_deg_pred = G.out_degree(pred)
                if out_deg_pred == 1:
                    nodes_to_combine.append(node)
        # continue until there is no qualified nodes
        if len(nodes_to_combine) == 0:
            break
        for curr in nodes_to_combine:
            # predecessor of current node
            pred = G.predecessors(curr)[0]
            # incoming edges to the combined node
            predspreds = G.predecessors(pred)
            # outgoing edges from the combined node
            succs = G.successors(curr)
            # create a new node
            new_node = combine_two_nodes(pred, curr, G)
            # make links between the new node and other nodes
            for p in predspreds:
                old_edge_data = G[p][pred]
                G.add_edge(p, new_node, old_edge_data)
            for s in succs:
                old_edge_data = G[curr][s]
                G.add_edge(new_node, s, old_edge_data)
            # remove used nodes
            G.remove_node(curr)
            G.remove_node(pred)
            if curr in nodes_to_combine:
                nodes_to_combine.remove(curr)
            if pred in nodes_to_combine:
                nodes_to_combine.remove(pred)

def combine_two_nodes(node1, node2, G):
    combined_basename = G.graph['combined_basename']
    combined_index = G.graph['combined_index']
    name = '%s%d' % (combined_basename, combined_index)
    G.graph['combined_index'] += 1
    overlap = G[node1][node2]['overlap']
    seq = G.node[node1]['seq'][:-overlap] + G.node[node2]['seq'] 
    length = len(seq)
    coverage = (G.node[node1]['coverage'] * G.node[node1]['length'] + 
                G.node[node2]['coverage'] * G.node[node2]['length']) / length 
    members = G.node[node1]['members'] | G.node[node2]['members']
    tags = G.node[node1]['tags'] + G.node[node2]['tags']
    new_node_data = dict(seq=seq, length=length, coverage=coverage, 
                         members=members, tags=tags)     
    assert not G.has_node(name)
    G.add_node(name, new_node_data)
    return name

def get_all_graph_tips(G, coverage_thres):
    tips = set()
    for node in G.nodes_iter():
        if is_out_tip(node, G, coverage_thres) or \
           is_in_tip(node, G, coverage_thres):
            tips.add(node)
    #for tip in tips:
    #    print tip, G.node[tip]['descs']
    #exit()
    return tips

# recursively remove tips.
def remove_graph_tips(G):
    COVERAGE_THRES = 2.0
    tips = get_all_graph_tips(G, COVERAGE_THRES)
    while tips:
        for node in tips:
            #print node, G.node[node]
            G.remove_node(node)
        tips = get_all_graph_tips(G, COVERAGE_THRES)

def is_out_tip(node, G, coverage_thres):
    if G.out_degree(node) == 0:
        if G.in_degree(node) > 1:
            return False
        if G.node[node]['coverage'] < coverage_thres: 
            return True
    return False

def is_in_tip(node, G, coverage_thres):
    if G.in_degree(node) == 0:
        if G.out_degree(node) > 1:
            return False
        if G.node[node]['coverage'] < coverage_thres: 
            return True
    return False

def is_inferior_node(node, G):
    assert G.out_degree(node) == 0 and G.in_degree(node) < 2
    if not G.predecessors(node):
        return False
    coverage = G.node[node]['coverage']
    pred = G.predecessors(node)[0]
    for succ in G.successors_iter(pred):
        if succ != node and G.node[succ]['coverage'] > coverage:
            return True
    return False  

# remove rectangles. a step in tour bus.
def remove_rectangles(rectangles, G):
    for begin, end in rectangles:
        paths = list(nx.all_simple_paths(G, begin, end))
        path_num = len(paths)
        if path_num <= 1:
            continue
        #assert path_num > 1
        weights = [get_path_weight(path, G) for path in paths]
        removed_indices = set()
        for i in xrange(path_num-1):
            for j in xrange(i+1, path_num):
                path1 = paths[i]
                weight1 = weights[i]
                path2 = paths[j]
                weight2 = weights[j]
                if path_is_redundant(path1, path2, G):
                    if weight1 >= weight2:
                        #print path1, path2, weight1, weight2
                        removed_indices.add(j)
                    else:
                        #print path1, path2, weight1, weight2
                        removed_indices.add(i)
        for index in removed_indices:
            remove_edges_from_path(paths[index], G)

def path_is_redundant(path1, path2, G):
    seq1 = get_seq_from_path(path1, G) 
    seq2 = get_seq_from_path(path2, G) 
    if len(seq1) != len(seq2): 
        return False
    avg_length = G.graph['avg_length']
    HAMMING_THRES = len(seq1) / avg_length * 2
    hamming = get_hamming_distance(seq1, seq2)
    return hamming <= HAMMING_THRES

def get_seq_from_path(path, G):
    assert path
    node = path[0]
    path_seq = G.node[node]['seq']
    for i in range(1, len(path)):
        node1 = path[i-1]
        node2 = path[i] 
        overlap = G[node1][node2]['overlap']
        seq = G.node[node2]['seq']
        path_seq += seq[overlap:]
    return path_seq

def remove_edges_from_path(path, G):
    for i in range(len(path)-1):
        if G.has_edge(path[i], path[i+1]):
            G.remove_edge(path[i], path[i+1])

def is_singleton_node(node, G):
    return G.out_degree(node, G) == 0 and \
           G.in_degree(node, G) == 0

# the weight of a path without a root.
def get_path_weight(path, G):
    assert path
    node = path[0]
    size = G.node[node]['coverage'] * G.node[node]['length']
    length = G.node[node]['length'] 
    node_num = len(path)
    for i in range(1, node_num):
        node1 = path[i-1]
        node2 = path[i] 
        overlap = G[node1][node2]['overlap']
        size += G.node[node2]['coverage'] * G.node[node2]['length']
        length = length + G.node[node2]['length'] - overlap
    return size / length

def get_rectangles(source_node, G):
    # the ending node and a pred of an alternative path.
    # the first element is the revisited node and the second one is its pred.
    joints = []
    q = Queue.Queue()
    visited = {}
    q.put(source_node)
    visited[source_node] = True
    while not q.empty():
        curr = q.get()
        for succ in G.successors(curr):
            if succ in visited:
                joints.append((succ, curr))
            else:
                visited[succ] = True
                q.put(succ)
    # the begin and end nodes of rectangles.
    rectangles = set()
    for curr, pred in joints:
        for other_pred in G.predecessors_iter(curr):
            if other_pred != pred:
                lca = get_lca(other_pred, pred, G)
                if lca:
                    rectangles.add((lca, curr))
    return rectangles

def load_rectangles(source_node, G):
    dir = 'Tmp'
    if not os.path.exists(dir):
        os.makedirs(dir)
    cPickle_file_name = '%s/%s_rectangles.cPickle' % (dir, G.graph['id'])
    if os.path.exists(cPickle_file_name):
        f = open(cPickle_file_name, 'Ur')
        return cPickle.load(f)
    else:
        rectangles = get_rectangles(source_node, G)
        f = open(cPickle_file_name, 'wb')
        cPickle.dump(rectangles, f, protocol=2)
        f.close()
        return rectangles

def get_all_preds(node_set, G):
    preds = set()
    for node in node_set:
        preds |= set(G.predecessors(node))
    return preds

# get all lca of two nodes.
# do not include these two nodes.
# return empty set if none exits.
def get_lca(node1, node2, G):
    lca = None
    all_preds1 = set()
    all_preds2 = set()
    all_preds1.add(node1)
    all_preds2.add(node2)
    preds1 = set(G.predecessors(node1))
    preds2 = set(G.predecessors(node2))
    while preds1 or preds2:
        all_preds1 |= preds1
        all_preds2 |= preds2
        ca_set = all_preds1 & all_preds2
        if ca_set:
            sorted_nodes_list = nx.topological_sort(G)
            ca_list = list(ca_set)
            lca = max(ca_list, key=lambda ca: sorted_nodes_list.index(ca))
            break
        else:
            preds1 = get_all_preds(preds1, G)
            preds2 = get_all_preds(preds2, G)
    return lca

# remove alternative paths that have the same sequence in rectangles.
def tour_bus(G):
    source_nodes = get_all_source_nodes(G)
    for source_node in source_nodes:
        rectangles = get_rectangles(source_node, G)
        if rectangles:
            #pprint(rectangles)
            remove_rectangles(rectangles, G)

def show_all_rectangles(G):
    source_nodes = get_source_nodes(G)
    for source_node in source_nodes:
        rectangles = get_rectangles(source_node, G)
        for begin, end in rectangles:
            paths = list(nx.all_simple_paths(G, begin, end))
            if len(paths) > 1:
                print '>', begin, end
                for path in paths:
                    print ' '.join(path)
                    for node in path:
                        cds_set = G.node[node]['cds_set']
                        print ' '.join(cds_set),
                    print

def set_all_nodes_cds_set(G):
    for node in G.nodes_iter():
        descs = G.node[node]['descs']
        cds_set = get_cds_set_from_descs(descs)
        G.node[node]['cds_set'] = cds_set
            

def show_all_node_cds_set(G):
    for node in G.nodes_iter():
        cds_set = G.node[node]['cds_set']
        print node, ' '.join(cds_set)

def show_all_edge_cds_set(G): 
    for node1, node2 in G.edges_iter():
        cds_set1 = G.node[node1]['cds_set']
        cds_set2 = G.node[node2]['cds_set']
        print '>', node1, node2
        print ' '.join(cds_set1)
        print ' '.join(cds_set2)

def show_all_contig_cds_set(contigs, G):
    for contig in contigs:
        print '>', contig.name, contig.cds_set, len(contig.tags), contig.avg_overlap
        continue
        print ' '.join(contig.nodes)
        for node in contig.nodes:
            cds_set = G.node[node]['cds_set']
            first = get_cds_from_desc(G.node[node]['descs'][0])
            last = get_cds_from_desc(G.node[node]['descs'][-1])
            print '(%s %s %s)' % (first, node, last),
        print

def show_region_contigs(contigs, G, cds_dict):
    #for contig in contigs:
    #    print contig.name, contig.cds_set, ' '.join(contig.nodes)
    region_cds_avg_overlap = get_region_cds_avg_overlap(contigs)
    compare_avg_overlap(region_cds_avg_overlap, cds_dict)

# predict cds avg overlap.
def get_region_cds_avg_overlap(contigs):
    region_cds_avg_overlap_dict = {}
    count_dict = {}
    READ_NUM_THRES = 5
    for contig in contigs:
        name = contig.name
        cds_set = contig.cds_set
        tag_num = len(contig.tags)
        avg_overlap = contig.avg_overlap
        #print name, tag_num, cds_set, avg_overlap
        if len(cds_set) < 3 and tag_num >= READ_NUM_THRES:
            for composite_cds in cds_set:
                for cds_name in composite_cds.split(','):
                    count_dict.setdefault(cds_name, 0)
                    count_dict[cds_name] += 1
                    region_cds_avg_overlap_dict.setdefault(cds_name, 0)
                    region_cds_avg_overlap_dict[cds_name] += avg_overlap
    for cds_name in region_cds_avg_overlap_dict:
        region_cds_avg_overlap_dict[cds_name] /= count_dict[cds_name]
    return region_cds_avg_overlap_dict

def compare_avg_overlap(region_cds_avg_overlap_dict, cds_dict):
    for cds_name in region_cds_avg_overlap_dict:
        if cds_name in cds_dict:
            region_avg_overlap = region_cds_avg_overlap_dict[cds_name]
            cds = cds_dict[cds_name]
            print cds_name, len(cds.profile), len(cds.aligned_profile), \
                  np.mean(cds.overlap_dist), region_avg_overlap

def show_all_contigs_node_features(contigs, G):
    for contig in contigs:
        show_contig_node_features(contig, G)

def show_contig_node_features(contig, G):
    print '>', contig.name, contig.cds_set, contig.shared_cds_set
    for node in contig.nodes:
        print '%s %s %d %d' % (node, G.node[node]['cds_set'], G.in_degree(node), G.out_degree(node))
    node_num = len(contig.nodes)
    for i in xrange(node_num-1):
        node1 = contig.nodes[i]
        for j in xrange(i+1, node_num):
            node2 = contig.nodes[j]
            pair_set = get_node_pair_set(node1, node2, G)
            if pair_set:
                print node1, node2, len(pair_set)
             
# get the pairs that are in two separate nodes.
def get_node_pair_set(node1, node2, G):
    members1 = G.node[node1]['members']    
    members2 = G.node[node2]['members']    
    pair_set = get_pairs_from_two_sets(members1, members2)
    return pair_set

def show_all_nodes_pair_num(G):
    all_nodes = G.nodes()
    node_num = len(all_nodes)
    for i in xrange(node_num-1):
        node1 = all_nodes[i]
        for j in xrange(i+1, node_num):
            node2 = all_nodes[j]
            pair_set = get_node_pair_set(node1, node2, G)
            if pair_set:
                print '>', node1, node2, len(pair_set)

def get_all_nodes_pair_num_dict(G):
    pair_num_dict = {}
    all_nodes = G.nodes()
    node_num = len(all_nodes)
    for i in xrange(node_num-1):
        node1 = all_nodes[i]
        for j in xrange(i+1, node_num):
            node2 = all_nodes[j]
            pair_set = get_node_pair_set(node1, node2, G)
            if pair_set:
                pair_num_dict[(node1, node2)] = len(pair_set)
    return pair_num_dict

def check_all_nodes_pair_set(G):
    all_nodes = G.nodes()
    node_num = len(all_nodes)
    for i in xrange(node_num-1):
        node1 = all_nodes[i]
        cds_set1 = G.node[node1]['cds_set']
        shared_set1 = get_shared_cds_set(cds_set1)
        for j in xrange(i+1, node_num):
            node2 = all_nodes[j]
            cds_set2 = G.node[node2]['cds_set']
            shared_set2 = get_shared_cds_set(cds_set2)
            pair_set = get_node_pair_set(node1, node2, G)
            if pair_set:
                if shared_set1 & shared_set2:
                    print '>', 'good', node1, cds_set1, shared_set1, node2, cds_set2, shared_set2
                else:
                    print '>', 'bad', node1, cds_set1, shared_set1, node2, cds_set2, shared_set2

# print out the names of the nodes that begin the cds conflict.
def get_contig_conflict_nodes(contig, G):
    conflict_nodes = ()
    node_num = len(contig.nodes)
    if node_num <= 1: 
        return conflict_nodes
    for i in xrange(1, node_num):
        node1 = contig.nodes[i]
        cds_set1 = G.node[node1]['cds_set']
        for j in xrange(0, i-1): # do not check neighboring nodes.
            node2 = contig.nodes[j]
            cds_set2 = G.node[node2]['cds_set']
            if not has_shared_set(cds_set1, cds_set2):
                return node2, node1
    return conflict_nodes

def has_shared_set(cds_set1, cds_set2):
    for cds1 in cds_set1:
        set1 = set(cds1.split('&'))
        for cds2 in cds_set2:
            set2 = set(cds2.split('&'))
            if set1 & set2:
                return True
    return False        

def show_all_contigs_conflict_nodes(contigs, G):
    for contig in contigs:
        print '>', contig.name, contig.descs[0], contig.descs[-1]
        print ' '.join(contig.nodes)
        conflict_nodes = get_contig_conflict_nodes(contig, G)
        print conflict_nodes

def show_contig_ranges(contigs, G):
    for contig in contigs:
        conflict_nodes = get_contig_conflict_nodes(contig, G)
        print '>', contig.name, contig.nodes, conflict_nodes
        begin_set = set()
        end_set = set()
        for desc in contig.descs[0].split('&'):
            cds = desc.split(',')[0]
            begin_pos = desc.split(',')[2]
            begin_set.add((cds, begin_pos))

        for desc in contig.descs[-1].split('&'):
            cds = desc.split(',')[0]
            end_pos = desc.split(',')[2]
            end_set.add((cds, end_pos))
        for cds, begin_pos in sorted(list(begin_set)):
            print '%s:%s' % (cds, begin_pos),
        print
        for cds, end_pos in sorted(list(end_set)):
            print '%s:%s' % (cds, end_pos),
        print

def show_covered_regions(contigs, G):
    for contig in contigs:
        print '> %s %s' % (contig.name, G.graph['domain'])
        begin_set = set()
        end_set = set()
        begin_cds_set = set()
        end_cds_set = set()         

        for desc in contig.descs[0].split('&'):
            cds = desc.split(',')[0]
            begin_pos = desc.split(',')[2]
            begin_set.add((cds, begin_pos))
            begin_cds_set.add(cds)

        for desc in contig.descs[-1].split('&'):
            cds = desc.split(',')[0]
            end_pos = desc.split(',')[2]
            end_set.add((cds, end_pos))
            end_cds_set.add(cds)
  
        begin_list = sorted(list(begin_set))
        end_list = sorted(list(end_set))
        for begin_cds, begin_pos in begin_list:
            if begin_cds not in end_cds_set:
                print begin_cds, begin_pos, '$', '$'
            else:
                for end_cds, end_pos in end_list:
                    if end_cds == begin_cds:
                        print begin_cds, begin_pos, end_cds, end_pos

# get the cds regions in the assembled contigs.
def get_contigs_cds_read_dict(contigs):
    contigs_cds_read_dict = {}
    for contig in contigs:
        name = contig.name
        contigs_cds_read_dict.setdefault(name, {})
        # reads that are not from cds.
        contigs_cds_read_dict[name]['random'] = set() 
        for desc in contig.descs:
            if ',' not in desc:
                for item in desc.split('&'):
                    contigs_cds_read_dict[name]['random'].add(item)
            else:
                for item in desc.split('&'):
                    info = item.split(',')
                    cds = info[0]
                    strand = info[1]
                    begin_pos = int(info[2])
                    end_pos = int(info[3])
                    contigs_cds_read_dict[name].setdefault(cds, set())
                    contigs_cds_read_dict[name][cds].add((begin_pos, end_pos))  
    return contigs_cds_read_dict
        
def get_contigs_covered_cds_dict(contigs_cds_read_dict):
    contigs_covered_cds_dict = {}
    for name in contigs_cds_read_dict:
        contigs_covered_cds_dict[name] = {}
        for cds in contigs_cds_read_dict[name]:
            if cds != 'random':
                pos_list = sorted(list(contigs_cds_read_dict[name][cds]))
                cds_begin = pos_list[0][0]
                cds_end = pos_list[-1][1]
                contigs_covered_cds_dict[name][cds] = (cds_begin, cds_end)
    return contigs_covered_cds_dict
    
def output_contigs_covered_cds_dict(contigs_covered_cds_dict, contigs):
    for contig in contigs:
        name = contig.name
        nodes = contig.nodes
        print '>', name, ','.join(nodes)
        for cds in contigs_covered_cds_dict[name]:
            assert cds != 'random'
            cds_begin, cds_end = contigs_covered_cds_dict[name][cds]
            print cds, cds_begin, cds_end 

# get cds covered by all contigs
def get_cds_covered_dict(contigs_covered_cds_dict):
    cds_covered_dict = {}
    for name in contigs_covered_cds_dict:
        for cds in contigs_covered_cds_dict[name]:
            cds_covered_dict.setdefault(cds, set())
            cds_covered_dict[cds].add(contigs_covered_cds_dict[name][cds])
    for cds in cds_covered_dict:
        cds_covered_dict[cds] = sorted(list(cds_covered_dict[cds]))
    return cds_covered_dict

def output_cds_covered_dict(cds_covered_dict):
    for cds in cds_covered_dict:
        print cds,
        for begin_pos, end_pos in cds_covered_dict[cds]:
            print '(%d,%d)' % (begin_pos, end_pos),
        print

# merge covered cds regions if they have overlap.
def get_merged_cds_covered_dict(cds_covered_dict):
    # random reads will not be included.
    merged_cds_covered_dict = {}
    for cds in cds_covered_dict:
        merged_regions = get_merged_regions(cds_covered_dict[cds])
        merged_cds_covered_dict[cds] = merged_regions
    return merged_cds_covered_dict

def output_merged_cds_covered_dict(merged_cds_covered_dict):
    for cds in merged_cds_covered_dict:
        print cds,
        for begin_pos, end_pos in merged_cds_covered_dict[cds]:
            print '(%d,%d)' % (begin_pos, end_pos),
        print
            
def get_merged_regions(regions):
    if len(regions) == 1:
        return regions[:]
    OVERLAP_THRES = 1
    merged_regions = []
    sorted_regions = sorted(regions)
    curr_begin = regions[0][0]
    curr_end = regions[0][1] 
    for i in xrange(1, len(sorted_regions)):
        begin_pos = sorted_regions[i][0]
        end_pos = sorted_regions[i][1]
        overlap = get_pos_overlap_length(curr_begin, curr_end, begin_pos, end_pos)
        if overlap >= OVERLAP_THRES:
            curr_end = max(curr_end, end_pos)
        else:
            merged_regions.append((curr_begin, curr_end))
            curr_begin = begin_pos
            curr_end = end_pos
    merged_regions.append((curr_begin, curr_end))
    return merged_regions  

# get the cds shared by all members in the same cds set.
def get_shared_cds_set(cds_set):
    shared_set = set()
    if cds_set:
        cds_list = list(cds_set)
        shared_set = set(cds_list[0].split('&'))
        for i in range(len(cds_list)):
            shared_set &= set(cds_list[i].split('&'))
    return shared_set
   
def draw_subgraphs(G):
    domain = G.graph['domain']
    folder_name = 'Graph/'
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    i = 1
    for subgraph in nx.weakly_connected_component_subgraphs(G):
        file_name = '%s/%s_subgraph%d.png' % (folder_name, domain, i)
        i += 1
        plt.figure(figsize=(16, 12))
        pos = nx.spring_layout(subgraph)
        nx.draw_networkx_nodes(subgraph, pos, alpha=0.5, node_color='r')
        plt.axis('off')
        plt.savefig(file_name)

# test whether path1 is subpath of path2.
def is_subpath(path1, path2):
    length1 = len(path1)
    length2 = len(path2)
    if length1 > length2:
        return False
    for i in xrange(length2-length1+1):
        if path1 == path2[i:i+length1]:
            return True
    return False

'''
# get the nodes that characterize alternative paths.
def get_target_nodes(path, rectangles, G):
    target_nodes = []
    for begin, end in rectangles:
        for connected_path in nx.all_simple_paths(G, begin, end):
            if is_subpath(connected_path, path):
                target_nodes.append(tuple(connected_path[1:-1]))
    sink = path[-1]
    for pred in G.predecessors_iter(sink):
      if G.out_degree(pred) > 1:
          target_nodes.append((sink,))
    target_nodes.sort(key=lambda nodes: path.index(nodes[0]))
    return target_nodes
'''

# connect target nodes that have pairs into a single graph g.
def get_target_nodes_graph(target_nodes, pair_num_dict, pair_num_thres):
    g = nx.Graph()
    target_num = len(target_nodes)
    if target_num == 0:
        return g
    else:
        for target_node in target_nodes:
            g.add_node(target_node)

    for i in xrange(target_num-1):
        nodes1 = target_nodes[i]
        for j in xrange(i+1, target_num):
            nodes2 = target_nodes[j]
            pair_num = 0
            for node1 in nodes1:
                for node2 in nodes2:
                    pair_num += pair_num_dict.get((node1, node2), 0)
            if pair_num >= pair_num_thres:
                #print nodes1, nodes2, pair_num
                g.add_edge(nodes1, nodes2)
    return g

# decide whether a path is supported by pairs.
def is_path_supported_by_pairs(target_nodes, G, pair_num_dict):
    g = get_target_nodes_graph(target_nodes, pair_num_dict, 1)
    if g.number_of_nodes() == len(target_nodes):
        return True
    else:
        return False

# select paths based on pairs.
def select_paths(paths, rectangles, G, pair_num_dict):
    all_paths_target_nodes = [(path, get_target_nodes(path, rectangles, G)) 
                               for path in paths]
    selected_paths = []
    unpaired_paths = []
    for path, target_nodes in all_paths_target_nodes:
        if len(target_nodes) < 2:
            selected_paths.append(path)
        if is_path_supported_by_pairs(target_nodes, G, pair_num_dict):
            selected_paths.append(path)
        else:
            unpaired_paths.append(path)
    return selected_paths

def get_contig_target_node_graphs(contigs, rectangles, G, pair_num_dict):
    contig_target_node_graphs = {}
    for contig in contigs:
        name = contig.name
        path = contig.path
        target_nodes = get_target_nodes(path, rectangles, G)
        g = get_target_nodes_graph(target_nodes, pair_num_dict, 1)
        contig_target_node_graphs[name] = g
    return contig_target_node_graphs

def get_repr_contigs(contig_fasta_file, contigs):
    repr_contigs = []
    repr_contig_names = []
    cdhit_file_name = contig_fasta_file + '.cdhit.fna'
    cmd = 'cd-hit-est -i %s -o %s' % (contig_fasta_file,
                                      cdhit_file_name)
    dummy_file = open('/dev/null')
    subprocess.call(cmd, stdout=dummy_file, shell=True)
    with open(cdhit_file_name, 'Ur') as f:
        for record in SeqIO.parse(f, 'fasta'):
            contig_name = '_'.join(record.id.split('_')[:2])   
            repr_contig_names.append(contig_name)
    for contig in contigs:
        if contig.name in repr_contig_names:
            repr_contigs.append(contig)
    return repr_contigs

def get_cds_domain_region_mapping_dict(in_file_name, target_domain):
    cds_domain_region_mapping_dict = {}
    with open(in_file_name, 'Ur') as f:
        for line in f:
            if line.strip():
                items = line.rstrip().split()
                cds = items[0]
                domain = items[1]
                if domain != target_domain:
                    continue
                begin_pos = int(items[2])
                end_pos = int(items[3])
                mapping_coverage = float(items[4])
                nscore = float(items[5])
                cds_domain_region_mapping_dict.setdefault(cds, {})    
                cds_domain_region_mapping_dict[cds].setdefault(domain, {})    
                cds_domain_region_mapping_dict[cds][domain][(begin_pos, end_pos)] = (mapping_coverage, nscore)    
    return cds_domain_region_mapping_dict

# this dic only has non-zero coverage regions.
def get_cds_domain_coverage_dict(cds_domain_mapping_region_dict, 
                                 merged_cds_covered_dict, 
                                 target_domain, nscore_thres):
    cds_domain_coverage_dict = {}
    for cds in cds_domain_mapping_region_dict:
        for domain in cds_domain_mapping_region_dict[cds]:
            if domain != target_domain:
                continue
            for region in cds_domain_mapping_region_dict[cds][domain]:
                mapping_coverage = cds_domain_mapping_region_dict[cds][domain][region][0]
                nscore = cds_domain_mapping_region_dict[cds][domain][region][1]
                if nscore >= nscore_thres:
                    cds_domain_coverage_dict.setdefault(cds, {})
                    cds_domain_coverage_dict[cds].setdefault(domain, {})
                    coverage = 0.0
                    if cds in merged_cds_covered_dict:
                        covered_regions = merged_cds_covered_dict[cds]            
                        coverage = get_region_coverage(region, covered_regions)
                    cds_domain_coverage_dict[cds][domain][region] = coverage
    return cds_domain_coverage_dict
                
def get_region_coverage(region, covered_regions):
    region_begin = region[0]
    region_end = region[1]
    region_length = region_end - region_begin + 1
    flags = [0] * region_length 
    for begin_pos, end_pos in covered_regions:
        overlap = get_pos_overlap_length(region_begin, region_end,
                                         begin_pos, end_pos)                     
        if overlap > 0:
            for i in xrange(begin_pos, end_pos+1):
                rel_pos = i - region_begin
                if rel_pos >= 0 and rel_pos < region_length and flags[rel_pos] == 0:
                    flags[rel_pos] = 1
    covered_length = sum(flags)
    coverage = float(covered_length) / region_length
    #print region, covered_regions, coverage
    return coverage

def output_cds_domain_coverage_dict(cds_domain_coverage_dict, 
                                    cds_domain_region_mapping_dict,
                                    merged_cds_covered_dict):
    for cds in cds_domain_coverage_dict:
        for domain in cds_domain_coverage_dict[cds]:
            for region in cds_domain_coverage_dict[cds][domain]:
                coverage = cds_domain_coverage_dict[cds][domain][region]
                mapping_coverage, nscore = cds_domain_region_mapping_dict[cds][domain][region]
                covered_regions = merged_cds_covered_dict.get(cds, [])
                covered_region_num = len(covered_regions)
                print cds, domain, region[0], region[1], \
                      mapping_coverage, nscore, covered_region_num, coverage

# get the portion of different cds in contigs.
def get_contig_cds_portion(contigs, read_desc_dict):
    contig_cds_count = {}
    contig_cds_portion = {}
    for contig in contigs:
        name = contig.name
        total_desc_num = 0 # only consider reads from cds regions.
        contig_cds_count = {}
        for desc in contig.descs:
            # do not consider random reads.
            if ',' in desc:
                total_desc_num += 1
                cds_set = set(get_cds_from_desc(desc).split('&'))
                for cds in cds_set:
                    contig_cds_count.setdefault(cds, 0)
                    contig_cds_count[cds] += 1
        count_list = contig_cds_count.items()
        if count_list:
            portion = 0.0  
            max_cds, max_count = max(count_list, key=operator.itemgetter(1))
            total_read_num = len(contig.descs)
            for desc in contig.descs:
                if ',' in desc:
                    cds_set = set(get_cds_from_desc(desc).split('&'))
                    if max_cds in cds_set:
                        portion += 1             
            contig_cds_portion[name] = (total_desc_num, portion / total_read_num, max_cds)
        else:
            contig_cds_portion[name] = (0, 0.0, 'random') 
    return contig_cds_portion

def output_contigs_cds_portion(contig_cds_portion, contig_dict):
    for name in contig_cds_portion:
        contig = contig_dict[name]
        descs = contig.descs
	print '>', name, contig_cds_portion[name][0], contig_cds_portion[name][1], \
              contig_cds_portion[name][2]

def get_contig_pair_num_dict(contigs):
    contig_pair_num_dict = {}
    contig_num = len(contigs)
    for i in xrange(contig_num-1):
        for j in xrange(i+1, contig_num):
            contig1 = contigs[i]        
            contig2 = contigs[j]        
            set1 = contig1.members
            set2 = contig2.members
            pairs = get_pairs_from_two_sets(set1, set2)
            pair_num = len(pairs)
            if pair_num > 0:
                contig_pair_num_dict[(contig1.name, contig2.name)] = pair_num
                contig_pair_num_dict[(contig2.name, contig1.name)] = pair_num
    return contig_pair_num_dict

def show_all_contig_descs(contigs):
    for contig in contigs:
        name = contig.name
        cds_set = contig.cds_set  
        descs = contig.descs
        print '>', name, cds_set
        pprint(contig.nodes)
        pprint(descs)

def get_avg_overlap(compressed_read_dict):
    avg_overlap = 0
    if compressed_read_dict:
        key = compressed_read_dict.keys()[0]
        read_length = compressed_read_dict[key].length
        read_num = len(compressed_read_dict)
        total_size = read_num * read_length
        min_begin_state = compressed_read_dict[key].begin_state
        max_end_state = compressed_read_dict[key].end_state
        for tag in compressed_read_dict:
            begin_state = compressed_read_dict[tag].begin_state            
            end_state = compressed_read_dict[tag].end_state            
            min_begin_state = min(min_begin_state, begin_state)
            max_end_state = max(max_end_state, end_state)
        domain_length = (max_end_state - min_begin_state) * 3
        #print min_begin_state, max_end_state, domain_length
        avg_overlap = int((total_size - domain_length) / (read_num - 1))
    return avg_overlap

def get_avg_read_length(compressed_read_dict):
    avg_read_length = 0.0
    for tag in compressed_read_dict:
        avg_read_length += compressed_read_dict[tag].length
    avg_read_length /= len(compressed_read_dict)
    return avg_read_length

def get_cds_dict(in_file_name, cds_length_dict, hits):
    cds_read_profile = {}
    cds_dict = {}
    with open(in_file_name, 'Ur') as f:
        for line in f:
            if line.strip():
                items = line.rstrip().split()
                read = items[0]
                cds_name = items[2]
                begin = int(items[3]) + 1
                length = len(items[4])
                end = begin + length - 1
                cds_read_profile.setdefault(cds_name, [])
                cds_read_profile[cds_name].append((read, begin, end))
        for cds_name in cds_read_profile:
            cds_read_profile[cds_name].sort(key=operator.itemgetter(1))
            profile = cds_read_profile[cds_name]
            cds = CDS(cds_name)
            cds_length = cds_length_dict[cds_name]
            cds.set_length(cds_length)
            cds.set_profile(profile)
            for read_name, begin, end in profile:
                cds.add_read(read_name)
            cds.set_coverage()
            cds.set_overlap_dist()
            cds.set_aligned_profile(hits)
            cds_dict[cds_name] = cds
    return cds_dict

def get_cds_length_dict(in_file_name):
    cds_length_dict = {}
    with open(in_file_name, 'Ur') as f:
        for record in SeqIO.parse(f, 'fasta'):
            cds = record.id
            cds_length = len(record.seq)
            cds_length_dict[cds] = cds_length
    return cds_length_dict

# get aigned tags in an aligned region.
def get_region_reads(target_read_dict, region):
    region_reads = {}
    for read_name in target_read_dict:
        read = target_read_dict[read_name]
        begin_state = read.begin_state
        end_state = read.end_state
        if begin_state >= region[0] and \
           end_state <= region[1]:
           region_reads[read_name] = read
    return region_reads

def get_compressed_region_reads(region_reads):
    prefix = 'regiontag'
    seq_read_dict = {}
    for read_name in region_reads:
        seq = region_reads[read_name].seq
        seq_read_dict.setdefault(seq, [])
        seq_read_dict[seq].append(read_name)
    compressed_region_reads = {}
    index = 1
    for seq in seq_read_dict:
        tag_name = '%s%d' % (prefix, index)
        index += 1
        # use the first read that has the same sequence.
        first_read_name = seq_read_dict[seq][0]
        first_read = region_reads[first_read_name]
        tag = Read(tag_name)
        tag.set_seq(first_read.seq)
        tag.set_begin_state(first_read.begin_state)
        tag.set_end_state(first_read.end_state)
        tag.set_score(first_read.score)
        # keep the information
        tag_size = len(seq)
        tag.set_size(tag_size)
        for read_name in seq_read_dict[seq]:
            read_descs = region_reads[read_name].descs
            tag.add_member(read_name)
            for desc in read_descs:
                tag.add_desc(desc)
        compressed_region_reads[tag_name] = tag
    return compressed_region_reads

def show_compressed_region_reads_cds(compressed_region_reads, repr_cds_dict):
    cds_read_dict = {}
    for read_name in compressed_region_reads:
        read = compressed_region_reads[read_name]
        for desc in read.descs:
            cds_name = desc.split(',')[0]
            cds_read_dict.setdefault(cds_name, set())     
            cds_read_dict[cds_name].add(read_name)
    for cds_name in cds_read_dict:
        if cds_name in repr_cds_dict:
            print cds_name, len(cds_read_dict[cds_name])

# all read names have trailing "$+" or "$-".
def get_other_read_name(read_name):
    assert read_name[-1] == '+' or read_name[-1] == '-'
    assert read_name[-3] == '1' or read_name[-3] == '2'
    assert read_name[-4] == '.'
    base_name = read_name[:-4]
    if read_name[-1] == '+':
        strand = '-'
    else:
        strand = '+'
    if read_name[-3] == '1':
        end = '2'
    else:
        end = '1'
    return base_name + '.%s$%s' % (end, strand)

# remove reads that are not present in the fasta file.
def remove_empty_other_end(other_target_read_dict):
    for read_name in other_target_read_dict.keys():
        if not other_target_read_dict[read_name].seq:
            del other_target_read_dict[read_name] 


def get_repr_cds(cds_dict):
    g = nx.DiGraph()
    cds_num = len(cds_dict)
    cds_list = cds_dict.keys()
    for cds_name in cds_list:
        g.add_node(cds_name)

    for i in xrange(cds_num-1):
        cds_name1 = cds_list[i]
        set1 = cds_dict[cds_name1].read_set
        for j in xrange(i+1, cds_num):
            cds_name2 = cds_list[j]
            set2 = cds_dict[cds_name2].read_set
            if set1 <= set2:
                g.add_edge(cds_name1, cds_name2)
            elif set1 > set2:
                g.add_edge(cds_name2, cds_name1)
    cds_num = g.number_of_nodes()
    sinks = get_sink_nodes(g)
    return sinks

# keep repr cds or non-mapped reads.
def keep_repr_descs(compressed_region_reads, repr_cds):
    for tag_name in compressed_region_reads:
        tag = compressed_region_reads[tag_name]
        descs = copy.deepcopy(tag.descs)
        for desc in descs:
            if ',' not in desc:
                continue
            cds_name = desc.split(',')[0]
            if cds_name not in repr_cds:
                tag.descs.remove(desc)

def load_overlap_graph(compressed_read_dict, target_domain, 
                       overlap_thres, relative_diff_thres):
    return G

def load_other_overlap_graph(compressed_read_dict, target_domain, overlap_thres):
    G = build_other_overlap_graph(compressed_read_dict,
                                  target_domain, overlap_thres)
    return G

# calculate relative overlap difference.
def get_overlap_diff(state_overlap, overlap):
    diff = abs(float(state_overlap-overlap)) / state_overlap   

# consider state information.
def build_simple_region_overlap_graph(compressed_region_reads, target_domain, 
                                      overlap_thres, relative_diff_thres):
    region_read_num = len(compressed_region_reads)
    region_read_list = compressed_region_reads.values()
    region_read_list.sort(key=lambda read: (read.begin_state, read.end_state))
    # reads that have close state begin positions will be tested for both direction.
    STATE_THRES = 3
    G = nx.DiGraph()
    for read in region_read_list:
        G.add_node(read.name,
                   begin_state=read.begin_state,
                   end_state=read.end_state, 
                   score=read.score,
                   seq=read.seq,
                   length=read.length,
                   coverage=float(read.size)/read.length,
                   members=copy.deepcopy(read.members), #original read names.
                   tags=[read.name]) # only consider single node or tags.
    for i in xrange(region_read_num-1):
        read1 = region_read_list[i]
        read1_length = len(read1.seq)
        for j in xrange(i+1, region_read_num):
            read2 = region_read_list[j]
            read2_length = len(read2.seq)
            begin_state1, end_state1 = read1.begin_state, read1.end_state
            begin_state2, end_state2 = read2.begin_state, read2.end_state
            # state overlap is in bp.
            state_overlap = get_pos_overlap_length(begin_state1, end_state1,
                                                   begin_state2, end_state2) * 3
            if state_overlap <= overlap_thres:
                break
            # handle the case that two reads have similar starting state.
            if abs(begin_state1 - begin_state2) <= STATE_THRES:
                seq_overlap_length1, h1 = get_seq_overlap_length(read1.seq, read2.seq)
                seq_overlap_length2, h2 = get_seq_overlap_length(read2.seq, read1.seq)
                if max(seq_overlap_length1, seq_overlap_length2) >= overlap_thres:
                    if seq_overlap_length1 >= seq_overlap_length2:  
                        overlap = seq_overlap_length1
                        diff = get_overlap_diff(state_overlap, overlap) 
                        # weight is the overlap.
                        if overlap >= overlap_thres and diff <= relative_diff_thres:
                            G.add_edge(read1.name, read2.name, 
                                       overlap=overlap, hamming=h1, state_overlap=state_overlap)
                    else:
                        overlap = seq_overlap_length2
                        diff = get_overlap_diff(state_overlap, overlap) 
                        # weight is the overlap.
                        if diff <= relative_diff_thres:
                            G.add_edge(read2.name, read1.name, 
                                       overlap=overlap, hamming=h2, state_overlap=state_overlap)
                # add an edge between read1 and read2 if they have significant overlap.
                # weight of the edge will be score of read2, which is the child.
            else:
                seq_overlap_length, h = get_seq_overlap_length(read1.seq, read2.seq)
                if seq_overlap_length >= overlap_thres:
                    overlap = seq_overlap_length
                    diff = get_overlap_diff(state_overlap, overlap)
                    if diff <= relative_diff_thres:
                        G.add_edge(read1.name, read2.name, 
                                   overlap=overlap,
                                   hamming=h, state_overlap=overlap)
    return G

def get_transitive_edges(G):
    transitive_edges = set()
    edges = G.edges()
    for begin_node, end_node in edges:
        current_edge_data = G.edge[begin_node][end_node]
        G.remove_edge(begin_node, end_node)
        if nx.has_path(G, begin_node, end_node):
            transitive_edges.add((begin_node, end_node))
            print begin_node, end_node
        else: 
            G.add_edge(begin_node, end_node, current_edge_data)

def process_overlap_graph(target_domain, G):
    set_overlap_graph_parameters(target_domain, G)
    remove_cycles(G)
    combine_single_paths(G)
    remove_graph_tips(G)
    #tour_bus(G)
    combine_single_paths(G)

def get_transitive_dict(transitive_nodes, G):
    transitive_dict = {}
    nodes = G.nodes()
    node_num = len(nodes)
    for i in xrange(node_num-1):
        node1 = nodes[i]
        for j in xrange(i+1, node_num):
            node2 = nodes[j]
            if has_transitive_edge(G, node1, node2, transitive_nodes):
                transitive_dict[(node1, node2)] = True
    return transitive_dict

def has_transitive_edge(G, node1, node2, transitive_nodes):
    single_nodes1 = G.node[node1]['tags']
    single_nodes2 = G.node[node2]['tags']
    for single_node1 in single_nodes1:
        for single_node2 in single_nodes1:
            if (single_node1, single_node2) in transitive_nodes or \
               (single_node2, single_node1) in transitive_nodes: 
                return True
    return False 

def get_supports(transitive_dict, pair_num_dict):
    supports = {}
    for node1, node2 in transitive_dict:
        supports[(node1, node2)] = True
    for node1, node2 in pair_num_dict:
        supports[(node1, node2)] = True
    return supports

def get_target_node_set(target_nodes):
    target_node_set = set()
    for nodes in target_nodes:
        for node in nodes:
            target_node_set.add(node)
    return target_node_set

def has_support(curr_node, node, supports):
    return (curr_node, node) in supports or \
           (node, curr_node) in supports

def get_target_nodes(rectangles, G):
    target_nodes = []
    for begin, end in rectangles:
        for connected_path in nx.all_simple_paths(G, begin, end):
            target_nodes.append(tuple(connected_path[1:-1]))
    sinks = get_sink_nodes(G)
    for sink in sinks:
        for pred in G.predecessors_iter(sink):
          if G.out_degree(pred) > 1 and G.in_degree(pred) > 1:
              target_nodes.append((sink,))
    return target_nodes

def get_target_node_set(target_nodes):
    target_node_set = set()
    for nodes in target_nodes:
        for node in nodes:
            target_node_set.add(node)
    return target_node_set

def is_valid_path_node(curr_node, path, target_nodes,
                       target_node_set, supports):
    if curr_node not in target_node_set:
        return True
    elif not set(path) & target_node_set:
        return True
    else:
        # if any of its previous nodes in the same rectangle is added.
        for nodes in target_nodes:
            if curr_node in nodes and nodes[0] != curr_node:
                return True
        for nodes in target_nodes:
            if curr_node not in nodes:
                for node in nodes:
                    if node in path and has_support(curr_node, node, supports):
                        return True
    return False

def dfs(G, curr_node, path, paths,
        target_nodes, target_node_set, supports):
    if not is_valid_path_node(curr_node, path, target_nodes,
                              target_node_set, supports):
        return
    #G.node[curr_node]['color'] += 1
    path.append(curr_node)
    if G.out_degree(curr_node) == 0:
        paths.append(path[:])
    else:
        for succ in G.successors_iter(curr_node):
            dfs(G, succ, path, paths,
                target_nodes, target_node_set, supports)
    path.pop()

def show_supports(supports):
    for node1, node2 in supports:
        print node1, node2, supports[(node1, node2)]

def get_other_repr_contigs_from_overlap_graph(G, target_domain):
    transitive_nodes = remove_redundant_edges(G)
    process_overlap_graph(target_domain, G)
    transitive_dict = get_transitive_dict(transitive_nodes, G)
    add_root_to_subgraph(G)
    paths = get_all_paths('root', G)
    contigs = get_contigs_from_paths(paths, G)
    return contigs

def get_repr_contigs_from_overlap_graph(G, target_domain):
    transitive_nodes = remove_redundant_edges(G)
    process_overlap_graph(target_domain, G)
    transitive_dict = get_transitive_dict(transitive_nodes, G)
    all_nodes_pair_num_dict = get_all_nodes_pair_num_dict(G)
    supports = get_supports(transitive_dict, all_nodes_pair_num_dict)
    keep_single_out_edge(G, supports)
    keep_single_in_edge(G, supports)
    add_root_to_subgraph(G)
    rectangles = get_rectangles('root', G)
    target_nodes = get_target_nodes(rectangles, G)
    target_node_set = get_target_node_set(target_nodes)
    path = []
    paths = []
    curr_node = 'root'
    dfs(G, curr_node, path, paths, target_nodes, target_node_set, supports)
    contigs = get_contigs_from_paths(paths, G)
    return contigs
   
def dump_contigs(contig_dict, out_file_name):
    contigs = contig_dict.items()
    contigs.sort(key=lambda item: int(item[0].split("contig")[1]))
    with open(out_file_name, 'wb') as f:
        for name, contig in contigs:
            record = SeqRecord(Seq(contig.seq), id=contig.name, description='')
            SeqIO.write(record, f, 'fasta')

def get_renamed_contig_dict(contig_dict):
    index = 1
    renamed_contig_dict = {}
    for name in contig_dict:
        contig = contig_dict[name]
        new_name = 'contig' + str(index)
        contig.set_name(new_name)
        renamed_contig_dict[new_name] = contig
        index += 1
    return renamed_contig_dict

def get_scaffold_list(contig_dict):
    scaffold_list = [] # find regions that are connected by pairs.
    G = nx.Graph()
    name_list = sorted(contig_dict.keys())
    contig_num = len(name_list)
    for i in xrange(contig_num-1):
        name1 = name_list[i]
        contig1 = contig_dict[name1]
        members1 = contig1.members
        reads1 = get_original_set(members1)
        for j in xrange(i+1, contig_num-1):
            name2 = name_list[j]
            contig2 = contig_dict[name2]
            members2 = contig2.members
            reads2 = get_original_set(members2)
            pairs = get_pairs_from_two_original_sets(reads1, reads2)
            if len(pairs) > 0:
                G.add_edge(name1, name2, pair_num=len(pairs))
    for nodes in nx.connected_components(G):
        scaffold_list.append(nodes)
    return scaffold_list

def get_original_set(members):
    sample = list(members)[0]
    if sample[-1] != '+' and sample[-1] != '-':
        return copy.copy(members)
    else:
        reads = set([member[:-2] for member in members])
        return reads

# the read names end with '+' or '-'
def get_pairs_from_two_original_sets(read_set1, read_set2):
    pairs = set()
    for read in read_set1:
        if read == 'root':
            continue
        name = read
        assert name[-1] == '1' or name[-1] == '2'
        if name[-1] == '1':
            other = name[:-1] + '2'
            if other in read_set2:
                pairs.add(name[:-2])
        else:
            other = name[:-1] + '1'
            if other in read_set2:
                pairs.add(name[:-2])
    return pairs

def output_scaffolds(scaffold_list, out_file_name):
    with open(out_file_name, 'wb') as f:
        i = 1
        for nodes in scaffold_list:
            f.write('>%s%d\n' % ('scaffold', i))
            f.write(' '.join(nodes)+'\n')
            i += 1

 
def main(): 
    if len(sys.argv) != 7:
        print >> sys.stderr, 'Usage: <alignment file> ' \
                             '<read fasta file> <target family> ' \
                             '<alignment thres> <relative diff thres> <output folder>' 
        sys.exit(2)
    alignment_file_name = sys.argv[1]
    fasta_file_name = sys.argv[2]
    target_domain = sys.argv[3]
    overlap_thres = int(sys.argv[4])
    diff_thres = float(sys.argv[5])
    out_dir = sys.argv[6]
    NSCORE_THRES = 0
    EVALUE_THRES = 1000
    aligned_read_dict = get_hmmer_aligned_read_dict(alignment_file_name, 
                                                    target_domain,
                                                    EVALUE_THRES)
    target_read_dict = aligned_read_dict
    if not aligned_read_dict:
        sys.stderr('No read belongs to the gene family! No contig is generated!')
        sys.exit(2)
    set_read_seq(fasta_file_name, target_read_dict)
    compressed_read_dict = get_compressed_read_dict(target_read_dict, 'tag')
    G = build_simple_region_overlap_graph(compressed_read_dict,
                                          target_domain, 
                                          overlap_thres,
                                          diff_thres)
    G.graph['id'] = 'graph'
    G.graph['contig_basename'] = 'contig'
    G.graph['out_dir'] = out_dir
    # handle the other end.
    other_target_read_dict = get_other_target_read_dict(target_read_dict)
    set_read_seq(fasta_file_name, other_target_read_dict)
    remove_empty_other_end(other_target_read_dict)
    other_compressed_read_dict = get_compressed_read_dict(other_target_read_dict,
                                                          'othertag')
    other_G =  build_other_overlap_graph(other_compressed_read_dict, 
                                         target_domain, 
                                         overlap_thres)
    other_G.graph['id'] = 'othergraph'
    other_G.graph['contig_basename'] = 'othercontig'
    other_G.graph['out_dir'] = out_dir
    # get both repr contigs.
    repr_contigs = []
    other_repr_contigs = []
    if G.number_of_nodes() > 0:
        repr_contigs = get_repr_contigs_from_overlap_graph(G, target_domain)
    if other_G.number_of_nodes() > 0:
        other_repr_contigs = get_other_repr_contigs_from_overlap_graph(other_G, target_domain) 
    both_repr_contigs = repr_contigs + other_repr_contigs
    if not both_repr_contigs:
        return
    both_repr_contig_dict = dict([contig.name, contig] for contig in both_repr_contigs)
    renamed_contig_dict = get_renamed_contig_dict(both_repr_contig_dict)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_file_name = '%s/%s_contigs.fa' % (out_dir, target_domain)
    dump_contigs(renamed_contig_dict, out_file_name)
    scaffold_list = get_scaffold_list(renamed_contig_dict)
    scaffold_file_name = '%s/%s_scaffolds.txt' % (out_dir, target_domain)
    output_scaffolds(scaffold_list, scaffold_file_name)
  

if __name__ == '__main__':
    main()
