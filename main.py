import pandas as pd
import networkx as nx
from tqdm import tqdm

def rec_list2df(filename):
    '''
    Args:
        filename (str): name of a rec_list file whose rows follow the format
                        data_num bc1 bc2 bc3 bc4 bc5 count flag group_num

    Returns:
        df (dataframe): table (with headings from above) populated with data
                        from rec_list file
    '''

    df = pd.read_csv(filename, sep = '\t')
    df.columns = ['data_num', 'bc1', 'bc2', 'bc3', 'bc4', 'bc5', 'count',
                              'flag', 'group_num']
    return df

def seq2int(seq):
    '''
    Args:
        seq (str): DNA sequence to be converted to int

    Returns:
        integer representation of seq

    TODO:
        make this base 4 (or 5...)
    '''
    d = {"A":"1", "C":"2", "G":"3", "T":"4" }
    decstring = ""
    for i in range(len(seq)):
        decstring += d[seq[i]]
    return int(decstring)

def int2seq(n):
    '''
    Args:
        n (int): integer representation of seq

    Returns:
        seq (string): DNA sequence as string
    '''
    d = {"1":"A", "2":"C", "3":"G", "4":"T"}
    decstring = str(n)
    seq = ""
    for i in range(0, len(decstring)):
        seq += d[decstring[i]]
    return seq


def add_to_graph_from_recs(G, rec_df, num_bcs, as_int=False):
    '''
    Args:
        G (networkx graph object): graph, either empty or already populated
        rec_df (dataframe): dataframe with the records that hold the
                            adjacency information we would like to add
        num_bcs (int): the number of barcodes that will be found in rec_df

    Returns:
        None
    '''
    if as_int:
        bclist = ['bc1', 'bc2', 'bc3', 'bc4', 'bc5']
        print("Adding paths to graph...")
        for index, row in tqdm(rec_df.iterrows()):
            for i in range(num_bcs - 1):
                G.add_edge(seq2int(row[bclist[i]]), seq2int(row[bclist[i+1]]))

    bclist = ['bc1', 'bc2', 'bc3', 'bc4', 'bc5']
    print("Adding paths to graph...")
    for index, row in tqdm(rec_df.iterrows()):
        for i in range(num_bcs - 1):
            G.add_edge(row[bclist[i]], row[bclist[i+1]])



def save_graph(G, filename, fileformat='adjlist'):
    '''
    Args:
        G (networkx graph object): graph to be saved
        filename (string): destination filename
        fileformat (string): adjlist for adjacency, anything else for pickle

    Returns:
        None
    '''

    if fileformat == 'adjlist':
        nx.write_adjlist(G, filename)
    else:
        nx.write_gpickle(G, filename)

reclist = input("enter rec list filename: ")

bc4recs = rec_list2df(reclist)

graph4 = nx.Graph()

add_to_graph_from_recs(graph4, bc4recs, 4, as_int=True)

save_graph(graph4, 'g4.pickle', fileformat='pickle')
