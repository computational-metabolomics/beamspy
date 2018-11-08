#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3
from beams import statistics
import networkx as nx


def group_features(df, db_out, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="pearson", block=5000, ncpus=None):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS groups")

    cursor.execute("""CREATE TABLE groups (
                   group_id int(11) DEFAULT NULL,
                   peak_id_a int(11) DEFAULT NULL,
                   peak_id_b int(11) DEFAULT NULL,
                   degree_a int(11) DEFAULT NULL,
                   degree_b int(11) DEFAULT NULL,
                   r_value decimal(12,2) DEFAULT NULL,
                   p_value decimal(12,2) DEFAULT NULL,
                   rt_diff decimal(12,2) DEFAULT NULL,
                   mz_diff decimal(12,2) DEFAULT NULL,                 
                   PRIMARY KEY (peak_id_a, peak_id_b));""")

    df_coeffs = statistics.correlation_coefficients(df, max_rt_diff, coeff_thres, pvalue_thres, method, block, ncpus)
    graph = statistics.correlation_graphs(df_coeffs, df)
    sub_graphs = list(graph.subgraph(c) for c in nx.weakly_connected_components(graph))

    for i in range(len(sub_graphs)):
        sub_graphs[i].graph["group_id"] = i + 1 # not stored in output - place holder
        sub_graph_edges = []
        for edge in sub_graphs[i].edges(data=True):
            sub_graph_edges.append((i+1,
                                    str(edge[0]), str(edge[1]),
                                    sub_graphs[i].degree(edge[0]), sub_graphs[i].degree(edge[1]),
                                    round(float(edge[2]["rvalue"]), 2), float(edge[2]["pvalue"]),
                                    float(edge[2]["rtdiff"]), float(edge[2]["mzdiff"])))
        cursor.executemany("""insert into groups (group_id, peak_id_a, peak_id_b, degree_a, degree_b,
                              r_value, p_value, rt_diff, mz_diff) values (?,?,?,?,?,?,?,?,?)""", sub_graph_edges)
    conn.commit()
    conn.close()
    return graph
