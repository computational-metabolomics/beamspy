#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sqlite3
from beamspy import statistics
import networkx as nx


def group_features(df, db_out, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="pearson", positive=True, block=5000, ncpus=None):

    conn = sqlite3.connect(db_out)
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS groups")

    cursor.execute("""CREATE TABLE groups (
                   group_id INTEGER DEFAULT NULL,
                   peak_id_a TEXT DEFAULT NULL,
                   peak_id_b TEXT DEFAULT NULL,
                   degree_a INTEGER DEFAULT NULL,
                   degree_b INTEGER DEFAULT NULL,
                   r_value REAL DEFAULT NULL,
                   p_value REAL DEFAULT NULL,
                   rt_diff REAL DEFAULT NULL,
                   mz_diff REAL DEFAULT NULL,                 
                   PRIMARY KEY (peak_id_a, peak_id_b));""")

    df_coeffs = statistics.correlation_coefficients(df, max_rt_diff, coeff_thres, pvalue_thres, method, positive, block, ncpus)
    graph = statistics.correlation_graphs(df_coeffs, df)
    sub_graphs = list(graph.subgraph(c) for c in nx.weakly_connected_components(graph))
    for i in range(len(sub_graphs)):
        sub_graphs[i].graph["groupid"] = i + 1 # not stored in output - place holder
        sub_graph_edges = []
        # sort edges
        edges = sorted(sub_graphs[i].edges(data=True), key=lambda e: (e[0], e[1]))
        for edge in edges:
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
