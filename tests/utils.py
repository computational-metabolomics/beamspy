#!/usr/bin/env python
#  -*- coding: utf-8 -*-
import os
import sqlite3


def to_test_data(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_data", *args)

def to_test_results(*args):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_results", *args)

def sqlite_records(db, table):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    cursor.execute("select * from {}".format(table))
    records = cursor.fetchall()
    conn.close()
    return records

def sqlite_count(db, table):
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    cursor.execute("select count(*) from {}".format(table))
    records = cursor.fetchone()[0]
    conn.close()
    return records