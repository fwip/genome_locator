# hello_world.py

from search import read_table_from, match_dna

from flask import Flask, request
app = Flask(__name__)


table = read_table_from("GRCh38_no_alts.index.pickle")


@app.route('/')
def hello_world():
    return 'Hello World!'


@app.route('/search/')
def search():
    query = "GTAATCTTAGCACTTTGGGAGGCGGAGACGGATGTATCGCTTGAGCTCAGGAGTTGAAGACCAGCCTGGGCAACATACTGAGACTCCGTCTTGTATAATTTAATTAAAATTTAAAAAAAGAAGAGAAAAAGACCTGTGTT"
    if 'query' in request.args:
        query = request.args['query']
    return str(match_dna(table, query))
