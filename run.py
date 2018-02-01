# hello_world.py

from search import read_table_from, match_dna

from flask import Flask, request, jsonify
app = Flask(__name__)


table = read_table_from("GRCh38_no_alts.2bit.M10.Q10.index")


@app.route('/')
def hello_world():
    return 'Hello World!'


@app.route('/search/')
def search():
    query = "GTAATCTTAGCACTTTGGGAGGCGGAGACGGATGTATCGCTTGAGCTCAGGAGTTGAAGACCAGCCTGGGCAACATACTGAGACTCCGTCTTGTATAATTTAATTAAAATTTAAAAAAAGAAGAGAAAAAGACCTGTGTT"
    if 'query' in request.args:
        query = request.args['query']
    return jsonify([{
        "chromosome": match[0],
        "start": match[1],
        "end": match[1] + len(query),
    } for match in match_dna(table, query)])
