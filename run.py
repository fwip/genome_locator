# hello_world.py
from flask import Flask, request, jsonify
from twobitreader import TwoBitFile

import search
app = Flask(__name__)

search.reference = TwoBitFile("GRCh38_no_alts.2bit")

index = "GRCh38_no_alts.2bit.M10.Q10.index.hdf5"


@app.route('/search/')
def search_for_needle():
    query = ""
    if 'query' in request.args:
        query = request.args['query']
    return jsonify([{
        "chromosome": match[0],
        "start": match[1],
        "end": match[1] + len(query),
    } for match in search.match_file(index, query)])
