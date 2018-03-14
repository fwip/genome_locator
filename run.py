#!/usr/bin/env python3
from flask import Flask, request, jsonify
from twobitreader import TwoBitFile

from search import GenomeSearcher

app = Flask(__name__)

search.reference = TwoBitFile("GRCh38_no_alts.2bit")
index = "GRCh38_no_alts.2bit.M10.Q10.index.hdf5"


@app.route('/search/')
def search_for_needle():
    query = ""
    if 'query' in request.args:
        query = request.args['query']
    search = GenomeSearcher(reference_file="GRCh38_no_alts.2bit")
    return jsonify([{
        "chromosome": match[0],
        "start": match[1],
        "end": match[1] + len(query),
    } for match in search.match_file(index, query)])


if __name__ == "__main__":
    app.run(host="0.0.0.0")
