# hello_world.py
from flask import Flask, request, jsonify
from twobitreader import TwoBitFile

import h5py

import search
app = Flask(__name__)


search.reference = TwoBitFile("GRCh38_no_alts.2bit")

index = "GRCh38_no_alts.2bit.M10.Q10.index.hdf5"

table = search.read_table_from_h5(index)

#table = {
#    chrom: search.read_table_from("{}.M10.Q10.index.gz".format(chrom))
#    for chrom in search.reference.keys()
#}


@app.route('/search/')
def search_for_needle():
    query = ""
    if 'query' in request.args:
        query = request.args['query']
    return jsonify([{
        "chromosome": match[0],
        "start": match[1],
        "end": match[1] + len(query),
    } for match in search.match_dna(table, query)])
