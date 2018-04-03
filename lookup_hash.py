#!/usr/bin/env python3

import numpy


class LookupHash(object):

    positions = []
    offsets = []
    counts = []

    def __init__(self, table=None, h5_group=None):
        if table is not None:
            self.import_table(table)
        if h5_group is not None:
            self.positions = h5_group['positions']
            self.offsets = h5_group['offsets']
            self.counts = h5_group['counts']

    def keys(self):
        return [i for i in range(len(self.counts)) if self.counts[i] != 0]

    def __contains__(self, key):
        if type(self.counts) is dict:
            return key in self.counts
        return key < len(self.counts) and self.counts[key] != 0

    def __getitem__(self, key):
        count = self.counts[key]
        if count == 0:
            return []
        offset = self.offsets[key]

        return [int(x) for x in self.positions[offset:offset+count]]

    @staticmethod
    def init_array(length, maximum):
        typ = None
        if maximum < 2**8:
            typ = numpy.uint8
        elif maximum < 2**16:
            typ = numpy.uint16
        elif maximum < 2**32:
            typ = numpy.uint32
        elif maximum < 2**64:
            typ = numpy.uint64
        return numpy.zeros((length,), dtype=typ)

    def import_table(self, table):

        key_set = set([row[0] for row in table])
        # Sort
        table = table[table.argsort(axis=0)][:, 0, :]

        key_set2 = set([row[0] for row in table])
        assert key_set == key_set2
        last_key = -1
        for row in table:
            key = row[0]
            assert last_key <= key
            last_key = key

        list_len = table.shape[0]
        array_len = int(max(key_set) + 1)
        len_max = 2 ** 60

        self.positions = table[:, 1]

        self.offsets = LookupHash.init_array(array_len, list_len)
        self.counts = LookupHash.init_array(array_len, len_max)

        offset = 0
        last_key = 0
        count = 0
        for row in table:
            key = row[0]
            if False:
                pos = row[1]
                assert self.positions[offset] == pos
            if key != last_key:
                self.offsets[last_key] = offset - count
                self.counts[last_key] = count
                last_key = key
                count = 0
            count += 1
            offset += 1
        if count:
            self.offsets[key] = offset - count
            self.counts[key] = count
