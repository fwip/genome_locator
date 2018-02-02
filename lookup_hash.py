#!/usr/bin/env python3

import numpy
import struct


class LookupHash(object):

    positions = []

    offsets = []
    counts = []

    def __init__(self, table):
        self.import_table(table)

    def keys(self):
        return [i for i in range(len(self.counts)) if self.counts[i] != 0]

    def __contains__(self, key):
        if type(self.counts) is dict:
            return key in self.counts
        return key < len(self.counts) and self.counts[key] != 0

    def __getitem__(self, key):
        offset = self.offsets[key]
        count = self.counts[key]

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
        offset = 0
        list_len = sum(len(v) for v in table.values())
        max_list = max(max(v) for v in table.values())
        array_len = max(table.keys()) + 1

        self.positions = LookupHash.init_array(list_len, max_list)

        # If we've got very sparse buckets, dicts will be more size-optimal
        occupancy = 1.0 * len(table.keys()) / array_len
        print("Occupancy is", occupancy)
        if occupancy < 0.25:
            self.offsets = {}
            self.counts = {}
        else:
            self.offsets = LookupHash.init_array(array_len, array_len)
            self.counts = LookupHash.init_array(array_len, array_len)

        for (key, value) in table.items():
            entry_count = len(value)
            self.offsets[key] = offset
            self.counts[key] = entry_count
            self.positions[offset:offset+entry_count] = value
            offset += entry_count


#    @staticmethod
#    def pack(subtable):
#        new_dict = {}
#        offset = 0
#        list_len = sum(len(v) for v in subtable.values())
#        array_len = max(subtable.keys()) + 1
#        offsets = numpy.zeros((array_len,), dtype=numpy.uint32)
#        counts = numpy.zeros((array_len,), dtype=numpy.uint32)
#
#        # Try smaller bins
#        # max_count = max(len(v) for v in subtable.values())
#        # if max_count < 2**16:
#        #     counts = numpy.zeros((array_len,), dtype=numpy.uint16)
#        # if max_count < 2**8:
#        #     counts = numpy.zeros((array_len,), dtype=numpy.uint8)
#
#        # print(list_len)
#        location_list = bytearray(4 * list_len)
#        for (key, value) in subtable.items():
#            entry_count = len(value)
#            offsets[key] = offset
#            counts[key] = entry_count
#            new_dict[key] = struct.pack('2I', offset, entry_count)
#            struct.pack_into('{}I'.format(entry_count),
#                             location_list, offset, *value)
#            offset += entry_count * 4
#
#        # key_count = len(new_dict)
#        # print("{buckets} buckets, {positions} positions: {avg} Avg positions, {size}MB Approx size, {bucketperc}% buckets filled".format(
#        #     buckets=key_count,
#        #     positions=list_len,
#        #     avg=list_len / key_count,
#        #     size=((24 + 24 + 16) * key_count + (4 * list_len)) / 1000000,
#        #     bucketperc=100 * key_count / max(new_dict.keys()),
#        # ))
#
#        return (new_dict, location_list, offsets, counts)
