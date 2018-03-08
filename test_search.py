#!/usr/bin/env python3

from twobitreader import TwoBitFile
import search

search.reference = TwoBitFile("sample.2bit")

index = "sample.2bit.M10.Q10.index"


def test_convert():
    tmpfile = "scratch"
    search.write_tables_from_2bit("sample.2bit", outname=tmpfile)

    query = "TGGATTTCTGGTGCTTTCTACTCTGCCATATTCTCTGAATCCTCCTCTCTGGTTAAATATTTTTAAGGAATGCAGAGTCTGCAAACACCAACATTGCTGTGAACTGAGGTGTTGCTTTTTTTTTTTTTTTTAAGGAAAAAGGAAAAAAAAGATTAACATAACCACCTGTTTCTACTTGGAAGAAACTCTAGAGCGCAAATGCATTTAAAT"

    for i in range(0, 12):
        matches = search.match_file(tmpfile, query[i:])
        assert len(matches) == 1
        assert matches[0] == ("chr3", 71 + i)


#def test_chrUn_lookup():
#    query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#    matches = search.match_file(index, query)
#    assert len(matches) == 1
#    assert matches[0][1] == 139301
#    assert matches[0][0] == "chrUn_GL000218v1"


# def test_lookup_bench_M10_Q10(benchmark):
#     index = "GRCh38_no_alts.2bit.M10.Q10.index.hdf5"
#     query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#     benchmark(search.match_file, index, query)
# 
# def test_lookup_bench_M10_Q10_nozip(benchmark):
#     index = "GRCh38_no_alts.2bit.M10.Q10.index.hdf5.nocompress"
#     query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#     benchmark(search.match_file, index, query)
# 
# def test_lookup_bench_M10_Q10_lzf(benchmark):
#     index = "GRCh38_no_alts.2bit.M10.Q10.index.hdf5.lzf"
#     query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#     benchmark(search.match_file, index, query)
# 
# def test_lookup_bench_M09_Q11(benchmark):
#     index = "GRCh38_no_alts.2bit.M9.Q11.index.hdf5"
#     query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#     benchmark(search.match_file, index, query)
# 
# def test_lookup_bench_M09_Q11_nozip(benchmark):
#     index = "GRCh38_no_alts.2bit.M9.Q11.index.hdf5.nocompress"
#     query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#     benchmark(search.match_file, index, query)
# 
# def test_lookup_bench_M09_Q11_lzf(benchmark):
#     index = "GRCh38_no_alts.2bit.M9.Q11.index.hdf5.lzf"
#     query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#     benchmark(search.match_file, index, query)
# 
# 
# def test_lookup_bench_M08_Q12(benchmark):
#     index = "GRCh38_no_alts.2bit.M8.Q12.index.hdf5"
#     query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#     benchmark(search.match_file, index, query)
# 
# def test_lookup_bench_M07_Q14(benchmark):
#     index = "GRCh38_no_alts.2bit.M7.Q14.index.hdf5"
#     query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#     benchmark(search.match_file, index, query)
# 
# def test_lookup_bench_M09_Q11_C10K_lzf(benchmark):
#     index = "GRCh38_no_alts.2bit.M10.Q10.C10K.index.hdf5.lzf"
#     query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#     benchmark(search.match_file, index, query)
# 
# def test_lookup_bench_M10_Q10_C10K_lzf(benchmark):
#     index = "GRCh38_no_alts.2bit.M10.Q10.C10K.index.hdf5.lzf"
#     query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#     benchmark(search.match_file, index, query)



# def test_lookup_bench_M06_Q16_C10K_lzf(benchmark):
#     index = "GRCh38_no_alts.2bit.M6.Q16.C10K.index.hdf5.lzf"
#     query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#     benchmark(search.match_file, index, query)
