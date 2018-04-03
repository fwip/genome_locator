#!/usr/bin/env python3

from twobitreader import TwoBitFile
from search import GenomeSearcher
import random
import subprocess
import tempfile
from hypothesis import example, given, strategies, settings


index = "sample.2bit.M10.Q10.index.hdf5"


def test_convert():
    tmpfile = "scratch"
    search = GenomeSearcher(reference_file="sample.2bit")

    search.write_tables_from_2bit("sample.2bit", outname=tmpfile)

    query = "TGGATTTCTGGTGCTTTCTACTCTGCCATATTCTCTGAATCCTCCTCTCTGGTTAAATATTTTTAAGGAATGCAGAGTCTGCAAACACCAACATTGCTGTGAACTGAGGTGTTGCTTTTTTTTTTTTTTTTAAGGAAAAAGGAAAAAAAAGATTAACATAACCACCTGTTTCTACTTGGAAGAAACTCTAGAGCGCAAATGCATTTAAAT"

    for i in range(0, 12):
        matches = search.match_file(tmpfile, query[i:])
        assert len(matches) == 1
        assert matches[0] == ("chr3", 71 + i)

def check_region(ref, chrom, pos, length):
    tmpfile = "scratch"
    search = GenomeSearcher(reference_file=ref)

    search.write_tables_from_2bit(ref, outname=tmpfile)

    tbf = TwoBitFile("sample.2bit")
    query = tbf[chrom][pos:pos+length]
    matches = search.match_file(tmpfile, query)
    assert (chrom, pos+1) in matches

def test_problem():
    check_region("sample.2bit", "chr1", 303467, 116)


def test_alternate_params():
    tmpfile = "scratch"
    search = GenomeSearcher(reference_file="sample.2bit", M=3, Q=12)

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


#def test_lookup_bench_M10_Q10(benchmark):
#    index = "GRCh38_no_alts.2bit.M10.Q10.index.hdf5"
#    query = "TGTATGTTTTTCTATCTCCACACACTCCTGAACATAGAAAGACCAAGTAACATCCCTGGTTAAGATGTGTACAGGTTACAAGACATGTCTAAATATATTCACCAAGAGGTTTATT"
#    benchmark(search.match_file, index, query)


def test_random_lookup():
    index = "scratch_random"
    search = GenomeSearcher(reference_file="sample.2bit")
    search.write_tables_from_2bit("sample.2bit", outname=index)
    random.seed(9001)
    tbf = TwoBitFile("sample.2bit")
    min_length = 100
    max_length = 120
    sizes = tbf.sequence_sizes()
    for i in range(0, 200):
        chrom = random.choice(list(sizes.keys()))
        length = random.randint(min_length, max_length)
        pos = random.randrange(sizes[chrom]-length-1)
        query = tbf[chrom][pos:pos+length]
        if "N" in query:
            continue
        matches = search.match_file(index, query)
        print(chrom, pos+1, length, matches)
        assert (chrom, pos+1) in matches


def write_dna_to_2bit(dna, bitfile=None):
    if bitfile is None:
        f, bitfile = tempfile.mkstemp()
    fastafile, fastaname = tempfile.mkstemp()
    with open(fastafile, 'w') as fasta:
        fasta.write('>chr1\n')
        fasta.write(dna)
        fasta.write('\n')
        print(fastaname)

    errcode = subprocess.call(["faToTwoBit", fastaname, bitfile])

    return bitfile


@given(dna=strategies.text(alphabet='ACTGN', min_size=100),
       pos=strategies.integers(min_value=0),
       length=strategies.integers(min_value=0),
       M=strategies.integers(min_value=2),
       Q=strategies.integers(min_value=3),
       )
@settings(max_examples=200,
          deadline=None)
@example(dna="CCCCTTTTGGAAAGG",
         pos=2,
         length=5,
         M=2,
         Q=2,
         )
def test_simple_alphabet(dna, pos, length, M, Q):
    if len(dna) < M * Q:
        return
    if len(dna) < pos + length:
        return
    if length < M * Q:
        return

    query = dna[pos:pos+length]
    chrom = "chr1"
    index = "test.index"
    bitfile = write_dna_to_2bit(dna)

    search = GenomeSearcher(reference_file=bitfile, M=M, Q=Q)
    search.write_tables_from_2bit(bitfile, index)

    matches = search.match_file(index, query)

    assert (chrom, pos+1) in matches


@given(chrom=strategies.sampled_from(['chr1', 'chr2', 'chr3', 'chr4']),
       pos=strategies.integers(min_value=0),
       length=strategies.integers(min_value=100, max_value=120))
@example('chr1', 359416, 100)
def test_sample_hypothesis(chrom, pos, length):
    index = "scratch_random"
    search = GenomeSearcher(reference_file="sample.2bit")
    tbf = TwoBitFile("sample.2bit")
    sizes = tbf.sequence_sizes()
    if pos + length >= sizes[chrom]:
        return

    query = tbf[chrom][pos:pos+length]
    if "N" in query:
        return

    print(pos, chrom, query)
    matches = search.match_file(index, query)
    assert (chrom, pos+1) in matches

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
