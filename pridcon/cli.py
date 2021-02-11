import click
import os
import pandas as pd

from assembly_predictor.digraph import deBruijnGraph
from assembly_predictor.utils import read_fasta, reads_generator, write_fastq, read_fastq, write_fasta
from assembly_predictor.startup import CACHE_PATH
from assembly_predictor.contig import Contig
from assembly_predictor.api_processor import ApiProcessor


@click.group()
def main():
    """Entry method"""
    pass


@main.command()
@click.argument('input_fasta', type=click.Path(exists=True))
@click.argument('output_fastq')
@click.option('-k', '--kmer', default=30, help='length of k-mer')
@click.option('-r', '--read', default=100, help='length of reads')
def createQ(input_fasta: str, output_fastq: str, kmer: int, read: int):
    """
    Simulates perfect sequencing. Reads are written into fastq file format.
    """
    reads_list = reads_generator(input_fasta, read, kmer)
    write_fastq(reads_list, output_fastq)


@main.command()
@click.argument('input_fasta', type=click.Path(exists=True))
@click.argument('output_path', type=click.Path())
@click.option('-k', '--kmer', default=30, help='length of k-mer')
@click.option('-c', '--contig', is_flag=True, default=True, help='Write fasta file of assembled sequences')
@click.option('-r', '--rna', is_flag=True, default=True, help='Write fasta file of mrna sequences from assembly')
@click.option('-a', '--aa', is_flag=True, default=True, help='Write fasta file of amino acid sequences from assembly')
@click.option('-c', '--curated', default=True, is_flag=True, help='Choose between manually/automatic curated results')
@click.option('-h', '--hits', default=5, help='Number of hits to be displayed per AA sequence')
@click.option('-r', '--ranked', default=True, is_flag=True, help='Choose if the results are ranked as per E-value')
@click.option('-b', '--blast', default=True, is_flag=True, help='Write csv file of predicted protein information')
@click.option('-u', '--uniprot', default=False, is_flag=True, help='Write csv file of Uniprot information')
@click.option('-e', '--ensembl', default=False, is_flag=True, help='Write csv file of Ensembl information')
def predictQ(input_fasta: str, output_path: str, kmer: int, contig: bool, rna: bool, aa: bool,
             curated: bool, hits: int, ranked: bool, blast: bool, uniprot: bool, ensembl: bool):
    """
    Sequence assembly and protein prediction from sequence reads. Fastq and multifasta accepted.
    Output files are written into output directory (which must be already created).

    The BLAST query is tuned to search in both UNIPROT KB and TrEMBL databases.
    The information is parsed and converted to a dataframe.
    """
    reads_dict = read_fastq(input_fasta)
    debruijn = deBruijnGraph(reads_dict, k=kmer)
    contigID, sequence = next(iter(debruijn.get_contigs().items()))
    prediction = Contig(contigID, sequence)

    if contig:
        write_fasta(os.path.join(output_path, 'contig_{}{}'.format(contigID, '.fasta')), prediction.contig_sequence())
    if rna:
        write_fasta(os.path.join(output_path, 'mrna_{}{}'.format(contigID, '.fasta')), prediction.mrna_sequence())
    if aa:
        write_fasta(os.path.join(output_path, 'orfs_{}{}'.format(contigID, '.fasta')), prediction.found_orfs())

    blast_request = ApiProcessor(prediction.found_orfs(), curated=curated, hits=hits)
    blast_request.api_request()

    if ranked:
        blast_request.ranking()

    blast_request.api_request_to_uniprot_id()
    blast_request.uniprot_id_to_ensembl_id()
    blast_request.wrapper()

    if blast:
        query_df = blast_request.QUERY_DF
        query_df.to_csv(os.path.join(output_path, "blast_dataframe.csv"), encoding='utf-8')
    if uniprot:
        uniprot_df = blast_request.UNIPROT_DF
        uniprot_df.to_csv(os.path.join(output_path, "uniprot_dataframe.csv"), encoding='utf-8')
    if ensembl:
        ensembl_df = blast_request.ENSEMBL_DF
        ensembl_df.to_csv(os.path.join(output_path, "ensembl_dataframe.csv"), encoding='utf-8')


if __name__ == "__main__":
    main()
