import urllib.parse
import urllib.request
import requests
import json
import xmltodict
import pprint
import collections
import pandas as pd
import copy
from IPython.display import HTML

class ApiAccess:
    """
    ApiProcessor class is responsible for provision of available information from
    https://web.expasy.org - Swiss Institute of Bioinformatics Vital-IT Center  for high-performance computing

    Protein BLAST is done via API get requests.
    The query is tuned to look across both UNIPROT KB and TrEMBL databases.
    The information is parsed and converted to dictionary data structure.
    """

    def __init__(self, query, curated=True, hits=5):
        
        self._request_dict = self._api_request(query, hits, curated)

    def get_request(self):
        """Returns BLAST request"""
        return self._request_dict


    def _api_request(self, query: dict, hits: int = 5, curated: bool = True) -> dict:
        """
        API access to SIB with respective protein query.

        Arguments
        ---------
            query: dictionary containing amino acid sequences
        Returns
        -------
            query_dict: dictionary containing retrieved API information

        """
        query_dict = {}
        curate = "curated=on&" if curated else ''
        for seq_number, aa_sequence in query.items():
            request_query = f'https://web.expasy.org/cgi-bin/blast/blast.pl?seq={aa_sequence}&prot_db1=UniProtKB&{curate}ethr=10&showsc={hits}&showal={hits}&format=xml'
            request = self.decode_request(request_query)
            if request: 
                for hit in range(len(request)):
                    query_index = request[hit]
                    uniprot_id = query_index['Hit_def'].split("|")[1]
                    protein_name = query_index['Hit_def'].split("|")[2]
                    hit_length = query_index['Hit_len']
                    bit_score = query_index['Hit_hsps']['Hsp']['Hsp_bit-score']
                    identity = query_index['Hit_hsps']['Hsp']['Hsp_identity']
                    e_value = query_index['Hit_hsps']['Hsp']["Hsp_evalue"]
                    hit_sequence = query_index['Hit_hsps']['Hsp']['Hsp_hseq']

                    query_dict['{}.{}'.format(seq_number, hit)] = { 'hit number ' : hit,
                                                                    'uniprot id': uniprot_id,
                                                                    'protein name': protein_name,
                                                                    'bitscore': bit_score,
                                                                    'E-value': e_value,
                                                                    'Identity' : identity
                                                                    }
        return query_dict
    
    def decode_request(self, request_query: str) -> dict:
        """Decodes API request"""
        
        response = requests.get(request_query)
        if response.status_code == requests.codes.ok:
            decoded = response.content 
            xml_to_dict = xmltodict.parse(decoded)
            json_format = json.dumps(xml_to_dict)
            try: 
                request = json.loads(json_format)['BlastOutput']['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit']
                return request
            except: raise Exception("No Hits found")  
                
    def request_to_pandas(self)-> pd.DataFrame:
        """request dictionary to pandas DataFrame"""
        request_df = pd.DataFrame.from_dict(self._request_to_uniprot_id(), orient='index')
        request_df.reset_index(level=0, inplace=True)
        request_df.rename(columns = {"index": "Sequence ID"},inplace = True) 
        request_df['Sequence ID'] = [seqID.split('.')[0] for seqID in request_df['Sequence ID']]
        
        return HTML(request_df.to_html(escape=False))
    
    def _request_to_uniprot_id(self) -> dict:
        """Converts UniProt ID to href to database"""
        request_copy = copy.deepcopy(self._request_dict)
        for sequence_ID, value in request_copy.items():
            uniprot_id = value['uniprot id']
            request_copy[sequence_ID]['uniprot id'] = '<a href="https://www.uniprot.org/uniprot/{}">{}</a>'.format(str(uniprot_id), uniprot_id) 
        return request_copy
        