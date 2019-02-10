from Bio import Entrez
from Bio import SeqIO
import io
import pandas as pd
import re

__author__ = 'Natalia Quinones-Olvera'
#__email__ = "nquinones@g.harvard.edu"
#Entrez.email = __email__

# .............................FUNCTIONS................................


def parse_loc(loc, ref_start):
    'Function to correctly parse location object'
    # start
    if str(type(loc.start)) == "<class 'Bio.SeqFeature.ExactPosition'>":
        start = loc.start + ref_start
        start_str = start
    elif str(type(loc.start)) == "<class 'Bio.SeqFeature.BeforePosition'>":
        start = loc.start + ref_start
        start_str = '<' + str(start)
    elif str(type(loc.start)) == "<class 'Bio.SeqFeature.AfterPosition'>":
        start = loc.start + ref_start
        start_str = '>' + str(start)
    else:
        start_str = 'unknown position type'

    # end
    if str(type(loc.end)) == "<class 'Bio.SeqFeature.ExactPosition'>":
        end = loc.end + ref_start - 1
        end_str = end
    elif str(type(loc.end)) == "<class 'Bio.SeqFeature.BeforePosition'>":
        end = loc.end + ref_start - 1 
        end_str = '<' + str(end)
    elif str(type(loc.end)) == "<class 'Bio.SeqFeature.AfterPosition'>":
        end = loc.end + ref_start -1
        end_str = '>' + str(end)
    else:
        end_str = 'unknown position type'

    if loc.strand == 1:
        strand = '+'
    else:
        strand = '-'

    return [start, end, start_str, end_str, strand]


def get_prot_coords(protein_acc):
    '''
    Helper function to find protein coordinates from accession.
    Has limited functionality
    '''
    # check if it is non-redundant
    if protein_acc.find('WP') == 0:
        print('Warning: Seems that you\'re providing a non-redundant protein accession.')

    # send request to NCBI protein database
    epost = Entrez.epost('protein', id=protein_acc)
    request = Entrez.read(epost)

    response = Entrez.efetch(db='protein',
                             webenv=request['WebEnv'],
                             query_key=request['QueryKey'],
                             rettype="gb",
                             retmode="text")

    response_io = io.StringIO(response.read())

    feat_list = []

    for gb_record in SeqIO.parse(response_io, "genbank"):
        for feat in gb_record.features:
            feat_list.append(feat.type)
                # source feature
            if feat.type == 'CDS':
                coded_str = feat.qualifiers['coded_by'][0].split(':')
                acc = coded_str[0]
                coords = coded_str[1].split('..')
                start = coords[0]
                end = coords[1]

    # if no CDS feature found, return error message
    if 'CDS' not in feat_list:
        print('Error: Couldn\'t find CDS feature. Unable to find coordinates.')
        return None

    return acc, (start, end)


def gene_viking(acc, start, end, thresh, output):
    '''
    Main functionality.
    '''
    start_2 = start - thresh
    end_2 = end + thresh

    if start_2 < 0:
        start_2 = 0

    epost = Entrez.epost('nuccore', id=acc)
    request = Entrez.read(epost)
    response = Entrez.efetch(db='nuccore',
                             webenv=request['WebEnv'],
                             query_key=request['QueryKey'],
                             rettype="gb",
                             retmode="text",
                             seq_start=start_2,
                             seq_stop=end_2)

    response_io = io.StringIO(response.read())

    query_range = set(range(start, end))

    all_l = []

    for record in SeqIO.parse(response_io, 'genbank'):
        for feat in record.features:
            if feat.type == 'CDS':
                if 'pseudo' not in feat.qualifiers:
                    fstart, fend, fstart_str, fend_str, fstrand = parse_loc(feat.location, start_2)
                    prot_id = feat.qualifiers['protein_id'][0]
                    product = feat.qualifiers['product'][0]

                    feat_range = set(range(fstart, fend))

                    if len(query_range.intersection(feat_range)) == 0:
                        status = ''
                    else:
                        status = '!'

                    feats = [prot_id, product, fstart_str, fend_str, fstrand, status]

                    all_l.append(feats)

    df = pd.DataFrame(all_l)

    df.columns = ['prot_acc', 'product', 'start', 'end', 'strand', 'query_overlap']

    if output is not None:
        df.to_csv(output, sep='\t', index=None)

    return df
