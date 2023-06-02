"""Description of what the script does"""

# Imports

import argparse

import textwrap as _textwrap
import urllib
import os
import sys
from collections import defaultdict
import pandas as pd

from Bio import SeqIO, Seq, SeqFeature,SeqRecord
from Bio.SeqFeature import BeforePosition, AfterPosition

# Class definitions

# This reclasses the argparse.HelpFormatter object to have newlines in the help text for paragraphs
class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent,
                                                 subsequent_indent=indent
                                                 ) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text


# Global variables

# Function definitions

def loadnamevariants(report=False):
    if report:
        print("The following are the gene name variants currently known.\n")
    output = {}
    fullparse = {}
    alltypes = set()
    url = "https://raw.githubusercontent.com/tjcreedy/constants/master/gene_name_variants.txt"
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8').strip()
        description, variants = line.split(":")
        name, annotype, fullname = description.split(";")
        variants = variants.split(',')
        variants.extend([name, fullname.upper()])
        fullvariants = []
        for v in variants:
            for g in ['', ' ']:
                v = v.replace(g, '')
                for s in ['', ' GENE', ' '+annotype.upper()]:
                    var = v+s
                    fullvariants.append(var)
                    output[v+s] = name
        alltypes.add(annotype)
        fullparse[name] = {'type': annotype, 'variants': fullvariants}
        if report:
            fullvariants = [v.replace(' ', '\u00A0') if len(v) < 12 else v for v in fullvariants]
            fullvariants = _textwrap.fill(', '.join(fullvariants), width=80,
                                          initial_indent='\t', subsequent_indent='\t')
            print(f"Standard name = {name}, type = {annotype}, full name = {fullname}:\n"
                  f"{fullvariants}")
    if not report:
        return output, alltypes, fullparse

def parse_specs(nameconvert):
    # Parse the mitocorrect master specifications file

    specs = dict()

    url = "https://raw.githubusercontent.com/tjcreedy/mitocorrect/bc42c7926947c0f35436d84d2ba7e76fd52d7cf5/specifications/mitocorrect_specifications_coleoptera_2022-02-24.tsv"
    sh = urllib.request.urlopen(url)

    # Extract the header row and split to list
    header = sh.readline().strip().split('\t')

    if len(set(header)) < len(header):
        sys.exit(f"Error: duplicate headers in {url}")

    ln = 1
    for line in sh:
        # line = sh.readline()

        ln += 1

        # Extract the values of the row and split to list
        items = line.strip().split('\t')

        # Get correct gene name
        name = None
        if items[0].upper() in nameconvert:
            name = nameconvert[items[0].upper()]

        if name is None:
            sys.exit(f"Error: gene name {items[0]} in first column of line {ln} is not recognised")

        # Generate holding dict for this line's specs
        hold = dict()

        # Work through specifications
        for spec, value in zip(header[1:], items[1:]):
            # spec, value = [header[2], items[2]]
            # Do not record if value is empty
            if value == '':
                continue
            # If a list specification, generate list
            if '/' in value:
                value = value.split('/')
            # If a dict specification, generate dict
            if (type(value) is list and ',' in value[0] or
                    type(value) is str and ',' in value):

                if type(value) is list:
                    value = [v.split(',') for v in value]
                else:
                    value = [value.split(',')]

                value = {k: v for k, v in value}

            # Add to specs dict
            hold[spec] = value

        # Do input checking

        for spec in hold.keys():
            # spec = 'overlap'
            value = hold[spec]
            err = f"Error: value {value} for {spec} on line {ln}"

            if spec == 'length':
                # should be > 0 integer
                if not str_is_int(value):
                    sys.exit(f"{err} is not an integer")
                value = int(value)
                if value < 1:
                    sys.exit(f"{err} is not greater than 0")
                hold[spec] = value

        del hold['end']
        
        # Initialise subdict if not already
        if name not in specs:
            specs[name] = dict()

        # Extract length and percent variation
        for i in ['length', 'lengthvariation']:
            if i in hold:
                specs[name][i] = hold.pop(i)

    sh.close()

    return specs




def get_feat_name(feat):
    featname = "unknown"
    nametags = ['gene', 'product', 'label', 'standard_name']
    if any(t in feat.qualifiers.keys() for t in nametags):
        for t in nametags:
            if t in feat.qualifiers.keys():
                featname = feat.qualifiers[t][0].upper()
                break
    return featname


def getcliargs(arglist=None, knowngenes=None, knowntypes=None):
    # Initialise the parser object and give a description
    parser = argparse.ArgumentParser(description="""
This script finds and extracts sequences corresponding to CDS, rRNA and/or tRNA annotations 
    from a set of sequences from one or more genbank-format flat files. Specify which regions to 
    extract using the -r/--regiontypes argument (default is CDS)
    |n
    The region name is identified using the 'standard_name' qualifier, and if that is not 
    available, the 'gene', 'product', or 'label' qualifiers. If none of these qualifiers are 
    found, the region is labeled as 'unknown'. The extracted sequences are written to separate 
    fasta-format files for each region type.
    """, formatter_class=MultilineFormatter)

    # Add arguments
    parser.add_argument("-g", "--genbank", type=str, metavar='PATH', required=True, nargs='+',
                        help="the path(s) to one or more genbank files from which gene sequences "
                             "should be extracted")
    parser.add_argument("-o", "--output", type=str, metavar='PATH', required=True,
                        help="the path to write the output csv")
    parser.add_argument("-r", "--genetypes", type=str, metavar='TYPE', nargs='+',
                        help="one or more gene types to be output",
                        choices=knowntypes, default=knowntypes)
    parser.add_argument("-m", "--lenerrormult", type=float, metavar='N', default = 5,
                        help=" multiple of the length variation percentage below which to denote "
                        "annotation as truncated")
    parser.add_argument("-s", "--showgenes", action='store_true',
                        help="print the known gene name variants")

    # Parse arguments
    args = parser.parse_args(arglist) if arglist else parser.parse_args()
    # Do some checking of the inputs

    # If the arguments are all OK, output them
    return args

# Start the actual script
if __name__ == "__main__":

        # Get the gene name variants
        nameconvert, annotypes, namevariants = loadnamevariants()

        # Get the length specification
        specs = parse_specs(nameconvert)

        # Get the arguments
        args = getcliargs(None, nameconvert.keys(), annotypes)
        # args = getcliargs('-g /home/thomas/work/iBioGen_postdoc/MMGdatabase/newdata_2022-04-15_syrphid/testnewsequences.gb -o testextract -w -k -p testpresence.txt'.split(' '),  nameconvert.keys(), annotypes)  # Read from a string, good for testing

        # Print gene name variants
        if args.showgenes:
            loadnamevariants(report=True)
            exit()

        # Set up for unrecognised genes
        unrecgenes = defaultdict(list)

        # *Create a DataFrame to store the results
        df = pd.DataFrame(
            columns=['NAME', 'SEQLENGTH', 'PRESENT GENES', 'PRESENT&COMPLETE GENES', 'ND2', 'COX1', 'COX2', 'ATP8',
                     'ATP6', 'COX3', 'ND3', 'ND5', 'ND4', 'ND4L', 'ND6',
                     'CYTB', 'ND1'])

        # Work through genbank files
        outfh = {}
        
        for gbpath in args.genbank:
            # gbpath = args.genbank[0]
            for seqrecord in SeqIO.parse(gbpath, "genbank"):
                # gb = SeqIO.parse(gbpath, "genbank")
                # seqrecord = next(gb)
                
                # Extract names
                seqname = seqrecord.name
                outname = seqname
                if args.organism and seqrecord.organism:
                    outname = seqrecord.organism.replace(' ', '_')

                # *Extract sequence length 序列长度
                seqlength = len(seqrecord.seq)

                # *Create a row to store the information
                row = pd.Series()
                row["NAME"] = outname
                row['SEQLENGTH'] = seqlength

                # *Set up containers
                prenottrunc = defaultdict(list)
                preandtrunc = defaultdict(list)

                # Work through features

                for feat in seqrecord.features:
                    # feat = seqrecord.features[3]
                    # Skip if not requested type
                    if feat.type not in args.genetypes:
                        continue

                    # Convert feat name
                    name = get_feat_name(feat)
                    if name in nameconvert:
                        stdname = nameconvert[name]
                    else:
                        unrecgenes[name].append(seqname)
                        if args.writeunknowns:
                            stdname = name
                        else:
                            continue
                    
                    # Get sequence
                    featsequence = feat.extract(seqrecord.seq)

                    # Check if length is too short
                    maxerr = specs[stdname]['length'] * (specs[stdname]['length']/100) * args.lenerrormult
                    short = len(featsequence) < (specs[stdname]['length'] - maxerr)
                    
                    # *store genes present and not truncated
                    if isinstance(feat.location.start, BeforePosition) or isinstance(feat.location.end, AfterPosition) or short:
                        preandtrunc[stdname].append(featsequence)
                    else:
                        prenottrunc[stdname].append(featsequence)

                # *write the presence to the row
                row['PRESENT GENES'] = len(prenottrunc) + len(preandtrunc)
                row['PRESENT&COMPLETE GENES'] = len(prenottrunc)
                for i in prenottrunc.keys():
                    if i in df.columns:
                        row[i] = "P"
                for i in preandtrunc.keys():
                    if i in df.columns:
                        row[i] = 'T'
                    # *Add the presence of genes in the current sequence to the dataframe将当前序列的基因存在情况添加到数据表中
                df = df.append(row, ignore_index=True)

        # *Output the dataframe as a CSV file
        if args.presence:
            df.to_csv(args.output, index=False)

        # Report unrecognised genes
        if len(unrecgenes) > 0:
            sys.stderr.write(f"Warning: unrecognised genes present in the following entries.\n")
            for gene, entries in unrecgenes.items():
                sys.stderr.write(f"\t{gene} - {', '.join(sorted(list(set(entries))))}\n")
