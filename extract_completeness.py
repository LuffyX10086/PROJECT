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
                        help="the path to a directory in which to place a fasta file for each "
                             "gene, this will be created if it doesn't exist")
    parser.add_argument("-m", "--mingenes", type=int, metavar='N', default=0,
                        help="the minimum number of genes that must be present in an entry for "
                             "any genes to be output")
    parser.add_argument("-q", "--reqgenes", type=str, metavar='GENE', nargs='*',
                        help="a list of genes that must be present in an entry for any gene "
                             "sequences to be output",
                        choices=knowngenes)
    parser.add_argument("-r", "--genetypes", type=str, metavar='TYPE', nargs='+',
                        help="one or more gene types to be output",
                        choices=knowntypes, default=knowntypes)
    parser.add_argument("-f", "--filter", type=str, metavar='PATH',
                        help="path to a file listing LOCUS names, one per line; if supplied, only "
                             "these entries will be output")
    parser.add_argument("-n", "--organism", action='store_true',
                        help="use the organism field instead of the locus field for sequence "
                             "headers in the output files")
    parser.add_argument("-w", "--writeunknowns", action='store_true',
                        help="write fastas for unidentifiable annotations anyway")
    parser.add_argument("-k", "--keepframe", action='store_true',
                        help="remove excess out-of-frame bases at the beginning of truncated "
                             "annotations")
    parser.add_argument("-p", "--presence", type=str, metavar='PATH',
                        help="if desired, the path to a text file to which will be written a list "
                             "of genes present for each input entry")
    parser.add_argument("-s", "--showgenes", action='store_true',
                        help="print the known gene name variants")

    # Parse arguments
    args = parser.parse_args(arglist) if arglist else parser.parse_args()
    # Do some checking of the inputs
    if args.mingenes < 0:
        parser.error(f"{args.mingenes} is not a possible minimum number of genes")

    # If the arguments are all OK, output them
    return args

# Start the actual script
if __name__ == "__main__":

        # Get the gene name variants
        nameconvert, annotypes, namevariants = loadnamevariants()

        # Get the arguments
        args = getcliargs(None, nameconvert.keys(), annotypes)
        # args = getcliargs('-g /home/thomas/work/iBioGen_postdoc/MMGdatabase/newdata_2022-04-15_syrphid/testnewsequences.gb -o testextract -w -k -p testpresence.txt'.split(' '),  nameconvert.keys(), annotypes)  # Read from a string, good for testing

        # Print gene name variants
        if args.showgenes:
            loadnamevariants(report=True)
            exit()

        # Make output directory if necessary
        if not os.path.exists(args.output):
            os.makedirs(args.output)

        # Set up for unrecognised genes
        unrecgenes = defaultdict(list)

        # Set up gene presence handle
        if args.presence:
            args.presence = open(args.presence, 'w')

        # Read filter list
        filter = []
        if args.filter:
            with open(args.filter) as fh:
                for line in fh:
                    filter.append(line.rstrip())

        # Work through genbank files
        nrejected = 0
        outfh = {}
        # *Create a DataFrame to store the results
        df = pd.DataFrame(
            columns=['NAME', 'SEQLENGTH', 'PRESENT GENES', 'PRESENT&COMPLETE GENES', 'ND2', 'COX1', 'COX2', 'ATP8',
                     'ATP6', 'COX3', 'ND3', 'ND5', 'ND4', 'ND4L', 'ND6',
                     'CYTB', 'ND1'])
        for gbpath in args.genbank:
            # gbpath = args.genbank[0]
            for seqrecord in SeqIO.parse(gbpath, "genbank"):
                # gb = SeqIO.parse(gbpath, "genbank")
                # seqrecord = next(gb)
                # If filtering, check if this entry should be used
                if args.filter and seqrecord.name not in filter:
                    nrejected += 1
                    continue

                # Extract names
                seqname = seqrecord.name
                outname = seqname
                if args.organism and seqrecord.organism:
                    outname = seqrecord.organism.replace(' ', '_')

                # *Extract sequence length 序列长度
                seqlength = len(seqrecord.seq)

                # *Create a row to store the information
                row = pd.Series()
                row["NAME"]=outname
                row['SEQLENGTH']=seqlength





                # *Set up containers
                foundgenes = defaultdict(list)
                prenottrunc=defaultdict(list)
                preandtrunc=defaultdict(list)

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



                    # Store sequence for writing
                    featsequence = feat.extract(seqrecord.seq)
                    if args.keepframe and 'codon_start' in feat.qualifiers:
                        featsequence = featsequence[(int(feat.qualifiers['codon_start'][0]) - 1):]
                    foundgenes[stdname].append(featsequence)

                    # *store genes present and not truncated
                    if isinstance(feat.location.start, BeforePosition) or isinstance(feat.location.end, AfterPosition):
                        preandtrunc[stdname].append(featsequence)
                    else:
                        prenottrunc[stdname].append(featsequence)

                # *write the presence to the row
                row['PRESENT GENES'] = len(foundgenes)
                row['PRESENT&COMPLETE GENES'] = len(prenottrunc)
                for i in prenottrunc:
                    if i in df.columns:
                        row[i] = "P"
                for i in preandtrunc:
                    if i in df.columns:
                        row[i] = 'T'
                    # *Add the presence of genes in the current sequence to the dataframe将当前序列的基因存在情况添加到数据表中
                df = df.append(row, ignore_index=True)

                # Check to see if there are enough found features and all required genes are present,
                # if so write sequences out
                if len(foundgenes) >= args.mingenes:
                    if not args.reqgenes or all([g in foundgenes for g in args.reqgenes]):
                        for gene, seqs in foundgenes.items():
                            # Check the number of sequences and warn if there are > 1
                            if len(seqs) > 1:
                                sys.stderr.write(f"Warning: entry {seqname} has multiple annotations "
                                                 f"of the target type ({'/'.join(args.genetypes)}) "
                                                 f"matching the standard name \"{gene}\"\n")
                            # Check to see if there's already a file handle for this gene, if not open
                            if gene not in outfh:
                                outfh[gene] = open(os.path.join(args.output, f"{gene}.fasta"), 'w')
                            # Write out sequences
                            for seq in seqs:
                                outfh[gene].write(f">{outname}\n{seq}\n")

                        #if args.presence:
                            #args.presence.write(f"{outname},{seqlength},{len(foundgenes)},{len(prenottrunc)},{','.join(foundgenes.keys())},{','.join(prenottrunc.keys())}\n")

                    else:
                        sys.stderr.write(f"Warning: entry {seqname} does not include all required "
                                         f"genes, it will be skipped")
                else:
                    sys.stderr.write(f"Warning: entry {seqname} has {len(foundgenes)} known "
                                     f"annotations with sufficient information, it will be skipped")

        # *Output the dataframe as a CSV file
        if args.presence:
            df.to_csv('completeness.csv', index=False)

        # Close file handles
        for gene, fh in outfh.items():
            fh.close()
        if args.presence:
            args.presence.close()

        # Report unrecognised genes
        if len(unrecgenes) > 0:
            sys.stderr.write(f"Warning: unrecognised genes present in the following entries.\n")
            for gene, entries in unrecgenes.items():
                sys.stderr.write(f"\t{gene} - {', '.join(sorted(list(set(entries))))}\n")

        # Report rejected entries
        if nrejected > 0:
            sys.stderr.write(f"Warning: {nrejected} entries in genbank-format file(s) rejected for "
                             f"not matching filter list")
