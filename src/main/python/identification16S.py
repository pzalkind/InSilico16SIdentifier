#!/usr/bin/python3.4


# Needed modules
import sys
import subprocess
import gc
from operator import itemgetter
from itertools import chain


# Functions

def grep(pattern, word_list):
    return [elem for elem in word_list if pattern in elem]


# function that return a list of non-overlapping intervals from a list of (possibly) overlapping intervals
# (taken from StackOverflow)
def join_ranges(data, offset=0):
    data = sorted(chain.from_iterable(((start, 1), (stop + offset, -1)) for start, stop in data))
    c = 0
    for value, label in data:
        if c == 0:
            x = value
        c += label
        if c == 0:
            yield x, value - offset


def filterBlastOutput(blastOut, minAlignLength):
    hitsOutput = []
    if len(blastOut) > 0:
        hit_marker = 0
        hit = ''
        coords = []
        for line in blastOut.split('\n'):
            tmp = line.split('\t')
            if hit_marker == 0:
                hit_marker = 1
                pIdent = [float(tmp[2])]
                hit = tmp[0]
                contig = tmp[1]
                coords = [[int(tmp[4]), int(tmp[5])]]
                alignLength = 0
            elif tmp[0] != hit:
                coords = list(join_ranges(coords))
                for i in range(0, len(coords)):
                    alignLength = alignLength + abs(int(coords[i][1]) - (int(coords[i][0])-1))
                pIdent = float(sum(pIdent))/max(len(pIdent), 1)
                score = round(pIdent*alignLength)
                if alignLength >= minAlignLength:
                    name = hit.split('|')[1].replace(';', '')
                    if (name, contig, alignLength, coords, pIdent, score) not in hitsOutput:
                        hitsOutput.append((name, contig, alignLength, coords, pIdent, score))
                hit = tmp[0]
                contig = tmp[1]
                pIdent = [float(tmp[2])]
                coords = [[int(tmp[4]), int(tmp[5])]]
                alignLength = 0
            elif tmp[1] != contig:
                coords = list(join_ranges(coords))
                for i in range(0, len(coords)):
                    alignLength = alignLength + abs(int(coords[i][1]) - (int(coords[i][0])-1))
                coords = [[int(tmp[4]), int(tmp[5])]]
                contig = tmp[1]
                pIdent.append(float(tmp[2]))
            else:
                coords.append([int(tmp[4]), int(tmp[5])])
                pIdent.append(float(tmp[2]))
        coords = list(join_ranges(coords))
        for i in range(0, len(coords)):
            alignLength = alignLength + abs(int(coords[i][1]) - (int(coords[i][0])-1))
        pIdent = float(sum(pIdent))/max(len(pIdent), 1)
        score = round(pIdent*alignLength)
        if alignLength >= minAlignLength:
            name = hit.split('|')[1].replace(';', '')
            if (name, contig, alignLength, coords, pIdent, score) not in hitsOutput:
                hitsOutput.append((name, contig, alignLength, coords, pIdent, score))
    hitsOutput = sorted(hitsOutput, key=itemgetter(5, 4, 2), reverse=True)
    return hitsOutput


def blastContigsAgainst16SBank(contigsFilePath, database, minAlignLength, identityPerc):
    if isinstance(contigsFilePath, str):
        gc.collect()
        blast_proc = subprocess.Popen(' '.join(['blastn', '-query', contigsFilePath, '-db', database,
                                                '-max_target_seqs', '5', '-perc_identity', str(identityPerc),
                                                '-outfmt', '"6 sacc qseqid pident evalue qstart qend"']),
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        try:
            blastOut, blastErr = blast_proc.communicate()
            if blast_proc.returncode != 0:
                print('CRITICAL/nBlastN of the given file failed with following error message:\n{}\n'.format(blastErr.decode('utf-8')))
                sys.exit(1)
            blastOut = '\n'.join(sorted(blastOut.decode('utf-8').split('\n')[:-1]))
            hits = filterBlastOutput(blastOut, minAlignLength)
            return hits
        except subprocess.CalledProcessError:
            print('CRITICAL:\nBlastN of the given file failed with following error message:\n{}\n'.format(blastErr.decode('utf-8')))
            sys.exit(1)


def identification_16S(fastaFile, database):
    print('BEGIN 16S identification of {}\n'.format(fastaFile))
    print('Trying species-level identification, with min alignment = 1500 bases, and minimum identity = 97%\n')
    hits = blastContigsAgainst16SBank(fastaFile, database, 1500, 97.0)
    if hits == '':
        print('WARNING:\nNo BLAST hits found with sufficient identity and coverage for species-level identification.\n')
        print('Trying less stringent identification, with min alignment = 1200 bases, and minimum identity = 94%.\nResults will be genus-level confident\n')
        hits = blastContigsAgainst16SBank(fastaFile, database, 1200, 94.0)
        if hits == '':
            print('No BLAST hits found with sufficient identity and coverage for genus-level identification.')
    if len(hits) == 1:
        print('1 hit found')
    else:
        print('{} hits found'.format(len(hits)))
    i = 1
    for hit in hits:
        print('''Hit {i}:\n\tOrganism: {org}
        \n\tContig with 16S sequence: {seq}
        \n\tCoordinates on contig: {coord}
        \n\tIdentity with Query: {id}
        \n\t(alignment length: {len}, score: {score})\n
        '''.format(i=i, org=hit[0], seq=hit[1], coord=hit[3], id=hit[4], len=hit[2], score=hit[5]))
        i += 1
    print('END 16S identification\n')
    return 0


# Commands to process if called as a main program

if __name__ == '__main__':
    fastaFile = sys.argv[1]
    database = sys.argv[2]
    identification_16S(fastaFile, database)

# End of file #
