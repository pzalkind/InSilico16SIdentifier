#!/usr/bin/python3.4


# Needed modules
import sys
import subprocess
import logging
import gc
import datetime
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
                logging.critical('BlastN of the given file failed with following error message:\n{}'.format(blastErr.decode('utf-8')))
                sys.exit(1)
            blastOut = '\n'.join(sorted(blastOut.decode('utf-8').split('\n')[:-1]))
            hits = filterBlastOutput(blastOut, minAlignLength)
            return hits
        except subprocess.CalledProcessError:
            logging.critical('BlastN of the given file failed with following error message:\n{}'.format(blastErr.decode('utf-8')))
            sys.exit(1)


def findOrganisms(fastaFile, database):
    logging.info('BEGIN 16S identification of {}'.format(fastaFile))
    logging.info('Trying species-level identification, with min alignment = 1500 bases, and minimum identity = 97%')
    hits = blastContigsAgainst16SBank(fastaFile, database, 1500, 97.0)
    if hits == '':
        logging.warning('No BLAST hits found with sufficient identity and coverage for species-level identification.')
        logging.info('Trying less stringent identification, with min alignment = 1200 bases, and minimum identity = 94%.\nRsults will be genus-level confident')
        hits = blastContigsAgainst16SBank(fastaFile, database, 1200, 94.0)
        if hits == '':
            logging.warning('No BLAST hits found with sufficient identity and coverage for genus-level identification.')
    logging.info('{} hits found'.format(len(hits)))
    i = 1
    for hit in hits:
        results = '''Hit n°{i}:\n\tOrganism: {org}
                     \n\tContig with 16S sequence: {seq}
                     \n\tCoordinates on contig: {coord}
                     \n\tIdentity with Query: {id}
                     \n\t(alignment length: {len}, score: {score})\n
                     '''.format(i=i, org=hit[0], seq=hit[1], coord=hit[3], id=hit[4], len=hit[2], score=hit[5])
        print(results)
        logging.info(results)
        i += 1
    logging.info('END 16S identification of assembled genome\n')
    return 0


# Main function

def identification_16S(fastaFile, database, logFile):
    logging.basicConfig(filename=logFile, format='%(asctime)s - %(levelname)s : %(message)s', level=logging.INFO)
    logging.FileHandler(logFile, mode='w')
    logging.info('BEGIN Pipeline (16S identification, reference search)\n')
    # 16S identification, blastn of the assembly against the specified 16S database
    findOrganisms(fastaFile, database)
    logging.info('END Pipeline')


# Commands to process if called as a main program

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Not enough argument, expected 2, found: {}'.format(sys.argv[1:]))
    else:
        fastaFile = sys.argv[1]
        database = sys.argv[2]
        logFile = 'log_{}.txt'.format(datetime.datetime.today().strftime('%Y_%m_%d_%H_%M_%S_%f'))
        identification_16S(fastaFile, database, logFile)

# End of file #
