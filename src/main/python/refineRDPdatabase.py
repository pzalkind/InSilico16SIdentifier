#!/usr/bin/python3.4


import os
import sys


def refineToNamed(dbOut, *dbFastaIn):
    dbBuffer = []
    for dbFasta in dbFastaIn:
        with open(dbFasta, 'r') as db:
            dbBuffer = db.read().strip('>').rstrip('\n').split('\n>')
    buffer2 = []
    for elem in dbBuffer:
        # take the first line (header) and remove the part after the tab, which corresponds to the lineage in RDP files
        header = elem.split('\n')[0].split('\t')[0]
        sequence = elem.split('\n')[1:]
        # first element of header corresponds to ID
        name = header.replace(';', '').split()[1:]
        if (len(name) >= 2) and (name[0][0].isupper()) and (name[1] != 'sp.') and (name[1] != 'bacterium'):
            # correct taxonomy suppose at least a Genus and a species, and Genus starts with an uppercase, and species must be known (not sp. nor bacterium)
            if ('Uncultivated' not in name) and ('Uncultured' not in name) and ('Unidentified' not in name) and ('Candidatus' not in name):
                # must corresponds to a known and cultivated bacterium, and validated (no Candidatus)
                new_header = '>{0}|{1}'.format(header.split()[0], '_'.join(name))
                buffer2.append('{0}\n{1}'.format(new_header, '\n'.join(sequence)))
    with open(dbOut, 'w') as newDb:
        newDb.write('\n'.join(buffer2))


if __name__ == '__main__':
    # go to the input directory, refine the RDP database and save new database with only named organisms, in the same directory
    os.chdir(os.path.dirname(os.path.abspath(sys.argv[1])))
    refineToNamed('16S_named_RDP.fasta', sys.argv[1:])

# EOF #
