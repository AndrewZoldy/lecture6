import pymysql
import sequence_worker as sw
import argparse

def create_table(conn_features, table_name, table_fields):
    conn = pymysql.connect(*conn_features)
    with conn.cursor() as cur:
        cur.execute('create table ' + table_name + ' (' + ', '.join(table_fields) + ') ENGINE=InnoDB')
    conn.commit()
    return 0


def add_to_table(conn_features, table_name, columns, data_to_add):
    conn = pymysql.connect(*conn_features)
    with conn.cursor() as cur:
        try:
            cur.execute('insert into ' + table_name + ' ('+', '.join(columns) + ') VALUES ('
                    + ', '.join(data_to_add) +')')
        except pymysql.err.IntegrityError:
            return 1
    conn.commit()
    return 0

def acces_table(conn_features, table_name, columns, condition=None):
    conn = pymysql.connect(*conn_features)
    with conn.cursor() as cur:
        if not condition:
            cur.execute('select '+', '.join(columns)+' from %s' % table_name)
        else:
            cur.execute('select ' + ', '.join(columns) + ' from %s' % table_name+
                        ' WHERE %s' % condition)
        return cur.fetchall()

def check_table(conn_features, name):
    conn = pymysql.connect(*conn_features)
    with conn.cursor() as cur:
        tables_number = cur.execute('show tables')
        tables = [table_name for (table_name,) in cur]
        if name in tables:
            return True
        else:
            return False


def make_index(conn_features, table_name, col_name, index_name, ind_arg=''):
    conn = pymysql.connect(*conn_features)
    with conn.cursor() as cur:
        cur.execute('CREATE '+ind_arg+'INDEX '+index_name+' ON ' + table_name+' ('+col_name+')')
    conn.commit()


def make_orf_table_names(name):
    dna_name = name + '_orfs_dna'
    rna_name = name + '_orfs_rna'
    aa_name = name + '_orfs_aa'
    return [dna_name, rna_name, aa_name]

def extraction(Seq, type):
    if type == 'dna':
        dna = Seq.as_dna()
        return [str(dna), dna.weight()]
    if type == 'rna':
        rna = Seq.as_rna()
        return [str(rna), rna.weight()]
    if type == 'protein':
        prt = Seq.as_protein()
        return [str(prt), prt.weight()]



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='tool for working around biological sequence data')
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('-db', '--database', type=str, default='hw_db', help='specify db to work with')
    parser.add_argument('-u', '--user', type=str, default='hw_user', help='specify user to work with db')
    parser.add_argument('-p', '--password', type=str, default='qwerty', help='password to login user')
    parser.add_argument('--host', type=str, default='localhost', help='host, where DB is present')
    parser.add_argument('-i', '--file', type=str, default='hw_db', help='specify file to work with')
    parser.add_argument('-s', '--start', type=int, default=0, help='set index from which you want to get sequence')
    parser.add_argument('-e', '--end', type=int, default=100, help='set index until which you want to get sequence')
    parser.add_argument('-f', '--field', type=str, default='Seq_name', help='specify table field to find sequence '
                                                                            '(currently only "Seq_name" for main sequence'
                                                                            ' and "ORF_id" for ORFs are allowed; ORF_id'
                                                                            'is parent`s {Seq_name}_{number of orf},'
                                                                            ' i.e. gene1_1).')
    parser.add_argument('-c', '--condition', type=str, default='gene1', help='from which sequence do you want to'
                                                                             ' retrieve subsequence')
    parser.add_argument('-x', '--extract', type=str, default='dna', help='type of molecule to extract')
    group.add_argument('-a', '--add_entries', action='store_true', help='add new entries in table '
                                                                        '(you have to specify input file)')
    group.add_argument('-g', '--get_sequence', action='store_true', help='get sequence and weight from table')
    args = parser.parse_args()
    database = args.database
    user = args.user
    password = args.password
    host = args.host
    file = args.file
    start = args.start
    end = args.end
    field = args.field
    condition = args.condition
    extract = args.extract
    add_entries = args.add_entries
    get_sequence = args.get_sequence

    connection_features = [host, user, password, database]
    with open('lec6_hw.fsa', 'r') as inp:
        data = inp.read()

    seq_table_name = file.split('.')[0] + '_sequences'
    orf_table_name = file.split('.')[0] + '_orfs'
    if add_entries:
        if not file:
            raise Exception('You have to specify input file for this mode!')
        index_name = 'fk_orf_to_sequence'
        colname = 'Seq_name'

        if not check_table(connection_features,seq_table_name):
            create_table(connection_features, seq_table_name,
                         ['Seq_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY','Seq_name VARCHAR(300)',
                                                                   'Seq TEXT', 'Seq_identity VARCHAR(200)',
                                                                   'Weight VARCHAR(100)'])
            make_index(connection_features, seq_table_name, 'Seq_identity', 'unique_' + seq_table_name, 'UNIQUE ')

        if not check_table(connection_features, orf_table_name):
            make_index(connection_features, seq_table_name, colname, index_name + orf_table_name)
            create_table(connection_features, orf_table_name, ['id INT NOT NULL AUTO_INCREMENT PRIMARY KEY',
                                                               'Seq_name VARCHAR(200)',
                                                               'ORF_id VARCHAR(250)',
                                                               'ORF_dna TEXT', 'ORF_rna TEXT', 'ORF_protein TEXT',
                                                               'Weight_dna VARCHAR(100)',
                                                               'Weight_rna VARCHAR(100)',
                                                               'Weight_protein VARCHAR(100)',
                                                               'CONSTRAINT `'+index_name+
                                                               orf_table_name+'` FOREIGN KEY (Seq_name) '
                                                               'REFERENCES '+seq_table_name+' (Seq_name) '
                                                                'ON DELETE CASCADE ON UPDATE RESTRICT'])
            make_index(connection_features, orf_table_name, 'ORF_id', 'unique_'+orf_table_name, 'UNIQUE ')

        activation_code = True
        ORF_dict = dict()
        for seq in data.split('>')[1:]:
            ORF_dict.update({seq.split('\n')[0]: sw.create_orf_objs(''.join(seq.split('\n')[1:]))})
            seq_classified = sw.create_obj(''.join(seq.split('\n')[1:]))
            dna = seq_classified.as_dna()
            add_code = add_to_table(connection_features, seq_table_name,
                                    ['Seq_name', 'Seq', 'Weight', 'Seq_identity'],
                                    ['"%s"' % seq.split('\n')[0], '"%s"' % str(dna), '"%s"' %dna.weight(),
                                     '"%s"' % seq.split('\n')[0]])
            if add_code == 1:
                print('Entry %s is already in database, skipping..' % seq.split('\n')[0])
                activation_code = False

        if activation_code:
            for seq in ORF_dict:
                orf_counter = 0
                for orf in ORF_dict[seq]:
                    orf_counter += 1
                    molecules = [orf.as_dna(), orf.as_rna(), orf.as_protein()]
                    columns = ['Seq_name', 'ORF_id']
                    values = ['"%s"' % seq.split('\n')[0], '"{}_{}"'.format(seq.split('\n')[0], orf_counter)]
                    for molecule in molecules:
                        columns += ['ORF_'+molecule.show_type(), 'Weight_'+molecule.show_type()]
                        values += ['"%s"' % str(molecule), '"%s"' % molecule.weight()]
                    add_code = add_to_table(connection_features, orf_table_name,
                                 columns, values)
                    if add_code == 1:
                        print('Entry {}_{} is already in database, skipping..'.format(seq.split('\n')[0], orf_counter))

    if get_sequence:
        if not (start or end or field or condition):
            raise Exception('You have to specify start position, end position, word to search (condition) and the table'
                            ' field by which you want to choose sequence (possible fields: Seq_name,'
                            ' ORF_id)')
        if field == 'Seq_name':
            data = acces_table(connection_features, seq_table_name, ['Seq'], 'Seq_name = ' + '"%s"' % condition)
            molecule = sw.create_obj(data[0][0][start:end])
        elif field == 'ORF_id':
            data = acces_table(connection_features, orf_table_name, ['ORF_'+extract], 'ORF_id = ' + '"%s"' % condition)
            molecule = sw.create_obj(data[0][0][start:end])
        else:
            raise Exception('Not supported field')
        extracted_seq = extraction(molecule, extract)
        print('Molecule ({}): {}; Weight: {}'.format(extract, extracted_seq[0], extracted_seq[1]))



