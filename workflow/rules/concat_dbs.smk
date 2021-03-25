# Original source: https://github.com/laelbarlow/jackhmmerer
# Licence notice: MIT License, copyright (c) 2021 Lael D. Barlow

rule concat_dbs:
    """
    Concatenate database files to make one big FASTA with all the sequences to
    be searched.
    """
    input:
        db_files = db_fasta_files

    output:
        db_file = 'results/concat_dbs/database.faa'

    run:
        with open(output.db_file, 'w') as o:
            for f in input.db_files:
                with open(f) as infh:
                    for i in infh:
                        # Format EukProt sequence headers differently.
                        if i.startswith('>EP'):
                            o.write(i.split(' ')[0] + '\n')
                        # By default, insert filenames in headers.
                        elif i.startswith('>'):
                            o.write('>' + os.path.basename(f).rsplit('.', \
                            1)[0].replace(' ', '_') + '__' + i[1:].split(' ')[0] + \
                            '\n')
                        else:
                            o.write(i)


