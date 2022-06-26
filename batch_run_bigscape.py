import os
import subprocess
import sys
from Bio import SeqIO


# bis_out = '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bis_gbk_out/'
# out = '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/'


def batch_run_big(sel_gbk_dir, output):
    if len(os.listdir(sel_gbk_dir)) == 0:
        print('Pfam no of interest is not found in the downstream or upstream of short peptide.')
        print('0 gbk file fetched')
        sys.exit()
    files = os.listdir(sel_gbk_dir)
    pre_enz_out = '%s/pre_enz.txt' % output
    with open(pre_enz_out, 'w') as f:
        f.write(
            '%s%s%s%s%s' % ('gbk_file\t','nucleotide_id\t','seq_pre\t', 'adjacent protein\t', 'pf number\n'))
        for each_file in files:
            id_to_seq = {}
            pre_enz_set = []
            enz_id = each_file.replace('-', '=').split('__')
            nuc_id = each_file.split('_pre__')[0]
            for i in range(3, len(enz_id) - 1, 2):
                pre_enz_pai = [enz_id[1], enz_id[i], enz_id[i + 1]]
                pre_enz_set.append(pre_enz_pai)

            for record in SeqIO.parse('%s%s' % (sel_gbk_dir, each_file), 'genbank'):
                for record_feature in record.features:
                    seq_id = record_feature.qualifiers['note'][0].split(';')[0]
                    seq = record_feature.qualifiers['translation'][0]
                    id_to_seq[seq_id] = seq
            for sub_pre_enz in pre_enz_set:
                new_sub_pre = ['%s\t%s\t' %(str(each_file),nuc_id)]
                new_sub_pre.extend(
                    ['%s\t' % id_to_seq[i] if i in id_to_seq else '%s\n' % i for i in sub_pre_enz]
                )
                f.writelines(new_sub_pre)
    f.close()

    big_scape_run = ['bigscape', '-i', sel_gbk_dir, '--outputdir', '%s/bigscape_out' % output, '--hybrids-off',
                     '--include_singletons',
                     '-c', '48', '--clan_cutoff', '0.5', '0.8', '--mode', 'auto']
    subprocess.run(big_scape_run)

#
# batch_run_big('/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bis_gbk_out/',
#               '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/')
