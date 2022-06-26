import os
import shutil
import subprocess
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run_hmm_fet(pfam_no, hmm_file_dir):
    pfam_acc_file = '/data/bbhe/Hmmscan/pfam34.0/pfam_34.0_acc.txt'  # 固定
    hmm_file_path = '/data/bbhe/Hmmscan/pfam34.0/Pfam-A.hmm'  # 固定脚本文件所在路径

    pfam_34_id = {}
    for line in open(pfam_acc_file, 'r'):
        pfam_34_id[line.split('\t')[0]] = line.split('\t')[1].replace('\n', '')

    with open(hmm_file_dir + '/fet_hmm_id.txt', 'w') as f:
        if ',' in pfam_no:
            for pfam_id in pfam_no.split(','):
                f.write(pfam_34_id[pfam_id] + '\n')
        else:
            f.write(pfam_34_id[pfam_no])
    f.close()

    fet_hmm_id = hmm_file_dir + '/fet_hmm_id.txt'
    fet_hmm_file = hmm_file_dir + '/fet_hmm_id.hmm'

    hmmfetch_run = ['hmmfetch', '-o', fet_hmm_file, '-f', hmm_file_path, fet_hmm_id]
    subprocess.run(hmmfetch_run)


run_hmm_fet('PF04055,PF01478', '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/hmm_file/')
