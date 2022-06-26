import os
import shutil
import subprocess
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def batch_run_hmm(pfam_no,faa_in_path, fet_hmm_file,bgc_hit_dir):
    hmm_tblout_path = faa_in_path.replace('.faa', '_hmm_out.txt')
    if os.path.getsize(faa_in_path) != 0:
        hmmsearch_run = ['hmmsearch', '-E', '1e-5', '--cpu', '24','--domtblout', hmm_tblout_path, fet_hmm_file, faa_in_path]  # domtblout
        subprocess.run(hmmsearch_run)

    pfam_no_list = pfam_no.split(',')
    pfam_out_acc = []
    enzyme_id = []  # 修改gbk文件时region取pre_id 和 enzyme_id 中距离最大的值
    bio_syn_id = []
    bio_syn_region = []
    for result in SearchIO.parse(hmm_tblout_path, 'hmmsearch3-domtab'):
        # if result.accession.split('.')[0] == str(pfam):
        pfam_out_acc.append(result.accession.split('.')[0])
        enzyme_id.append(result.hit_keys[0].split(';')[-1])
        enzyme_id.append(result.accession.split('.')[0])
        bio_syn_id.append(result.hit_keys[0].split('__')[1].replace('-', '='))
        bio_syn_id.append(result.hit_keys[0].split(';')[-1].replace('-', '='))
        uni_syn_id = list(set(bio_syn_id))  # set()后变成集合，集合没有索引取值操作，再转换成list就可以索引取值
        print(f'search for {uni_syn_id}')

        sta_id = min([int(uni_syn_id[i].split('_')[-1]) for i in range(0, len(uni_syn_id))])  # 最小orf
        end_id = max([int(uni_syn_id[i].split('_')[-1]) for i in range(0, len(uni_syn_id))])  # 最大orf

        for reg in range(sta_id - 1, end_id + 2):  # 选定pre-修饰酶向外各扩增1个orf来做bigscape分析
            bio_syn_region.append('ID=1_' + str(reg))

    if set(pfam_no_list) == set(pfam_out_acc):  # 会有多个orf是input的pfam number，因此需要set()去重复
        gbk_src = faa_in_path.replace('.faa', '.gbk')
        bio_syn_gbk = open(gbk_src, 'r')
        new_file_name = faa_in_path.split('/')[-1].replace('.faa', '__' + '__'.join(enzyme_id) + '__region_.gbk')
        new_path = os.path.join(bgc_hit_dir, new_file_name)
        bgc_hit = open(new_path, 'w')
        for seq_record in SeqIO.parse(bio_syn_gbk, 'genbank'):
            bgc_hit_record = seq_record
            for i in range(0, len(seq_record.features)):
                orf_id = seq_record.features[i].qualifiers['note'][0].split(';')[0]
                if orf_id in bio_syn_region:
                    bgc_hit_record.features[i].qualifiers['gene_kind'] = "biosynthetic"
                    # print(bgc_hit_record.features[i])
            SeqIO.write(bgc_hit_record, bgc_hit, 'genbank')


# def batch_run_hmm(pfam_no, faa_in_path, fet_hmm_file):
# batch_run_hmm( '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_8/bgc_8__ID'
#                                  '-1_16__1462-21561bp_sub_bio_GBk.faa',
#               '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/hmm_file/fet_hmm_id.hmm',
#               '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_hits/')
