import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from batch_run_hmm_func import batch_run_hmm

'''
fetch up and down by user-defined bp, default = 10000 bp
ext_sta_loc = max(0, up)
ext_end_loc = min(seq_feature.location.end + 10000, len(str(seq_record.seq)))
'''


def fetch_on_bp(pfam_on_bp,pre_in_min,pre_in_max,bp_up,bp_down, parent_gbk_path, ext_gbk_dir):
    # gbk_file = open('C:/Coding/domain_scanner/test_data/test_data_out/bgc_20/bgc_20_bio_GBk.gbk', 'r')
    gbk_file = open(parent_gbk_path, 'r')
    for seq_record in SeqIO.parse(gbk_file, 'genbank'):
        # print(seq_record.features[0].qualifiers['translation'])
        # print(type(seq_record.features[0].qualifiers))
        for seq_feature in seq_record.features:
            pre_seq = seq_feature.qualifiers['translation'][0]
            pre_len = len(seq_feature.qualifiers['translation'][0])
            # rbs_type = str(seq_feature.qualifiers['note'][0]).split(';')[3]
            start_cond = str(str(seq_feature.qualifiers['note'][0]).split(';')[2]).split('=')[1]
            #  user defines pre len
            if int(pre_in_min) <= pre_len <= int(pre_in_max) and start_cond in ['ATG', 'TTG', 'GTG']: #  and rbs_type != 'rbs_motif=None'
                pre_seq = seq_feature.qualifiers['translation']
                pre_id = seq_feature.qualifiers['note'][0].split(';')[0].replace('=', '-')
                # print(pre_id)
                # up = seq_feature.location.start - 10000  用 str.isdigit 是因为start和end会带 '>' 符号
                up = int(''.join(filter(str.isdigit, str(seq_feature.location.start)))) - int(bp_up)
                # down = seq_feature.location.end + 10000 - len(str(seq_record.seq))
                ext_sta_loc = max(0, up)
                ext_end_loc = min(int(''.join(filter(str.isdigit, str(seq_feature.location.end)))) +
                                  int(bp_down), len(str(seq_record.seq)))
                sub_record = seq_record[ext_sta_loc:ext_end_loc]
                sub_record.id = seq_record.id + '_' + str(ext_sta_loc) + '-' + str(ext_end_loc) + 'bp'
                sub_record.name = seq_record.id + '_partial'
                sub_record.description = 'extracted from_' + seq_record.id
                sub_record.annotations = {"molecule_type": "DNA"}

                file_name_prefix = seq_record.id + '__' + pre_id + '__' + str(ext_sta_loc) + '-' + str(
                    ext_end_loc) + 'bp' + '_sub_bio_GBk'  # pre_region记录目标pre的

                sub_gbk = open(ext_gbk_dir + '/' + file_name_prefix + '.gbk', 'w')
                sub_faa_handle = open(ext_gbk_dir + '/' + file_name_prefix + '.faa', 'w')
                for sub_sel_seq_feature in sub_record.features:
                    sub_faa_handle.write('>' + file_name_prefix + ';' +
                                         sub_sel_seq_feature.qualifiers['note'][0].split(';')[0].replace('=',
                                                                                                         '-') + '\n' +
                                         sub_sel_seq_feature.qualifiers['translation'][0] + '\n')
                SeqIO.write(sub_record, sub_gbk, 'genbank')
                sub_gbk.close()
                sub_faa_handle.close()
                # pfam number user input
                hmm_file_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(ext_gbk_dir))), 'hmm_info/fet_hmm_id.hmm')
                bis_gbk_out = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(ext_gbk_dir))), 'bis_gbk_out/')
                if not os.path.exists(bis_gbk_out):
                    os.makedirs(bis_gbk_out)
                pfam_to_hmm = pfam_on_bp
                batch_run_hmm(pfam_to_hmm,ext_gbk_dir + '/' + file_name_prefix + '.faa', hmm_file_dir, bis_gbk_out)
                # pfam_no
                # faa_in_path, fet_hmm_file, bgc_hit_dir
                # fet_hmm_file = ext_gbk_dir 上级目录/hmm_info/fet_hmm_id.hmm
                # fet_hmm_file = make ext_gbk_dir 上级目录/bgc_hits
    gbk_file.close()

# fetch_on_bp('/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_8/bgc_8_bio_GBk.gbk',
#             '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_8/')
