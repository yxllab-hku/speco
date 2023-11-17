import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from batch_run_hmm_func import batch_run_hmm


def fetch_on_orf(pfam_on_orf, pre_in_min,pre_in_max,orf_up, orf_down, parent_gbk_path, ext_gbk_dir):
    id2index = {}  # key = cds_id, value = 'CDS_start','CDS_end'
    gbk_file = open(parent_gbk_path, 'r')
    # gbk_file = open('C:/Coding/domain_scanner/test_data/test_data_out/bgc_31/bgc_31_bio_GBk.gbk', 'r')
    for seq_record in SeqIO.parse(gbk_file, 'genbank'):
        for i in range(0, len(seq_record.features)):
            id2index_sub = {
                str(i + 1): [''.join(filter(str.isdigit, str(seq_record.features[i].location.start))),
                             ''.join(filter(str.isdigit, str(seq_record.features[i].location.end)))]}
            id2index.update(id2index_sub)

        for seq_feature in seq_record.features:
            pre_len = len(seq_feature.qualifiers['translation'][0])
            # rbs_type = str(seq_feature.qualifiers['note'][0]).split(';')[3]
            start_cond = str(str(seq_feature.qualifiers['note'][0]).split(';')[2]).split('=')[1]
            #  user defines pre len
            if int(pre_in_min) <= pre_len <= int(pre_in_max) and start_cond in ['ATG', 'TTG', 'GTG']: #  and rbs_type != 'rbs_motif=None'
                pre_id = seq_feature.qualifiers['note'][0].split(';')[0].replace('=', '-')
                pre_loc = int(pre_id.split('_')[-1])
                # num_orf = 10
                num_orf_up = int(orf_up)
                num_orf_down = int(orf_down)
                cds_sta = max(1, pre_loc - num_orf_up)
                cds_end = min(pre_loc + num_orf_down, len(seq_record.features))
                ext_sta_loc = id2index[str(cds_sta)][0]
                ext_end_loc = id2index[str(cds_end)][1]
                sub_record = seq_record[int(ext_sta_loc):int(ext_end_loc)]
                sub_record.id = 'pre_' + str(pre_loc)
                sub_record.name = seq_record.id + '_partial'
                sub_record.description = 'extracted from_' + seq_record.id
                sub_record.annotations = {"molecule_type": "DNA"}

                file_name_prefix = seq_record.id + '_pre__' + pre_id + '__up_down_' + str(num_orf_up) + '_' + str(
                    num_orf_down) + '_sub_bio_GBk'
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
                hmm_file_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(ext_gbk_dir))),
                                            'hmm_info/fet_hmm_id.hmm')
                bis_gbk_out = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(ext_gbk_dir))),
                                           'bis_gbk_out/')
                if not os.path.exists(bis_gbk_out):
                    os.makedirs(bis_gbk_out)
                pfam_to_hmm = pfam_on_orf
                batch_run_hmm(pfam_to_hmm, ext_gbk_dir + '/' + file_name_prefix + '.faa', hmm_file_dir, bis_gbk_out)
                # pfam_no

    gbk_file.close()

# fetch_on_orf('/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_8/bgc_8_bio_GBk.gbk',
#              '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_8/')
