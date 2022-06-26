import os
import shutil
import subprocess
from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# from batch_run_prodigal import batch_run_prod

# fna_in_path, out_dir
# batch_run_prod('/home/bbhe/Python_Script/domain_scanner/test_data/bgc_32_33.fna',
#                '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out')

# fetch_on_orf('/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_8/bgc_8_bio_GBk.gbk',
#              '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_8/')

# def batch_run_hmm(pfam_no,faa_in_path,bgc_hit_dir)
    # # hmm_file_path = '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/hmm_file/rSAM.hmm'
hmm_file_path = '/data/bbhe/Hmmscan/pfam34.0/Pfam-A.hmm'  # 固定脚本文件所在路径
fet_hmm_id = '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/hmm_file/fet_hmm_id.txt'
# out_dir + /hmm_file/fet_hmm_id.txt
fet_hmm_file = '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/hmm_file/fet_hmm_id.hmm'
# out_dir + /hmm_file/fet_hmm_id.hmm
pfam_acc_file = '/data/bbhe/Hmmscan/pfam34.0/pfam_34.0_acc.txt'  # 固定
pfam_no = 'PF04055,PF01478'  # 参数，user input, seperated by comma

pfam_34_id = {}
pfam_no_list = pfam_no.split(',')

for line in open(pfam_acc_file,'r'):
    pfam_34_id[line.split('\t')[0]] = line.split('\t')[1].replace('\n','')

with open(fet_hmm_id, 'w') as f:
    if ',' in pfam_no:
        for pfam_id in pfam_no.split(','):
            f.write(pfam_34_id[pfam_id] + '\n')
    else:
        f.write(pfam_34_id[pfam_no])
f.close()

hmmfetch_run = ['hmmfetch', '-o', fet_hmm_file, '-f', hmm_file_path, fet_hmm_id]
subprocess.run(hmmfetch_run)
faa_in_path = '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_8/bgc_8__ID-1_8__0-17177bp_sub_bio_GBk.faa'
# 参数 prodigal_faa_path
hmm_tblout_path = '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/hmm_out/bgc_8_ID=1_8_0-17177bp_sub_bio_GBk_hmm_out.txt'
# prodigal_faa_path.replace('.faa','_hmm_out.txt')
bgc_hit_dir = '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_hits'
# 参数
hmmsearch_run = ['hmmsearch', '-E', '1e-5', '--domtblout', hmm_tblout_path, fet_hmm_file, faa_in_path]  # domtblout
subprocess.run(hmmsearch_run)

# hmm_handle = open(hmm_tblout_path,'r')
pfam_out_acc = []
enzyme_id = []  # 修改gbk文件时region取pre_id 和 enzyme_id 中距离最大的值
bio_syn_id = []
bio_syn_region = []
for result in SearchIO.parse(hmm_tblout_path,'hmmsearch3-domtab'):
    # if result.accession.split('.')[0] == str(pfam):
    pfam_out_acc.append(result.accession.split('.')[0])
    enzyme_id.append(result.hit_keys[0].split(';')[-1])
    enzyme_id.append(result.accession.split('.')[0])
    print(result.hit_keys, result.id, result.accession.split('.')[0])  # output this
    bio_syn_id.append(result.hit_keys[0].split('__')[1].replace('-','='))
    bio_syn_id.append(result.hit_keys[0].split(';')[-1].replace('-','='))
    uni_syn_id = list(set(bio_syn_id))  # set()后变成集合，集合没有索引取值操作，再转换成list就可以索引取值

    sta_id = min([int(uni_syn_id[i].split('_')[-1]) for i in range(0,len(uni_syn_id))])  # 最小orf
    end_id = max([int(uni_syn_id[i].split('_')[-1]) for i in range(0,len(uni_syn_id))])  # 最大orf

    for reg in range(sta_id-1,end_id+2):  # 选定pre-修饰酶向外各扩增1个orf来做bigscape分析
        bio_syn_region.append('ID=1_' + str(reg))
    print(bio_syn_region)

if set(pfam_no_list) == set(pfam_out_acc):  # 会有多个orf是input的pfam number，因此需要set()去重复
    gbk_src = faa_in_path.replace('.faa', '.gbk')
    bio_syn_gbk = open(gbk_src,'r')
    new_file_name = faa_in_path.split('/')[-1].replace('.faa','__' + '__'.join(enzyme_id)+'__region_.gbk')
    new_path = os.path.join(bgc_hit_dir,new_file_name)
    bgc_hit = open(new_path,'w')
    for seq_record in SeqIO.parse(bio_syn_gbk,'genbank'):
        bgc_hit_record = seq_record
        for i in range(0,len(seq_record.features)):
            orf_id = seq_record.features[i].qualifiers['note'][0].split(';')[0]
            if orf_id in bio_syn_region:
                bgc_hit_record.features[i].qualifiers['gene_kind'] = "biosynthetic"
                print(bgc_hit_record.features[i])
        SeqIO.write(bgc_hit_record,bgc_hit,'genbank')

# big_scape_out = '/home/bbhe/TxtProcess/test_batch_mod_bigscape_out_11'
# big_scape_run = ['bigscape', '-i',bgc_hit_dir, '--outputdir', big_scape_out, '--hybrids-off', '--include_singletons', '--cutoffs', '1.0', '--clan_cutoff', '1.0', '1.0']
# subprocess.run(big_scape_run)




# fna_in_path, out_dir
# batch_run_prod('/home/bbhe/Python_Script/domain_scanner/test_data/bgc_32_33.fna',
#                '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out')

# fetch_on_orf('/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_8/bgc_8_bio_GBk.gbk',
#              '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out/bgc_8/')

    # shutil.copy(gbk_src,new_path)
# print(pfam_out_acc)
#     print(result.hit_keys, result.id, result.seq_len, result.description, result.accession, result.hits)

    # print(result.hits)
    # print(dir(result))
    # print(result.hit_keys)
    # num_hits = len(result.hits)
    # if num_hits > 0:
    #     for i in range(0,num_hits):
    #         print(result.hits[i].accession, result.hits[i].id,
    #               result.hits[i].description,result.hits[i].description.query_id)
            # print(result.hits[i].accession)
# df = pd.read_csv(pfam_acc_file) #,sep = '; ', header = None, names=['input','id'])
# print(df.head(2))

# df = pd.read_excel(pfam_acc_file,converters={'input': str, 'no': str})
# # print(df.head(2))
# # print(df.shape)
# # print(df.dtypes)
# # print(df.loc[0])
# # print(df.loc[1:3])
# myseries = pd.Series(df.input)
# # print(myseries[myseries == 'PF04055'])
# pf_idx = myseries[myseries == 'PF04055'].index
# # print(Name(myseries.get_loc(PF04055)))
# # print(df.loc[[15028,15029]])
# # print(df.loc[pf_idx])
# print(df.loc[pf_idx],['no'])
# # print(df.input)