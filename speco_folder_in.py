import os
import subprocess
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

parser = argparse.ArgumentParser(
    description='Short Peptide-Enzyme Co-Occurrence analysis. speco requires prodigal-short, hmmsearch, bigscape and fasttree. If <20000 bp sequence exists, multiple non-coding nucleotide "N" will be concatenated in the 3-end accordingly in order to be processed by prodigal-short. Usage example:python3 speco.py -p PF00067 -l 20 -m 100 -t orf -u 2 -d 2 -i /data/lab/ -o /data/lab/output')
parser.add_argument('-p', '--pfam_no', type=str, metavar='', required=True,
                    help='Pfam number of enzyme adjacent to short peptide, e.g., "PF04055". For multiple pfam numbers, each of them should be seperated by comma')
parser.add_argument('-l', '--pre_min', type=int, metavar='', required=True,
                    help='Minimal length of short peptide, input type is integer. The lower limit of here is "19" which is defined by default output of prodigal-short')
parser.add_argument('-m', '--pre_max', type=int, metavar='', required=True,
                    help='Maximal length of short peptide, input type is integer')
parser.add_argument('-t', '--sele_type', type=str, metavar='', required=True,
                    help='Fetching type should be "bp" (nucleotide base pair) or "orf" (open reading frame)')
parser.add_argument('-u', '--up', type=int, metavar='', required=True,
                    help='Length (or number) of bp (orf) at the upstream side of short peptide. If input length (number) exceed the bondary, fetching will be started from the 1st position of the sequence. E.g., if "sele_type" is "bp" and "up" is 10000, upstream 10000 bp will be fetched')
parser.add_argument('-d', '--down', type=int, metavar='', required=True,
                    help='Length(or number) of bp(orf) at the downstream side of short peptide. If input length (number) exceed the bondary, fetching will be ended in the last position of the sequence. E.g., if "sele_type" is "orf" and "down" is 5, downstream 5 orfs will be fetched')
parser.add_argument('-i', '--fna_in_path', type=str, required=True, metavar='',
                    help='Directory that contains all the nucleotide sequences to be analyzed, which must be in fasta '
                         'format and sequence id should not contain illegal characters like *, &, ^, $, _ |, / ? '
                         'etc. E.g., /home/lab/')
parser.add_argument('-o', '--out_dir', type=str, required=True, metavar='',
                    help='Output directory where you want to save the output data. E.g., /home/lab/test')
args = parser.parse_args()


def run_prod_hmmf(pfam_no, pre_min, pre_max, sele_type, up, down, fna_in_path, out_dir):
    print('Fetching %s.hmm file' % pfam_no)
    pfam_acc_file = '/data/bbhe/Hmmscan/pfam34.0/pfam_34.0_acc.txt'  # 固定
    hmm_file_path = '/data/bbhe/Hmmscan/pfam34.0/Pfam-A.hmm'  # 固定脚本文件所在路径

    pfam_34_id = {}
    for line in open(pfam_acc_file, 'r'):
        pfam_34_id[line.split('\t')[0]] = line.split('\t')[1].replace('\n', '')

    hmm_file_dir = out_dir + '/' + 'hmm_info'
    if not os.path.exists(hmm_file_dir):
        os.makedirs(hmm_file_dir)
    with open(hmm_file_dir + '/fet_hmm_id.txt', 'w') as f:
        if ',' in pfam_no:
            for pfam_id in pfam_no.split(','):
                if pfam_id in pfam_34_id.keys():
                    f.write(pfam_34_id[pfam_id] + '\n')
                if pfam_id not in pfam_34_id.keys():
                    print(
                        'Error: %s is not a valid pfam number, you may refer to Pfam database: http://pfam.xfam.org/' % pfam_id)
                    sys.exit()
        else:
            if str(pfam_no) in pfam_34_id.keys():
                f.write(pfam_34_id[pfam_no])
            if str(pfam_no) not in pfam_34_id.keys():
                print(
                    'Error: %s is not a valid pfam number, you may refer to Pfam database: http://pfam.xfam.org/' % pfam_no)
                sys.exit()
    f.close()

    with open(hmm_file_dir + '/input_hmm_id.txt', 'w') as f_in_id:
        f_in_id.write(pfam_no)
    f_in_id.close()

    fet_hmm_id = hmm_file_dir + '/fet_hmm_id.txt'
    fet_hmm_file = hmm_file_dir + '/fet_hmm_id.hmm'

    hmmfetch_run = ['hmmfetch', '-o', fet_hmm_file, '-f', hmm_file_path, fet_hmm_id]
    subprocess.run(hmmfetch_run)
    print('Running prodigal-short, processing input sequences')
    files = os.listdir(fna_in_path)
    for each_nuc_file in files:
        print('processing %s' % str(each_nuc_file))
        nuc_seq = open(os.path.join(fna_in_path, each_nuc_file), 'r')
        for seq_record in SeqIO.parse(nuc_seq, 'fasta'):
            each_fna_in = out_dir + '/' + 'prodigal_out' + '/' + seq_record.id
            if not os.path.exists(each_fna_in):
                os.makedirs(each_fna_in)
            fna_out_dir = each_fna_in + '/'
            with open(fna_out_dir + seq_record.id + '.fna', 'w') as f:
                if len(seq_record) >= 20000:
                    f.write('>' + seq_record.id + '\n')
                    f.write(str(seq_record.seq))
                else:
                    print(
                        'Note: %s is shorter than 20k bp. Concatenating %d non-coding nucleotide "N" for prodigal-short '
                        'processing' % (seq_record.id, 20001 - len(seq_record)))
                    f.write('>' + seq_record.id + '\n')
                    f.write(str(seq_record.seq) + 'N' * (20001 - len(seq_record)))
            f.close()

            prodigal_in = fna_out_dir + seq_record.id + '.fna'
            prodigal_faa_path = fna_out_dir + seq_record.id + '.faa'
            prodigal_gbk_path = fna_out_dir + seq_record.id + '_gbk.txt'
            prodigal_run = ['prodigal-short', '-i', prodigal_in, '-a', prodigal_faa_path, '-o', prodigal_gbk_path]
            p1 = subprocess.run(prodigal_run, capture_output=True)

            '''
            write prodigal_out gbk to biopython readable gbk
            '''
            line = []
            with open(prodigal_gbk_path, 'r+') as rewrite_gbk:
                for lines in rewrite_gbk:
                    if "DEFINITION  " in lines:
                        line.append(
                            (lines.split('  ')[
                                 0] + '\n'))  # 原始prodigal输出的gbk文件‘DEFINITION’过长，会报错‘Annotation %r too long"
                    else:
                        line.append(lines)
            rewrite_gbk.close()
            line.insert(0, 'LOCUS       ' + seq_record.id + '                 '
                        + str(len(seq_record.seq)) +
                        ' bp    DNA              UNK 01-JAN-1980\n')
            line.insert(-1, 'ORIGIN\n')
            s = ''.join(line)
            with open(fna_out_dir + seq_record.id + '_readable_GBK.gbk', 'w+') as rewrite_gbk_out:
                rewrite_gbk_out.write(s)
            rewrite_gbk_out.close()

            '''
            write prodigal .faa into readable .gbk
            1. build a dic to store protein and desc information:{"desc","prot_seq"}
            '''
            prot_file = open(prodigal_faa_path, 'r')
            prot_dic = {}
            for prot_record in SeqIO.parse(prot_file, 'fasta'):
                # print(prot_record.description.split(' ')[-1])  # 如果读取id，会被‘#’打断，读取desc则不会
                # print(prot_record.seq)
                prot_dic[prot_record.description.split(' # ')[-1]] = str(prot_record.seq).replace('*', '')
            prot_file.close()

            with open(fna_out_dir + seq_record.id + '_readable_GBK.gbk', 'r') as rewrite_gbk_out:
                for gbk_record in SeqIO.parse(rewrite_gbk_out, 'genbank'):
                    for gbk_feature in gbk_record.features:
                        prot_id = (gbk_feature.qualifiers['note'][0].split(";conf=")[0])
                        # print(prot_id)
                        gbk_feature.qualifiers['translation'] = prot_dic[prot_id]
                writeAA2gbk = open(fna_out_dir + seq_record.id + '_writeAA2gbk.gbk', 'w')
                SeqIO.write(gbk_record, writeAA2gbk, 'genbank')
                writeAA2gbk.close()
            rewrite_gbk_out.close()

            '''
            fetch the whole record_feature of writeAA2gbk into each_fna file to give bio.gbk
            '''
            writeAA2gbk = open(fna_out_dir + seq_record.id + '_writeAA2gbk.gbk', 'r')
            each_fna = open(fna_out_dir + seq_record.id + '.fna', 'r')
            for seq_record_each in SeqIO.parse(each_fna, 'fasta'):
                seq_object = seq_record_each.seq
                # print(seq_object)
                record = SeqRecord(seq_object, id=seq_record_each.id, name=seq_record_each.id,
                                   description=seq_record_each.id, annotations={'molecule_type': 'DNA'})
                record.features = [gbk_record_feature for gbk_record in
                                   SeqIO.parse(writeAA2gbk, 'genbank')
                                   for gbk_record_feature in
                                   gbk_record.features]
                # print(record.features)
                bio_gbk = open(fna_out_dir + seq_record_each.id + '_bio_GBk.gbk', 'w')
                SeqIO.write(record, bio_gbk, 'genbank')
                bio_gbk.close()

                pfam_id_to_sele = pfam_no
                bp_up = int(up)
                bp_down = int(down)
                num_orf_up = int(up)
                num_orf_down = int(down)
                pre_in_min = pre_min
                pre_in_max = pre_max
                if sele_type == 'bp':
                    from pre_sele_bp import fetch_on_bp
                    fetch_on_bp(pfam_id_to_sele, pre_in_min, pre_in_max, bp_up, bp_down,
                                fna_out_dir + seq_record_each.id + '_bio_GBk.gbk', fna_out_dir)
                if sele_type == 'orf':
                    from pre_sele_orf import fetch_on_orf
                    fetch_on_orf(pfam_id_to_sele, pre_in_min, pre_in_max, num_orf_up, num_orf_down,
                                 fna_out_dir + seq_record_each.id + '_bio_GBk.gbk', fna_out_dir)
                if sele_type not in ['orf', 'bp']:
                    print('Invalid selection type, should be "bp" or "orf".')
            each_fna.close()
        nuc_seq.close()
    print('Prodigal-short completed')
    print('GBK file fetching completed')
    print('Start running bigscape')

    from batch_run_bigscape import batch_run_big
    batch_run_big('%s/bis_gbk_out/' % out_dir, out_dir)


#
#
# run_prod_hmmf('PF04055', 10, 101, 'orf', 8, 8,
#               '/home/bbhe/Python_Script/domain_scanner/test_data/bgc_8.fna',
#               '/home/bbhe/Python_Script/domain_scanner/test_data/test_data_out')
if __name__ == '__main__':
    run_prod_hmmf(args.pfam_no.upper(), args.pre_min, args.pre_max, args.sele_type, args.up, args.down,
                  args.fna_in_path,
                  args.out_dir)
