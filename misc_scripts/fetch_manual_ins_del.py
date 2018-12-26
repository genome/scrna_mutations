import pysam
import subprocess as sp
import argparse

npm1 = ['5','171410538', '171410546', '-', 'TCTG', 'frame_shift_ins']
cebpa = ['19', '33301989', '33301990', 'CC', '-', 'CEBPA', 'frame_shift_del']


def find_barcodes_ins(good_barcodes,bam_file1, gene, upn):
    """

    :param good_barcodes:
    :param bam_file1:
    :param gene:
    :param upn:
    :return:
    """
    v = []
    a = []
    barc = {}
    with open(good_barcodes, 'r') as barcodes:
        for line in barcodes:
            lines = line.strip()
            barc[lines] = 1

    with pysam.AlignmentFile(bam_file1, 'rb') as pile, \
        open(upn + '_' + gene + '_all', 'w+') as _all, \
            open(upn + '_' + gene + '_var', 'w+') as var:
        for pileupcolumn in pile.pileup(reference=npm1[0], start=int(npm1[1]) - 1, end=int(npm1[2]),
                                        truncate=True, stepper="nofilter",
                                        max_depth=100000000):
                for read in pileupcolumn.pileups:
                    if read.alignment.has_tag('CB') and read.alignment.has_tag('UB'):
                        if read.alignment.get_tag('CB') not in barc:  # filter reads with only good barcodes
                            continue
                        if read.is_refskip or read.is_del:
                            continue
                        CB = read.alignment.get_tag('CB')
                        UB = read.alignment.get_tag('UB')
                        qstring_all = CB + ":" + UB
                        a.append('{}'.format(qstring_all))
                        _all.write('{}\n'.format(qstring_all))

                        if '4I' in read.alignment.cigarstring:
                            CB = read.alignment.get_tag('CB')
                            UB = read.alignment.get_tag('UB')
                            alt_qstring = CB + ":" + UB
                            v.append('{}'.format(alt_qstring))
                            var.write('{}\n'.format(alt_qstring))
    return a, v


def find_barcodes_del(good_barcodes,bam_file, gene, upn):
    """

    :param bam_file:
    :param gene:
    :param upn:
    :return:
    """
    _all = open(upn + '_' + gene + '_all', 'w+')
    _var = open(upn + '_' + gene + '_var', 'w+')
    a = []
    with pysam.AlignmentFile(bam_file, 'rb') as pile:
        for pileupcolumn in pile.pileup(reference='19', start=int(33301989) - 1, end=int(33301990),
                                        truncate=True, stepper="nofilter",
                                        max_depth=100000000):
            for read in pileupcolumn.pileups:
                if not read.is_refskip:
                    if read.alignment.has_tag('CB') and read.alignment.has_tag('UB'):
                        CB = read.alignment.get_tag('CB')
                        UB = read.alignment.get_tag('UB')
                        a.append('{}:{}'.format(CB, UB))
                        _all.write('{}:{}'.format(CB, UB))

    v = []
    with pysam.AlignmentFile(bam_file, 'rb') as pile:
        for pileupcolumn in pile.pileup(reference='19', start=int(33301989)- 1, end=int(33301990),
                                        truncate=True, stepper="nofilter",
                                        max_depth=100000000):
            for read in pileupcolumn.pileups:
                if not read.is_refskip:
                    if '2D' in read.alignment.cigarstring and read.alignment.has_tag('CB') and \
                               read.alignment.has_tag('UB'):
                        CB = read.alignment.get_tag('CB')
                        UB = read.alignment.get_tag('UB')
                        v.append('{}:{}'.format(CB, UB))
                        _var.write('{}:{}'.format(CB, UB))
    bar = []
    with open(good_barcodes, 'r') as barcodes:
        for line in barcodes:
            lines = line.strip()
            bar.append(lines)
    #
    tbar = []
    for i in bar:
        for j in a:
            if i == j.split(':')[0]:
                tbar.append(j)

    vbar = []
    for i in bar:
        for j in v:
            if i == j.split(':')[0]:
                vbar.append(j)

    return tbar, vbar


def consensus(all_file, var, outfile, gene):
    """

    :param all_file:
    :param var:
    :param outfile:
    :param gene:
    :return:
    """
    uniq_raw_barcodes_ins = set()
    for tupr in all_file:
        # if tup not in uniq_barcodes_ins:
        uniq_raw_barcodes_ins.add(tupr)

    uniq_barcodes_ins = set()
    for tup in var:
        uniq_barcodes_ins.add(tup)

    indelbarU_count = {'ref': [], 'alt': []}

    for all_counts in all_file:
        if all_counts not in var:
            indelbarU_count['ref'].append(all_counts)
        else:
            indelbarU_count['alt'].append(all_counts)
    di = {}
    iUBuniq_barcodes_raw = []
    if gene == 'NPM1':
        head = "5\t171410539\t171410540\t-\tTCTG\tNPM1\tframe_shift_ins"
    else:
        head = "19\t33301989\t33301990\tCC\t-\tCEBPA\tframe_shift_del"
    with open(outfile + '_counts_CB.tsv', 'w+') as w, \
            open(outfile + '_counts_UB.tsv', 'w+') as sai:  # sai == UB , w == CB
        for iutags in uniq_raw_barcodes_ins:
            if iutags in indelbarU_count['ref'] and iutags in indelbarU_count['ref']:
                # print 'both: ', utags)
                iUnref = indelbarU_count['ref'].count(iutags)
                iUnalt = indelbarU_count['alt'].count(iutags)
                iutotl = iUnref + iUnalt
                iUn1ref = indelbarU_count['ref'].count(iutags)
                iUn1alt = indelbarU_count['alt'].count(iutags)
                iut1totl = iUn1ref + iUn1alt
                if iUnref / iutotl < float(0.75):
                    if iUnalt / iutotl < float(0.75):
                        continue
                    else:
                        iUn1alt = 1
                        iUn1ref = 0
                        iu1totl = iUn1ref + iUn1alt
                else:
                    if iUnalt / iutotl < float(0.75):
                        iUn1alt = 0
                        iUn1ref = 1
                        iu1totl = iUn1ref + iUn1alt
                    else:
                        iUn1alt = 1
                        iUn1ref = 1
                        iu1totl = iUn1ref + iUn1alt
            else:
                iUn1ref = list(set(indelbarU_count['ref'])).count(iutags)
                iUn1alt = list(set(indelbarU_count['alt'])).count(iutags)
                iu1totl = iUn1ref + iUn1alt
                iUnref = indelbarU_count['ref'].count(iutags)
                iUnalt = indelbarU_count['alt'].count(iutags)
                iutotl = iUnref + iUnalt

            if iUn1alt:
                iUBuniq_barcodes_raw.append(iutags)

            iUBtager = '{chrm}\t{bar}\t{r}\t{v}\t{tot}\n'.format(
                chrm=head,
                bar=iutags,
                r=iUnref,
                v=iUnalt, tot=iutotl)

            sai.write(iUBtager)

            if iutags.split(':')[0] not in di:
                di[iutags.split(':')[0]] = {'ref': iUn1ref, 'alt': iUn1alt}
            else:
                di[iutags.split(':')[0]]['ref'] += iUn1ref
                di[iutags.split(':')[0]]['alt'] += iUn1alt

        iUuni_alt = []
        itUuni_alt = []
        iCBuniq_barcodes_raw = []
        for ii, vv in di.items():
            iCBuniq_barcodes_raw.append(ii)
            iUuni_alt.append(vv['alt'])
            itot1 = vv['alt'] + vv['ref']
            itUuni_alt.append(itot1)
            iCBtager = '{chrm}\t{bar}\t{ref1}\t{alt1}\t{tot}\n'.format(
                chrm=head, bar=ii,
                alt1=vv['alt'], ref1=vv['ref'], tot=itot1)
            w.write(iCBtager)
        iUuni_alt_c = sum(iUuni_alt)
        itUuni_alt_c = sum(itUuni_alt)
        iUBfin_umi = ','.join(iUBuniq_barcodes_raw)
        iCBfin_umi = ','.join(iCBuniq_barcodes_raw)
        try:
            # _indelU = round(float(insmutU) / float(instotalU), 2)  # raw-vaf
            _uni_indelU = round(float(iUuni_alt_c) / float(itUuni_alt_c), 2)  # uni-vaf
        except ZeroDivisionError:
            _indelU = 0
            _uni_indelU = 0

        if gene == 'CEBPA':
            head = "19\t33301989\t33301990\tCC\t-\tCEBPA\tframe_shift_del"
            print('{}\t{}\t{}\t{}\t{}\t{}'.format(head, itUuni_alt_c, iUuni_alt_c, _uni_indelU, iCBfin_umi, iUBfin_umi))
        elif gene == 'NPM1':
            head = "5\t171410539\t171410540\t-\tTCTG\tNPM1\tframe_shift_ins"
            print('{}\t{}\t{}\t{}\t{}\t{}'.format(head, itUuni_alt_c, iUuni_alt_c, _uni_indelU, iCBfin_umi, iUBfin_umi))

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse CB barcodes from Single cell rna seq data')
    parser.add_argument('bam_file', help='BAM file')
    parser.add_argument('gene', help='NPM1 or CEBPA')
    parser.add_argument('barcodes', help='list of good barcodes file')
    parser.add_argument('upn', help='upn/sample name: will be used as prefix for out_file')
    args = parser.parse_args()

    _bam_file = args.bam_file
    _gene = args.gene
    barcodes_good = args.barcodes
    _outfile = args.upn

    if _gene == 'NPM1':
        _a, _v = find_barcodes_ins(barcodes_good,_bam_file, _gene, _outfile)
        consensus(_a, _v, _outfile, _gene)
    elif _gene == 'CEBPA':
        _a, _v = find_barcodes_del(barcodes_good,_bam_file, _gene, _outfile)
        consensus(_a, _v, _outfile, _gene)



