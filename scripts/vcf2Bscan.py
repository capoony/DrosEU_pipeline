#!/usr/bin/env python2.7
import vcf
import sys
import argparse
import os
import csv
import operator


def SNPprinter(vcf_reader,cont,Q):
    print '\t'.join(["CHROM","POS","AVG_DP"])
    if cont == "all":
        for record in vcf_reader:
            #if the record is a SNP
            if record.is_snp:
                dp_l = [sample['DP'] for sample in record.samples]
                printline = [record.CHROM,str(record.POS), str(sum(dp_l)/len(dp_l))]
                print '\t'.join(printline)
    else:
        for record in vcf_reader.fetch(cont):
            #if the record is a SNP
            if record.is_snp:
                dp_l = [sample['DP'] for sample in record.samples]
                printline = [record.CHROM,str(record.POS), str(sum(dp_l)/len(dp_l))]
                print '\t'.join(printline)

def SNPdensity(vcf_reader,cont,W_l,dp):
    print '\t'.join(["CHROM","WIND","START","END","N_SNPs"])
    if cont == "all":
        for contig in vcf_reader.contigs.keys():
            #always start at POS 1
            cont_l = vcf_reader.contigs[contig].length
            start = 0
            end = start + W_l
            wind = 1
            while end < cont_l:
                snps = 0
                for record in vcf_reader.fetch(contig,start,end):
                    dp_l = [sample['DP'] for sample in record.samples]
                    if any([d is None for d in dp_l]) == True:
                        #one sample has no depth for SNP: Skip it!
                        #print record, dp_l
                        continue
                    m_dp = float(sum(dp_l))/float(len(dp_l))
                    #if the record is a SNP and mean depth > dp
                    if record.is_snp and m_dp > dp:
                        snps = snps + 1
                LINE = [contig,str(wind),str(start),str(end),str(snps)]
                print '\t'.join(LINE)
                start = end
                end = start + W_l
                if end > cont_l:
                    end = cont_l
                wind = wind + 1
    else:
        #print vcf_reader.contigs[cont] #script tester line
        cont_l = vcf_reader.contigs[cont].length
        start = 0
        end = start + W_l
        wind = 1
        while end < cont_l:
            snps = 0
            for record in vcf_reader.fetch(cont,start,end):
                dp_l = [sample['DP'] for sample in record.samples]
                if any([d is None for d in dp_l]) == True:
                    #one sample has no depth for SNP: Skip it!
                    #print record, dp_l
                    continue
                m_dp = float(sum(dp_l))/float(len(dp_l))
                #if the record is a SNP and mean depth > dp
                if record.is_snp and m_dp > dp:
                    snps = snps + 1
            LINE = [cont,str(wind),str(start),str(end),str(snps)]
            print '\t'.join(LINE)
            start = end
            end = start + W_l
            if end > cont_l:
                end = cont_l
            wind = wind + 1
        

def Bscandict(vcf_reader,cont,prefix,outdir):
    if cont == "all":
        loci=0
        loci_t = '[loci]='
        pops = len(vcf_reader.samples)
        pops_t = '[populations]='+str(pops)
        snps = {}
        #for contig in vcf_reader.contigs.keys():
        snps={}
        snp = 1
        #also write snp id and snp position to separate file
        if outdir == 'none':
            snpfildir = os.getcwd()
        else:
            snpfildir = outdir
        snpfil = open(snpfildir+'/'+prefix+'_'+cont+'_snps.txt','w')
        snpfilwriter = csv.writer(snpfil,delimiter = "\t")
        for record in vcf_reader:
            if record.is_snp:
                dp_l = [sample['DP'] for sample in record.samples]
                if any([d is None for d in dp_l]) == True:
                    #one sample has no depth for SNP: Skip it!
                    #print record, dp_l
                    continue
                snpfilwriter.writerow([snp,record.CHROM,
                                       record.POS,record.REF])
                for sample in vcf_reader.samples:
                    data = record.genotype(sample).data
                
                    if sample in snps.keys():
                        snps[sample][str(snp)]=[data.DP,
                                                        '2',
                                                        data.AD,
                                                        data.RD]
                    else:
                        snps[sample]={}
                        snps[sample][str(snp)]=[data.DP,
                                                        '2',
                                                        data.AD,
                                                        data.RD]
                snp=snp+1
                loci=loci+1
        loci_t=loci_t+str(loci)
        snpfil.close()
        return [snps,loci_t,pops_t,loci,pops]
    else:
        loci=0
        loci_t = '[loci]='
        pops = len(vcf_reader.samples)
        pops_t = '[populations]='+str(pops)
        snps = {}
        #for contig in vcf_reader.contigs.keys():
        snps={}
        snp = 1
        #also write snp id and snp position to separate file
        if outdir == 'none':
            snpfildir = os.getcwd()
        else:
            snpfildir = outdir
        snpfil = open(snpfildir+'/'+prefix+'_'+cont+'_snps.txt','w')
        snpfilwriter = csv.writer(snpfil,delimiter = "\t")
        for record in vcf_reader.fetch(cont):
            if record.is_snp:
                dp_l = [sample['DP'] for sample in record.samples]
                if any([d is None for d in dp_l]) == True:
                    #one sample has no depth for SNP: Skip it!
                    #print record, dp_l
                    continue
                snpfilwriter.writerow([snp,record.CHROM,
                                       record.POS,record.REF])
                for sample in vcf_reader.samples:
                    data = record.genotype(sample).data
                
                    if sample in snps.keys():
                        snps[sample][str(snp)]=[data.DP,
                                                        '2',
                                                        data.AD,
                                                        data.RD]
                    else:
                        snps[sample]={}
                        snps[sample][str(snp)]=[data.DP,
                                                        '2',
                                                        data.AD,
                                                        data.RD]
                snp=snp+1
                loci=loci+1
        loci_t=loci_t+str(loci)
        snpfil.close()
        return [snps,loci_t,pops_t,loci,pops]

def printBscan(snps_data,cont,prefix,outdir):
    pops = 1
    #also write pop id and pop nr to separate file:
    if outdir == 'none':
	outfildir = os.getcwd()
    else:
	outfildir = outdir
    popfil = open(outfildir+'/'+prefix+'_'+cont+'_pops.txt','w')
    bscanfil = open(outfildir+'/'+prefix+'_'+cont+'_bscan.in','w')
    bscanfilwriter = csv.writer(bscanfil,delimiter = " ")
    popfilwriter = csv.writer(popfil,delimiter = "\t")
    #print snps_data[1]
    bscanfilwriter.writerow([snps_data[1]])
    #print snps_data[2]
    bscanfilwriter.writerow([snps_data[2]])
    for pop in snps_data[0].keys():
        #print '[pop]='+str(pops)
        bscanfilwriter.writerow(['[pop]='+str(pops)])
        popfilwriter.writerow([pops,pop])
	for snp in range(1,snps_data[3]+1):
		printlist = [str(snp)]
		for dat in snps_data[0][pop][str(snp)]:
			printlist.append(str(dat))
		bscanfilwriter.writerow(printlist)
	pops=pops+1
    popfil.close()
    bscanfil.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=globals()['__doc__'])

    #MAIN CONTROLS
    parser.add_argument('-filename',
                       help = 'tabix indexed .vcf file name')

    #FILTER CONTROLS
    parser.add_argument('-qual', default = 30,
                        help = 'Minimum SNP Phred score (default: 30).')

    parser.add_argument('-depth', default = 5,
                        help = 'Minimum SNP depth (default: 5).')

    parser.add_argument('-task', default = 'd',
                        help = 'Which task to perform:\n'+\
                        'd = snp density\n'+\
                        'dp = print average snp DP across samples\n'+\
                        'bscan = print bscan output')
    
    parser.add_argument('-prefix', default = 'bscan',
                        help = 'Prefix for snp and pop files.')

    parser.add_argument('-outdir', default = 'none',
                        help = 'Output directory for snp and pop files.')
    
    parser.add_argument('-region',nargs='*',
                        help = 'Contig/Chromosome id for the regions'+\
                        ' to consider')

    parser.add_argument('-window', default = 1000,
                        help = 'Window length for window analyses (default 1000).')

    args = vars(parser.parse_args())

    try:
        if args.has_key('region'):
            #print type(args['region']) #script tester line
            if type(args['region']) == list:
                for cont in args['region']:
                    #print cont #script tester line
                    if args['task'] == 'd':
                        vcf_reader = vcf.Reader(open(args['filename'],'r'))
                        SNPdensity(vcf_reader,cont,int(args['window']),int(args['depth']))
                    if args['task'] == 'dp':
                        vcf_reader = vcf.Reader(open(args['filename'],'r'))
                        SNPprinter(vcf_reader,cont,args['qual'])
                    if args['task'] == 'bscan':
                        vcf_reader = vcf.Reader(open(args['filename'],'r'))
                        snps_data = Bscandict(vcf_reader,cont,args['prefix'],args['outdir'])
                        printBscan(snps_data,cont,args['prefix'],args['outdir'])
            elif type(args['region']) != str:
                cont = str(args['region'])
                if args['task'] == 'd':
                    vcf_reader = vcf.Reader(open(args['filename'],'r'))
                    SNPdensity(vcf_reader,cont,int(args['window']),int(args['depth']))
                if args['task'] == 'dp':
                    vcf_reader = vcf.Reader(open(args['filename'],'r'))
                    SNPprinter(vcf_reader,cont,args['qual'])
                if args['task'] == 'bscan':
                    vcf_reader = vcf.Reader(open(args['filename'],'r'))
                    snps_data = Bscandict(vcf_reader,cont,args['prefix'],args['outdir'])
                    printBscan(snps_data,cont,args['prefix'],args['outdir'])
            else:
                cont = args['region']
                if args['task'] == 'd':
                    vcf_reader = vcf.Reader(open(args['filename'],'r'))
                    SNPdensity(vcf_reader,cont,int(args['window']),int(args['depth']))
                if args['task'] == 'dp':
                    vcf_reader = vcf.Reader(open(args['filename'],'r'))
                    SNPprinter(vcf_reader,cont,args['qual'])
                if args['task'] == 'bscan':
                    vcf_reader = vcf.Reader(open(args['filename'],'r'))
                    snps_data = Bscandict(vcf_reader,cont,args['prefix'],args['outdir'])
                    printBscan(snps_data,cont,args['prefix'],args['outdir'])

                      
        else:
            sys.stderr.write("Warning: No region given, performing task for"+\
                             "entire file.")
            cont = "all"
            if args['task'] == 'd':
                vcf_reader = vcf.Reader(open(args['filename'],'r'))
                SNPdensity(vcf_reader,cont,int(args['window']),int(args['depth']))
            if args['task'] == 'dp':
                vcf_reader = vcf.Reader(open(args['filename'],'r'))
                SNPprinter(vcf_reader,cont,args['qual'])
            if args['task'] == 'bscan':
                vcf_reader = vcf.Reader(open(args['filename'],'r'))          
                snps_data = Bscandict(vcf_reader,cont,args['prefix'],args['outdir'])
                printBscan(snps_data,cont,args['prefix'],args['outdir'])

            
    except IOError:
        sys.stderr.write("stdout is closed")
        #stdout is closed, no point continuing
        #close explicitly to prevent cleanup problems
        try: sys.stdout.close()
        except IOError: pass
        try: sys.stdout.close()
        except IOError: pass
