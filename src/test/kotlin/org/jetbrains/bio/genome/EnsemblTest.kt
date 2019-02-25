package org.jetbrains.bio.genome

import org.junit.Test
import kotlin.test.assertEquals

/**
 * @author Roman.Chernyatchik
 */
class GtfReaderTest {
    val ENST00000541371 = """# ENST00000541371, hg19: -chr11:60666013-60674041
# [no stop codonl no UTR3]
1	havana	transcript	60666013	60674041	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	60673836	60674041	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "1"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; exon_id "ENSE00002273630"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	60673836	60673854	.	-	0	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "1"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; protein_id "ENSP00000440266"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	start_codon	60673852	60673854	.	-	0	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "1"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	60671184	60671333	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "2"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; exon_id "ENSE00003588469"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	60671184	60671333	.	-	2	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "2"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; protein_id "ENSP00000440266"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	60670931	60671007	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "3"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; exon_id "ENSE00003654388"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	60670931	60671007	.	-	2	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "3"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; protein_id "ENSP00000440266"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	60670212	60670353	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "4"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; exon_id "ENSE00003585619"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	60670212	60670353	.	-	0	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "4"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; protein_id "ENSP00000440266"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	60670055	60670128	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "5"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; exon_id "ENSE00000721018"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	60670055	60670128	.	-	2	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "5"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; protein_id "ENSP00000440266"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	60669875	60669937	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "6"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; exon_id "ENSE00001167760"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	60669875	60669937	.	-	0	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "6"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; protein_id "ENSP00000440266"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	60668971	60669012	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "7"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; exon_id "ENSE00003599721"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	60668971	60669012	.	-	0	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "7"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; protein_id "ENSP00000440266"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	60668767	60668841	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "8"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; exon_id "ENSE00001167737"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	60668767	60668841	.	-	0	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "8"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; protein_id "ENSP00000440266"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	60668326	60668401	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "9"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; exon_id "ENSE00001167731"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	60668326	60668401	.	-	0	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "9"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; protein_id "ENSP00000440266"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	60666013	60666050	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "10"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; exon_id "ENSE00002229042"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	60666013	60666050	.	-	2	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; exon_number "10"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; protein_id "ENSP00000440266"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	five_prime_utr	60673855	60674041	.	-	.	gene_id "ENSG00000110107"; gene_version "4"; transcript_id "ENST00000541371"; transcript_version "1"; gene_name "PRPF19"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "PRPF19-003"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000396339"; havana_transcript_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
"""

    val ENST00000541351 = """# ENST00000541351, hg19: -chr12:10853874-10863313
# [no start codon; several UTR3; no UTR5]
1	havana	transcript	10853874	10863313	.	-	.	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	exon	10863253	10863313	.	-	.	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "1"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; exon_id "ENSE00002310575"; exon_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	CDS	10863253	10863313	.	-	1	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "1"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; protein_id "ENSP00000441447"; protein_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	exon	10862507	10862713	.	-	.	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "2"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; exon_id "ENSE00003601328"; exon_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	CDS	10862507	10862713	.	-	0	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "2"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; protein_id "ENSP00000441447"; protein_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	exon	10856650	10856747	.	-	.	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "3"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; exon_id "ENSE00003572809"; exon_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	CDS	10856650	10856747	.	-	0	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "3"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; protein_id "ENSP00000441447"; protein_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	exon	10856122	10856218	.	-	.	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "4"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; exon_id "ENSE00002300734"; exon_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	CDS	10856122	10856218	.	-	1	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "4"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; protein_id "ENSP00000441447"; protein_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	exon	10854559	10854733	.	-	.	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "5"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; exon_id "ENSE00003668113"; exon_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	CDS	10854689	10854733	.	-	0	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "5"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; protein_id "ENSP00000441447"; protein_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	stop_codon	10854686	10854688	.	-	0	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "5"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	exon	10853874	10853952	.	-	.	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; exon_number "6"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; exon_id "ENSE00002220183"; exon_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	three_prime_utr	10854559	10854685	.	-	.	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
1	havana	three_prime_utr	10853874	10853952	.	-	.	gene_id "ENSG00000060138"; gene_version "8"; transcript_id "ENST00000541351"; transcript_version "1"; gene_name "YBX3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "YBX3-014"; transcript_source "havana"; transcript_biotype "nonsense_mediated_decay"; havana_transcript "OTTHUMT00000399636"; havana_transcript_version "1"; tag "cds_start_NF"; tag "mRNA_start_NF";
"""

    val ENST00000311042="""# ENST00000311042, hg19: +chr3:39509171-39557492
# [several UTR5; STOP codon split in 2 exons (several stop codon records)]
1	havana	transcript	39509171	39557492	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; tag "basic";
1	havana	exon	39509171	39509231	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; exon_number "1"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; exon_id "ENSE00001773788"; exon_version "1"; tag "basic";
1	havana	exon	39521531	39521614	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; exon_number "2"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; exon_id "ENSE00001738261"; exon_version "1"; tag "basic";
1	havana	exon	39543557	39543766	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; exon_number "3"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; exon_id "ENSE00001735986"; exon_version "1"; tag "basic";
1	havana	CDS	39543561	39543766	.	+	0	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; exon_number "3"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; protein_id "ENSP00000312293"; protein_version "6"; tag "basic";
1	havana	start_codon	39543561	39543563	.	+	0	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; exon_number "3"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; tag "basic";
1	havana	exon	39543954	39544367	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; exon_number "4"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; exon_id "ENSE00001717989"; exon_version "1"; tag "basic";
1	havana	CDS	39543954	39544365	.	+	1	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; exon_number "4"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; protein_id "ENSP00000312293"; protein_version "6"; tag "basic";
1	havana	exon	39554874	39557492	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; exon_number "5"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; exon_id "ENSE00001759995"; exon_version "1"; tag "basic";
1	havana	stop_codon	39544366	39544367	.	+	0	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; exon_number "4"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; tag "basic";
1	havana	stop_codon	39554874	39554874	.	+	1	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; exon_number "5"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; tag "basic";
1	havana	five_prime_utr	39509171	39509231	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; tag "basic";
1	havana	five_prime_utr	39521531	39521614	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; tag "basic";
1	havana	five_prime_utr	39543557	39543560	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; tag "basic";
1	havana	three_prime_utr	39554875	39557492	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000311042"; transcript_version "6"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-015"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS63598"; havana_transcript "OTTHUMT00000343653"; havana_transcript_version "1"; tag "basic";
"""

    val ENST00000441980 = """# ENST00000441980 hg19: +chr3:39509171-39544518
# several UTR5, one UTR3, start & stop codons
1	havana	transcript	39509171	39544518	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; tag "basic";
1	havana	exon	39509171	39509231	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; exon_number "1"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; exon_id "ENSE00001773788"; exon_version "1"; tag "basic";
1	havana	exon	39540958	39541023	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; exon_number "2"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; exon_id "ENSE00001654736"; exon_version "1"; tag "basic";
1	havana	exon	39543557	39543766	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; exon_number "3"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; exon_id "ENSE00001735986"; exon_version "1"; tag "basic";
1	havana	CDS	39543561	39543766	.	+	0	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; exon_number "3"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; protein_id "ENSP00000388827"; protein_version "2"; tag "basic";
1	havana	start_codon	39543561	39543563	.	+	0	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; exon_number "3"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; tag "basic";
1	havana	exon	39544026	39544518	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; exon_number "4"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; exon_id "ENSE00001693771"; exon_version "2"; tag "basic";
1	havana	CDS	39544026	39544368	.	+	1	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; exon_number "4"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; protein_id "ENSP00000388827"; protein_version "2"; tag "basic";
1	havana	stop_codon	39544369	39544371	.	+	0	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; exon_number "4"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; tag "basic";
1	havana	five_prime_utr	39509171	39509231	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; tag "basic";
1	havana	five_prime_utr	39540958	39541023	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; tag "basic";
1	havana	five_prime_utr	39543557	39543560	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; tag "basic";
1	havana	three_prime_utr	39544372	39544518	.	+	.	gene_id "ENSG00000168314"; gene_version "13"; transcript_id "ENST00000441980"; transcript_version "2"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000343719"; havana_transcript_version "2"; tag "basic";
"""

    val ENST00000441980_legacy = """# ENST00000441980 (v.75), hg19, +chr3:39509171-39544518
# start, stop codons
1	protein_coding	transcript	39509171	39544518	.	+	.	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana";
1	protein_coding	exon	39509171	39509231	.	+	.	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; exon_number "1"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; exon_id "ENSE00001773788";
1	protein_coding	exon	39540958	39541023	.	+	.	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; exon_number "2"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; exon_id "ENSE00001654736";
1	protein_coding	exon	39543557	39543766	.	+	.	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; exon_number "3"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; exon_id "ENSE00001735986";
1	protein_coding	CDS	39543561	39543766	.	+	0	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; exon_number "3"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; protein_id "ENSP00000388827";
1	protein_coding	start_codon	39543561	39543563	.	+	0	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; exon_number "3"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana";
1	protein_coding	exon	39544026	39544518	.	+	.	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; exon_number "4"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; exon_id "ENSE00001693771";
1	protein_coding	CDS	39544026	39544368	.	+	1	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; exon_number "4"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana"; protein_id "ENSP00000388827";
1	protein_coding	stop_codon	39544369	39544371	.	+	0	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; exon_number "4"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana";
1	protein_coding	UTR	39509171	39509231	.	+	.	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana";
1	protein_coding	UTR	39540958	39541023	.	+	.	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana";
1	protein_coding	UTR	39543557	39543560	.	+	.	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana";
1	protein_coding	UTR	39544372	39544518	.	+	.	gene_id "ENSG00000168314"; transcript_id "ENST00000441980"; gene_name "MOBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "MOBP-012"; transcript_source "havana";
"""

    val ENSG00000215788 = """# ENSG00000215788, hg19 : -chr1:6522697-6525690
# split START codon (several start_codons records), CDS of size 1
1	havana	transcript	6522697	6525690	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	6525500	6525690	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "1"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; exon_id "ENSE00001823169"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	6525148	6525282	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "2"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; exon_id "ENSE00003582174"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	6524612	6524779	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "3"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; exon_id "ENSE00003497646"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	6524435	6524513	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "4"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; exon_id "ENSE00003529010"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	6523132	6523187	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "5"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; exon_id "ENSE00003657231"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	6523132	6523132	.	-	0	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "5"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; protein_id "ENSP00000465381"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	start_codon	6523132	6523132	.	-	0	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "5"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	start_codon	6523029	6523030	.	-	2	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "6"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	6522923	6523030	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "6"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; exon_id "ENSE00003684931"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	6522923	6523030	.	-	2	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "6"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; protein_id "ENSP00000465381"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	exon	6522697	6522723	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "7"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; exon_id "ENSE00001831041"; exon_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	CDS	6522697	6522723	.	-	2	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; exon_number "7"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; protein_id "ENSP00000465381"; protein_version "1"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	five_prime_utr	6525500	6525690	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	five_prime_utr	6525148	6525282	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	five_prime_utr	6524612	6524779	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	five_prime_utr	6524435	6524513	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; tag "cds_end_NF"; tag "mRNA_end_NF";
1	havana	five_prime_utr	6523133	6523187	.	-	.	gene_id "ENSG00000215788"; gene_version "5"; transcript_id "ENST00000481401"; transcript_version "1"; gene_name "TNFRSF25"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TNFRSF25-012"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000099349"; havana_transcript_version "2"; tag "cds_end_NF"; tag "mRNA_end_NF";
"""

    val ENST00000589266 = """# ENST00000589266, hg19: +chr19:1266802-1272995
# No stop codon, but [CDS.3end + 1 .. CDS.3end + 3] is part of UTR3 although supposed to be 3 stop codon nucleotides
1	havana	transcript	1266802	1272995	.	+	.	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; tag "mRNA_end_NF";
1	havana	exon	1266802	1267164	.	+	.	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; exon_number "1"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; exon_id "ENSE00002972777"; exon_version "1"; tag "mRNA_end_NF";
1	havana	exon	1270927	1271035	.	+	.	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; exon_number "2"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; exon_id "ENSE00003501315"; exon_version "1"; tag "mRNA_end_NF";
1	havana	CDS	1270933	1271035	.	+	0	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; exon_number "2"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; protein_id "ENSP00000467138"; protein_version "1"; tag "mRNA_end_NF";
1	havana	start_codon	1270933	1270935	.	+	0	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; exon_number "2"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; tag "mRNA_end_NF";
1	havana	exon	1271139	1271245	.	+	.	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; exon_number "3"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; exon_id "ENSE00003619232"; exon_version "1"; tag "mRNA_end_NF";
1	havana	CDS	1271139	1271245	.	+	2	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; exon_number "3"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; protein_id "ENSP00000467138"; protein_version "1"; tag "mRNA_end_NF";
1	havana	exon	1271328	1271360	.	+	.	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; exon_number "4"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; exon_id "ENSE00002852443"; exon_version "1"; tag "mRNA_end_NF";
1	havana	CDS	1271328	1271360	.	+	0	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; exon_number "4"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; protein_id "ENSP00000467138"; protein_version "1"; tag "mRNA_end_NF";
1	havana	exon	1272968	1272995	.	+	.	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; exon_number "5"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; exon_id "ENSE00002944318"; exon_version "1"; tag "mRNA_end_NF";
1	havana	CDS	1272968	1272992	.	+	0	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; exon_number "5"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; protein_id "ENSP00000467138"; protein_version "1"; tag "mRNA_end_NF";
1	havana	five_prime_utr	1266802	1267164	.	+	.	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; tag "mRNA_end_NF";
1	havana	five_prime_utr	1270927	1270932	.	+	.	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; tag "mRNA_end_NF";
1	havana	three_prime_utr	1272993	1272995	.	+	.	gene_id "ENSG00000099622"; gene_version "9"; transcript_id "ENST00000589266"; transcript_version "1"; gene_name "CIRBP"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "CIRBP-017"; transcript_source "havana"; transcript_biotype "protein_coding"; havana_transcript "OTTHUMT00000449976"; havana_transcript_version "2"; tag "mRNA_end_NF";
"""

    val ENST00000457222 = """# ENST00000457222, hg19: +chrY:9236030-9238826
# Stop codon is part of CDS
1	ensembl_havana	gene	9236030	9307357	.	+	.	gene_id "ENSG00000228927"; gene_version "5"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding";
1	ensembl_havana	transcript	9236030	9238826	.	+	.	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; tag "basic";
1	ensembl_havana	exon	9236030	9236561	.	+	.	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "1"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; exon_id "ENSE00001774191"; exon_version "1"; tag "basic";
1	ensembl_havana	CDS	9236076	9236561	.	+	0	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "1"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; protein_id "ENSP00000398163"; protein_version "2"; tag "basic";
1	ensembl_havana	start_codon	9236076	9236078	.	+	0	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "1"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; tag "basic";
1	ensembl_havana	exon	9237169	9237246	.	+	.	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "2"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; exon_id "ENSE00003692666"; exon_version "1"; tag "basic";
1	ensembl_havana	CDS	9237169	9237246	.	+	0	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "2"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; protein_id "ENSP00000398163"; protein_version "2"; tag "basic";
1	ensembl_havana	exon	9237375	9237486	.	+	.	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "3"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; exon_id "ENSE00001747539"; exon_version "1"; tag "basic";
1	ensembl_havana	CDS	9237375	9237486	.	+	0	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "3"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; protein_id "ENSP00000398163"; protein_version "2"; tag "basic";
1	ensembl_havana	exon	9237588	9237733	.	+	.	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "4"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; exon_id "ENSE00001719444"; exon_version "1"; tag "basic";
1	ensembl_havana	CDS	9237588	9237733	.	+	2	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "4"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; protein_id "ENSP00000398163"; protein_version "2"; tag "basic";
1	ensembl_havana	exon	9237840	9237921	.	+	.	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "5"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; exon_id "ENSE00001781353"; exon_version "1"; tag "basic";
1	ensembl_havana	CDS	9237840	9237921	.	+	0	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "5"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; protein_id "ENSP00000398163"; protein_version "2"; tag "basic";
1	ensembl_havana	exon	9238616	9238826	.	+	.	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "6"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; exon_id "ENSE00001954357"; exon_version "1"; tag "basic";
1	ensembl_havana	CDS	9238616	9238638	.	+	2	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "6"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; protein_id "ENSP00000398163"; protein_version "2"; tag "basic";
1	ensembl_havana	stop_codon	9238636	9238638	.	+	0	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; exon_number "6"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; tag "basic";
1	ensembl_havana	five_prime_utr	9236030	9236075	.	+	.	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; tag "basic";
1	ensembl_havana	three_prime_utr	9238639	9238826	.	+	.	gene_id "ENSG00000228927"; gene_version "5"; transcript_id "ENST00000457222"; transcript_version "2"; gene_name "TSPY3"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "TSPY3-001"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS48204"; havana_transcript "OTTHUMT00000099505"; havana_transcript_version "1"; tag "basic";
"""

    val ENST00000423104 =  """# ENST00000423104, hg19, -chr2:204306801-204360076
# short 3bp exon with stop codon
2	protein_coding	transcript	204306801	204360076	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana";
2	protein_coding	exon	204359957	204360076	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "1"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00003567360";
2	protein_coding	CDS	204359957	204360076	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "1"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	start_codon	204360074	204360076	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "1"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana";
2	protein_coding	exon	204355937	204356042	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "2"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001199380";
2	protein_coding	CDS	204355937	204356042	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "2"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204354307	204354812	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "3"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001199377";
2	protein_coding	CDS	204354307	204354812	.	-	2	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "3"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204341799	204341879	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "4"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001199327";
2	protein_coding	CDS	204341799	204341879	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "4"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204326571	204326648	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "5"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001199371";
2	protein_coding	CDS	204326571	204326648	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "5"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204325972	204326131	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "6"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001199366";
2	protein_coding	CDS	204325972	204326131	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "6"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204324630	204324751	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "7"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001199364";
2	protein_coding	CDS	204324630	204324751	.	-	2	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "7"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204322253	204322318	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "8"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001199362";
2	protein_coding	CDS	204322253	204322318	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "8"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204320160	204320303	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "9"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001199357";
2	protein_coding	CDS	204320160	204320303	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "9"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204319153	204319263	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "10"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001199353";
2	protein_coding	CDS	204319153	204319263	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "10"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204313461	204313559	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "11"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001199349";
2	protein_coding	CDS	204313461	204313559	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "11"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204312682	204312802	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "12"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001199342";
2	protein_coding	CDS	204312682	204312802	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "12"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204309591	204309733	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "13"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00001344391";
2	protein_coding	CDS	204309591	204309733	.	-	2	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "13"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; protein_id "ENSP00000397751";
2	protein_coding	exon	204306801	204306803	.	-	.	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "14"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana"; exon_id "ENSE00002483457";
2	protein_coding	stop_codon	204306801	204306803	.	-	0	gene_id "ENSG00000173166"; transcript_id "ENST00000423104"; exon_number "14"; gene_name "RAPH1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "RAPH1-015"; transcript_source "havana";
"""

    val ENST00000402757_legacy = """# ENST00000402757, hg18, -NT_113923:82145-100313
1	protein_coding	exon	100125	100313	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "1"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	exon	99411	99581	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "2"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	99411	99574	.	-	0	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "2"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	start_codon	99572	99574	.	-	0	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "2"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	exon	96180	96310	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "3"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	96180	96310	.	-	1	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "3"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	exon	95959	96045	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "4"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	95959	96045	.	-	2	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "4"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	exon	90645	90837	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "5"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	90645	90837	.	-	2	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "5"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	exon	90077	90233	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "6"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	90077	90233	.	-	1	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "6"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	exon	87896	88049	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "7"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	87896	88049	.	-	0	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "7"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	exon	87235	87371	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "8"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	87235	87371	.	-	2	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "8"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	exon	83497	83684	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "9"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	83497	83684	.	-	0	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "9"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	exon	83216	83264	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "10"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	83216	83264	.	-	1	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "10"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	exon	83138	83213	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "11"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	83138	83213	.	-	0	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "11"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	exon	82924	83042	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "12"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	82924	83042	.	-	2	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "12"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	exon	82729	82842	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "13"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	82729	82842	.	-	0	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "13"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	exon	82145	82424	.	-	.	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "14"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
1	protein_coding	CDS	82281	82424	.	-	0	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "14"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201"; protein_id "ENSP00000385238";
1	protein_coding	stop_codon	82278	82280	.	-	0	 gene_id "ENSG00000216522"; transcript_id "ENST00000402757"; exon_number "14"; gene_name "AL356585.7-5"; transcript_name "AL356585.7-201";
"""

    val ENSMUST00000168714 = """# ENSMUST00000168714, mm9, +chr11:120534397-120570386
# 2 bp exons
1	protein_coding	exon	120534397	120534484	.	+	.	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "1"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018";
1	protein_coding	exon	120535354	120535409	.	+	.	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "2"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018";
1	protein_coding	exon	120538825	120538939	.	+	.	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "3"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018";
1	protein_coding	CDS	120538898	120538939	.	+	0	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "3"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018"; protein_id "ENSMUSP00000129462";
1	protein_coding	start_codon	120538898	120538900	.	+	0	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "3"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018";
1	protein_coding	exon	120539734	120539834	.	+	.	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "4"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018";
1	protein_coding	CDS	120539734	120539834	.	+	0	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "4"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018"; protein_id "ENSMUSP00000129462";
1	protein_coding	exon	120547958	120548015	.	+	.	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "5"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018";
1	protein_coding	CDS	120547958	120548015	.	+	1	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "5"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018"; protein_id "ENSMUSP00000129462";
1	protein_coding	exon	120569972	120570023	.	+	.	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "6"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018";
1	protein_coding	CDS	120569972	120570023	.	+	0	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "6"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018"; protein_id "ENSMUSP00000129462";
1	protein_coding	exon	120570244	120570313	.	+	.	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "7"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018";
1	protein_coding	CDS	120570244	120570313	.	+	2	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "7"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018"; protein_id "ENSMUSP00000129462";
1	protein_coding	exon	120570385	120570386	.	+	.	 gene_id "ENSMUSG00000025142"; transcript_id "ENSMUST00000168714"; exon_number "8"; gene_name "Aspscr1"; gene_biotype "protein_coding"; transcript_name "Aspscr1-018";
"""

    @Test
    fun noStartCodon() {
        val transcript = readTranscript(ENST00000541351)
        assertEquals(Location(10853873, 10863313, Chromosome("to1", "chr1"),
                              Strand.MINUS),
                     transcript.location)
        assertEquals(Range(10854688, 10863313), transcript.cdsRange)
    }

    @Test
    fun noUTR5() {
        val transcript = readTranscript(ENST00000541351)
        assertEquals(0, transcript.utr5.size)
    }

    @Test
    fun noUTR3() {
        val transcript = readTranscript(ENST00000541371)
        assertEquals(0, transcript.utr3.size)
    }

    @Test
    fun singleUTR5() {
        val transcript = readTranscript(ENST00000541371)
        assertEquals(1, transcript.utr5.size)
        assertEquals(Range(60673854, 60674041), transcript.utr5.map(Location::toRange).first())
    }

    @Test
    fun singleUTR3() {
        val transcript = readTranscript(ENST00000441980)
        assertEquals(1, transcript.utr3.size)
        assertEquals(Range(39544371, 39544518), transcript.utr3.map(Location::toRange).first())
    }

    @Test
    fun splitStartCodon() {
        val transcript = readTranscript(ENSG00000215788)

        assertEquals(Range(6522696, 6523132), transcript.cdsRange)
        assertEquals(Range(6523131, 6523132), transcript.cds.map(Location::toRange).last())

        assertEquals(5, transcript.utr5.size)
        assertEquals(listOf(Range(6523132, 6523187), Range(6524434, 6524513), Range(6524611, 6524779),
                            Range(6525147, 6525282), Range(6525499, 6525690)),
                     transcript.utr5.map(Location::toRange))
    }

    @Test
    fun splitStopCodon() {
        val transcript = readTranscript(ENST00000311042)

        assertEquals(Range(39543953, 39544365), transcript.cds.map(Location::toRange).last())
        assertEquals(Range(39543560, 39544365), transcript.cdsRange)

        assertEquals(1, transcript.utr3.size)
        assertEquals(listOf(Range(39554874, 39557492)), transcript.utr3.map(Location::toRange))
    }

    @Test
    fun stopCodonInNextExon() {
        val transcript = readTranscript(ENST00000423104)

        assertEquals(Range(204306800, 204306803), transcript.exons.sorted()[0].toRange())
        assertEquals(Range(204309590, 204309733), transcript.exons.sorted()[1].toRange())
        assertEquals(Range(204309590, 204309733), transcript.cds.map(Location::toRange).first())
        assertEquals(204309590, transcript.cdsRange!!.startOffset)
        // only stop codon in last exon, no space for UTR3"
        assertEquals(listOf(), transcript.utr3.map(Location::toRange))
    }

    @Test
    fun multipleUTR5() {
        val transcript = readTranscript(ENST00000441980)
        assertEquals(3, transcript.utr5.size)
        assertEquals(listOf(Range(39509170, 39509231), Range(39540957, 39541023), Range(39543556,39543560)),
                     transcript.utr5.map(Location::toRange))
    }

    @Test
    fun multipleUTR3() {
        val transcript = readTranscript(ENST00000541351)

        assertEquals(2, transcript.utr3.size)
        assertEquals(listOf(Range(10853873, 10853952), Range(10854558, 10854685)),
                     transcript.utr3.map(Location::toRange))
    }


    @Test
    fun legacyFormat1() {
        // With 'UTR' instead of 'three_prime_utr', 'five_prime_utr'
        val transcript = readTranscript(ENST00000441980_legacy)

        // Transcript
        assertEquals(Range(39509170, 39544518), transcript.location.toRange())

        // Exons
        assertEquals(4, transcript.exons.size)
        assertEquals(listOf(Range(39509170, 39509231), Range(39540957, 39541023), Range(39543556, 39543766),
                            Range(39544025, 39544518)),
                     transcript.exons.map { it.toRange() })

        // CDS
        assertEquals(2, transcript.cds.size)
        assertEquals(listOf(Range(39543560, 39543766), Range(39544025, 39544368)), transcript.cds.map(Location::toRange))

        // CDS Bounds
        assertEquals(Range(39543560, 39544368), transcript.cdsRange)

        // UTR5
        assertEquals(3, transcript.utr5.size)
        assertEquals(listOf(Range(39509170, 39509231), Range(39540957, 39541023), Range(39543556, 39543560)),
                     transcript.utr5.map(Location::toRange))

        // UTR3
        assertEquals(1, transcript.utr3.size)
        assertEquals(listOf(Range(39544371, 39544518)), transcript.utr3.map(Location::toRange))
    }

    @Test
    fun legacyFormat2() {
        // Without 'transcript' and utrs
        val transcript = readTranscript(ENST00000402757_legacy)

        // Transcript
        assertEquals(Range(82144, 100313), transcript.location.toRange())

        // Exons
        assertEquals(14, transcript.exons.size)
        val exons = transcript.exons.map { it.toRange() }.sorted()
        assertEquals(Range(82144, 82424), exons.first())
        assertEquals(Range(83496, 83684), exons[5])
        assertEquals(Range(100124, 100313), exons.last())

        // CDS
        assertEquals(13, transcript.cds.size)
        assertEquals(Range(82280, 82424), transcript.cds.map(Location::toRange).first())
        assertEquals(Range(83496, 83684), transcript.cds.map(Location::toRange)[5])
        assertEquals(Range(99410, 99574), transcript.cds.map(Location::toRange).last())

        // CDS Bounds
        assertEquals(Range(82280, 99574), transcript.cdsRange)

        // UTR5
        assertEquals(2, transcript.utr5.size)
        assertEquals(listOf(Range(99574, 99581), Range(100124, 100313)),
                     transcript.utr5.map(Location::toRange))

        // UTR3
        assertEquals(1, transcript.utr3.size)
        assertEquals(listOf(Range(82144, 82277)), transcript.utr3.map(Location::toRange))
    }

    @Test
    fun shortCDS() {
        val transcript = readTranscript(ENSG00000215788)
        // cds
        assertEquals(3, transcript.cds.size)
        assertEquals(1, transcript.cds.map(Location::toRange).last().length())
        assertEquals(Range(6523131, 6523132), transcript.cds.map(Location::toRange).last())
    }

    @Test
    fun shortExon() {
        val transcript = readTranscript(ENSMUST00000168714)

        // utr3
        assertEquals(listOf(), transcript.utr3)
//        assertEquals(3, transcript.cds.size)
//        assertEquals(1, transcript.cds.map(Location::toRange).last().length())
//        assertEquals(Range(6523131, 6523132), transcript.cds.map(Location::toRange).last())
    }

    @Test
    fun multipleCDS() {
        val transcript = readTranscript(ENST00000541371)
        // cds
        assertEquals(10, transcript.cds.size)
        assertEquals(Range(60666012, 60666050), transcript.cds.map(Location::toRange).first())
        assertEquals(Range(60670054, 60670128), transcript.cds.map(Location::toRange)[5])
        assertEquals(Range(60673835, 60673854), transcript.cds.map(Location::toRange).last())
    }
    
    @Test
    fun multipleExons() {
        val transcript = readTranscript(ENST00000541371)
        // exons
        assertEquals(10, transcript.exons.size)
        assertEquals(Range(60666012, 60666050), transcript.exons.first().toRange())
        assertEquals(Range(60669874, 60669937), transcript.exons[4].toRange())
        assertEquals(Range(60673835, 60674041), transcript.exons.last().toRange())
    }
    
    @Test
    fun startCodon() {
        val transcript = readTranscript(ENST00000541371)
        assertEquals(Range(60666012, 60673854), transcript.cdsRange)
    }

    @Test
    fun stopCodon() {
        val transcript = readTranscript(ENST00000541351)
        assertEquals(10854688, transcript.cdsRange!!.startOffset)
    }

    @Test
    fun noStopCodon() {
        val transcript = readTranscript(ENST00000541371)

        assertEquals(Range(60666012, 60674041), transcript.location.toRange())
        assertEquals(Range(60666012, 60666050), transcript.cds.map(Location::toRange).first())
        assertEquals(Range(60666012, 60673854), transcript.cdsRange)
        assertEquals(0, transcript.utr3.size)
    }

    @Test
    fun noStopCodonWithUTR3() {
        val transcript = readTranscript(ENST00000589266)

        assertEquals(Range(1266801, 1272995), transcript.location.toRange())
        assertEquals(Range(1272967, 1272992), transcript.cds.map(Location::toRange).last())
        assertEquals(Range(1270932, 1272992), transcript.cdsRange)

        // here stop codon should be [1272992, 1272994] => actually no UTR3, but
        // according to GTF:
        assertEquals(listOf(Range(1272992, 1272995)), transcript.utr3.map(Location::toRange))
    }

    @Test
    fun stopCodonInCDS() {
        val transcript = readTranscript(ENST00000457222)

        assertEquals(Range(9236029, 9238826), transcript.location.toRange())
        assertEquals(Range(9238615, 9238826), transcript.exons.last().toRange())
        assertEquals(Range(9238615, 9238638), transcript.cds.map(Location::toRange).last())
        assertEquals(9238638, transcript.cdsRange!!.endOffset)

        // Stop codon should be [9238638..9238640] => UTR3 is [9238640, 9238826),
        // But in GTF 'stop_codon' [9238635..9238637] is in CDS, and UTR3 is:
        assertEquals(listOf(Range(9238638, 9238826)), transcript.utr3.map(Location::toRange))
    }

    @Test
    fun determineUTR3End5() {
        val chr = Chromosome("to1", "chr1")

        // basic:
        assertEquals(83, GtfReader.determineUTR3End5(
                Location(20, 80, chr, Strand.PLUS),
                listOf(Range(10, 90)),
                "ENCT"
        ))

        assertEquals(16, GtfReader.determineUTR3End5(
                Location(20, 80, chr, Strand.MINUS),
                listOf(Range(10, 90)),
                "ENCT"
        ))

        // Next exon:
        assertEquals(93, GtfReader.determineUTR3End5(
                Location(20, 80, chr, Strand.PLUS),
                listOf(Range(10, 80), Range(90, 94)),
                "ENCT"
        ))
        assertEquals(10, GtfReader.determineUTR3End5(
                Location(20, 80, chr, Strand.MINUS),
                listOf(Range(10, 14), Range(20, 90)),
                "ENCT"
        ))

        // Split:
        assertEquals(92, GtfReader.determineUTR3End5(
                Location(20, 80, chr, Strand.PLUS),
                listOf(Range(10, 81), Range(90, 93)),
                "ENCT"
        ))
        assertEquals(11, GtfReader.determineUTR3End5(
                Location(20, 80, chr, Strand.MINUS),
                listOf(Range(11, 14), Range(19, 90)),
                "ENCT"
        ))

        // No space:
        assertEquals(-1, GtfReader.determineUTR3End5(
                Location(20, 80, chr, Strand.PLUS),
                listOf(Range(10, 80), Range(90, 93)),
                "ENCT"
        ))

        assertEquals(-1, GtfReader.determineUTR3End5(
                Location(20, 80, chr, Strand.MINUS),
                listOf(Range(10, 13), Range(20, 90)),
                "ENCT"
        ))
    }

    private fun readTranscript(content: String): Transcript {
        val reader = GtfReader(content.reader().buffered(), Genome["to1"])
        val transcripts = reader.readTranscripts()
        assertEquals(1, transcripts.size)
        return transcripts[0]
    }
}