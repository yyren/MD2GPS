import math
import os
import random
import re
import time
import sys
import argparse
import copy
import networkx as nx
import numpy as np
import grandalf
import pandas as pd
import queue
import threading
import heapq
import statsmodels.stats.multitest as smm
from matplotlib import pylab as plt
from grandalf.layouts import SugiyamaLayout
from grandalf.graphs import Vertex,Edge,Graph,graph_core
from scipy import stats
from multiprocessing import Process



############### 初始化HPO树状图 ###################
def SubgraphInit(HPOoboFileWay: str):
    global Node_Num
    G = nx.DiGraph()
    with open(HPOoboFileWay, 'r') as f:
        omim_obo = re.split("\[Term\]", f.read())
        omim_obo = omim_obo[1:]
        for o in omim_obo:
            o = o.strip()
            osplit = re.split("\n", o)
            id = ''
            parents = []
            for l in osplit:
                m = re.match('id:\s(HP:\d{7})', l)
                if m is not None:
                    id = m.group(1)
                    if id not in G:
                        G.add_node(id)
                m2 = re.match('is_a:\s(HP:\d{7})\s\!', l)
                if m2 is not None:
                    parents.append(m2.group(1))
                    if m2.group(1) not in G:
                        G.add_node(m2.group(1))
            for p in parents:
                G.add_edge(p, id)
    Node_Num = len(G.nodes())
    return G

def getVariantGeneHPO(disease_name, disease_gene, variant_genes, gene_relation_data, gene_hpo_data):
    #disease related direct and indirect variant genes
    global disease_inheritance_df
    global gene_genotype
    global aachage
    global variant_data
    
    direct_genes=[]
    disease_related_pathway_genes=[]
    variant_gene_HPO=[]
    kegg_pathway_genes=gene_relation_data.columns.to_list()
    hpo_genes=gene_hpo_data['geneSymbol'].to_list()
    
    variant_gene=copy.deepcopy(variant_genes)
    disease_inheritance=disease_inheritance_df.loc[disease_inheritance_df['diseaseID']==disease_name]['inheritance'].to_list()
    # if len(disease_inheritance)>0:
        # #disease_inheritances=re.split(";",disease_inheritance[0])
        # #print(genotype['KMT2C'])
        # for gene in variant_genes:
            # if (not re.findall('dominant', disease_inheritance[0])) and ('Hom' not in gene_genotype[gene]):
                # variant_gene.remove(gene)
    for gene in disease_gene:
        if gene in variant_gene:
            #remove the gene against the inheritance of disease
            if len(disease_inheritance)>0:
                if (not re.findall('dominant', disease_inheritance[0])) and ('Hom' not in gene_genotype[gene]):
                    #candidate compound heterozygous
                    if len(variant_data.loc[variant_data['Gene_Name']==gene])>1:
                        disease_related_pathway_genes.append(gene)
                        direct_genes.append(gene)
                else:
                    disease_related_pathway_genes.append(gene)
                    direct_genes.append(gene)
            else:
                disease_related_pathway_genes.append(gene)
                direct_genes.append(gene)
            #get indirected gene
            if gene in kegg_pathway_genes:
                #rej, p_adj = smm.multipletests(gene_relation_data[gene].to_list(), alpha=0.05, method='fdr_i')[:2]
                rej, p_adj = smm.multipletests(gene_relation_data[gene].to_list(), alpha=0.05, method='bonferroni')[:2]
                p_adj_df=pd.DataFrame(p_adj)
                p_adj_df.columns=['pvalue']
                pathway_genes_sig=gene_relation_data.columns[p_adj_df.index[p_adj_df['pvalue']<0.05].to_list()].to_list()
                #indirected variant gene
                for path_gene in pathway_genes_sig:
                    if path_gene in variant_gene:
                        disease_related_pathway_genes.append(path_gene)
    disease_related_pathway_genes=list(set(disease_related_pathway_genes))
    #get variant genes related HPOs
    for gene in disease_related_pathway_genes:
        if gene in hpo_genes:
            hpo=gene_hpo_data.loc[gene_hpo_data['geneSymbol']==gene]['HPO_ID'].to_list()
            hpos=re.split(" ", hpo[0])
            for hpo_temp in hpos:
                variant_gene_HPO.append(hpo_temp)
    variant_gene_HPO=list(set(variant_gene_HPO))
    variant_gene_HPO.sort()
    if variant_gene_HPO != '':
        variant_gene_HPO.append('NA')
    direct_genes_str=''
    amino_acid_change=''
    if len(direct_genes) != 0:
        direct_genes.sort()
        for gene in direct_genes:
            aa_rec=aachange[gene]
            str_split='|'
            new_gene='true'
            for aa_change in aa_rec:
                if amino_acid_change=="":
                    amino_acid_change=aa_change
                else:
                    if new_gene == 'true':
                        amino_acid_change=str(amino_acid_change)+"|"+str(aa_change)
                    else:
                        amino_acid_change=str(amino_acid_change)+";"+str(aa_change)
                new_gene='false'
            if direct_genes_str == "":
                direct_genes_str=str(gene)
            else:
                direct_genes_str=str(direct_genes_str)+';'+str(gene)
    else:
        direct_genes_str='NA'
        amino_acid_change='NA'
    disease_related_pathway_genes_str=''
    if len(disease_related_pathway_genes) != 0:
        disease_related_pathway_genes.sort()
        for gene in disease_related_pathway_genes:
            if disease_related_pathway_genes_str == "":
                disease_related_pathway_genes_str=str(gene)
            else:
                disease_related_pathway_genes_str=str(disease_related_pathway_genes_str)+';'+str(gene)
    else:
        disease_related_pathway_genes_str='NA'
    return (variant_gene_HPO,direct_genes_str,amino_acid_change,disease_related_pathway_genes_str)

# get the parent HPOs of the specified disease
def get_diagnosis_results(threadID: str, sid: str, out_dir: str, start_idx: int, end_idx: int):
    global G
    global variant_data
    global variant_genes
    global gene_relation_df
    global disease_gene_df
    global gene_hpo_df
    global disease_inheritance_df
    global query_hpos
    global query_topology_hpos
    global disease_hpo_file
    
    enrich_results=[]
    HPO_pvalues=[]
    geneHPO_pvalues=[]
    hit_hpos=[]
    hit_genehpos=[]
    hit_key_hpos=[]
    genehpos_hit_key_query=[]
    idx=0
    with open(disease_hpo_file, 'r') as f:
        topology_hpos=set()
        start_time=time.asctime(time.localtime(time.time()))
        for row in f.readlines():
            if int(idx)>=start_idx and int(idx)<=end_idx:
                #print('tread '+str(threadID)+' | '+str(start_idx)+'_'+str(end_idx)+': '+str(idx))
                recs = row.strip().split('\t')
                if recs[0] != 'diseaseID':
                    disease_name=recs[0]
                    inheritances=disease_inheritance_df.loc[disease_inheritance_df['diseaseID']==disease_name]['inheritance'].to_list()
                    inheritances=list(set(inheritances))
                    inheritance=""
                    for inheri in inheritances:
                        if inheritance == "":
                            inheritance=inheri
                        else:
                            inheritance=str(inheritance)+";"+str(inheri)
                    if inheritance == "":
                        inheritance='NA'
                    
                    target_disease_topology_hpos=recs[1].split()
                    target_disease_topology_hpos.sort()
                    disease_edges=[]
                    #only retain the edges with both disease HPO terms
                    for ed in G.edges():
                        tag='retain'
                        for i in range(2):
                            if not (ed[i] in target_disease_topology_hpos):
                                tag='rm'
                        if tag=='retain':
                            disease_edges.append(ed)
                    disease_Sub=nx.DiGraph()
                    disease_Sub.add_edges_from(disease_edges)
                    hit_hpo_numb=0
                    hit_hpo=""
                    hit_key_hpo=""
                    hit_hpo_freq=[]
                    hit_obohpo=[]
                    hit_query_hpos=[]
                    for x in disease_Sub.nodes():
                        if x in query_hpos:
                            hit_query_hpos.append(x)
                            if hit_key_hpo=="":
                                hit_key_hpo=str(x)
                            else:
                                hit_key_hpo=str(hit_key_hpo)+';'+str(x)
                        if x in query_topology_hpos:
                            hpo_freq_temp=disease_hpo_freq_df.loc[disease_hpo_freq_df['diseaseID']==disease_name][x].item()
                            hit_hpo_freq.append(hpo_freq_temp)
                            hit_obohpo.append(x)
                            hit_hpo_numb += 1
                            if hit_hpo=="":
                                hit_hpo=str(x)
                            else:
                                hit_hpo=str(hit_hpo)+';'+str(x)
                    # calculate disease enrichment score
                    hit_hpo_freq_min=1.0
                    hit_obohpo_min="NA"
                    if len(hit_hpo_freq)>0:
                        hit_hpo_freq_min=min(hit_hpo_freq)
                        index_temp=hit_hpo_freq.index(hit_hpo_freq_min)
                        hit_obohpo_min=hit_obohpo[index_temp]
                        if len(hit_query_hpos)>0:
                            for x in hit_query_hpos:
                                if hit_hpo_freq_min==disease_hpo_freq_df.loc[disease_hpo_freq_df['diseaseID']==disease_name][x].item():
                                    hit_obohpo_min=x
                    hit_hpo_numb_statistic=hit_hpo_numb-1
                    all_disease_hpo_numb=8818
                    disease_hpo_numb=len(disease_Sub.nodes)
                    query_hpo_numb=len(query_topology_hpos)
                    enrich_pvalue_hpo=stats.hypergeom.sf(hit_hpo_numb_statistic,all_disease_hpo_numb,disease_hpo_numb,query_hpo_numb)
                    enrich_pvalue_hpo=enrich_pvalue_hpo*hit_hpo_freq_min
                    ##enrichment using the gene related HPO (directed or indirected)
                    disease_gene_str=disease_gene_df.loc[disease_gene_df['diseaseID']==disease_name]['geneSymbol'].to_list()
                    disease_gene=re.split(" ",disease_gene_str[0])
                    genehpos,hit_direct_gene,AAchange,pathwayGenes=getVariantGeneHPO(disease_name, disease_gene, variant_genes, gene_relation_df, gene_hpo_df)
                    #hit_direct_genes.append(hit_direct_gene)
                    #aachanges.append(AAchange)
                    #predict_pathwayGenes.append(pathwayGenes)
                    
                    # if hit_direct_gene == 'NA':
                        # genotypes.append('NA')
                        # gene_frequency.append('NA')
                        # gene_effect_impact.append('NA')
                        # gene_effect.append('NA')
                    # else:
                    genotype=""
                    frequency_list=[]
                    effect_impact=""
                    effect=""
                    if hit_direct_gene == 'NA':
                        genotype="NA"
                        frequency_list.append(0)
                        effect_impact="NA"
                        effect="NA"
                    else:
                        hit_genes=re.split(";",hit_direct_gene)
                        for gene in hit_genes:
                            if genotype == "":
                                genotype_temp=variant_data.loc[variant_data['Gene_Name']==gene]['Genome_Type'].to_list()
                                frequency_temp=variant_data.loc[variant_data['Gene_Name']==gene]['VarFreq_precent_max'].to_list()
                                effect_impact_temp=variant_data.loc[variant_data['Gene_Name']==gene]['Effect_Impact'].to_list()
                                effect=variant_data.loc[variant_data['Gene_Name']==gene]['Effect'].to_list()
                                genotype=";".join(genotype_temp)
                                effect=str(effect[0])
                                effect_impact=str(effect_impact_temp[0])
                                for temp in frequency_temp:
                                    frequency_list.append(temp)
                                for temp in effect_impact_temp:
                                    if re.findall('HIGH',str(temp)):
                                        effect_impact='HIGH'
                            else:
                                genotype_temp=variant_data.loc[variant_data['Gene_Name']==gene]['Genome_Type'].to_list()
                                genotype_temp_str=";".join(genotype_temp)
                                frequency_temp=variant_data.loc[variant_data['Gene_Name']==gene]['VarFreq_precent_max'].to_list()
                                effect_impact_temp=variant_data.loc[variant_data['Gene_Name']==gene]['Effect_Impact'].to_list()
                                effect_temp=variant_data.loc[variant_data['Gene_Name']==gene]['Effect'].to_list()
                                for temp in frequency_temp:
                                    frequency_list.append(temp)
                                genotype=str(genotype)+'|'+str(genotype_temp_str)
                                #frequency=str(frequency[0])+'|'+str(frequency_temp)
                                effect=str(effect_temp[0])
                                effect_impact=str(effect_impact_temp[0])
                                for temp in effect_impact_temp:
                                    if re.findall('HIGH',str(temp)):
                                        effect_impact='HIGH'
                    
                    # genotypes.append(genotype)
                    # gene_frequency.append(frequency)
                    # gene_effect_impact.append(effect_impact)
                    # gene_effect.append(effect)
                    frequency=max(frequency_list)
                    query_genehpo_numb=0
                    genehpo_hit_numb_statistic=0
                    genehpo_hit_numb=0
                    genehpo_hit=""
                    genehpo_hit_min=""
                    hits_genehpo_diseasefreq=[]
                    hits_genehpos_temp=[]
                    #get the hit hpos located on the downstream of the queried HPO and mapped to the disease HPO tree
                    if genehpos != 'NA':
                        query_genehpos=[]
                        query_diseasehpos_sub=[]
                        for x in genehpos:
                            if x not in query_topology_hpos:
                                query_genehpos.append(x)
                        
                        for x in disease_Sub.nodes():
                            if x not in query_topology_hpos:
                                query_diseasehpos_sub.append(x)
                        
                        for x in query_diseasehpos_sub:
                            if x in query_genehpos:
                                freq_temp=disease_hpo_freq_df.loc[disease_hpo_freq_df['diseaseID']==disease_name][x].item()
                                hits_genehpo_diseasefreq.append(freq_temp)
                                hits_genehpos_temp.append(x)
                                genehpo_hit_numb += 1
                                if genehpo_hit=="":
                                    genehpo_hit=str(x)
                                else:
                                    genehpo_hit=str(genehpo_hit)+';'+str(x)
                                # if x in query_hpos:
                                    # if genehpo_hit_query=="":
                                        # genehpo_hit_query=str(x)
                                    # else:
                                        # genehpo_hit_query=str(genehpo_hit_query)+';'+str(x)
                        
                        # for x in disease_Sub.nodes():
                            # if x in genehpos:
                                # freq_temp=disease_hpo_freq_df.loc[disease_hpo_freq_df['diseaseID']==disease_name][x].item()
                                # hits_genehpos_temp.append(x)
                                # hits_genehpo_diseasefreq.append(freq_temp)
                        # hits_genehpo_diseasefreq2=copy.deepcopy(hits_genehpo_diseasefreq)
                        # hits_genehpos_numb=len(hits_genehpos_temp)
                        # if genehpo_hit>query_hpo_numb:
                            # hits_genehpos_numb=query_hpo_numb
                        # if hits_genehpos_numb>0:
                            # for _ in range(hits_genehpos_numb):
                                # number=min(hits_genehpo_diseasefreq2)
                                # index_temp=hits_genehpo_diseasefreq2.index(number)
                                # genehpo_min_freqs.append(number)
                                # hits_genehpos.append(hits_genehpos_temp[index_temp])
                    
                    query_genehpo_numb=len(query_genehpos)
                    genehpo_min_freq=1.0
                    genehpo_hit_arr=[]
                    hits_genehpo_diseasefreq2=copy.deepcopy(hits_genehpo_diseasefreq)
                    if len(hits_genehpo_diseasefreq)>0:
                        genehpo_min_freq=min(hits_genehpo_diseasefreq)
                        for _ in range(len(hits_genehpo_diseasefreq)):
                            if len(genehpo_hit_arr)<5:
                                min_freq=min(hits_genehpo_diseasefreq2)
                                min_idx=hits_genehpo_diseasefreq2.index(min_freq)
                                genehpo_hit_arr.append(hits_genehpos_temp[min_idx])
                                hits_genehpo_diseasefreq2[min_idx]=1
                                if genehpo_hit_min=="":
                                    genehpo_hit_min=str(hits_genehpos_temp[min_idx])
                                else:
                                    genehpo_hit_min=str(genehpo_hit_min)+';'+str(hits_genehpos_temp[min_idx])
                    hits_genehpo_diseasefreq2=[]
                    genehpo_hit_numb_statistic=genehpo_hit_numb-1
                    enrich_pvalue_genehpo=stats.hypergeom.sf(genehpo_hit_numb_statistic,all_disease_hpo_numb,disease_hpo_numb,query_genehpo_numb)
                    enrich_pvalue_genehpo=enrich_pvalue_genehpo*genehpo_min_freq
                    variant_data_colname=variant_data.columns.to_list()
                    
                    pop_freq=1
                    pop_freq_placeholder=1
                    if hit_direct_gene != 'NA':
                        genes=re.split(";",hit_direct_gene)
                        for gene in genes:
                            pop_freq_value=0
                            for col_name in variant_data_colname:
                                if re.findall('gnomAD', str(col_name)):
                                    pop_freq_temp=variant_data.loc[variant_data['Gene_Name']==gene][col_name].to_list()
                                    for pop_f in pop_freq_temp:
                                        if pop_f == 'NA' or pop_f == '.':
                                            pop_freq_temp=0.0
                                        else:
                                            pop_freq_temp=float(pop_f)
                                        if pop_freq_temp>pop_freq_value:
                                            pop_freq_value=pop_freq_temp
                            if pop_freq_value == 0:
                                pop_freq_placeholder=0
                            else:
                                pop_freq=pop_freq_value
                        if pop_freq_placeholder == 0:
                            pop_freq=0
                    #print(str(pop_freq)+'\n')
                    # tree count, query count, hit count, pvalue
                    result_temp=(disease_name, hit_obohpo_min, hit_hpo_freq_min, hit_direct_gene, inheritance, pop_freq, genotype, AAchange, effect_impact, effect, frequency, pathwayGenes, disease_hpo_numb, query_hpo_numb, hit_hpo_numb, query_genehpo_numb,genehpo_hit_numb, enrich_pvalue_hpo, enrich_pvalue_genehpo)
                    enrich_results.append(result_temp)
                    HPO_pvalues.append(enrich_pvalue_hpo)
                    geneHPO_pvalues.append(enrich_pvalue_genehpo)
                    if hit_hpo != "":
                        hit_hpos.append(hit_hpo)
                    else:
                        hit_hpos.append('NA')
                    if hit_key_hpo == "":
                        hit_key_hpos.append('NA')
                    else:
                        hit_key_hpos.append(hit_key_hpo)
                    if genehpo_hit != "":
                        hit_genehpos.append(genehpo_hit)
                    else:
                        hit_genehpos.append('NA')
                    if genehpo_hit_min == "":
                        genehpos_hit_key_query.append('NA')
                    else:
                        genehpos_hit_key_query.append(genehpo_hit_min)
            idx = idx + 1
        end_time=time.asctime(time.localtime(time.time()))
        print('Thread'+str(threadID)+' Start: '+str(start_time)+'\n'+'Thread '+str(threadID)+'Finish: '+str(end_time)+'\n')
    
    f.close()
    
    # rej, hpo_p_adj = smm.multipletests(HPO_pvalues, alpha=0.05, method='bonferroni')[:2]
    # rej, genehpo_p_adj = smm.multipletests(geneHPO_pvalues, alpha=0.05, method='bonferroni')[:2]
    # HPO_p_adj_df=pd.DataFrame(hpo_p_adj)
    # genehpo_p_adj_df=pd.DataFrame(genehpo_p_adj)
    enrich_result_df=pd.DataFrame(enrich_results)
    hit_key_hpos_df=pd.DataFrame(hit_key_hpos)
    hit_key_genehpos_df=pd.DataFrame(genehpos_hit_key_query)
    hit_hpos_df=pd.DataFrame(hit_hpos)
    hit_genehpos_df=pd.DataFrame(hit_genehpos)
    out_recs=pd.concat([enrich_result_df,hit_key_hpos_df,hit_key_genehpos_df,hit_hpos_df,hit_genehpos_df],axis=1)
    out_recs.columns=['DiseaseName','Hit_Minimum_HPO', 'Hit_Minimum_HPO_Freq','Pathogenic_Gene','Disease_Inheritance','Population_Freq','Genotype','Amino_Change','Effect_Impact','Effect','VarFreq','Interaction_variantGenes','Disease_AncestorHPO_Count','Query_observedAncestorHPO_Count','Hit_observedHPO_Count','Query_Ancestor_GeneHPO_Count','Hit_GeneHPO_Count','Enrich_ObservedHPO_pvalue','Enrich_GeneHPO_pvalue','Hit_observedHPO','Hit_geneHPO','Hit_observedAncestorHPO','Hit_geneAncestorHPO']
    
    out_recs=out_recs.sort_values(by='Enrich_ObservedHPO_pvalue',ascending=True)
    #out_recs_result=out_recs.loc[(out_recs['GeneHPO_pvalue_Bonf']<0.05) & (out_recs['GeneHPO_pvalue_Bonf']<0.05)]
    out_file=str(out_dir)+'/'+str(sid)+'_'+str(threadID)+'.txt'
    out_recs.to_csv(out_file,sep='\t',index=False)
    

parser = argparse.ArgumentParser(description='get parent HPOs for the queried HPOs, output the their parent HPOs' )
parser.add_argument('--hpo', type=str, default="query_HPO.txt", help="file including the patient hpos")
parser.add_argument('--input', type=str, default="gene_variants.txt", help="input file including the gene variants")
parser.add_argument('--disease_hpo', type=str, default="Disease_HPO_tree.txt", help="file including the disease related hpos")
parser.add_argument('--disease_gene', type=str, default="Disease_Gene_2col.txt", help="file including the disease related genes")
parser.add_argument('--gene_relation', type=str, default="KEGG_Gene_Relation_pvalue_symbolID.txt", help="file including the gene relation statistic results")
parser.add_argument('--gene_hpo', type=str, default="Gene_HPO_symbol_2col.txt", help="file including the gene-hpos data")
parser.add_argument('--hpo_obo', type=str, default="hp.obo", help="the hpo tree in obo format")
parser.add_argument('--inheritance', type=str, default="Disease_inheritance.txt", help="the disease-inheritance data")
parser.add_argument('--threads', type=int, default=4, help="the threads")
parser.add_argument('--sid', type=str, default='test', help="the sample id")
parser.add_argument('--freq', type=str, default='Disease_HPO_freq.txt', help="the propotion of disease containing the HPO in the perticular disease")
parser.add_argument('--out_all', type=str, default="diagnosis_results.txt", help="output file including the diagnosis results")
parser.add_argument('--out_filter', type=str, default="diagnosis_results_filtered.txt", help="output file including the diagnosis results after filtering")
parser.add_argument('--out_filter2', type=str, default="diagnosis_results_filtered2.txt", help="output file including the diagnosis results after filtering")


args = parser.parse_args()
query_hpo_file=args.hpo
query_gene_file=args.input
disease_hpo_file=args.disease_hpo
disease_gene_file=args.disease_gene
gene_relation_file=args.gene_relation
gene_hpo_file=args.gene_hpo
HPO_obo_file=args.hpo_obo
inheritance_file=args.inheritance
threads=args.threads
sid=args.sid
freq_file=args.freq
out_all=args.out_all
out_filter=args.out_filter
out_filter2=args.out_filter2

# query_hpo_file='/export/home/renyongyong/project/PHD/rare_disease/results/CN-2100561/target_hpo.txt'
# query_gene_file='/export/home/renyongyong/project/PHD/rare_disease/results/CN-2100561/results_all_filtered.txt'
# disease_hpo_file='/export/home/renyongyong/project/PHD/rare_disease/database/Disease_HPO_tree.txt'
# disease_gene_file='/export/home/renyongyong/project/PHD/rare_disease/database/Disease_Gene_2col.txt'
# gene_relation_file='/export/home/renyongyong/project/PHD/rare_disease/database/KEGG_Gene_Relation_pvalue_symbolID.txt'
# gene_hpo_file='/export/home/renyongyong/project/PHD/rare_disease/database/Gene_HPO_symbol_2col.txt'
# HPO_obo_file='/export/home/renyongyong/project/PHD/rare_disease/database/hp.obo'
# inheritance_file='/export/home/renyongyong/project/PHD/rare_disease/database/Disease_inheritance.txt'
# freq_file='/export/home/renyongyong/project/PHD/rare_disease/database/Disease_HPO_frequency.txt'
# out_all='/export/home/renyongyong/project/PHD/rare_disease/results/CN-2100561/test_diagnosis_all.txt'
# out_filter='/export/home/renyongyong/project/PHD/rare_disease/results/CN-2100561/test_diagnosis_all_filter.txt'
# out_filter2='/export/home/renyongyong/project/PHD/rare_disease/results/CN-2100561/test_diagnosis_all_filter2.txt'
# sid='CN-2100561'
# threads=30



query_hpos=[]
query_topology_hpos=[]
G = nx.DiGraph()
G = SubgraphInit(HPOoboFileWay=HPO_obo_file)

# get variant gene
variant_data=pd.read_csv(query_gene_file, header=0, sep='\t', low_memory=False, keep_default_na=False)
variant_genes=variant_data['Gene_Name']
variant_genes=variant_genes.drop_duplicates()
variant_genes=variant_genes.tolist()

# load gene-gene relation
gene_relation_df=pd.read_csv(gene_relation_file, header=0, sep=' ', low_memory=False)

# load disease-gene data
disease_gene_df=pd.read_csv(disease_gene_file, header=0, sep='\t', low_memory=False)

#load gene-hpo data
gene_hpo_df=pd.read_csv(gene_hpo_file, header=0, sep='\t', low_memory=False)

#load disease-inheritance data
disease_inheritance_df=pd.read_csv(inheritance_file, header=0, sep='\t', low_memory=False)

#load disease-hpo frequency data
disease_hpo_freq_df=pd.read_csv(freq_file, header=0, sep=' ', low_memory=False)

# get genotype of each gene
gene_genotype={}
aachange={}
with open(query_gene_file,'r') as f:
    idx=0
    for row in f.readlines():
        idx = idx +1
        if idx > 1:
            recs=row.strip().split('\t')
            bases=recs[13].split('/')
            geno='Hom'
            if bases[0] != bases[1]:
                geno='Het'
            if recs[20] not in gene_genotype.keys():
                gene_genotype[recs[20]]=[geno]
            else:
                if geno not in gene_genotype[recs[20]]:
                    gene_genotype[recs[20]].append(geno)
            if recs[20] not in aachange.keys():
                aachange[recs[20]]=[recs[18]]
            else:
                if recs[18] not in aachange[recs[20]]:
                    aachange[recs[20]].append(recs[18])

f.close()
# get the parent HPOs of the queried HPOs
with open(query_hpo_file, 'r') as f:
    topology_hpos=set()
    for row in f.readlines():
        hpos = row.strip().split()
        for hpo in hpos:
            if hpo not in query_hpos:
                query_hpos.append(hpo)
            if hpo in G.nodes:
                ancestor_hpos=nx.ancestors(G, hpo)
                for ancensor_hpo in ancestor_hpos:
                    topology_hpos.add(ancensor_hpo)
            topology_hpos.add(hpo)
    query_topology_hpos=list(topology_hpos)
    query_topology_hpos.sort()

f.close()

#run with multiple process
def myThread(threadID, sid, out_dir, start_idx: int, end_idx: int):
    print ("Start " + str(sid) + ':' + str(start_idx) + '_' + str(end_idx))
    get_diagnosis_results(threadID, sid, out_dir, start_idx, end_idx)

threads_list = []
threadID = 1

# creat and start multiple threads
out_dir=os.path.abspath(os.path.dirname(out_all))
for i in range(threads):
    splice_count=int(5685/int(threads))#5685 diseases in OMIM database
    start_idx=i*splice_count
    end_idx=start_idx+splice_count-1
    if i==0:
        start_idx=1
    elif i == (threads-1):
        end_idx=5685
    #thread = myThread(threadID, sid, out_dir, start_idx, end_idx)
    thread = Process(target=myThread,args=(threadID, sid, out_dir, start_idx, end_idx,))
    thread.start()
    threads_list.append(thread)
    threadID += 1

for t in threads_list:
    t.join()

print ("Start merge results from multiple threads")
out_results_df=pd.DataFrame()
for i in range(threads):
    idx=i+1
    temp_out_file=str(out_dir)+'/'+str(sid)+'_'+str(idx)+'.txt'
    temp_df=pd.read_csv(temp_out_file, header=0, sep='\t', low_memory=False, keep_default_na=False)
    out_results_df=pd.concat([out_results_df,temp_df],axis=0,ignore_index=True)
    os.remove(temp_out_file)

rej, hpo_p_adj = smm.multipletests(out_results_df['Enrich_ObservedHPO_pvalue'].to_list(), alpha=0.05, method='bonferroni')[:2]
rej, genehpo_p_adj = smm.multipletests(out_results_df['Enrich_GeneHPO_pvalue'].to_list(), alpha=0.05, method='bonferroni')[:2]
HPO_p_adj_df=pd.DataFrame(hpo_p_adj)
genehpo_p_adj_df=pd.DataFrame(genehpo_p_adj)


out_recs=pd.concat([out_results_df,HPO_p_adj_df,genehpo_p_adj_df],axis=1)
out_recs.columns=['DiseaseName','Hit_Minimum_HPO','Hit_Minimum_HPO_Freq','Pathogenic_Gene','Disease_Inheritance','Population_Freq','Genotype','Amino_Change','Effect_Impact','Effect','VarFreq','Interaction_variantGenes','Disease_AncestorHPO_Count','Query_observedAncestorHPO_Count','Hit_observedHPO_Count','Query_Ancestor_GeneHPO_Count','Hit_GeneHPO_Count','Enrich_ObservedHPO_pvalue','Enrich_GeneHPO_pvalue','Hit_observedHPO','Hit_geneHPO','Hit_observedAncestorHPO','Hit_geneAncestorHPO','Enrich_ObservedHPO_pvalue_Bonf','Enrich_GeneHPO_pvalue_Bonf']
out_recs['VarFreq_Group']=np.where(out_recs['VarFreq'] <30,0,1)

temp_df=pd.DataFrame()
min_pvalue_1=out_recs.loc[out_recs['Enrich_ObservedHPO_pvalue_Bonf']!=0,'Enrich_ObservedHPO_pvalue_Bonf'].min()
temp_df['Enrich_ObservedHPO_pvalue_Bonf']=out_recs['Enrich_ObservedHPO_pvalue_Bonf']+min_pvalue_1
p1=(temp_df['Enrich_ObservedHPO_pvalue_Bonf']).map(lambda x : -np.log10(x))

min_pvalue_2=out_recs.loc[out_recs['Enrich_GeneHPO_pvalue_Bonf']!=0,'Enrich_GeneHPO_pvalue_Bonf'].min()
temp_df['Enrich_GeneHPO_pvalue_Bonf']=out_recs['Enrich_GeneHPO_pvalue_Bonf']+min_pvalue_2
p2=(temp_df['Enrich_GeneHPO_pvalue_Bonf']).map(lambda x : -np.log10(x))

max_value_1=max(p1)
max_value_2 = max(p2)
min_value_2 = min(p2)
p2_normalize=(temp_df['Enrich_GeneHPO_pvalue_Bonf']).map(lambda x : ((-np.log10(x)-min_value_2)/(max_value_2-min_value_2)*max_value_1))

#p2=(out_recs['Enrich_GeneHPO_pvalue_Bonf']+(1E-500)).map(lambda x : -np.log10(x))
out_recs['Enrich_observeHPO_score']=p1
out_recs['Enrich_GeneHPO_score']=p2_normalize
out_recs['Score_max'] =out_recs[['Enrich_observeHPO_score','Enrich_GeneHPO_score']].max(axis=1).map(lambda x : x)
out_recs['Score_total'] =p1+p2_normalize
#out_recs['Pvalue_log_max'] = temp_df[['Enrich_ObservedHPO_pvalue_Bonf','Enrich_GeneHPO_pvalue_Bonf']].min(axis=1).map(lambda x : -np.log10(x))

p3=np.where((out_recs['Effect_Impact']=='HIGH'),1,0)
p4=np.where((out_recs['Hit_observedHPO']!='NA'),1,0)
out_recs['HPO_geneImpact_group']=p3+p4

col_order=['DiseaseName','Hit_Minimum_HPO','Hit_Minimum_HPO_Freq','Enrich_ObservedHPO_pvalue_Bonf','Enrich_GeneHPO_pvalue_Bonf','Enrich_observeHPO_score','Enrich_GeneHPO_score','Score_max','Score_total','HPO_geneImpact_group','VarFreq_Group','Population_Freq','Effect_Impact','Disease_Inheritance','Genotype','Pathogenic_Gene','Amino_Change','Effect','VarFreq','Interaction_variantGenes','Disease_AncestorHPO_Count','Query_observedAncestorHPO_Count','Hit_observedHPO_Count','Query_Ancestor_GeneHPO_Count','Hit_GeneHPO_Count','Enrich_ObservedHPO_pvalue','Enrich_GeneHPO_pvalue','Hit_observedHPO','Hit_geneHPO','Hit_observedAncestorHPO','Hit_geneAncestorHPO']
out_recs=out_recs[col_order]
out_recs=out_recs.sort_values(by='Score_max',ascending=False)
out_recs.to_csv(out_all,sep='\t',index=False)

# filter and sort
print('Start output filtered results: '+str(out_filter))
impact_effect_order=['HIGH','MODERATE','MODIFIER','NA']
#out_recs_filtered=out_recs.loc[(out_recs['Enrich_ObservedHPO_pvalue_Bonf']<0.05) & (out_recs['Enrich_GeneHPO_pvalue_Bonf']<0.05) & (out_recs['Hit_observedHPO'] !='NA') & (out_recs['Pathogenic_Gene'] !='NA')]
#out_recs_filtered=out_recs.loc[(out_recs['VarFreq']>=30) & (out_recs['Enrich_ObservedHPO_pvalue_Bonf']<0.05) & (out_recs['Enrich_GeneHPO_pvalue_Bonf']<0.05) & (out_recs['Pathogenic_Gene'] !='NA') & (not out_recs['Effect'].str.contains('inframe_deletion'))]
#out_recs_filtered=out_recs.loc[(out_recs['VarFreq']>=30) & (~ out_recs['Effect'].str.contains('inframe_deletion'))]
out_recs_filtered=out_recs.loc[(~ out_recs['Effect'].str.contains('inframe_deletion'))]
out_recs_filtered['Effect_Impact']=out_recs_filtered['Effect_Impact'].astype('category')
out_recs_filtered['Effect_Impact'].cat.set_categories(impact_effect_order, inplace=True)

out_recs_filtered_sort=out_recs_filtered_sort.sort_values(by=['HPO_geneImpact_group','Hit_Minimum_HPO_Freq','Score_total'], ascending=(False,True,False))
out_recs_filtered_sort.to_csv(out_filter,sep='\t',index=False)

out_recs_filtered=out_recs.loc[(out_recs['VarFreq']>=30) & (~ out_recs['Effect'].str.contains('inframe_deletion'))]
out_recs_filtered['Effect_Impact']=out_recs_filtered['Effect_Impact'].astype('category')
out_recs_filtered['Effect_Impact'].cat.set_categories(impact_effect_order, inplace=True)

out_recs_filtered_sort2=out_recs_filtered_sort.sort_values(by=['HPO_geneImpact_group','Hit_Minimum_HPO_Freq','Score_total'], ascending=(False,True,False))
out_recs_filtered_sort2.to_csv(out_filter2,sep='\t',index=False)

