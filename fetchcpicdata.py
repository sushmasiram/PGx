import os
import json
import pandas as pd
import ast
import re
import requests
import xml.etree.ElementTree as ET

caution_terms =['Shorted feasible course',
                'use drug cautiously',
                'Initiate, increase, absence of adverse',
                'start with reduced starting doses',
                'Start with reduced doses',
                'reduce starting dose',
                'lowest recommended starting dose',
                'consider combination therapy',
                'possible increased risk for',
                'titrating to a higher',
                '%% of starting dose',
                'standard doses for the shortest feasible course with therapeutic dose monitoring',
                'Use drug cautiously',
                'reduction of recommended starting dose',
                'Start with drastically reduced dose',
                'therapeutic drug monitoring to guide dose adjustments',
                'cautiously consider',
                'slower titration', 
                "less than typical maintenance dose",
                "%% reduction",
                '%% of normal dose',
                'possible increased risk for, disease-specific guidelines',
                'Increase starting daily dose by 100%',
                '^Initiate, reduction',
                'decreased dose',
                'with caution and',
                'Increase starting dose',
                'increased risk to',
                'Daily dose may be given in divided doses',
                'Consider initiating .* with decreased dose',
                'Prescribe .*mg as a starting dose, combination therapy']

avoid_terms = ['naive',
                'nave',
                '^avoid',
                'do not use',
                '^Choose an alternative agent',
                '^Select alternative drug',
                '^Prescribe an alternative',
                'alternative','is contraindicated',
                'is not recommended',
                'if no response, consider non-codeine opioid',
                'If no response and opioid use is warranted']

normal_terms = ['standard dosing guidelines',
                'No reason to avoid',
                'with recommended starting dose',
                'with recommended standard of care dosing',
                'label recommended age- or weight-specific dosing',
                'desired starting dose',
                'start with normal starting dose',
                'Use, per standard dosing guidelines',
                '^Prescribe desired starting dose',
                'possible increased risk for',
                'standard recommended dose',
                'No adjustments needed',
                'For first dose, use typical initial or loading dose',
                'no need to avoid prescribing',
                'use at standard dose',
                '^Initiate .* with standard dosing',
                '^Use, according to the product label',
                '^Initiate therapy with standard recommended dose',
                'Use at standard doses with caution and with close monitoring for anemia']

gtf_file = '../GenCode/Homo_sapiens.GRCh38.112.chr.gtf'
allele_loc = pd.read_csv('../DBs/CPIC/allele_location.csv')
allele_def = pd.read_csv('../DBs/CPIC/allele_definition.csv')
seq_loc = pd.read_csv('../DBs/CPIC/sequence_location.csv')
genes = pd.read_csv('../DBs/CPIC/gene.csv',encoding='latin')
allele = pd.read_csv('../DBs/CPIC/allele.csv')

population = pd.read_csv('../DBs/CPIC/population.csv')
publication = pd.read_csv('../DBs/CPIC/publication.csv')
guideline = pd.read_csv('../DBs/CPIC/guideline.csv')

pair = pd.read_csv('../DBs/CPIC/pair.csv')

test_alert = pd.read_csv('../DBs/CPIC/test_alert.csv')

drug = pd.read_csv('../DBs/CPIC/drug.csv')
recommendation = pd.read_csv('../DBs/CPIC/recommendation.csv')

gene_result_diplotype = pd.read_csv('../DBs/CPIC/gene_result_diplotype.csv')
gene_result_lookup = pd.read_csv('../DBs/CPIC/gene_result_lookup.csv')
gene_result = pd.read_csv('../DBs/CPIC/gene_result.csv')

allele_freq = pd.read_csv('../DBs/CPIC/allele_frequency.csv')

def find_df_name(df):
    name = [name for name, obj in globals().items() if id(obj) == id(df)]
    return name[0] if name else None


def clean_data():
    allele_def.rename(columns={'id':'alleledefinitionid'},inplace=True)
    seq_loc.rename(columns={'id':'locationid'},inplace=True)
    allele.rename(columns={'definitionid':'alleledefinitionid'},inplace=True)

    allele_def.rename(columns={'name':'name_alleledef'},inplace=True)
    seq_loc.rename(columns={'name':'name_seqloc'},inplace=True)
    allele.rename(columns={'name':'name_allele'},inplace=True)

    publication.rename(columns={'id':'publicationid'},inplace=True)
    population.rename(columns={'id':'populationid'},inplace=True)
    genes.rename(columns={'symbol':'genesymbol'},inplace=True)
    gene_result.rename(columns={'id':'phenotypeid'},inplace=True)
    gene_result_lookup.rename(columns={'id':'functionphenotypeid'},inplace=True)
    gene_result_diplotype.rename(columns={'id':'gene_result_diplotype_id'},inplace=True)
    allele.rename(columns={'id':'alleleid'},inplace=True)
    recommendation.rename(columns={'id':'recommendationid'},inplace=True)
    test_alert.rename(columns={'id':'test_alertid'},inplace=True)
    guideline.rename(columns={'id':'guidelineid'},inplace=True)


def add_suffix():

    list_of_dfs = [allele,allele_def,allele_loc,genes,gene_result,gene_result_diplotype,gene_result_lookup,population,
                publication,guideline,test_alert,drug,recommendation,seq_loc,allele_freq]  # Your list of DataFrames
    suffix = '_version'

    for df in list_of_dfs:
        if 'version' in df.columns:
            df_name = [name for name, frame in globals().items() if frame is df][0]
            new_column_name = df_name + suffix
            df.rename(columns={'version': new_column_name}, inplace=True)

        for col in df.columns:
            df[col] = df[col].astype(str).str.replace(r'[^\x00-\x7F]+', '', regex=True)
            
    clean_data()
    cols_exclude = ['alleledefinitionid','locationid','guidelineid','test_alertid','recommendationid',
    'alleleid','functionphenotypeid','phenotypeid','genesymbol', 'populationid','diplotypekey','diplotype',
    'publicationid','name_allele', 'name_seqloc','name_alleledef','pairid','drugid','result','gene_result_diplotype_id']
    for df in list_of_dfs:
        # Iterate over each column in the DataFrame
        for column in df.columns:
            # Check if the column is not in the excluded list x
            if column not in cols_exclude:
                # Rename the column as column_dfname_columnname_dataframe_name
                df_name = find_df_name(df)
                new_column_name = f'{column}_{df_name}'
                df.rename(columns={column: new_column_name}, inplace=True)

# Function to convert JSON string to dictionary if necessary
def str_to_dict(json_str):
    if isinstance(json_str, str):
        return json.loads(json_str)
    return json_str

def check_match(row, query):  
    if (list(row["diplotypekey"].keys())[0] == query[0]):
        try:
            value1 = (row["diplotypekey"][query[0]][query[1]])
            value2 = (row["diplotypekey"][query[0]][query[2]])
            if value1==1 and value2==1:
                return True
        except:
            return False
    else:
        return False
    
# Function to check if the query matches the values in diplotypekey
def check_diplo_match(row,query):
    for key,value in query.items():        
        if list(row["diplotypekey"].keys())[0] == key:
            
            if (value == row['diplotype']):
                return True
            else:
                return False
        else:
            return False

def get_gene_diplo_info(diplo_query):
    dict={}
    for key,value in diplo_query.items():
        for i in value:
            dict[key]= i

    # Apply the function to the 'diplotypekey' column
    gene_result_diplotype["diplotypekey"] = gene_result_diplotype["diplotypekey"].apply(str_to_dict)

    num = 0 
    filtered_df_diplo = pd.DataFrame() 
    for key,value in diplo_query.items():
        for i in value:        
            num = num+1   
            df = gene_result_diplotype[gene_result_diplotype.apply(lambda row: check_diplo_match(row, {key:i}), axis=1)]
            filtered_df_diplo = pd.concat([filtered_df_diplo,df],ignore_index=True)

    merged_diplo_lookup = pd.merge(filtered_df_diplo,gene_result_lookup,how='left')
    merged_lookup_result = pd.merge(merged_diplo_lookup,gene_result,on=['phenotypeid'],how='left')
    df_genes = pd.merge(merged_lookup_result,genes,on=['genesymbol'],how='left')
    return df_genes

def make_lookupdict(df_genes):
    lookup_key={}
    for index, row in df_genes[['genesymbol','totalactivityscore_gene_result_lookup','lookupmethod_genes','result']].drop_duplicates().iterrows():
        row['lookup_key_genes'] = {}
        if row['lookupmethod_genes']=='ACTIVITY_SCORE':
            row['lookup_key_genes'] = {row['genesymbol']:row['totalactivityscore_gene_result_lookup']}
            lookup_key.update({row['genesymbol']:row['totalactivityscore_gene_result_lookup']})
        else:
            row['lookup_key_genes'] = {row['genesymbol']:row['result']}
            lookup_key.update({row['genesymbol']:row['result']})
    return lookup_key, df

def match_row(row, match_dict):
    row_dict = json.loads(row)
    return all(item in match_dict.items() for item in row_dict.items())

def get_reco(lookup_key, relevant_rows):    
    recommendation['matches'] = recommendation['lookupkey_recommendation'].apply(lambda x: match_row(x, lookup_key))
    matching_rows = recommendation[recommendation['matches']]
    matching_rows['guidelineid'] = matching_rows['guidelineid'].astype(str)
    drug['guidelineid'] = drug['guidelineid'].astype(str)
    drug['guidelineid'] = drug['guidelineid'].replace('.0','')
    recom_req = drug[['drugid', 'pharmgkbid_drug', 'name_drug']].merge(matching_rows, on='drugid', how='right')
    recom_req = recom_req.merge(relevant_rows, how='left', on=['genesymbol',])
    return recom_req

def get_reco1(lookup_key):    
    recommendation['matches'] = recommendation['lookupkey_recommendation'].apply(lambda x: match_row(x, lookup_key))
    matching_rows = recommendation[recommendation['matches']]
    matching_rows['guidelineid'] = matching_rows['guidelineid'].astype(str)
    drug['guidelineid'] = drug['guidelineid'].astype(str)
    drug['guidelineid'] = drug['guidelineid'].replace('.0','' )   
    recom_req = drug[['drugid','pharmgkbid_drug','name_drug']].merge(matching_rows,on='drugid',how='right')
    return recom_req

def str_to_dict(json_str):
    if isinstance(json_str, str):
        return json.loads(json_str)
    return json_str

def filter_reportablereco(CPIC_rec):
    CPIC_pair = pair.copy()
    CPIC_pair_cl = CPIC_pair[CPIC_pair['cpiclevel'].str.contains('A|B|A/B')]
    CPIC_pair_pl = CPIC_pair[CPIC_pair['pgkbcalevel'].str.contains('1A|1B|2A|2B',na=False)]
    CPIC_pair_FDA = CPIC_pair[CPIC_pair['pgxtesting'].str.contains('Actionable PGx|Testing Recommended|Informative PGx|Testing Required', na=False)]
    CPIC_pair_filtered = pd.concat([CPIC_pair_cl,CPIC_pair_pl,CPIC_pair_FDA],ignore_index=True)
    CPIC_pair_filtered_nodups = CPIC_pair_filtered.drop_duplicates()

    Gene_drugPair = CPIC_pair_filtered_nodups[['drugid','genesymbol','guidelineid','pgxtesting']].drop_duplicates()
    Gene_drugPair['guidelineid'] = Gene_drugPair['guidelineid'].astype(str).replace(r'\.0$', '', regex=True)

    CPIC_rec['Genes'] = CPIC_rec['lookupkey_recommendation'].apply(lambda x: list(str_to_dict(x).keys()))
    CPIC_rec_filt = pd.DataFrame()
    
    for index,row in Gene_drugPair.iterrows():
        filtered = CPIC_rec[(CPIC_rec['drugid'] == row['drugid'])&(CPIC_rec['guidelineid'] == row['guidelineid']) & (CPIC_rec['Genes'].apply(lambda x: row['genesymbol'] in x))]
        filtered['pgxtesting'] = row['pgxtesting']
        CPIC_rec_filt = pd.concat([CPIC_rec_filt, filtered], ignore_index=True)
    CPIC_rec_filt.drop(columns='Genes',inplace=True)
    CPIC_rec_filt_nodup = CPIC_rec_filt.drop_duplicates()
        
    CPIC_rec_filt_nodup_class=CPIC_rec_filt_nodup[CPIC_rec_filt_nodup['classification_recommendation'].isin(['Strong', 'Moderate'])]
    
    return CPIC_rec_filt_nodup_class

def check_terms(row, terms_list, sec_term):    
    if row['category'] is None:
        for term in terms_list:
            if ',' in term:  # If term is comma-separated
                subterms = term.split(', ')
    
                if all(re.search(f'\\b{subterm}\\b', row['drugrecommendation_recommendation'], re.IGNORECASE) for subterm in subterms):
                    return sec_term
            else:
                if re.search(term, row['drugrecommendation_recommendation'], re.IGNORECASE):
                    return sec_term
    else:
        return row['category']
    
def add_category(recom_req):
    rec_fil = recom_req.copy()
    rec_fil['category'] = None
    rec_fil['category'] = rec_fil.apply(lambda row: check_terms(row, avoid_terms,'Use Alternative Drug'), axis=1)
    rec_fil['category'] = rec_fil.apply(lambda row: check_terms(row, caution_terms,'Use with Caution'), axis=1)
    rec_fil['category'] = rec_fil.apply(lambda row: check_terms(row, normal_terms,'Normal'), axis=1)
    return rec_fil

def add_reference(df):
    publication[df['guidelineid'] == publication['guidelineid']]
    
    
def get_citation_for_pmid(pmid):
    
    base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
    params = {
        'db': 'pubmed',
        'id': pmid,
        'retmode': 'xml'
    }
    response = requests.get(base_url, params=params)
    if response.status_code == 200:
        root = ET.fromstring(response.content)
        docsum = root.find('DocSum')
        if docsum is not None:
            title = docsum.findtext("Item[@Name='Title']")
            authors = docsum.findtext("Item[@Name='AuthorList']")
            source = docsum.findtext("Item[@Name='Source']")
            citation = f"{authors}. {title}. {source}. PMID: {pmid}"
            return citation
        else:
            return ""
    else:
        return ""

def add_citation(rec_df):
    rec_df_fil = rec_df[rec_df['guidelineid'].notna()]
    publication['guidelineid'] = publication['guidelineid'].str.replace('.0','')
    rec_guid = pd.merge(rec_df_fil, guideline, on='guidelineid', how='left')   
    rec_guid_pub = pd.merge(rec_guid, publication, on ='guidelineid', how='left')    
    
    rec_guid_pub['pmid_citation'] = rec_guid_pub['pmid_publication'].apply(lambda x: get_citation_for_pmid(x) if pd.notna(x) else "")
    
    return rec_guid_pub

def process_lookupkey(cell):
    key_value_pairs = []
    if isinstance(cell, str):
        entries = cell.split(';')  # Split the entries based on the delimiter
        for entry in entries:
            dict_entry = str_to_dict(entry)
            key_value_pairs.extend(dict_entry.items())  # Add all key-value pairs from the dictionary to the list
    elif isinstance(cell, dict):
        key_value_pairs.extend(cell.items())
    return key_value_pairs


def reportable_alleles():    

    CPIC_pair_cl = pair[pair['cpiclevel'].str.contains('A|B|A/B')]
    CPIC_pair_pl = pair[pair['pgkbcalevel'].str.contains('1A|1B|2A|2B',na=False)]
    CPIC_pair_FDA = pair[pair['pgxtesting'].str.contains('Actionable PGx|Testing Recommended|Informative PGx|Testing Required', na=False)]
    CPIC_pair_filtered = pd.concat([CPIC_pair_cl,CPIC_pair_pl,CPIC_pair_FDA],ignore_index=True)
    CPIC_pair_filtered_nodups = CPIC_pair_filtered.drop_duplicates()

    Gene_drugPair = CPIC_pair_filtered_nodups[['drugid','genesymbol','guidelineid']].drop_duplicates()   
    Gene_drugPair['guidelineid'] = Gene_drugPair['guidelineid'].astype('str').str.replace('.0','')
    recommendation['Genes'] = recommendation['lookupkey_recommendation'].apply(lambda x: list(str_to_dict(x).keys()))
    CPIC_rec_filt = pd.DataFrame()

    for index,row in Gene_drugPair.iterrows():
        filtered = recommendation[(recommendation['drugid'] == row['drugid'])&(recommendation['guidelineid'] == row['guidelineid']) & (recommendation['Genes'].apply(lambda x: row['genesymbol'] in x))]
        CPIC_rec_filt = pd.concat([CPIC_rec_filt, filtered], ignore_index=True)
    
    CPIC_rec_filt.drop(columns='Genes',inplace=True)
    CPIC_rec_filt_nodup = CPIC_rec_filt.drop_duplicates()
    CPIC_rec_filt_nodup_class=CPIC_rec_filt_nodup[CPIC_rec_filt_nodup['classification_recommendation'].isin(['Strong', 'Moderate'])]
    key_value_pairs = CPIC_rec_filt_nodup_class['lookupkey_recommendation'].apply(process_lookupkey).to_dict()
    key_value_pairs_set=set()

    for i in key_value_pairs:
        for j in key_value_pairs[i]:
            key_value_pairs_set.add((j[0],j[1]))
    
    gr_filt = pd.DataFrame()
    gene_result['activityscore_gene_result'] = gene_result['activityscore_gene_result'].astype(str)
    
    for i in key_value_pairs_set:
        gene = i[0]
        result = i[1]

        new = gene_result.loc[
        (gene_result['genesymbol'] == gene) & 
        (
            ((gene_result['activityscore_gene_result'] == 'nan') & (gene_result['result'] == result)) |
            ((gene_result['activityscore_gene_result'] !='nan') & (gene_result['activityscore_gene_result'] == result))
        )
        ]
        gr_filt = pd.concat([gr_filt,new],ignore_index=True)

    grl_filt = pd.merge(gr_filt, gene_result_lookup,how='left')
    dip_filt = pd.merge(gene_result_diplotype, grl_filt,on='functionphenotypeid',how='inner')
    filtered_rows = dip_filt[dip_filt['diplotype'].str.startswith('*') | dip_filt['diplotype'].str.startswith('rs') | dip_filt['diplotype'].str.startswith('c.') | dip_filt['diplotype'].str.startswith('m.')]
    genes_to_haplotypes = {}

    # Group by the 'genes' column
    grouped = filtered_rows.groupby('genesymbol')['diplotype'].apply(list)

    # Iterate through the grouped data
    for gene, diplotypes_list in grouped.items():
        # Use a set to store unique haplotypes
        haplotypes = set()
        for diplotypes in diplotypes_list:
            haplotypes.update(diplotypes.split('/'))
        # Convert the set back to a list and sort it (optional)
        genes_to_haplotypes[gene] = sorted(haplotypes)
    
    genes_to_haplotypes_filt={}
    filtered_allele = pd.DataFrame()
    # Add   VKORC1	rs9923231 variant (T), VKORC1	rs9923231 reference (C)

    
    for i in genes_to_haplotypes:
        genes_to_haplotypes_filt[i] = []
        for j in genes_to_haplotypes[i]:
            found = allele[
                (allele['genesymbol'] == i) & 
                (allele['name_allele'] == j) & 
                (~allele['clinicalfunctionalstatus_allele'].isnull()) & 
                (allele['clinicalfunctionalstatus_allele']!='nan') &
                (~allele['clinicalfunctionalstatus_allele'].str.contains('Uncertain function|Unknown function', na=False)) & 
                (~allele['strength_allele'].isnull()) & 
                (allele['strength_allele'] != 'nan') & 
                (allele['strength_allele'].str.lower().str.strip() != 'limited')
            ]

            if found.empty:
                continue
            filtered_allele = pd.concat([filtered_allele,found],ignore_index=False)
            genes_to_haplotypes_filt[i].append(j)

    '''fil_allele_def = allele.merge(allele_def,on='alleledefinitionid',how='left')
    fil_allele_loc = fil_allele_def.merge(allele_loc,on='alleledefinitionid',how='left')
    fil_seq_loc = fil_allele_loc.merge(seq_loc, on = 'locationid', how ='left')'''
    fil_allele_def = filtered_allele.merge(allele_def,on='alleledefinitionid',how='left')
    fil_allele_loc = fil_allele_def.merge(allele_loc,on='alleledefinitionid',how='left')
    fil_seq_loc = fil_allele_loc.merge(seq_loc, on = 'locationid', how ='left')
    return fil_seq_loc


def read_gtf(gtf_file, gene_list):
    column_names = [
        'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'
    ]
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=column_names)
    exons = gtf[gtf['feature'] == 'exon']
    exons['gene_name'] = exons['attribute'].str.extract('gene_name "([^"]+)"')
    filtered_exons = exons[exons['gene_name'].isin(gene_list)]
    return filtered_exons


def check_location_in_exon(exons, position, gene_name):
    gene_exons = exons[(exons['gene_name'] == gene_name)]
    for _, exon in gene_exons.iterrows():
        if exon['start'] <= position <= exon['end']:
            return True
    return False

def repor_alleles_exonic():
    CPIC_rsids  = reportable_alleles()
    
    gene_list=CPIC_rsids['genesymbol'].unique().tolist()
    
    exons = read_gtf(gtf_file, gene_list)
    CPIC_rsids['IsExonic']= False
    CPIC_rsids['position_seq_loc'] = CPIC_rsids['position_seq_loc'].astype(float)
    for index, row in CPIC_rsids.iterrows():
        chr_pos = row['position_seq_loc']
        gene = row['genesymbol']
        is_exonic = check_location_in_exon(exons, chr_pos, gene)  # Added row['chr'] for chromosome
        CPIC_rsids.at[index, 'IsExonic'] = is_exonic
    
    CPIC_genes_filt = CPIC_rsids[CPIC_rsids['IsExonic']][['genesymbol', 'name_allele','dbsnpid_seq_loc','variantallele_allele_loc']]

    CPIC_genes_filt_nodups = CPIC_genes_filt.drop_duplicates()
    CPIC_rsids[CPIC_rsids['IsExonic']][['genesymbol', 'name_allele']].drop_duplicates().groupby('genesymbol')['name_allele'].apply(list)
    
    #grouped_data = CPIC_genes_filt_nodups.groupby('genesymbol')['name_allele'].apply(list)
    grouped_data=CPIC_rsids[CPIC_rsids['IsExonic']][['genesymbol', 'name_allele']].drop_duplicates().groupby('genesymbol')['name_allele'].apply(list)
    result_dict = grouped_data.to_dict()
    gene_dict={}
    for _, row in CPIC_genes_filt.iterrows():
        gene = row['genesymbol']
        dbsnpid = row['dbsnpid_seq_loc']
        name_allele = row['name_allele']
        alt_allele = row['variantallele_allele_loc']
        # If the gene is not already in the dictionary, add it with an empty dictionary as its value
        if gene not in gene_dict:
            gene_dict[gene] = []
        
    # Add the dbsnpid and name_allele to the gene's dictionary
        gene_dict[gene].append((name_allele,dbsnpid,alt_allele))

    return result_dict, gene_dict

def filter_genes(data_dict, gene_rsid_dict):
    filtered_dict = {}
    for gene, diplotypes in data_dict.items():
        if gene in gene_rsid_dict:
            valid_diplotypes = []
            for diplotype in diplotypes:
                alleles = diplotype.split('/')
                all_alleles_present = True
                for allele in alleles:                            
                    if not any(allele == entry for entry in gene_rsid_dict[gene]):
                        all_alleles_present = False
                        break
                if all_alleles_present:
                    valid_diplotypes.append(diplotype)
            if valid_diplotypes:
                filtered_dict[gene] = valid_diplotypes
    return filtered_dict

def main(diplo_query):
    clean_data()
    add_suffix()
    genewise_alleles_exonic_reportable, gene_rsid_allele_exonic = repor_alleles_exonic()
    filtered_dict = filter_genes(diplo_query, genewise_alleles_exonic_reportable)
    print(genewise_alleles_exonic_reportable)
    df_genes = get_gene_diplo_info(filtered_dict)
    '''lookup_key, df_genes_lookupkey = make_lookupdict(df_genes)  
    rec4query = get_reco(lookup_key, df_genes_lookupkey)
    fil_rec4query = filter_reportablereco(rec4query)
    rec4queryvthcategory = add_category(fil_rec4query)
    rec4queryvthcategory_publication = add_citation(rec4queryvthcategory)
    print(genewise_alleles_exonic_reportable)
    return rec4queryvthcategory_publication, df_genes, genewise_alleles_exonic_reportable'''

if __name__ == "__main__":
    data_dict = {
    'CYP2C19': ['*1/*2','*3/*8'],
    'CYP2D6': ['*1/*4'],
    'CYP3A4': ['*1/*1'],
    'CYP2C9': ['*1/*1'],
    'SLCO1B1': ['*14/*20']
    }
    main(data_dict)
    
