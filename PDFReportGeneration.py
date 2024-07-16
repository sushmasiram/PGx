import fetchcpicdata as fd
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Image
from reportlab.lib import colors
import json
import pandas as pd

data_dict_table1 = {
"Cardiovascular Diseases": {
    "Acute Coronary Syndrome": ["Clopidogrel"],
    "Arrhythmias, Cardiac": ["Propafenone"],
    "Atrial Fibrillation": ["Warfarin"],
    "Cardiotoxicity": ["Venlafaxine", "Tramadol"],
    "Coronary Artery Disease": ["Clopidogrel"],
    "Coronary Disease": ["Rosuvastatin"],
    "Heart Diseases": ["Carvedilol"],
    "Heart Failure": ["Metoprolol", "Carvedilol"],
    "Hypertension": ["Metoprolol"],
    "Myocardial Infarction": ["Clopidogrel"],
    "Tachycardia": ["Venlafaxine", "Propafenone"],
    "Thrombocytopenia": ["Thioguanine", "Mercaptopurine"],
    "Cardiovascular Diseases": ["Warfarin", "Clopidogrel"],
    "Heart valve replacement": ["Warfarin"],
    "Statin-related myopathy": ["Rosuvastatin"]
},
"Cancer and Neoplasms":{
    "Adenocarcinoma": ["Gefitinib"],
    "Breast Neoplasms": ["Tamoxifen"],
    "Carcinoma, Non-Small-Cell Lung": ["Gefitinib"],
    "Carcinoma, Renal Cell": ["Pazopanib"],
    "Colorectal Neoplasms": ["Capecitabine"],
    "Leukopenia": ["Mercaptopurine", "Azathioprine"],
    "Lymphoma": ["Irinotecan"],
    "Neoplasms": ["Capecitabine", "Irinotecan", "Fluorouracil"],
    "Precursor Cell Lymphoblastic Leukemia-Lymphoma": ["Thioguanine", "Mercaptopurine"]
},
"Mental Health Disorders": {
        "Agitation": ["Venlafaxine"],
        "Alzheimer Disease": ["Donepezil"],
        "Attention Deficit Disorder with Hyperactivity": ["Atomoxetine"],
        "Depression": ["Venlafaxine"],
        "Depressive Disorder": ["Venlafaxine", "Imipramine", "Citalopram", "Amitriptyline", "Fluvoxamine", "Nortriptyline", "Desipramine"],
        "Depressive Disorder, Major": ["Fluoxetine", "Escitalopram", "Venlafaxine", "Imipramine", "Citalopram", "Nortriptyline", "Desipramine", "Clomipramine"],
        "Obsessive-Compulsive Disorder": ["Venlafaxine", "Citalopram"],
        "Psychotic Disorders": ["Aripiprazole", "Risperidone"],
        "Schizoaffective Disorder": ["Aripiprazole"],
        "Schizophrenia": ["Aripiprazole", "Iloperidone", "Risperidone"],
        "Dysphoria": ["Venlafaxine"]
    },
    "Substance-Related Disorders": {
        "Alcohol-Related Disorders": ["Venlafaxine"],
        "Opioid-Related Disorders": ["Codeine"],
        "Ototoxicity": ["Gentamicin", "Amikacin", "Tobramycin", "Streptomycin"]
    },
    "Infectious Diseases": {
        "HIV Infections": ["Abacavir", "Efavirenz"]
    },
    

    "Metabolic Disorders": {
        "Hyperbilirubinemia": ["Nilotinib", "Capecitabine"],
        "Hypercholesterolemia": ["Rosuvastatin"],
        "Methemoglobinemia": ["Rasburicase", "Pegloticase"],
        "Protein Deficiency": ["Methylene Blue"]
    },
    "Pain and Symptom Management": {
        "Pain": ["Tramadol", "Codeine"],
        "Pain, Postoperative": ["Tramadol"],
        "Nausea": ["Venlafaxine"],
        "Vomiting": ["Ondansetron", "Venlafaxine"]
    },
    "Medical Procedures and Complications": {
        "Transplantation": ["Azathioprine"],
        "Postanesthesia apnea": ["Succinylcholine"],
        "Myelosuppression": ["Azathioprine"],
        "Malignant Hyperthermia": ["Sevoflurane", "Desflurane", "Isoflurane", "Succinylcholine"]
    },
    "Anticoagulation Therapy": {
        "Time in Therapeutic Range / Time to Therapeutic INR": ["Warfarin"],
        "Over-anticoagulation": ["Warfarin"]
    },
    "Hematological Disorders": {
        "Neutropenia": ["Thioguanine", "Mercaptopurine", "Irinotecan"]
    },
    "Diabetes": {
        "Diabetes mellitus 2": ["Glyburide"]
    },
    "Drug Reactions and Toxicities": {
        "Adverse Events": ["Sulfamethoxazole", "Trimethoprim"],
        "Drug Hypersensitivity": ["Abacavir", "Carbamazepine", "Allopurinol"],
        "Drug Reaction with Eosinophilia and Systemic Symptoms": ["Phenytoin", "Carbamazepine", "Allopurinol"],
        "Drug Toxicity": ["Phenytoin", "Fluorouracil", "Venlafaxine", "Gefitinib"],
        "Erythema Multiforme": ["Carbamazepine"],
        "Epidermal Necrolysis, Toxic": ["Phenytoin", "Carbamazepine", "Allopurinol"],
        "Maculopapular Exanthema": ["Carbamazepine", "Oxcarbazepine"],
        "Stevens-Johnson Syndrome": ["Phenytoin", "Carbamazepine", "Allopurinol", "Oxcarbazepine"],
        "Severe Cutaneous Adverse Reactions": ["Sulfamethoxazole", "Trimethoprim", "Allopurinol", "Carbamazepine", "Phenytoin"],
        "Dose reduction": ["Mercaptopurine", "Azathioprine"],
        "Exanthema": ["Gefitinib"]
    },
    "Digestive System Disorders": {
        "Crohn Disease": ["Thioguanine"],
        "Esophagitis": ["Omeprazole"],
        "Gastroesophageal Reflux": ["Esomeprazole", "Omeprazole", "Rabeprazole"],
        "Peptic Ulcer": ["Omeprazole", "Rabeprazole"],
        "Ulcer": ["Pantoprazole"],
        "Toxic liver disease": ["Valproic Acid", "Gefitinib"],
        "Inflammatory Bowel Diseases": ["Azathioprine"]
    },
    "Respiratory Disorders": {
        "Apnea": ["Succinylcholine"],
        "Cystic Fibrosis": ["Ivacaftor"],
        "Respiratory Insufficiency": ["Tramadol"],
        "Tuberculosis": ["Efavirenz"]
    },
    "Immune System Disorders": {
        "Anemia, Hemolytic": ["Sulfamethoxazole", "Trimethoprim", "Chloroquine", "Nitrofurantoin", "Methylene Blue"],
        "Hypersensitivity": ["Sulfamethoxazole", "Trimethoprim"],
        "Lupus Erythematosus, Systemic": ["Mercaptopurine"]
    },
    "Neurological Disorders": {
        "Dementia": ["Galantamine"],
        "Epilepsy": ["Phenytoin", "Clobazam", "Oxcarbazepine"],
        "Ischemic Attack, Transient": ["Clopidogrel"],
        "Stroke": ["Clopidogrel"],
        "Mental Disorders": ["Fluoxetine", "Escitalopram", "Imipramine", "Amitriptyline", "Trimipramine", "Desipramine", "Clomipramine"]
    }
}

tables_sect1={'Page1' : ['Cardiovascular Diseases', 'Cancer and Neoplasms'],
'Page2': ['Mental Health Disorders'],
'Page3': ['Substance-Related Disorders','Infectious Diseases','Metabolic Disorders','Pain and Symptom Management','Medical Procedures and Complications','Anticoagulation Therapy','Hematological Disorders', 'Diabetes'],
'Page4': ['Drug Reactions and Toxicities'],
'Page5': ['Digestive System Disorders','Respiratory Disorders'],
'Page6':['Immune System Disorders','Neurological Disorders']}
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, PageBreak
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet


def add_spans_for_repeated_words(data, column_index):
    spans = []
    current_start = None
    current_word = None
  
    for i, row in enumerate(data):
        word = row[column_index]
        
        if word == current_word:
            if current_start is None:
                current_start = i - 1
        else:
            if current_start is not None:
                spans.append((current_start, i - 1))
                current_start = None
            current_word = word
    
    if current_start is not None:
        spans.append((current_start, len(data) - 1))
    return spans

def generate_table_data(data, pageNo, drug_dict):
    styles = getSampleStyleSheet()    
    table_data_header = [
        Paragraph("CATEGORY", styles['Normal']),
        Paragraph("DISEASE/ PHENOTYPE", styles['Normal']),
        Paragraph("DRUGS", styles['Normal']),
        Paragraph("USE WITH CAUTION", styles['Normal']),
        Paragraph("CONSIDER ALTERNATIVES", styles['Normal'])
    ]
    
   
    use_with_caution_image_path = '../Report/images/caution.png'
    consider_alternatives_image_path = '../Report/images/alternate.png'
    
    use_with_caution = Image(use_with_caution_image_path)
    use_with_caution.drawHeight = 0.25 * inch
    use_with_caution.drawWidth = 0.25 * inch
    
    consider_alternatives = Image(consider_alternatives_image_path)
    consider_alternatives.drawHeight = 0.25 * inch
    consider_alternatives.drawWidth = 0.25 * inch
     
    '''table_data_header0 = [
        Paragraph("", styles['Normal']),
        Paragraph("", styles['Normal']),
        Paragraph("", styles['Normal']),
        Paragraph(use_with_caution, styles['Normal']),
        Paragraph(consider_alternatives, styles['Normal'])
    ]'''
    # Adjust the size of the images if necessary
  
    
    table_data = [table_data_header]
    table_data_span = [table_data_header]
    
    for major_condition in tables_sect1[pageNo]:
        diseases = data.get(major_condition, {})
        
        for disease_name, drugs in diseases.items():
            for drug_n in drugs:                
                row_1 = [major_condition, disease_name, drug_n]
                category = Paragraph(major_condition, styles['Normal'])
                disease = Paragraph(disease_name, styles['Normal'])
                drug = Paragraph(drug_n, styles['Normal'])
                use_with_caution = Paragraph('', styles['Normal'])
                consider_alternatives = Paragraph('', styles['Normal'])        
                if drug_n.lower() in drug_dict:                   
                    for i in drug_dict[drug_n.lower()]:                        
                        if i == 'Use Alternative Drug':                            
                            consider_alternatives = Image(consider_alternatives_image_path)
                            consider_alternatives.drawHeight = 0.25 * inch
                            consider_alternatives.drawWidth = 0.25 * inch   
                        elif i == 'Use with Caution':
                            use_with_caution  = Image(use_with_caution_image_path)
                            use_with_caution.drawHeight = 0.25 * inch
                            use_with_caution.drawWidth = 0.25 * inch
                
                
                row = [category, disease, drug, use_with_caution, consider_alternatives]
                table_data.append(row)
                table_data_span.append(row_1)    
    
    return table_data, table_data_span

def trade_names(drug_dict):
    #drug_df = pd.read_csv('../DBs/CPIC/drug.csv')
    chemicals_df = pd.read_table('../DBs/PharmGKB/chemicals.tsv')
    drugs = set(list(drug_dict.keys()))
    tradename_dict = {}
    for each in drugs:
        matches = chemicals_df['Trade Names'][(chemicals_df["Name"].str.lower() == each.lower()) & (chemicals_df['Trade Names'].notna())].to_list()
        
        if len(matches)>0:
            tradename_dict[each.upper()]=matches[0].replace('"','')
    return tradename_dict

def section_tradenames(tradename_dict):
    styles = getSampleStyleSheet() 
    tn_hd = ['Drug', 'Tradenames']
    data = [tn_hd] + [[Paragraph(drug,styles['Normal']),Paragraph(tradename_dict[drug],styles['Normal'])] for drug in tradename_dict]
    
    # Create the Table object with the prepared data
    table = Table(data)
    style = TableStyle([
    ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
    ('VALIGN', (0, 0), (-1, -1), 'TOP'),
    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
    ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
    ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
    ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ])
    table.setStyle(style)
    return table 

def rgb_to_hex(color):
    return "#{:02x}{:02x}{:02x}".format(int(color.red * 255), int(color.green * 255), int(color.blue * 255))

def create_recommendation_table(sec3_df):
    styles = getSampleStyleSheet()
    style_normal = styles['Normal']
    table_data = []

    for index, row in sec3_df.iterrows():
        if row['category'] == 'Use with Caution':
            image_var = Image('../Report/images/caution.png', width=15, height=15)
            text_color = colors.orange
        elif row['category'] == 'Use Alternative Drug':
            image_var = Image('../Report/images/alternate.png', width=15, height=15)
            text_color = colors.red
        else:
            image_var = None
            text_color = colors.black

        text_color_hex = rgb_to_hex(text_color)

        row_data = [
            image_var,
            Paragraph(f"<font color='{text_color_hex}'>{row['name_drug']}<br/><font size='8' color='black'>{row['population_recommendation']}</font></font>", style_normal),
            Paragraph(f"<font color='{text_color_hex}'> {row['phenotypes_recommendation1']} </font>", style_normal),
            Paragraph(f"<font color='{text_color_hex}'> {row['pgxtesting']} </font>", style_normal)    
        ]
        table_data.append(row_data)

        recommendation_row = [
            '',
            '',
            Paragraph(row['drugrecommendation_recommendation'], style_normal),
            ''
        ]
        table_data.append(recommendation_row)
        
        citation_row = [
            '',
            '',
            Paragraph(row['pmid_citation'], style_normal),
            ''
        ]

        table_data.append(citation_row)
    
    table = Table(table_data)
    style = TableStyle([
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
    ])
    table.setStyle(style)
    return table

def draw_image_and_heading(canvas, doc, image_path, heading_text):
    width, height = letter

    # Add image
    image = Image(image_path, width=2*inch, height=2*inch)  # Adjust size as needed
    image.drawOn(canvas, (width - 2*inch) / 2, height - 2*inch - inch)

    # Add heading
    styles = getSampleStyleSheet()
    heading = Paragraph(heading_text, styles['Title'])
    heading_width, heading_height = heading.wrap(width, height)
    heading.drawOn(canvas, (width - heading_width) / 2, height - 2*inch - inch - heading_height - 10)

    # Add red line
    canvas.setStrokeColor(colors.red)
    canvas.setLineWidth(2)
    canvas.line(inch, height - 2*inch - inch - heading_height - 20, width - inch, height - 2*inch - inch - heading_height - 20)


def create_genotype_table(df_genes):
    table_data = [['Gene Symbol', 'Diplotype', 'Phenotype']]
    table_data.extend(df_genes[['genesymbol', 'diplotype', 'result']].values.tolist())
    table = Table(table_data)
    style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
    ])
    table.setStyle(style)
    return table

def create_reportablealleles_table(df_genes, genes_alleles_dict):
    styles = getSampleStyleSheet()
    normal_style = styles['Normal']
    normal_style.fontSize = 8  # Decrease font size

    table_data = [['Gene Symbol', 'Diplotype', 'Phenotype', 'Alleles tested']]

    for gene, alleles in genes_alleles_dict.items():
    
        if gene in df_genes['genesymbol'].to_list():
            found = df_genes[['genesymbol', 'diplotype', 'result']][df_genes['genesymbol'] == gene]
            found = found.reset_index(drop=True)
            
            for index, row in found.iterrows():
                
                diplotype, result = row['diplotype'], row['result']
                if index == 0:
                    table_data.append([gene, diplotype, result, Paragraph(', '.join(alleles), normal_style)])
                else:
                    table_data.append([gene, diplotype, result, Paragraph('', normal_style)])
        else:
            diplotype, result = '', ''
            table_data.append([gene, diplotype, result, Paragraph(', '.join(alleles), normal_style)])

    table = Table(table_data)
    
    style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('VALIGN', (3, 0), (3, -1), 'TOP'),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 6),  # Adjust padding
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.black),  # Adjust grid thickness
        ('FONTSIZE', (0, 0), (-1, -1), 8)  # Set font size
    ])
    table.setStyle(style)
    table.splitByRow = True  # Allow table to split across pages
    return(table)

def create_genotype_table1(df_genes, genes_alleles_dict):
    table_data = [['Gene Symbol', 'Diplotype', 'Phenotype','Alleles tested']]
    table_data.extend(df_genes[['genesymbol', 'diplotype', 'result']].values.tolist())

    table = Table(table_data)
    style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
    ])
    table.setStyle(style)
    return table

def create_pdf(file_path, data, drug_dict,sec3_df, df_genes, genes_allele_dict):
    doc = SimpleDocTemplate(file_path, pagesize=letter)
    elements = []
    for pageNo in tables_sect1:
        table_data,table_data_span = generate_table_data(data, pageNo, drug_dict)   
        # Calculate row heights based on table data length    
        col_widths = [1.5 * inch, 2 * inch, 1.5 * inch, 1 * inch, 1.3 * inch]
        table = Table(table_data, colWidths=col_widths)
        spans_diseases = add_spans_for_repeated_words(table_data_span, 1)
        spans_major = add_spans_for_repeated_words(table_data_span, 0)
        style = TableStyle([
            ('BACKGROUND', (4, 0), (4, -1), colors.HexColor('#FFD5D5')),
            ('BACKGROUND', (3, 0), (3, -1), colors.HexColor('#FFF8E1')),
            ('BACKGROUND', (0, 0), (0, -1), colors.HexColor('#F2F2F2')),
            ('BACKGROUND', (1, 0), (1, -1), colors.HexColor('#E7E6E6')),
            ('BACKGROUND', (2, 0), (2, -1), colors.HexColor('#D9D9D9')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.black),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('GRID', (0, 0), (-1, -1), 0.05, colors.black),
            ('GRID', (0, 0), (-1, -1), 1, colors.HexColor('#A5A5A5')), 
        ])
        
        for start, end in spans_diseases:
            style.add('SPAN', (1, start), (1, end))
        for start, end in spans_major:
            style.add('SPAN', (0, start), (0, end))

        table.setStyle(style)
        elements.append(table)        
        
        elements.append(PageBreak())
    
    elements.append(create_genotype_table(df_genes))
    elements.append(PageBreak())
    
    tradenames_dict = trade_names(drug_dict)
    
    table_tradenames = section_tradenames(tradenames_dict)
    elements.append(table_tradenames)
    elements.append(PageBreak())  
    
    elements.append(create_recommendation_table(sec3_df))    
    elements.append(PageBreak())
    elements.append(create_reportablealleles_table(df_genes, genes_allele_dict))
    doc.build(elements)

def get_data_cpic(data_dict):
    rec_data, gene_data, genewise_alleles_exonic_reportable = fd.main(data_dict)
    sec3_df = rec_data[['name_drug', 'drugrecommendation_recommendation', 'phenotypes_recommendation', 'category', 'population_recommendation','pgxtesting','pmid_citation']][rec_data['category'].isin(['Use with Caution', 'Use Alternative Drug'])]
    sec3_df['phenotypes_recommendation1'] = sec3_df['phenotypes_recommendation'].apply(make_string)
    return sec3_df, gene_data, genewise_alleles_exonic_reportable

def make_string(row):
    row_dict = json.loads(row)
    return ' and '.join([f"{key}: {value}" for key, value in row_dict.items()])

def main():
    #rec_fil = {"category": ["Use with Caution", "Use Alternative Drug"]}  # Example filter, adjust as needed
    data_dict = {
    'CYP2C19': ['*1/*2','*3/*8'],
    'CYP2D6': ['*1/*4'],
    'CYP3A4': ['*1/*1'],
    'CYP2C9': ['*1/*1'],
    'SLCO1B1': ['*14/*20']
    }
    sec3_df, df_genes, genes_alleles_dict = get_data_cpic(data_dict)
    
    sec1_df = sec3_df[['name_drug','category']].drop_duplicates()
    grouped_df = sec1_df.groupby('name_drug')['category'].apply(lambda x: list(set(x))).reset_index()
    drug_category_dict = dict(zip(grouped_df['name_drug'], grouped_df['category']))
    create_pdf("output.pdf", data_dict_table1, drug_category_dict,sec3_df, df_genes, genes_alleles_dict)

if __name__ == "__main__":
    main()