import streamlit as st
import pandas as pd
import numpy as np
import requests
import base64

st.title('U-BRITE Tech Demo')
st.markdown('*Jelai Wang, Zongliang Yue, Abakash Samal, Dale Johnson, Patrick Dezenzio, Matt Wyatt, Christian Stackhouse, Lara Ianov, Chris Willey, Jake Chen*')

# Return GBM PDX clinical data as a data frame.
def load_clinical_data():
	# Note these data are results from a UWS API query performed by Abakash for GBM cohort demographic data. See Nov 19, 2019 e-mail for further detail.
	df = pd.read_csv('getalli2b2demographics-rdalej-27676.csv')
	# Remove age field due to possible confusion created by -1 values, see https://gitlab.rc.uab.edu/jelaiw/infrastructure-development/issues/146#note_18590 for further detail and context.
	return df.drop(columns=['Age(in years)'])

# Return DEG results for JX12T pairwise comparison as data frame.
@st.cache
def load_deg_results():
	# See https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_table.html.
	# See https://gitlab.rc.uab.edu/jelaiw/infrastructure-development/issues/146#note_19095.
	df = pd.read_csv('JX12T_sig_DE_Results.txt', sep='\t')
	# Add sample name field to clarify that these data belong to sample JX12T, see https://gitlab.rc.uab.edu/jelaiw/infrastructure-development/issues/146.
	df['Sample Name'] = 'JX12T'
	# Reorder columns so that sample name field is first.
	# See https://stackoverflow.com/questions/13148429/how-to-change-the-order-of-dataframe-columns.
	df = df[['Sample Name', 'Unnamed: 0', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'symbol', 'ensembl', 'external_gene', 'gene_biotype', 'description', 'chromosome_name', 'start_position', 'end_position', 'strand']]
	return df

# Call PAGER REST API to perform hypergeometric test and return enriched PAGs associated with given list of genes as a data frame.
# See pathFun() in PAGER R SDK at https://uab.app.box.com/file/529139337869.
def run_pager(genes, sources, fdr):
	# Set up the call parameters as a dict.
	params = {}
	# Work around PAGER API form encode issue.
	params['genes'] = '%20'.join(genes)
	params['source'] = '%20'.join(sources)
	params['type'] = 'All'
	params['sim'] = '0.01'
	params['olap'] = '1'
	params['organism'] = 'All'
	params['cohesion'] = '0'
	params['pvalue'] = 0.05
	params['FDR'] = fdr
	params['ge'] = 1
	params['le'] = 2000

	response = requests.post('http://discovery.informatics.uab.edu/PAGER/index.php/geneset/pagerapi', data=params)
#	print(response.request.body)
	return pd.DataFrame(response.json())

st.header('Query Clinical Data')
st.markdown("These data are read from U-BRITE's *Clinical Data Repository* programmatically and securely over the network via REST API call to the *Unified Web Services* (UWS) API at http://ubritedvapp1.hs.uab.edu:8080/UbriteServices/getalli2b2demographics?requestorid=rdalej&cohortid=27676&format=csv, which returns de-identified results.")

clinical_data_load_state = st.text('Loading query results ... ')
clinical_data = load_clinical_data()
clinical_data_load_state.text('Loading query results ... done!')
st.write(clinical_data)

st.header('Parse DEG Results')
st.markdown("This is from a differential gene expression analysis performed with a custom DESeq2-based pipeline (see source code for pipeline in access-controlled *Source Code Repository* at https://gitlab.rc.uab.edu/gbm-pdx/deseq2-rnaseq) on RNAseq data located in the *Omics Data Repository* (in GBM science team's dedicated project subdirectory).")
st.markdown("These particular results are for the JX12T pairwise comparison and read from the *Documentation Repository* location at https://uab.app.box.com/folder/63730723635, a shared Box folder where the GBM analysis team shares results, slides, and other works-in-progress.")
deg_results = load_deg_results()
if st.checkbox('Show DEG results table', value=True):
	st.write(deg_results)

# Remove nan from gene list.
genes = [x for x in deg_results['symbol'].tolist() if str(x) != 'nan']

st.header('Run PAGER Analysis')
st.markdown("The list of significantly differentially expressed genes is then passed to PAGER, which offers a network-accessible REST API for performing various gene-set, network, and pathway analyses.")

st.sidebar.subheader('Adjust PAGER Parameters')
sources = st.sidebar.multiselect('Available Data Sources',
	('KEGG', 'WikiPathway', 'BioCarta', 'MSigDB', 'Reactome', 'Spike'),
	('KEGG', 'WikiPathway')
)
fdr = st.sidebar.slider('FDR Cutoff', 0.0, 1.0, 0.05, 0.01)
pager_output = run_pager(genes, sources, fdr)

st.subheader('View/Filter Results')
# Convert GS_SIZE column from object to integer dtype.
# See https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.astype.html.
pager_output = pager_output.astype({'GS_SIZE': 'int32'})
gs_sizes = pager_output['GS_SIZE'].tolist()
# Figure out the min and max GS_SIZE within the PAGER output.
min_gs_size = min(gs_sizes)
max_gs_size = max(gs_sizes)
# Set up a range slider. Cool!
# See https://streamlit.io/docs/api.html#streamlit.slider.
user_min, user_max = st.slider('GS_SIZE Range', max_value=max_gs_size, value=(min_gs_size, max_gs_size))
filtered_output = pager_output[pager_output['GS_SIZE'].between(user_min, user_max)]
st.write(filtered_output)

# See "File Download Workaround" in Gallery at http://awesome-streamlit.org/ for background reading.
st.subheader('Download Results')
st.markdown('The filtered PAGER results can be now be downloaded for further review and post-processing.')
st.markdown('Note this is an unofficial file download workaround from the Awesome Streamlit Gallery at http://awesome-streamlit.org, implemented here to share a neat idea.')
# See https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_csv.html.
csv = filtered_output.to_csv(index=False)
b64 = base64.b64encode(csv.encode()).decode()
href = f'<a href="data:file/csv;base64,{b64}">Download CSV File</a> (right-click and save as &lt;some_name&gt;.csv)'
st.markdown(href, unsafe_allow_html=True)
