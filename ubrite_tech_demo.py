import streamlit as st
import pandas as pd
import requests

st.title('U-BRITE Tech Demo')

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
	return pd.read_table('JX12T_sig_DE_Results.txt')

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

st.header('Clinical data')
st.markdown("These data are read from U-BRITE's *Clinical Data Repository* programmatically from a secure call the *Unified Web Services* (UWS) API at http://ubritedvapp1.hs.uab.edu:8080/UbriteServices/getalli2b2demographics?requestorid=rdalej&cohortid=27676&format=csv.")

clinical_data_load_state = st.text('Loading data ... ')
clinical_data = load_clinical_data()
clinical_data_load_state.text('Loading data ... done!')
st.write(clinical_data)

st.header('DEG results')
st.markdown("This is from a differential gene expression analysis performed with a custom DESeq2-based pipeline (see source code for pipeline in access-controlled *Source Code Repository* at https://gitlab.rc.uab.edu/gbm-pdx/deseq2-rnaseq) on RNAseq data located in the *Omics Data Repository* (in GBM science team's dedicated project subdirectory).")
st.markdown("These particular results are for the JX12T pairwise comparison and read from the *Documentation Repository* location at https://uab.app.box.com/folder/63730723635, a shared Box folder where the GBM analysis team shares results, slides, and other works-in-progress.")
deg_results = load_deg_results()
if st.checkbox('Show DEG results table', value=True):
	st.write(deg_results)

# Remove nan from gene list.
genes = [x for x in deg_results['symbol'].tolist() if str(x) != 'nan']

st.header('PAGER results')
st.markdown("The list of significantly differentially expressed genes is then passed to PAGER, which offers a network-accessible REST API for performing various gene-set, network, and pathway analyses.")

pager_run_state = st.text('Calling PAGER REST API ... ')
sources = st.multiselect('Data sources',
	('KEGG', 'WikiPathway', 'BioCarta', 'MSigDB', 'Reactome', 'Spike'),
	('KEGG', 'WikiPathway')
)
fdr = st.slider('FDR cutoff', 0.0, 1.0, 0.05, 0.01)
pager_output = run_pager(genes, sources, fdr)
pager_run_state.text('Calling PAGER REST API ... done!')

st.write(pager_output)
