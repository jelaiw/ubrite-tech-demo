import streamlit as st
import pandas as pd
import numpy as np
import requests
import base64
from io import StringIO

st.title('U-BRITE Tech Demo')
st.markdown('*Jelai Wang, Zongliang Yue, Abakash Samal, Dale Johnson, Patrick Dezenzio, Matt Wyatt, Christian Stackhouse, Lara Ianov, Chris Willey, Jake Chen*')

# Return GBM PDX cohort demographic data as a data frame.
def load_demographic_data():
	# Set up request parameters for UWS API call, see https://ubrite.slack.com/files/UAVTLGHT7/FLDADUC72/unified_ws_client.py for example from Abakash.
	params = {}
	params['requestorid'] = "rdalej"
	params['cohortid'] = "27676"
	params['format'] = "csv"
	response = requests.get('https://ubritedvapp1.hs.uab.edu/UbriteServices/getalli2b2demographics', headers={'eppn': 'jelaiw@uab.edu'}, params=params)

	# Remove first two lines of CSV-formatted response, see https://ubrite.slack.com/archives/CJTQDGE30/p1579191182004100 for context.
	lines = response.text.split('\n')
	relevant_csv_text = '\n'.join(lines[2:])
	# See https://stackoverflow.com/questions/20696479/pandas-read-csv-from-string-or-package-data.
	df = pd.read_csv(StringIO(relevant_csv_text))

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

st.header('Overview')
st.image('team_ubrite_interaction_diagram.png', use_column_width=True)

st.header('Query Clinical Data')
st.markdown("These data are read from U-BRITE's *Clinical Data Repository* programmatically and securely over the network via REST API call to the *Unified Web Services* (UWS) API, which returns cohort-related clinical data sets.")
st.markdown("See https://ubritedvapp1.hs.uab.edu/UbriteServices/services.html for API specification.")

clinical_data_load_state = st.text('Loading query results ... ')
clinical_data = load_demographic_data()
clinical_data_load_state.text('Loading query results ... done!')
st.write(clinical_data)

st.header('Parse DEG Results')
st.markdown("These results are from a differential gene expression (DEG) analysis performed with a custom DESeq2-based pipeline on RNAseq data located in the *Omics Data Repository*.")
st.markdown("See source code for pipeline in *Source Code Repository* at https://gitlab.rc.uab.edu/gbm-pdx/deseq2-rnaseq.")
# Somewhat kludgy section to optionally show source code for st.write() call.
show_code = st.checkbox('Show source code')
if show_code:
	with st.echo():
		deg_results = load_deg_results()
		st.write(deg_results)
else:	
	deg_results = load_deg_results()
	st.write(deg_results)

# Remove nan from gene list.
genes = [x for x in deg_results['symbol'].tolist() if str(x) != 'nan']

st.header('Run PAGER Analysis')
st.markdown("The list of significantly differentially expressed genes (DEG) is then passed to PAGER, which offers a network-accessible REST API for performing various gene-set, network, and pathway analyses.")

st.sidebar.subheader('Set PAGER Parameters')
sources = st.sidebar.multiselect('Available Data Sources',
	('KEGG', 'WikiPathway', 'BioCarta', 'MSigDB', 'Reactome', 'Spike'),
	('KEGG', 'WikiPathway')
)
fdr = st.sidebar.slider('FDR Cutoff', 0.0, 1.0, 0.05, 0.01)
pager_output = run_pager(genes, sources, fdr)

st.subheader('Dynamically View/Filter Results')
st.markdown('We invite you to interactively, *in real-time*, filter these results below. Also, you can set PAGER parameters in the sidebar')
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
# Per Jake's feedback, show a row count to make it easier to tell that results table is being dynamically updated as filtering is performed.
row_count, _ = filtered_output.shape
st.markdown('Row count = **{0}**'.format(row_count))

# See "File Download Workaround" in Gallery at http://awesome-streamlit.org/ for background reading.
st.subheader('Download Interactive Results')
st.markdown('The interactively end user-filtered PAGER results can now be downloaded for further review and post-processing.')
# See https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_csv.html.
csv = filtered_output.to_csv(index=False)
b64 = base64.b64encode(csv.encode()).decode()
href = '<a href="data:file/csv;base64,{0}">Download CSV File</a> (right-click and save as &lt;some_name&gt;.csv)'.format(b64)
st.markdown(href, unsafe_allow_html=True)

# See Zenodo record at http://doi.org/10.5281/zenodo.3700076.
st.header('Please Cite As')
st.markdown('*Jelai Wang, Zongliang Yue, Abakash Samal, Dale Johnson, Patrick Dezenzio, Matt Wyatt, ... Jake Chen*. (2020, January 29). Technical Demo from U-BRITE 2.0 Launch Day. Zenodo. http://doi.org/10.5281/zenodo.3700076')
# See markdown image link kludge at https://discuss.streamlit.io/t/is-it-possible-to-add-link/1027/3.
st.markdown('[![DOI badge](https://zenodo.org/badge/DOI/10.5281/zenodo.3700076.svg)](http://doi.org/10.5281/zenodo.3700076)')

st.header('References')
with open('references.md', 'r') as ref_file:
	st.markdown(ref_file.read())
