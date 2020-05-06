# technology_data

sample and convert the technology data for pypsa-eur and pypsa-eur-sec

## structure
this repository has the following structure: 

* **inputs**: input data in format .csv or .xlsx
* **outputs**: technology data for technology saved as .csv format for defined years
* **docu**: reports, paper, additional information about the input data, format .pdf
* **scripts** : 
  * compile_cost_assumptions.py converts input data to wished output data
  * convert_pdf_frauenhofer_to_dataframe.py converts table in pdf to csv format for input data
  * retrieve_data_from_dea.py takes up to date data from DEA and saves it in input, as well as documentation of the                 data in docu
  
  ## sources
           
