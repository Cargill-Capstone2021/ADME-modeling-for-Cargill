## Use cases
1. Predict ADME/Tox of natural product derived drugs using their SMILE strings 
	- *User* - Researchers and pharmacists
	- *Input* - SMILE string of the natural product derived drugs 
	- *Function* - Pass the SMILE string to machine learning model and perform the estimates
	- *Output* - Facets of the input pharmaceutical lead in the ADME/Tox aspect

2. Detect if the input data is a chemical SMILE string
	- *User* - Researchers
	- *Input* - Chemical SMILE string
	- *Function* - Passes the input through the model and evaluate it
	- *Output* - Raises an error if the input is not SMILE string

3. Modify the model by adjusting the hyperparameters of the model  
	- *User* - Researchers, pharmacists, software engineers
	- *Input* - The model and new hyperparameters
	- *Function* - Change the hyperparameters of the model
	- *Output* - New model that meets the new expectations of the users
