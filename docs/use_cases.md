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

3. Detect if the input molecule is drug-like
	- *User* - Researchers and pharmacists
	- *Input* - Chemical SMILE string
	- *Function* - Passes the input through the model and evaluate to see if the molecule is drug-like
	- *Output* - Raise an error if the input molecule is not drug-like

4. Modify the model by adjusting the hyperparameters of the model  
	- *User* - Researchers, pharmacists, software engineers
	- *Input* - The model and new hyperparameters
	- *Function* - Change the hyperparameters of the model
	- *Output* - New model that meets the new expectations of the users

5. Predict similarity of natural product derived drugs to conventional drugs for a set application
    - *User* - Researchers
    - *Input* - List of SMILE strings of natural product derived drugs and the SMILE string of a target conventional drug
    - *Function* - Compare natural products to the conventional drug based on their AMDE/Tox and rank them based on prediction of fitness for application
    - *Output* - Ranking of the input natural product based on the suitability for the set application (which of the natural drugs is the best candidate)
