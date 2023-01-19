

## 1) General Information


### +++++ 1.1) Project info


### 1.2) Folders

* ./data folder contains the data that was used to train and test the ML models.
* ./notebooks contain the jupyter notebooks that were used to create the ML models. 
        - *Nanomet_processing_and_modeling.ipynb* : Contains the basic model (neural network) allong with the train testing.
        - *Nanomet_project_part2 (transfer learning).ipynb*: Contains a transfer learning approach. (this is not finished work)

* ./models cotnain the neural network keras models (forward and ivnerse). We have also the tflite versions of them.
* ./minmaxscalers contain the scripts for transforming the input (for forward) and output (for inverse) of the models.
#### 1.3) Training information
The forward model architect the was used:


```python
def get_model(n_inputs, n_outputs):
	net1Model=keras.models.Sequential()
	net1Model.add(Dense(n_inputs,input_shape=(4,)))
	net1Model.add(Dense(240))
	net1Model.add(LeakyReLU(alpha=0.05))
	net1Model.add(Dense(30))
	net1Model.add(LeakyReLU(alpha=0.05))
	net1Model.add(Dense(1024))
	net1Model.add(Dense(n_outputs))
	net1Model.compile(optimizer='Adam',loss='mse')
	return net1Model

model_forward.fit(X_train, y_train, verbose=0, epochs=900)
```
high number of epochs of course may lead to overfitting. We check in the end of the notebook if this is the case. (it looks like its not unless there is a leakage that we dont understand)


The inverse architect that was used:


```python
# get the model
def get_model(n_inputs, n_outputs):
    
	net1Model.add(Dense(n_inputs))
	net1Model.add(Dense(1024))
	net1Model.add(LeakyReLU(alpha=0.05))
	net1Model.add(Dense(240))
	net1Model.add(LeakyReLU(alpha=0.05))
	net1Model.add(Dense(30))
	net1Model.add(LeakyReLU(alpha=0.05))
    
	return net1Model

model.fit(X_train, y_train, verbose=0, epochs=900)
```

## 2) Instalation

requirements and specifications: python 3.8+ (these steps have been checked in win11 os)

Instalation steps:

    1)Download and install python 3.8+
    2)Clone this repo
    3)Inside the repo run `pip install -r requirments.txt`
    4)Use enviroment for notebooks (anaconda,vs studio etc.)
	5)Open Example.ipynb to run the examples

## 3) Running examples

Example.ipynb contains the examples for running examples for the forward and inverse case.

For the forward case it uses the keras neural network model_forward_v4 while for the inverse, it uses the model_inverse_v3. Both models exist in ./models folder.

Note that in Example.ipynb there is the step of minmaxscaling. This is required in order for the neural network models learn better the values.

