

## 1) General Information


#+++++ 1.1) Project info


#+++++ 1.2) Folders

./data folder contains the data that was used to train and test the ML models.
./notebooks contain the jupyter notebooks that were used to create the ML models. 
        - *Nanomet_processing_and_modeling.ipynb* : Contains the basic model (neural network) allong with the train testing.
        - *Nanomet_project_part2 (transfer learning).ipynb*: Contains a transfer learning approach. (this is not finished work)

#### 1.3) Training information
The forward model architect the was used:


```python
def get_model(n_inputs, n_outputs):
	net1Model=keras.models.Sequential()
	net1Model.add(Dense(n_inputs,input_shape=(4,)))
	net1Model.add(Dense(240))
	net1Model.add(Activation('relu'))
	net1Model.add(Dense(30))
	net1Model.add(Activation('relu'))
	net1Model.add(Dense(1024))
	net1Model.add(Dense(n_outputs))
	net1Model.compile(optimizer='Adam',loss='mse')
	return net1Model

model_forward.fit(X_train, y_train, verbose=0, epochs=500)
```
high number of epochs of course may lead to overfitting. We check in the end of the notebook if this is the case. (it looks like its not unless there is a leakage that we dont understand)


The inverse architect that was used:


```python
# get the model
def get_model(n_inputs, n_outputs):
    
	net1Model=keras.models.Sequential()
	net1Model.add(Dense(n_inputs))
	#net1Model.add(Activation('relu'))
	net1Model.add(Dense(240))
	net1Model.add(Activation('relu'))
	#net1Model.add(Dropout(0.5))
	net1Model.add(Dense(30))
	net1Model.add(Activation('relu'))
	net1Model.add(Dense(1024))
	net1Model.add(Dense(n_outputs))
	net1Model.compile(optimizer='Adam',loss='mse')
    
	return net1Model

model.fit(X_train, y_train, verbose=0, epochs=700)
```

## 2) Instalation

#### 2.1) Requirements

#### 2.2) How to install


## 3) Running examples

