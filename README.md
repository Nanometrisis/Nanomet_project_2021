The code is a bit messy. It is not for public to see hehe.


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

