import tensorflow as tf
from tensorflow.keras.models import load_model
from tensorflow.keras.callbacks import ModelCheckpoint
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, LeakyReLU, ReLU
import tensorflow as tf
import os

def shapeData(data):
    return data.reshape(len(data), len(data[0]) * len(data[0][0]))

class NN:
    def __init__(self, modelFile=""):
        self.model = None
        if modelFile != "":
            self.loadModel(modelFile)

    def createModel(self, type = "Laplace"):
        tf.keras.backend.clear_session()
        self.model = Sequential()
        if type == "Laplace":
            self.model.add(Dense(512, input_dim=500))
            self.model.add(LeakyReLU(alpha=0.3))
            self.model.add(Dense(256))
            self.model.add(LeakyReLU(alpha=0.3))
            self.model.add(Dense(128))
            self.model.add(LeakyReLU(alpha=0.3))
            self.model.add(Dense(32))
            self.model.add(LeakyReLU(alpha=0.3))
            self.model.add(Dense(2))
            self.model.add(ReLU())
        elif type == "Hooke":
            self.model.add(Dense(1024, input_dim=1000, activation="sigmoid"))
            self.model.add(Dense(512))
            self.model.add(LeakyReLU(alpha=0.3))
            self.model.add(Dense(256))
            self.model.add(LeakyReLU(alpha=0.3))
            self.model.add(Dense(16, activation="linear"))
            self.model.add(Dense(2, activation="sigmoid"))
        elif type == "HookeLarge":
            self.model.add(Dense(1024, input_dim=1000, activation="sigmoid"))
            self.model.add(Dense(512))
            self.model.add(LeakyReLU(alpha=0.3))
            self.model.add(Dense(256))
            self.model.add(LeakyReLU(alpha=0.3))
            self.model.add(Dense(16, activation="linear"))
            self.model.add(Dense(3, activation="sigmoid"))

        self.model.compile(optimizer='Adadelta', loss="mse")
        self.model.summary()
        print("Generated new model!")



    def loadModel(self, file):
        if os.path.isfile(file):
            tf.keras.backend.clear_session()
            self.model = load_model(file)
            return True
        else:
            print("Warning: Model file does not exist!")
            return False

    def predict(self, data):
        if len(data) == 0 or self.model is None:
            print("Error: Input missing or model not loaded!")
            exit(1)

        return self.model.predict(data)

    def __getCallbacks(self, saveFile, saveBest=True):
        save_callback = ModelCheckpoint(saveFile, monitor='val_loss',
                                                          verbose=0, save_best_only=saveBest,
                                                          save_weights_only=False, mode='auto', period=1)
        return [save_callback]

    def train(self, epochs, batch_size, saveFile, data, labels):
        data_train, labels_train, data_eval, labels_eval = data[:int(0.9*len(data))], labels[:int(0.9*len(labels))], data[int(0.9*len(data)):], labels[int(0.9*len(labels)):]

        self.model.fit(x=data_train, y=labels_train, epochs=epochs, batch_size=batch_size,
                        validation_data=(data_eval, labels_eval), callbacks=self.__getCallbacks(saveFile))
