import numpy
from utils.features import prepare_for_training
from utils.hypothesis import sigmoid, sigmoid_gradient
# from sklearn.neural_network import multilayer_preceptron

class MultilayerPerceptron:
    def __init__(self, data, labels, layers,normalize_data=False):
        data_processed = prepare_for_training(data, normalize_data)
        self.data = data_processed
        self.labels = labels
        self.layers = layers # 784 25 10
        self.normalize_data = normalize_data
        self.thetas = MultilayerPerceptron.thetas_init(layers)
    def thetas_init(layers):
        num_layers = len(layers)
        thetas = []
        for layer_index in range(num_layers - 1):
            """
            会执行两次，得到两组参数局矩阵：25*785, 10*26
            """
            in_count = layers[layer_index]
            out_count = layers[layer_index + 1]
def input():

def forward():
    pass