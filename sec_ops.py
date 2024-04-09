"""
mapping functions to transform one SecSeq to another
"""

class BlackBox:

    def __init__(self,input_shape,output_shape):
        self.input_shape = input_shape
        self.output_shape = output_shape
        return

    def pr_fit_function(self,input:list,output:list):
        return -1