from collections import defaultdict 

class SecSeq:

    def __init__(self,sequence,dependency_map=None,codependency_map=None,reference_map=None): 
        self.sequence = sequence
        self.dependency_map = dependency_map
        self.codependency_map = codependency_map
        self.reference_map = reference_map
        self.check_args()
        return

    def check_args(self):
        assert type(self.sequence) == list
        assert self.of_type() != -1

    def of_type(self):
        stat1 = type(self.dependency_map) == defaultdict 
        stat2 = type(self.codependency_map) == defaultdict
        stat3 = type(self.reference_map) == defaultdict
        
        stat4 = (stat1 and stat2) and not stat3
        stat5 = stat3 and (not stat1 and not stat2)

        if stat4:
            return "enc"
        if stat5:
            return "dec"
        return -1