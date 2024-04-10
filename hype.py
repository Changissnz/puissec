"""
The Hypothesis structure.

uses vector logs for 3rd-party scores
S with the `suspected_subbounds`, to 
produce an output
P(s in S) --> is s relevant?
"""
class HypStruct:

    def __init__(self,suspected_subbounds):
        self.suspected_subbounds = suspected_subbounds
        # every i'th element is vector of scores 
        # corresponding to the i'th element in 
        # the sequence `suspected_subbounds`
        self.subbound_scorevec = [None for i in range(len(self.suspected_subbounds))]
        return

    def mark_subbound(self,i):
        return -1

    def log_score(self):
        return -1

