class HypStruct:

    def __init__(self,suspected_subbounds):
        self.suspected_subbounds = suspected_subbounds
        self.subbound_scores = [None for i in range(len(self.suspected_subbounds))]
        return

