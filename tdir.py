"""
Was supposed to be named Traversal Directing [Wang Fong Qhong].
By and through the connection. 
By,through, and for the connection.
"""

class TDir:

    def __init__(self,loc):
        self.location = loc
        return

    """
    return: 
    - candidate destination nodeset
    - location of target iso-ring. 
    """
    def lookround(self,G):
        assert self.location in G

        return -1

    def __next__(self):
        return -1 
