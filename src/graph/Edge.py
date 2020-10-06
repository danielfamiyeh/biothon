class Edge:
    def __init__(self, endpoint, weight):
        self.endpoint = endpoint
        self.weight = weight

    def __repr__(self):
        return f"e:{self.weight}"
