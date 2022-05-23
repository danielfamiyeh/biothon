class RecurrenceRelation:
    def __init__(self, **kwargs):
        self.maps_from = kwargs.get("maps_from")
        self.i_offset = kwargs.get("i", 0)
        self.j_offset = kwargs.get("j", 0)
        self.score = kwargs.get("score", "x")

    def get_indices(self, i, j):
        if self.score == "x":
            return 0
        return i + self.i_offset, j + self.j_offset
