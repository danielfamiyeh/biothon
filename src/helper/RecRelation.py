class RecRelation:
    def __init__(self, m_from, m_to, offset_i, offset_j, offset_val):
        self.maps_from = m_from
        self.maps_to = m_to
        self.offset_i = offset_i
        self.offset_j = offset_j
        self.offset_val = offset_val

    def get_relation(self, i, j):
        return self.maps_to, i + self.offset_i, j + self.offset_j
