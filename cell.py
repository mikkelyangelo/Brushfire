from math import sqrt


class Cell:
    def __init__(self, x, y, end):
        self.x = x
        self.y = y

        self.distance = 0
        self.capability = 0

        self.evcl_distance = sqrt((x - end[0]) ** 2 + (y - end[1]) ** 2)

    def __repr__(self):
        return str((self.x, self.y))
        # return str(round(self.distance, 2))

    def __gt__(self, other):
        return self.distance > other.distance

    def __ge__(self, other):
        return self.distance >= other.distance

    def __lt__(self, other):
        if type(self) == type(other):
            return self.capability < other.capability
        else:
            return self.capability < other
