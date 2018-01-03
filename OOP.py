import random

class Cell:
    def __init__(self):
        pass # welche Parameter braucht jede Zelle? Größe? Größe könnte dann Schrittlänge bestimmen!
    def randomWalk2D(self, n):
        x, y = 0, 0
        for i in range(n):
            (dx, dy) = random.choice([(0,1),(0,-1),(1,0),(-1,0)])
            x += dx
            y += dy
        return (x, y)
    def randomWalk3D(self): # Geschwindigkeit? Schrittlänge? Dauer?
        pass
#    def binding(self):
#        pass

class Glycan: # Einteilung in Man- und Fuc-type?
    pass

class Lectin: # erstmal nur DC-SIGN
    pass

class EncoderCell(Cell): # can use all methods of class Cell
    def __init__(self):
        super().__init__()
        self.glycans=[]
        self.glycan_density=0
    def expressGlycans(self, glycan_string, glycan_density_int, glycan_receptor_list):
        self.glycans.append(Glycan(glycan_string, glycan_density_int, glycan_receptor_list))
        self.glycan_density += glycan_density_int
        if self.glycan_density > 1:
            print("ERROR!")
    def getAllGlycans(self):
        pass
    def showDensity(self):
        pass # Density aller Glycane

class DecoderCell(Cell):
    pass
    def binding(self, encoderCell): # Dauer abhängig von Glycan-Lectin-Paar?
        pass
    # Bindung, wenn Zellen genau in der selben Position sind? Erhöht WK des Re-bindings -- oder Abstoßung?
    # Bindung hängt von density der Glycane und Lectine ab -> WKVerteilung
    # gleichzeitige Bindung mehrerer Encoder and Decoder?
    def getAllLectins(self):
        pass
    def showDensity(self):
        pass
    def cytokineExpression(self): # Abhängig vom Glycan-Lectin-Paar
        pass

#if __name__ == ''__main__'' : # wofür is dat eijentlich?

### Creational Pattern: PROTOTYPE
# Encoder zu Decoder ca. 50:1

#if __name__ == "__main__":   # wofür is dat eijentlich?

test=EncoderCell()
print(test.randomWalk2D(25))