# git add -A .
# git commit -m ""
# git push -u origin master
#
## \brief Brief description.
#         Brief description continues.
#
# Detailed description starts here.
# @package HannesOOPProject
# Dokumentation for this module using doxygen
#
# Markdown is supported. Would be nice to add some at the end of the project.

import random

## Dokumentation for class Cell
#
# \param x,y spatial coordinates

class Sphere:
    def __init__(self, ID_string):
        self.id=ID_string

class Bead(Sphere):
    def __init__(self):
        super().__init__()
        self.glycans=[]
        self.glycans_ratio=0
    def attachGlycans(self, glycan_string, glycan_type_string):
        if glycan_string in self.glycans:
            print("ERROR") # Duplikate sind nicht zugelassen!
            return # funktioniert das so, wie ich denke??
        self.glycans.append(Glycan(glycan_string, glycan_type_string))
        self.glycans_ratio = 100/len(self.glycans)

class Glycan:
    def __init__(self, name_string, type_string): # Man-type vs. Fuc-type...
        self.name = name_string
        self.type = type_string

# DecoderCell und Lectin erstellen, wenn Bead und Glycan abgeschlossen sind.

class Cytokine:
    def __init__(self, name, x, y):
        self.cytokine_name=name
        self.coordinates=[x,y]

class Well:
    def __init__( self, x, y ):
        self.length=x
        self.height=y
        self.beads = {}
        self.decoderCells = {}
    def addBead( self, bead ):
        self.beads.update( {bead.id:bead} )
    def addDecoderCell( self, decoderCell ):
        self.decoderCells.update( {decoderCell.id:decoderCell} )

class Builder:
    def __init__(self):
        pass
    def buildModel( self, x, y ): #Fehlermeldung, wenn Well kleiner als Anzahl Objekte?
        self.well = Well(x, y)
    def buildBead( self, id ):
        self.well.addBead(id)
    def buildDecoderCell( self, id ):
        self.well.addDecoderCell(id)

class Simulation: ## enthält (sehr simples) Dictionary für Bindungsspezifität und Cytokinexpression
    def __init__(self, numberOfDecoders, numberOfEncoders):
        self.n_beads = numberOfEncoders
        self.n_decoder = numberOfDecoders
        self.cytokines=[]
        self.cytokine_dict={"Man":{"DC-SIGN": ("IL-6",), "Dectin-1": ("IL-6",)},
                            "Fuc":("DC-SIGN", "IL-27p28")}
    def createModel(self, builder):
        builder.buildModel()

        for i in range( self.n_beads ):
            builder.buildEncoderCell(i)
            builder.well.encoderCells.update( {i:{}} )

        for i in range( self.n_decoder ):
            builder.buildDecoderCell(i)
            builder.well.decoderCells.update( {i:{}} )
        return builder.well
    def simulate(self):
        randomWalk
        detectBinding # if type == type!
        # Bindung, wenn Zellen genau in der selben Position sind? Erhöht WK des Re-bindings -- oder Abstoßung?
        # Bindung hängt von density der Glycane und Lectine ab -> WKVerteilung
        # gleichzeitige Bindung mehrerer Encoder and Decoder?
    def detectBinding(self): # cytokines are produced and saved in container
        for glycan in glycan_list:
            for lectin in lectin_list:
                for product in self.cytokine_dict[glycan][lectin]:
                    self.cytokines.append(Cytokine(product, x, y))
    ## Very simple 2D Random walk.
    def randomWalk2D(self, n):
        x, y = 0, 0 # should better start at random points
        for i in range(n):
            (dx, dy) = random.choice([(0,1),(0,-1),(1,0),(-1,0)])
            x += dx
            y += dy
        return (x, y)
    def randomWalk3D(self): # Geschwindigkeit? Schrittlänge? Dauer?
        pass