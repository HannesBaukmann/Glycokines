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
    def addBead( self, bead):
        self.beads.update( {bead.id:bead} )
        bead.x = random.randint(0, self.length) # Beads obtain starting position upon placement into Well
        bead.y = random.randint(0, self.height)
    def addDecoderCell( self, decoderCell):
        self.decoderCells.update( {decoderCell.id:decoderCell} )
        decoderCell.x = random.randint(0, self.length)
        decoderCell.y = random.randint(0, self.height)

class Builder:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def buildModel( self): #Fehlermeldung, wenn Well kleiner als Anzahl Objekte?
        self.well = Well(self.x, self.y)
    def buildBead( self, id ):
        self.well.addBead(id)
    def buildDecoderCell( self, id ):
        self.well.addDecoderCell(id)

class Simulation: ## enthält (sehr simples) Dictionary für Bindungsspezifität und Cytokinexpression
    def __init__(self, numberOfDecoders, numberOfEncoders):
        self.n_beads = numberOfEncoders
        self.n_decoder = numberOfDecoders # Fehlermeldung, wenn eines davon 0 ist?
        self.cytokines=[]
        self.cytokine_dict={"Man":{"DC-SIGN": ("IL-6",), "Dectin-1": ("IL-6",)},
                            "Fuc":("DC-SIGN", "IL-27p28")}
    def createModel(self, builder):
        builder.buildWell()

        for i in range( self.n_beads ):
            builder.buildEncoderCell(i)
            builder.well.encoderCells.update( {i:{}} )

        for i in range( self.n_decoder ):
            builder.buildDecoderCell(i)
            builder.well.decoderCells.update( {i:{}} )
        return builder.well
    def simulate(self, well, n): # Kontrolle, ob überhaupt passende Paare im Well sind?
        for i in range(n):
            for j in range(well.beads): # erst bewegen sich alle Beads
                (dx, dy) = random.choice([(0,1),(0,-1),(1,0),(-1,0)]) # Grenzen festlegen...
                well.beads(j).x += dx
                well.beads(j).y += dy
            for k in range(well.encoderCells): # dann bewegen sich alle Zellen
                (dx, dy) = random.choice([(0,1),(0,-1),(1,0),(-1,0)]) # Grenzen festlegen...
                well.encoderCells(k).x += dx
                well.encoderCells(k).y += dy
                for l in range(well.beads): # für jeden Bead wird geprüft, ob an der Stelle Zelle ist. Beads können nacheinander von verschiedenen Zellen gebunden werden
                    if well.encoderCells(k).x == well.beads(l).x: # x==x?
                        if well.encoderCells(k).y == well.beads(l).y: #y==y?
                            if random.uniform(0, 1) <= len(well.beads(l).glycans) * len(well.encoderCells(k).lectins): #Wk der Interaktion
                                for m in range(well.encoderCells(k).lectins): # iteration über alle Lectine
                                    for n in range(well.beads(l).glycans): # Iteration über alle Glycane
                                         if well.encoderCells(k).lectins(m).type == well.beads(l).glycans(n).type:
                                             for product in self.cytokine_dict[well.beads(l).glycans(n).type][well.encoderCells(k).lectins(m).type]:
                                                 self.cytokines.append(Cytokine(product, well.encoderCell(k).x, well.encoderCell(k).y))
                                                print("Motherfukcing Bindung!")
        return self.cytokines
    