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
        self.ID=ID_string

class Bead(Sphere):
    def __init__(self, *args):
        super().__init__(*args)
        self.glycans=[]
        self.glycans_ratio=0
        self.x=None
        self.y=None
    def attachGlycans(self, glycan_string, glycan_type_string):
        self.glycans.append(Glycan(glycan_string, glycan_type_string))
        self.glycans_ratio = 100/len(self.glycans)

class Glycan:
    def __init__(self, name_string, type_string): # Man-type vs. Fuc-type...
        self.name = name_string
        self.type = type_string

class DecoderCell(Sphere):
    def __init__(self, *args):
        super().__init__(*args)
        self.lectins=[]
        self.lectins_ratio=0
        self.x=None
        self.y=None
    def expressLectins(self, lectin_string, lectin_type_string):
        self.lectins.append(Lectin(lectin_string, lectin_type_string))
        self.lectins_ratio = 100/len(self.lectins)

class Lectin:
    def __init__(self, name_string, type_string): # Man-type vs. Fuc-type...
        self.name = name_string
        self.type = type_string

class Cytokine:
    def __init__(self, name, x, y):
        self.cytokine_name=name
        self.coordinates=[x,y]

class Well:
    def __init__( self, x, y ):
        self.length=x
        self.height=y
        self.beads = [] # warum dict?
        self.decoderCells = []
    def addBead( self, bead, glyan_name_string, glycan_type_string):
        self.beads.append(bead)
        bead.x = random.randint(0, self.length) # Beads obtain starting position upon placement into Well
        bead.y = random.randint(0, self.height)
        bead.attachGlycans(glyan_name_string, glycan_type_string)
    def addDecoderCell( self, decoderCell, lectin_name_string, lectin_type_string):
        self.decoderCells.append(decoderCell)
        decoderCell.x = random.randint(0, self.length)
        decoderCell.y = random.randint(0, self.height)
        decoderCell.expressLectins(lectin_name_string, lectin_type_string)

class Builder:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def buildWell( self): #Fehlermeldung, wenn Well kleiner als Anzahl Objekte?
        self.well = Well(self.x, self.y)
    def buildBead( self, ID , glyan_name_string, glycan_type_string):
        self.well.addBead(Bead(ID), glyan_name_string, glycan_type_string)
    def buildDecoderCell( self, ID, lectin_name_string, lectin_type_string ):
        self.well.addDecoderCell(DecoderCell(ID), lectin_name_string, lectin_type_string)

class Simulation: ## enthält (sehr simples) Dictionary für Bindungsspezifität und Cytokinexpression
    def __init__(self, numberOfBeads, numberOfDecoderCells):
        self.n_beads = numberOfBeads
        self.n_decoder = numberOfDecoderCells # Fehlermeldung, wenn eines davon 0 ist?
        self.cytokines=[]
        self.cytokine_dict={"Man":{"DC-SIGN": ("IL-6",), "Dectin-1": ("IL-6",)},
                            "Fuc":("DC-SIGN", "IL-27p28")}
    def createModel(self, builder, glyan_name_string, glycan_type_string, lectin_name_string, lectin_type_string):
        builder.buildWell()

        for i in range( self.n_beads):
            ID = glyan_name_string+str(i)
            builder.buildBead(ID, glyan_name_string, glycan_type_string)
#            builder.well.beads.append()

        for i in range( self.n_decoder ):
            ID = lectin_name_string+str(i)
            builder.buildDecoderCell(ID, lectin_name_string, lectin_type_string)
#            builder.well.decoderCells.append( {i:{}} )
        return builder.well
    def simulate(self, well, n): # Kontrolle, ob überhaupt passende Paare im Well sind?
        for i in range(n):
            for j in range(len(well.beads)): # erst bewegen sich alle Beads
                (dx, dy) = random.choice([(0,1),(0,-1),(1,0),(-1,0)]) # Grenzen festlegen...
                well.beads[j].x += dx
                well.beads[j].y += dy
            for k in range(len(well.decoderCells)): # dann bewegen sich alle Zellen
                (dx, dy) = random.choice([(0,1),(0,-1),(1,0),(-1,0)]) # Grenzen festlegen...
                well.decoderCells[k].x += dx
                well.decoderCells[k].y += dy
                for l in range(len(well.beads)): # für jeden Bead wird geprüft, ob an der Stelle Zelle ist. Beads können nacheinander von verschiedenen Zellen gebunden werden
                    if well.decoderCells[k].x == well.beads[l].x: # x==x?
                        if well.decoderCells[k].y == well.beads[l].y: #y==y?
                            if random.uniform(0, 1) <= len(well.beads[l].glycans)*len(well.decoderCells[k].lectins): #Wk der Interaktion
                                for m in range(len(well.decoderCells[k].lectins)): # iteration über alle Lectine
                                    for n in range(len(well.beads[l].glycans)): # Iteration über alle Glycane
                                        if well.decoderCells[k].lectins[m].type == well.beads[l].glycans[n].type:
                                             for product in self.cytokine_dict[well.beads[l].glycans[n].type]:#[well.decoderCells[k].lectins[m].type]:
                                                 self.cytokines.append(Cytokine(product, well.decoderCells[k].x, well.decoderCells[k].y))
                                                 print("Motherfukcing Bindung!")
        for items in self.cytokines:
            print(items)
        return self.cytokines

if __name__ == "__main__":
#    well = Well(10,10)
#    well.addBead(Bead(b1))
#    well.addBead(Bead(b2))
#    well.addDecoderCell(DecoderCell(c1))

    model = Simulation(5, 2)
    builder = Builder(10, 10)
    well = model.createModel(builder, "Mannan", "Man", "DC-SIGN", "Man")
#    print(well.beads)
#    print(well.beads)
    model.simulate(well, 5)
#    print(well)


#    test=Bead("dings")
#   test.attachGlycans("Mannan", "Man")
#    print(test.glycans_ratio)