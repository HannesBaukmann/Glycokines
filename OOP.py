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
        self.ID = ID_string
        self.coordinates = []


class Bead(Sphere):
    def __init__(self, *args):
        super().__init__(*args)

    def attachGlycans(self, glycan_names_list, glycan_types_list, density_int):
        self.glycan = Glycan(glycan_names_list, glycan_types_list)
        self.glycan_density = density_int


class Glycan:
    def __init__(self, glycan_names_list, glycan_types_list):  # Fehlermeldung, wenn Listen nicht gleich lang sind!!
        r = random.randint(0, len(glycan_types_list)-1)
        self.name = glycan_names_list[r]
        self.type = glycan_types_list[r]


class DecoderCell(Sphere):
    def __init__(self, *args):
        super().__init__(*args)

    def expressLectins(self, lectins_list, density_int):
        self.lectin = Lectin(lectins_list)
        self.lectin_density = density_int


class Lectin:
    def __init__(self, lectins_list):
        self.name = random.choice(lectins_list)


class Cytokine:
    def __init__(self, name, coordinates_list):
        self.name = name
        self.coordinates = [coordinates_list]


class Well:
    def __init__(self, x, y, z):
        self.size = [x, y, z]
        self.beads = []
        self.decoderCells = []

    def addBead(self, bead, glyan_name_string, glycan_type_string, density_int):
        self.beads.append(bead)
        bead.coordinates = [random.randint(0, self.size[0]), random.randint(0, self.size[1]), random.randint(0, self.size[2])]
        bead.attachGlycans(glyan_name_string, glycan_type_string, density_int)

    def addDecoderCell(self, decoderCell, lectin_name_string, lectin_type_string):
        self.decoderCells.append(decoderCell)
        decoderCell.coordinates = [random.randint(0, self.size[0]), random.randint(0, self.size[1]), random.randint(0, self.size[2])]
        decoderCell.expressLectins(lectin_name_string, lectin_type_string)


class Builder:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def buildWell(self):  # Fehlermeldung, wenn Well kleiner als Anzahl Objekte?
        self.well = Well(self.x, self.y, self.z)

    def buildBead(self, ID, glyan_name_string, glycan_type_string, density_int):
        self.well.addBead(Bead(ID), glyan_name_string, glycan_type_string, density_int)

    def buildDecoderCell(self, ID, lectin_name_string, density_int):
        self.well.addDecoderCell(DecoderCell(ID), lectin_name_string, density_int)


class Simulation:  ## enthält (sehr simples) Dictionary für Bindungsspezifität und Cytokinexpression
    def __init__(self, numberOfBeads, numberOfDecoderCells):
        self.n_beads = numberOfBeads
        self.n_decoder = numberOfDecoderCells  # Fehlermeldung, wenn eines davon 0 ist?
        self.cytokines = []
        self.cytokine_dict = {"Man": {"DC-SIGN": ("IL-6"), "Dectin-1": ("IL-6")},
                              "Fuc": {"DC-SIGN": ("IL-27p28")}}

    def createModel(self, builder, glyan_names_list, glycan_types_list, glycan_density, lectin_name_string, lectin_density):
        builder.buildWell()

        for i in range(self.n_beads):
            ID = "Bead_" + str(i)
            builder.buildBead(ID, glyan_names_list, glycan_types_list, glycan_density)

        for i in range(self.n_decoder):
            ID = "DecoderCell_" + str(i)
            builder.buildDecoderCell(ID, lectin_name_string, lectin_density)
        return builder.well

    def simulate(self, well, n):  # Kontrolle, ob überhaupt passende Paare im Well sind?
        for i in range(n):
            for j in range(len(well.beads)):  # erst bewegen sich alle Beads
                (dx, dy, dz) = random.choice([(0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 0, -1), (0, -1, 0), (-1, 0, 0)])  # Grenzen festlegen...
                well.beads[j].coordinates = [sum(s) for s in zip(well.beads[j].coordinates, (dx, dy, dz))]
            for k in range(len(well.decoderCells)):  # dann bewegen sich alle Zellen
                (dx, dy, dz) = random.choice([(0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 0, -1), (0, -1, 0), (-1, 0, 0)])  # Grenzen festlegen...
                well.decoderCells[k].coordinates = [sum(s) for s in zip(well.decoderCells[k].coordinates, (dx, dy, dz))]
                for l in range(len(well.beads)):  # für jeden Bead wird geprüft, ob an der Stelle Zelle ist. Beads können nacheinander von verschiedenen Zellen gebunden werden
                    if well.decoderCells[k].coordinates == well.beads[l].coordinates:
                        if random.uniform(0, 10000) <= well.beads[l].glycan_density * well.decoderCells[k].lectin_density:
                            if well.decoderCells[k].lectin.name in self.cytokine_dict[well.beads[l].glycan.type]:
                                self.cytokines.append(Cytokine(self.cytokine_dict[well.beads[l].glycan.type][well.decoderCells[k].lectin.name], well.decoderCells[k].coordinates))# Fehlermeldung, wenn irgendwas nicht im dict
        #sowas wie Klassen Auswertung aufrufen?
        for i in range(len(self.cytokines)):
            print(self.cytokines[i].name, self.cytokines[i].coordinates)
        return self.cytokines


if __name__ == "__main__":
    model = Simulation(10, 22)  # number of beads an cells
    builder = Builder(11, 21, 5)  # Abmessungen des Wells
    well = model.createModel(builder, ["Mannan", "Lewis-X"], ["Man", "Fuc"], 100, ["DC-SIGN"], 100)
    model.simulate(well, 50)  # number of randomWalk steps