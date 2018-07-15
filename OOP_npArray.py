## Description of this project.
#
# Ueberschift!
# ============
#
# @package HannesOOPProject
# THP-1 cells are about 15 um in diameter, while beads typically are few um in size. Therefore, 1 pt in this program
# corresponds to 10 um in reality.
# Size of well:
# The (round) basis of a well has an area of 0.36 cm2. Here, we consider a quadratic basis, hence x=y=0.06 cm2. With
# 1pt = 10 um, x=y=600pt
# Typical volume is 25 ul. Hence, the height of the liquid in a well is \f$ 25 \mu l = \frac{25 \cdot 10^{-8} m^3}{0.36 cm^2} = 41.67 \cdot 10^{-4} m =
# 420 pt \f$.
# Typical experiment takes 30 minutes or 1800 s. How many steps does that correlate to? How fast does a particle of 10 um / 1 pt in diameter diffuse?
# Stokes-Einstein equation:
# \f$ x^2 = \frac{k_B \cdot T}{3r \pi \eta} t = 4.9 \cdot 10^{-13} \frac{m^2}{s} \cdot t \f$
# (D_H_2O (25°C) = 0.891 mPa \cdot s)
# At 1 s, the medium ... is 7 \cdot 10^{-7} m^2.
# Therefore, a typical experiment of 30 minutes corresponds to \f$ 1800 s \cdot 0.7 \mu m/s = 1260 \mu m \f$, or 126 steps.
# 25 grad

import random
from collections import Counter # um Elemente in cytokines zu zählen
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt # für Plot
import sys
import time
import operator


## \brief Sphere serves as superclass for the two spherical objects indroduced below.
#
# \param[in] ID Identifier
# \param coordinates Empty. Will be filled when added to Well.

class Sphere:
    def __init__(self, ID_string):
        self.ID = ID_string
        self.coordinates = []

## \brief Bead is a subclass of Sphere and represents beads (made of PMMA, or acrylic glass) loaded with glycan
# structures.
#

class Bead(Sphere):
    def __init__(self, *args):
        super().__init__(*args)

    ## Method to attach Glycan objects to Bead objects. Every Bead objects contains one type of Glycan with a given density.
    #  Parameters glycan_names_list and glycan_types_list are handed over to the constructor of Glycan.
    #
    # \param[in] glycan_names_list List of names of all Glycans to be added to the Bead. List may contain any number of
    # member, including 0
    # \param[in] see documentation for class Glycan.
    # \param[in] see documentation for class Glycan.
    def attachGlycans(self, glycan_names_list, glycan_types_list, density_percentage):
        if density_percentage > 100:
            print("Given density_percentage higher than 100%! Value was set to 100%")
            self.glycan_density = 100
        elif density_percentage < 0:
            print("Given density_percentage lower than zero! Value was set to 0%")
            self.glycan_density = 0
        else:
            self.glycan_density = density_percentage
        self.glycan = Glycan(glycan_names_list, glycan_types_list)


## \brief Objects of the class Glycan are attached to objects of the class Bead. They represent glycan structures which
# are divided into different types.
#
# \param[in] glycan_names_list List of Strings containing the "real" names of Glycan structures to be added to the Bead.
# List may contain any number of members, including 0
# \param[in] glycan_types_list List of Strings containing the respective types of Glycan structures. These types
# determine the outcome of the interaction between Bead and DecoderCell as specified in the cytokine dictionary. If the
# this list does not contain the same number of member as glycan_names_list, the program terminates.

class Glycan:
    def __init__(self, glycan_names_list, glycan_types_list):
        if len(glycan_names_list) != len(glycan_types_list):
            sys.exit("List of Glycan names and types have to have the same number of entries. Program terminated!")
        r = random.randint(0, len(glycan_types_list)-1)
        self.name = glycan_names_list[r]
        self.type = glycan_types_list[r]

## \brief DecoderCell is a subclass of Sphere and represents immune cells (such as THP-1 cells) that express Lectins
# which bind to Glycan structures and engage in immune response by secreting cytokines.
#
# | Right | Center | Left  |
# | ----: | :----: | :---- |
# | 10    | 10     |    10 |
# | 1000  |   1000 |  1000 |
#
# \f$ \sqrt{(x_2-x_1)^2+(y_2-y_1)^2} \f$

class DecoderCell(Sphere):
    def __init__(self, *args):
        super().__init__(*args)

    ## Method to attach Lectin objects to DecoderCell objects. Every DecoderCell objects contains one type of Lectin
    # with a given density. Parameter lectins_list is handed over to the constructor of Lectin.
    #
    # \param[in] lectin_list see documentation for class Lectin.
    def expressLectins(self, lectin_list, density_percentage):
        if density_percentage > 100:
            print("Given density_percentage higher than 100%! Value was set to 100%")
            self.lectin_density = 100
        elif density_percentage < 0:
            print("Given density_percentage lower than zero! Value was set to 0%")
            self.lectin_density = 0
        else:
            self.lectin_density = density_percentage
        self.lectin = Lectin(lectin_list)

## Objects of the class Lectin are attached to objects of the class DecoderCell. They represent lectin receptors
# which recognize certain sets of glycans on other cells, beads, viruses a.s.o.
#
# \param[in] lectin_list List of Strings containing the names of Lectin receptors to be added to the DecoderCell. These
# names appear in the dictionary. List may contain any number of members, including 0.

class Lectin:
    def __init__(self, lectin_list):
        self.name = random.choice(lectin_list)


## Objects of class Cytokine are produced if objects of the classes Bead and DecoderCell containing the correct pair of
# objects of the classes Glycan and Lectin, respectively, as described by the cytokine dictionary.
#
# \param[in] name Name of the cytokine.
# \param[in] coordinates_list List containing the three values for x, y, and z, where the Cytokine was produced.

class Cytokine:
    def __init__(self, name, coordinates_list):
        self.name = name
        self.coordinates = coordinates_list

## The class Well describes the sample, in which binding events take place. It serves a a superclass for two subclasses
# of which one uses NumPy Arrays as a container for onjects of the classes Bead and DecoderCell while the other one uses
# built-in Python lists for this purpose.
#
# \param[in] x, y, z Size of the well.

class Well:
    def __init__(self, x, y, z):
        if x <= 0 or y <= 0 or z <= 0:
            sys.exit("Illegal Well dimensions! Please enter values larger than zero!")
#        self.size = [x, y, z]

    ## The method borderControl ensures that objects' coordinates don't exceed the well's size. If a step in randomWalk
    # would lead to a forbidden value, the respective step value will be set to 0 and the object won't move this turn as
    # if repelled from the well's borders.
    #
    # \param[in] coordinates The object's coordinates before the move.
    # \param[in] dx, dy, dz The randomly chosen next steps.
    # \param[out] dx, dy, dz The revised next steps

    def borderControl(self, coordinates, dx, dy, dz):
        coordinates = [sum(s) for s in zip(coordinates, (dx, dy, dz))]
        if self.size[0] < abs(coordinates[0]) or 0 > coordinates[0]:
            dx = 0
        if self.size[1] < abs(coordinates[1]) or 0 > coordinates[1]:
            dy = 0
        if self.size[2] < abs(coordinates[2]) or 0 > coordinates[2]:
            dz = 0
        return dx, dy, dz

## Subclass of Well using built-in python lists as containers for objects of the classes Bead and DecoderCell.
#

class Well_list(Well):
    def __init__(self, x, y, z, n_beads, n_cells):
        super().__init__(x, y, z)
        self.size=[x, y, z]
        self.beads = []
        self.decoderCells = []

    ## Method to add objects of class Bead to the bead list. Coordinates of the Bead are randomly chosen between 0 and
    # size of the Well in the respective dimension.
    #
    # \param[in] i Not used here, but neccessary for convenient switching between container types.
    # \param[in] bead Object of class Bead to be added.
    # \param[in] 

    def addBead(self, i, bead, glyan_name_string, glycan_type_string, density_percentage):
        self.beads.append(bead)
        bead.coordinates = [random.randint(0, self.size[0]), random.randint(0, self.size[1]), random.randint(0, self.size[2])]
        bead.attachGlycans(glyan_name_string, glycan_type_string, density_percentage)

    def addDecoderCell(self, i, decoderCell, lectin_name_string, lectin_type_string):
        self.decoderCells.append(decoderCell)
        decoderCell.coordinates = [random.randint(0, self.size[0]), random.randint(0, self.size[1]), random.randint(0, self.size[2])]
        decoderCell.expressLectins(lectin_name_string, lectin_type_string)


class Well_npArray(Well):
    def __init__(self, x, y, z, n_beads, n_cells):
        super().__init__(x, y, z)
        self.size = np.array([x, y, z])
        self.beads = np.empty(n_beads, dtype=object)
        self.decoderCells = np.empty(n_cells, dtype=object)

    def addBead(self, i, bead, glyan_name_string, glycan_type_string, density_percentage):
        self.beads[i]=bead
        bead.coordinates = np.array([random.randint(0, self.size[0]), random.randint(0, self.size[1]), random.randint(0, self.size[2])])
        bead.attachGlycans(glyan_name_string, glycan_type_string, density_percentage)

    def addDecoderCell(self, i, decoderCell, lectin_name_string, lectin_type_string):
        self.decoderCells[i] = decoderCell
        decoderCell.coordinates = np.array([random.randint(0, self.size[0]), random.randint(0, self.size[1]), random.randint(0, self.size[2])])
        decoderCell.expressLectins(lectin_name_string, lectin_type_string)


class Builder:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def buildWell(self, container_type, n_beads, n_cells):
        self.well = container_type(self.x, self.y, self.z, n_beads, n_cells)

    def buildBead(self, i, ID, glyan_name_string, glycan_type_string, density_percentage):
        self.well.addBead(i, Bead(ID), glyan_name_string, glycan_type_string, density_percentage)

    def buildDecoderCell(self, i, ID, lectin_name_string, density_percentage):
        self.well.addDecoderCell(i, DecoderCell(ID), lectin_name_string, density_percentage)

## Dictionary is based on Geijtenbeek & Gringhuis C-type lectin receptors in the control of T helper cell differentiation Nature Reviews Immunology volume 16, pages 433–448 (2016)
class Simulation: ## enthält (sehr simples) Dictionary für Bindungsspezifität und Cytokinexpression
    def __init__(self, numberOfBeads, numberOfDecoderCells):
        self.n_beads = numberOfBeads
        self.n_decoder = numberOfDecoderCells
        self.cytokines = []
        self.cytokine_dict = {"Man": {"DC-SIGN": ("IL-6"), "Dectin-1": ("IL-6")},
                              "Fuc": {"DC-SIGN": ("IL-27p28")}}

    def createModel(self, builder, container_type, glycan_names_list, glycan_types_list, glycan_density, lectin_list, lectin_density):
        allLinD = []
        for i in range(len(glycan_types_list)):
            if glycan_types_list[i] not in self.cytokine_dict:
                sys.exit("Glycan Type not in Dictionary!")
        for lectins, cytokines in self.cytokine_dict.items():
            for l, c in cytokines.items():
                allLinD.append(l)
        allLinD=list(set(allLinD))
        for j in range(len(lectin_list)):
            if lectin_list[j] not in allLinD:
                sys.exit("Lectin not in Dictionary!")

        builder.buildWell(container_type, self.n_beads, self.n_decoder)

        for i in range(self.n_beads):
            ID = "Bead_" + str(i)
            builder.buildBead(i, ID, glycan_names_list, glycan_types_list, glycan_density)

        for i in range(self.n_decoder):
            ID = "DecoderCell_" + str(i)
            builder.buildDecoderCell(i, ID, lectin_list, lectin_density)
        return builder.well

#    def produceCytokines (self, well, coordinates, k, l):
#        if random.uniform(0, 10000) <= well.beads[l].glycan_density * well.decoderCells[k].lectin_density:
#            if well.decoderCells[k].lectin.name in self.cytokine_dict[well.beads[l].glycan.type]:
#                self.cytokines.append(Cytokine(self.cytokine_dict[well.beads[l].glycan.type][well.decoderCells[k].lectin.name],well.decoderCells[k].coordinates))

    def simulate(self, well, n):
        for i in range(n):
            for j in range(len(well.beads)):  # erst bewegen sich alle Beads
                (dx, dy, dz) = random.choice([(0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 0, -1), (0, -1, 0), (-1, 0, 0)])
#                well.beads[j].coordinates = [sum(s) for s in zip(well.beads[j].coordinates, well.borderControl(well.beads[j].coordinates, dx, dy, dz))]
                dc=well.borderControl(well.beads[j].coordinates, dx, dy, dz)
                for (i, coordinate) in enumerate(well.beads[j].coordinates):
                    well.beads[j].coordinates[i] += dc[i]
            for k in range(len(well.decoderCells)):  # dann bewegen sich alle Zellen
                (dx, dy, dz) = random.choice([(0, 0, 1), (0, 1, 0), (1, 0, 0), (0, 0, -1), (0, -1, 0), (-1, 0, 0)])
#                well.decoderCells[k].coordinates = [sum(s) for s in zip(well.decoderCells[k].coordinates, well.borderControl(well.decoderCells[k].coordinates, dx, dy, dz))]
                dc=well.borderControl(well.decoderCells[k].coordinates, dx, dy, dz)
                for (i, coordinate) in enumerate(well.decoderCells[k].coordinates):
                    well.decoderCells[k].coordinates[i] += dc[i]
                for l in range(len(well.beads)):  # für jeden Bead wird geprüft, ob an der Stelle Zelle ist. Beads können nacheinander von verschiedenen Zellen gebunden werden
                    same = True
                    for i in range(len(well.decoderCells[k].coordinates)):
                        same = (well.decoderCells[k].coordinates[i] == well.beads[l].coordinates[i])
#                        print(same)
                    if same:
#                    if (well.decoderCells[k].coordinates == well.beads[l].coordinates).all():
                        if random.uniform(0, 10000) <= well.beads[l].glycan_density * well.decoderCells[k].lectin_density:
                            if well.decoderCells[k].lectin.name in self.cytokine_dict[well.beads[l].glycan.type]:
                                self.cytokines.append(Cytokine(self.cytokine_dict[well.beads[l].glycan.type][well.decoderCells[k].lectin.name], well.decoderCells[k].coordinates))
        return self.cytokines


class Analysis:
    def __init__(self):
        pass

    def countCytokines(self, simulationResults):
        self.cytokineAmount = []
        names = []
        for i in range(len(simulationResults)):
            names.append(simulationResults[i].name)
        self.cytokineNames = list(set(names))
        for i in range(len(self.cytokineNames)):
            counter = 0
            for j in range(len(simulationResults)):
                if self.cytokineNames[i] == simulationResults[j].name:
                    counter += 1
            self.cytokineAmount.append(counter)
        print(self.cytokineNames)
        print(self.cytokineAmount)


    def plotCytokines(self, simulationResults):
        self.countCytokines(simulationResults)
        labels = self.cytokineNames
        values = self.cytokineAmount

        indexes = np.arange(len(labels))
        width = 0.75

        plt.bar(indexes, values, width)
        plt.xticks(indexes, labels)
        plt.show()

if __name__ == "__main__":
    comparison = input("Laufzeitvergleich machen (True/False)? ")
    if comparison == "True":
        t_list = []
        t_npArray = []
        n_range = [5, 10, 50, 100, 500]
        for i in range(len(n_range)):
            n_cells=n_range[i]
            n_beads = 5*n_cells
            n_total = n_cells + n_beads
            model = Simulation(n_beads, n_cells)  # number of beads an cells, realistisch 50.000 (oder nur 10.000 wegen Sedimentation); 1-10x (5x) so viele beads
            builder = Builder(6, 6, 42)  # Abmessungen des Wells; entspricht 2.5ul bei 1 pt=10 um
            start_time = time.time()
            well = model.createModel(builder, Well_list, ["Mannan", "Lewis-X"], ["Man", "Fuc"], 100, ["DC-SIGN"], 77)
            model.simulate(well, 1) #number of randomWalk steps
            t_list.append(time.time() - start_time)
            print(t_list)
            start_time = time.time()
            well = model.createModel(builder, Well_npArray, ["Mannan", "Lewis-X"], ["Man", "Fuc"], 100, ["DC-SIGN"], 77)
            model.simulate(well, 1) #number of randomWalk steps
            t_npArray.append(time.time() - start_time)
            print(t_npArray)



#        if isinstance(well.beads, np.ndarray):
#            print("si")
#        else:
#            print("nope")
#        print(type(well.beads[0].coordinates))
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.scatter(n_range, t_list, s=10, c='b', marker="s", label='list')
        ax1.scatter(n_range, t_npArray, s=10, c='r', marker="o", label='npArray')
        plt.legend(loc='upper left');
        plt.show()
    else:
        print("wibe")
#        auswertung = Analysis().plotCytokines(model.cytokines)