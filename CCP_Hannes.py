## @package Cytokines Hannes Baukmann's Course Completion Project
#
# This project aims to simulate the binding of immune cells to artificial beads decorated with molecules that induce a
# biological response of the cells. In nature, immene cells express receptors (so-called lectins) binding to
# carbohydrate structures, or glycans, on the surface of body's own cells or pathogens. Depending on the nature of the
# glycan ligand, cells express specific cytokines, e. g. the pro-inflammatory IL-6. The idea of the biochemical
# is to investigate the biological response of a monocytic cell line which recombinantly expresses a defined set of
# lectin receptors upon binding to glycan ligands that are chemically linked to acrylic glass beads in well-defined
# amounts.
#
# In the following, the conversion of real-world dimensions for this project will be explained.
#
# Monocytes are 5--20 um in diameter, while beads typically are few um in size.
# By definition, 1 pt in this program corresponds to 10 um in reality.
# The (round) basis of a well used in these experiments has an area of 0.36 qcm. Here, we consider a quadratic basis, hence
# \f$x=y=0.06 cm=600 pt\f$.
# Typical volume is 25 ul. Hence, the height of the liquid in a well is
# \f[ h = \frac{25 \mu l}{b} = \frac{25 \cdot 10^{-8} m^3}{0.36 cm^2} = 41.67 \cdot 10^{-4} m = 420 pt. \f]
# The number of cells usually used in these experiments is 50,000, incubated with a five-fold excess of beads.
# To make calculations feasible, the well dimension are reduced ten-fold in every dimension. In
# order to maintain equal density of particles, the number of cells and beads are consequently reduced to 50 and 250,
# respectively.
#
# The incubation time of the cells with beads amounts to 30 minutes, or 1800 s. In order to determine to many steps that
# correlates to, the Stokes-Einstein equation was utilized:
# \f[ x^2 = \frac{k_B \cdot T}{6r \pi \eta} t = \frac{1.381 \cdot 10^{-23} J/K \cdot 298 K}{6 \cdot 5 \mu m \cdot \pi \cdot 0.891 mPa \cdot s} = 9.8 \cdot 10^{-13} \frac{m^2}{s} \cdot t \f]
# At 1 s, a particle of 10 um in aqueous solution would diffuse (ignoring the influence of gravity on sinking particles):
# \f[ \sqrt{x^2} = \sqrt{4.9 \cdot 10^{-13} \frac{m^2}{s} \cdot 1s} = 7 \cdot 10^{-7} m.\f]
# Therefore, a typical experiment of 30 minutes corresponds to
# \f$ 1800 s \cdot 0.7 \mu m/s = 1260 \mu m = 126\mbox{ steps}. \f$
#
# In summary, the values used are:
# | Parameter       | Value (*in vitro*) | Value (*in silico*) | Reduced value |
# | :-------------- | :----------------- | :------------------ | :------------ |
# | Cell size       | 5-20 um            | 1 pt                | ---           |
# | Bead size       | few um             | 1 pt                | ---           |
# | Well basis      | 0.36 qcm           | 600 pt * 600 pt     | 60 pt, 60 pt  |
# | Volume          | 25 ul              | basis * 420 pt      | 42 pt         |
# | Number of cells | 50,000             | 50,000              | 50            |
# | Number of beads | 250,000            | 250,000             | 250           |
# | Incubation time | 30 min             | 126 steps           | ---           |
#
# The pairs of lectins, glycans, and cytokines implemented in the dictionary (vide infra) now (Aug 31, 2018) are:
#
# | Lectin   | Glycan  | Glycan type | Cytokine | Reference |
# | :------- | :------ | :---------- | :------- | :-------  |
# | DC-SIGN  | Mannan  | Mannose     | IL-6     | 1         |
# | DC-SIGN  | Lewis-Y | Fucose      | IL-27p28 | 1         |
# | Dectin-1 | Mannan  | Mannose     | IL-6     | 1         |
# 1) Geijtenbeek, T. B. H. and Gringhuis, S. I. (2016) C-type lectin receptors in the control of T helper cell differentiation. *Nat Rev Immunol* 16(7):433.

import random                   # random number
import numpy as np              # numpy arrays
import matplotlib.pyplot as plt # for plotting
import sys                      # system exit
import time                     # runtime meaurement

## Sphere serves as superclass for the two spherical objects indroduced below.
#
# \param[in] ID Identifier
# \param coordinates Empty. Will be filled when added to Well.

class Sphere:
    def __init__(self, ID_string):
        self.ID = ID_string
        self.coordinates = []

## Bead is a subclass of Sphere and represents beads (made of acrylic glass) loaded with glycan
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
        coordinates[0] = coordinates[0] + dx
        coordinates[1] = coordinates[1] + dy
        coordinates[2] = coordinates[2] + dz
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
                    if same:
#                    if (well.decoderCells[k].coordinates == well.beads[l].coordinates).all():
                        if random.uniform(0, 10000) <= well.beads[l].glycan_density * well.decoderCells[k].lectin_density:
                            if well.decoderCells[k].lectin.name in self.cytokine_dict[well.beads[l].glycan.type]:
                                self.cytokines.append(Cytokine(self.cytokine_dict[well.beads[l].glycan.type][well.decoderCells[k].lectin.name], well.decoderCells[k].coordinates))
        return self.cytokines


class Analysis:
    def __init__(self, simulationResults, cytokine_names=None):
        self.results = simulationResults
        if cytokine_names is None:
            c = []
            for i in range(len(simulationResults)):
                c.append(simulationResults[i].name)
            self.cytokineNames = list(set(c))
        else:
            self.cytokineNames=cytokine_names

    def countCytokines(self):
        self.cytokineAmount = []
        for i in range(len(self.cytokineNames)):
            counter = 0
            for j in range(len(self.results)):
                if self.cytokineNames[i] == self.results[j].name:
                    counter += 1
            self.cytokineAmount.append(counter)
        return [self.cytokineNames, self.cytokineAmount]

    def plotCytokines(self):
        self.countCytokines()
        labels = self.cytokineNames
        values = self.cytokineAmount

        indexes = np.arange(len(labels))
        width = 0.75

        plt.bar(indexes, values, width)
        plt.xticks(indexes, labels)
        plt.show()


if __name__ == "__main__":
    dependant = input("Show dependency of overall cytokine expression on (Lectins/Ratio_Cells_Beads/Density_Glycans/Runtime) ")
    if dependant == "Lectins":
        results = []
        lectins = [["DC-SIGN"], ["Dectin-1"], ["DC-SIGN", "Dectin-1"]]
        glycans = ["Mannan", "Lewis-Y"]
        glycans_types = ["Man", "Fuc"]
        cytokines = ["IL-6", "IL-27p28"]
        builder = Builder(60, 60, 42)
        for i in range(len(lectins)):
            model = Simulation(250, 50)
            well = model.createModel(builder, Well_list, glycans, glycans_types, 50, lectins[i], 50)
            model.simulate(well, 126) #number of randomWalk steps
            results.append(Analysis(model.cytokines, cytokines).countCytokines())

        plt.figure(1)

        p1 = plt.subplot(131)
        labels = cytokines
        values = results[0][1]
        indexes = np.arange(len(labels))
        plt.bar(indexes, values)
        plt.xticks(indexes, labels)
        plt.title(lectins[0])
        plt.ylabel("Cytokines Produced")

        p2 = plt.subplot(132, sharey=p1) # gleiche y achse wie erster subplot
        values = results[1][1]
        plt.bar(indexes, values)
        plt.xticks(indexes, labels)
        plt.title(lectins[1])
        plt.setp(p2.get_yticklabels(), visible=False) # y ticks wie oben festgelegt, keine beschrifttúng

        p3 = plt.subplot(133, sharey=p1) # gleiche y achse wie oben festgelegt
        values = results[2][1]
        plt.bar(indexes, values)
        plt.xticks(indexes, labels)
        plt.title(lectins[2])
        plt.setp(p3.get_yticklabels(), visible=False) # y ticks wie oben festgelegt, keine beschrifttúng

        plt.show()


    elif dependant == "Ratio_Cells_Beads":
        lectins = ["DC-SIGN", "Dectin-1"]
        glycans = ["Mannan", "Lewis-Y"]
        glycans_types = ["Man", "Fuc"]
        c = []
        r_range = [0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25]
        builder = Builder(60, 60, 42)
        n_total = 300
        for i in range(len(r_range)):
            n_cells = int(300/(r_range[i]+1)) # rundet auf ganze Zahl
            n_beads = n_total-n_cells
            model = Simulation(n_beads, n_cells) # number of beads an cells, realistisch 50.000 (oder nur 10.000 wegen Sedimentation); 1-10x (5x) so viele beads
            well = model.createModel(builder, Well_list, glycans, glycans_types, 50, lectins, 50)
            model.simulate(well, 126) #number of randomWalk steps
            c.append(len(model.cytokines))

        fig = plt.figure()
        plt.scatter(r_range, c, c='b', marker="s")
        plt.ylim(ymin=0)  # adjust the min leaving max unchanged
        plt.xlabel('Cell-to-Bead Ratio')
        plt.ylabel('Total Amount of Cytokines')
        plt.show()


    elif dependant == "Density_Glycans":
        lectins = ["DC-SIGN", "Dectin-1"]
        glycans = ["Mannan", "Lewis-Y"]
        glycans_types = ["Man", "Fuc"]
        c = []
        DG_range = [10, 25, 50, 75, 100]
        builder = Builder(60, 60, 42)
        model = Simulation(250, 50)
        for i in range(len(DG_range)):
            well = model.createModel(builder, Well_list, glycans, glycans_types, i, lectins, 50)
            model.simulate(well, 126) #number of randomWalk steps
            c.append(len(model.cytokines))

        fig = plt.figure()
        plt.scatter(DG_range, c, c='b', marker="s")
        plt.xlabel('Cell-to-Bead Ratio')
        plt.ylabel('Total Amount of Cytokines')
        plt.show()


    elif dependant == "Runtime":
        lectins = ["DC-SIGN", "Dectin-1"]
        glycans = ["Mannan", "Lewis-Y"]
        glycans_types = ["Man", "Fuc"]
        t_list = []
        t_npArray = []
        n_range = [1, 5, 10, 25, 50, 100, 250, 500]
        for i in range(len(n_range)):
            n_cells=n_range[i]
            n_beads = 5*n_cells
            n_total = n_cells + n_beads
            model = Simulation(n_beads, n_cells)  # number of beads an cells, realistisch 50.000 (oder nur 10.000 wegen Sedimentation); 1-10x (5x) so viele beads
            builder = Builder(60, 60, 42)  # Abmessungen des Wells; entspricht 2.5ul bei 1 pt=10 um
            start_time = time.time()
            well = model.createModel(builder, Well_list, glycans, glycans_types, 50, lectins, 50)
            model.simulate(well, 1) #number of randomWalk steps
            t_list.append(time.time() - start_time)
            start_time = time.time()
            well = model.createModel(builder, Well_npArray, glycans, glycans_types, 50, lectins, 50)
            model.simulate(well, 1) #number of randomWalk steps
            t_npArray.append(time.time() - start_time)

        fig = plt.figure()
        plt.semilogy(n_range, t_list, c='b', marker="s", label='list')
        plt.semilogy(n_range, t_npArray, c='r', marker="o", label='npArray')
        plt.legend(loc='upper left')
        plt.xlabel('Number of Cells')
        plt.ylabel('Runtime / s')
        plt.show()
    else:
        exit("Could not interpret key input. Input must be exactly as specified in the brackets in the promt statement.")