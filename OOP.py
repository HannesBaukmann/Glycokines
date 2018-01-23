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

class Cell:
    ## The constructor.
    def __init__(self, ID_string):
        self.ID=ID_string
    ## Very simple 2D Random walk.
    def randomWalk2D(self, n):
        x, y = 0, 0
        for i in range(n):
            (dx, dy) = random.choice([(0,1),(0,-1),(1,0),(-1,0)])
            x += dx
            y += dy
        return (x, y)
    def randomWalk3D(self): # Geschwindigkeit? Schrittlänge? Dauer?
        pass
    def binding(self):
        pass

class EncoderCell(Cell): # can use all methods of class Cell
    def __init__(self):
        super().__init__()
        self.glycans=[]
        self.glycan_density=0
    def expressGlycans(self, glycan_string, glycan_type_string, glycan_receptor_list, glycan_amount_int):
        self.glycans.append(Glycan(glycan_string, glycan_type_string, glycan_receptor_list, glycan_amount_int))
        self.glycan_density += glycan_amount_int
        if self.glycan_density > 1:
            print("ERROR!")
    def getAllGlycans(self):
        print(self.glycans)
        print()
    def showDensity(self):
        print(self.glycan_density)

class Glycan: # Einteilung in Man- und Fuc-type?
    def __init__(self, name_string, type_string, amount_int, cytokine_string): # Man-type vs. Fuc-type? # initialisieren mit Cytokinen?
        self.name = name_string
        self.type = type_string
        self.density = amount_int
        self.cytokine_evoked = cytokine_string
    def getSpecificity(self):
        print(self.receptors)

# für jeden Typ Glykan eine Klasse der Klasse Glycan!

class DecoderCell(Cell):
    def __init__(self, cytokine):
    def binding(self, encoderCell):
        pass
    # Bindung, wenn Zellen genau in der selben Position sind? Erhöht WK des Re-bindings -- oder Abstoßung?
    # Bindung hängt von density der Glycane und Lectine ab -> WKVerteilung
    # gleichzeitige Bindung mehrerer Encoder and Decoder?
    def getAllLectins(self):
        pass
    def showDensity(self):
        pass

class Cytokine:
    def __init__(self, name, x, y): # coordinaten
        self.cytokine_name=name

class Lectin: # erstmal nur DC-SIGN
    def __init__(self, name_string):
        self.name = name_string
    def setSpecificity(self, glycans_list_strings):
        self.ligands = glycans_list_strings

class Simulation: ## enhält Dictionary für Bindungsspezifität und Cytokinexpression
    def __init__(self, numberOfDecoders, numberOfEncoders):
        self.n = numberOfEncoders
        self.m = numberOfDecoders
    def createModel(self, builder):
        builder.buildModel()
    def __init__(self):
        self.cytokine_dict={"glycan1":{"DC-SIGN": ("IL6",), "Langerin": ("IL10",)},
                            "glycan2":("", "")}
        self.cytokines=[]
    def simulate(self):
        randomWalk
        detectBinding
    def detectBinding(self): # cytokines are produced and saved in container
        for glycan in glycan_list:
            for lectin in lectin_list:
                for product in self.cytokine_dict[glycan][lectin]:
                    self.cytokines.append(Cytokine(product, x, y))


# Creational Pattern: PROTOTYPE

test=EncoderCell()
print(test.randomWalk2D(25))

# um ein bestimmtes Verhältnis (z. B. von Glykanen) herzustellen:
# self.fraction initialisieren
# def buildNode(self, ID):
#   if random.random() < self.fraction:
#       self.graph.addNode(Node(ID))
#   else:
#       self.graph.addNode(EnergyTrap(ID))

# optionale Parameter: init...(bla, optional=10) optionale Parameter brauchen default-Wert!

if __name__ == "__main__":
    graph = TrapBarrierGraph()
    graph.addNode(Node(1))
    graph.addNode(Node(2))
    barrier = EnergyBarrier( 1, 2, 1)
    graph.addBarrier( barrier )
    print( graph.adjacency[1][2].height )

    model = Well( 10, 5 ) # Anzahl von Knoten und Mauern
    builder = Builder( .1,)  # meanHeight, meanDepth, fraction
    graph = model.createGraph( builder ) # createGraph aus TrapBarrierModel benutzt Methoden aus GraphBuilder
    print( graph.adjacency )

    walker = Walker( 0, graph )
    print( walker.position )
    for i in range( 20 ):
        walker.changeNode()
        print( walker.position )

class Builder: ## = GraphBuilder
    def __init__( self, fraction ):
        self.fraction = fraction # Hier fraction encoder / decoder?
    def createModel( self ):
        self.model = Well()
    def createCell( self, ID ):
        if random.random() < self.fraction:
            self.model.addCell( DecoderCell )
        else:
            self.model.addCell( EncoderCell )

class Well: ## = TrapBarrierGraph
    def __init__(self):

class TrapBarrierModel: ## = Simulation
    def __init__( self, numberOfNodes, numberOfBarriers ):
        if numberOfBarriers > ( numberOfNodes * ( numberOfNodes - 1 ) / 2 ):
            raise ValueError(" Too many Barriers! ")
        self.n = numberOfNodes
        self.m = numberOfBarriers
    def createGraph( self, builder ):
        builder.buildGraph()

        for i in range( self.n ):
            builder.buildNode(i)
            builder.graph.adjacency.update( {i:{}} )
        j = self.m
        while j > 0:
            nodes = random.sample( range( self.n ), 2 )
            if nodes[1] not in builder.graph.adjacency[ nodes[0] ]:
                builder.buildBarrier( nodes[0], nodes[1] )
                j -= 1

        return builder.graph

## checken, welche meiner Klassen Simons Klassen entsprechen und dann nochmal alles durchgehen.

## Klassendiagramm machen