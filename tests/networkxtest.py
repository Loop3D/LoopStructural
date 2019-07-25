import networkx as nx
import matplotlib.pyplot as plt

model = nx.Graph()

model.add_node('FaultFrame')
model.add_node('FaultPlaneField')
model.add_node('FaultSlipField')
model.add_node('FaultDepthField')

model.add_edge('FaultFrame','FaultPlaneField', assocaition='component')
model.add_edge('FaultFrame','FaultSlipField', assocaition='component')
model.add_edge('FaultFrame','FaultDepthField', assocaition='component')

model.add_node('FoldFrame')
model.add_node('AxialFoliation')
model.add_node('FoldAxisDirection')
model.add_node('FoldExtensionDirection')

model.add_edge('FaultFrame', 'FoldFrame')
model.add_edge('FoldFrame', 'AxialFoliation', assocaition='component')
model.add_edge('FoldFrame', 'FoldAxisDirection', assocaition='component')
model.add_edge('FoldFrame', 'FoldExtensionDirection', assocaition='component')


model.add_node('FoldAxisRotation')
model.add_node('FoldLimbRotation')

model.add_edge('FoldFrame', 'FoldLimbRotation' )
model.add_edge('FoldFrame', 'FoldAxisRotation' )


model.add_node('Stratigraphy')
model.add_edge('FoldFrame', 'Stratigraphy' )


nx.draw(model,with_labels=True)
plt.savefig("testgraph.png")

for n in model.nodes:
    print(model.edges(n))