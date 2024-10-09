from ..modelling.features import LambdaGeologicalFeature

X = LambdaGeologicalFeature(lambda pos: pos[:, 0], name="x")
Y = LambdaGeologicalFeature(lambda pos: pos[:, 1], name="y")
Z = LambdaGeologicalFeature(lambda pos: pos[:, 2], name="z")
