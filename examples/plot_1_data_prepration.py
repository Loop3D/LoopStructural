"""
============================
1. Preparing input datasets
============================
This tutorial will demonstrate how to setup a dataset for using with LoopStructural.

"""

#########################################################################################
# Overview of data types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LoopStructural uses two main types of input data
# 1. Orientation data
# 2. Value data
#
# Orientation data constrains the gradient of the interpolated scalar field either by controlling
# the components of the normal vector to the scalar field or by adding a constraint which enforces that
# the dot product between the gradient and the interpolated scalar field is 0. Value constraints control the value
# of the scalar field and could be thought of as the distance away from a reference horizon.
#
