# -*- coding: utf-8 -*-

###########################################################
# PyToxo
#
# A Python tool to calculate penetrance tables for 
# high-order epistasis models
#
# Copyright 2021 Borja González Seoane
#
# Contact: borja.gseoane@udc.es
###########################################################

from pytoxo.model import Model
from pytoxo.ptable import PTable

# Package PyToxo exports only the `Model` class
__all__ = ["Model", "PTable"]
