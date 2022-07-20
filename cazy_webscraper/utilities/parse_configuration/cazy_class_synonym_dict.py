#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# (c) University of St Andrews 2022
# (c) University of Strathclyde 2022
# (c) James Hutton Institute 2022
#
# Author:
# Emma E. M. Hobbs
#
# Contact
# eemh1@st-andrews.ac.uk
#
# Emma E. M. Hobbs,
# Biomolecular Sciences Building,
# University of St Andrews,
# North Haugh Campus,
# St Andrews,
# KY16 9ST
# Scotland,
# UK
#
# The MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""Accepted synonyms of CAZy class names"""


def cazy_synonym_dict():
    """Create a dictionary of accepted synonms for CAZy classes."""
    cazy_class_synonym_dict = {
        "Glycoside Hydrolases (GHs)": ["Glycoside-Hydrolases", "Glycoside-Hydrolases", "Glycoside_Hydrolases", "GlycosideHydrolases", "GLYCOSIDE-HYDROLASES", "GLYCOSIDE-HYDROLASES", "GLYCOSIDE_HYDROLASES", "GLYCOSIDEHYDROLASES", "glycoside-hydrolases", "glycoside-hydrolases", "glycoside_hydrolases", "glycosidehydrolases", "GH", "gh", "GHs", "ghs"],
        "GlycosylTransferases (GTs)": ["Glycosyl-Transferases", "GlycosylTransferases", "Glycosyl_Transferases", "Glycosyl Transferases", "GLYCOSYL-TRANSFERASES", "GLYCOSYLTRANSFERASES", "GLYCOSYL_TRANSFERASES", "GLYCOSYL TRANSFERASES", "glycosyl-transferases", "glycosyltransferases", "glycosyl_transferases", "glycosyl transferases", "GT", "gt", "GTs", "gts"],
        "Polysaccharide Lyases (PLs)": ["Polysaccharide Lyases", "Polysaccharide-Lyases", "Polysaccharide_Lyases", "PolysaccharideLyases", "POLYSACCHARIDE LYASES", "POLYSACCHARIDE-LYASES", "POLYSACCHARIDE_LYASES", "POLYSACCHARIDELYASES", "polysaccharide lyases", "polysaccharide-lyases", "polysaccharide_lyases", "polysaccharidelyases", "PL", "pl"],
        "Carbohydrate Esterases (CEs)": ["Carbohydrate Esterases", "Carbohydrate-Esterases", "Carbohydrate_Esterases", "CarbohydrateEsterases", "CARBOHYDRATE ESTERASES", "CARBOHYDRATE-ESTERASES", "CARBOHYDRATE_ESTERASES", "CARBOHYDRATEESTERASES", "carbohydrate esterases", "carbohydrate-esterases", "carbohydrate_esterases", "carbohydrateesterases", "CE", "ce", "CEs", "ces"],
        "Auxiliary Activities (AAs)": ["Auxiliary Activities", "Auxiliary-Activities", "Auxiliary_Activities", "AuxiliaryActivities", "AUXILIARY ACTIVITIES", "AUXILIARY-ACTIVITIES", "AUXILIARY_ACTIVITIES", "AUXILIARYACTIVITIES", "auxiliary activities", "auxiliary-activities", "auxiliary_activities", "auxiliaryactivities", "AA", "aa", "AAs", "aas"],
        "Carbohydrate-Binding Modules (CBMs)": ["Carbohydrate-Binding-Modules", "Carbohydrate_Binding_Modules", "Carbohydrate_Binding Modules", "CarbohydrateBindingModules", "CARBOHYDRATE-BINDING-MODULES", "CARBOHYDRATE_BINDING_MODULES", "CARBOHYDRATE_BINDING MODULES", "CARBOHYDRATEBINDINGMODULES", "carbohydrate-binding-modules", "carbohydrate_binding_modules", "carbohydrate_binding modules", "carbohydratebindingmodules", "CBMs", "CBM", "cbms", "cbm"]
    }
    return cazy_class_synonym_dict
