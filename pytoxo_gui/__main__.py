# -*- coding: utf-8 -*-

###########################################################
# PyToxo
#
# A Python library for calculating penetrance tables of any
# bivariate epistasis model.
#
# Copyright 2021 Borja González Seoane
#
# Contact: borja.gseoane@udc.es
###########################################################

"""Graphical user interface entry point."""

import PySimpleGUI as sg

import pytoxo
import pytoxo.errors


class PyToxoContext:
    """Support container class to maintain the PyToxo stuff state during GUI
    main loop"""

    model = None  # Loaded PyToxo model


# ####################### GUI DESIGN ######################
# General style settings
sg.theme("SystemDefaultForReal")
window_general_font = "Verdana, 13"

# Main menu
menu = sg.Menu(
    [
        ["File", ["Open model", "Exit"]],
        # ["Edit", ["Paste", "Undo"]],
        ["Help", "About PyToxo GUI"],
    ]
)

# Epistatic model table
headings = [
    " Genotype definition ",
    "Penetrance expression",
    "Calculated penetrance",
]
rows = [[col + row for col in range(len(headings))] for row in range(100)]

# Configuration frame
configuration_frame = sg.Frame(
    key="_CONFIG_FRAME_",
    title="Configuration",
    layout=[
        [
            sg.Combo(
                ("Heritability", "Prevalence"),
                key="prev_or_her_cb",
                size=(10, 1),
                readonly=True,
            ),
            sg.InputText(key="prev_or_her_in", size=(5, 1)),
            sg.Text("MAFs"),
            # sg.InputText(key="mafs_in", size=(5, 1)),
            sg.Text(
                "None model loaded",
                text_color="grey",
                tooltip="You need to have set the model before setting MAFs",
            ),
        ],
    ],
    element_justification="left",
)

# Model frame
model_frame = sg.Frame(
    key="_MODEL_FRAME_",
    title="Epistatic model",
    layout=[
        [
            sg.Text(
                "None model loaded",
                key="_MODEL_FRAME_TEXT_",
                text_color="grey",
                tooltip="Use the menu to load a model from a file or enter one manually",
                visible=True,  # Pending to be disabled when a model was loaded
            ),
            sg.Table(
                key="_MODEL_FRAME_TABLE_",
                headings=headings,
                values=rows,
                vertical_scroll_only=True,
                justification="center",
                visible=False,  # Pending to be enabled when a model was loaded
            ),
        ]
    ],
    element_justification="center",
)

# Layout composition
layout = [
    [menu],
    [model_frame],
    [configuration_frame],
    [
        sg.Button(
            "Calculate table",
            disabled=True,
            tooltip="You need to have set all configurations before calculating the table",
        )
    ],
]

# Window composition
window = sg.Window(
    "PyToxo GUI",
    layout,
    font=window_general_font,
    element_justification="center",
    size=(600, 600),
    finalize=True
    # background_color="#eeeeee",
)

# Window style patches
window["_MODEL_FRAME_"].expand(expand_x=True, expand_y=False)
# window["_CONFIG_FRAME_"].expand(expand_x=True, expand_y=False)
# #########################################################


# Create PyToxo context object
pytoxo_context = PyToxoContext()


# GUI event loop
def main():
    while True:
        event, values = window.read()

        if event == "Open model":
            filename = sg.popup_get_file(
                "Open model",
                no_window=True,  # To use a native approach
                file_types=(("Comma separated values", "*.csv"),),
            )
            try:
                pytoxo_context.model = pytoxo.Model(filename)
                # Print model in the GUI's table
                window["_MODEL_FRAME_TEXT_"].Update(visible=False)
                window["_MODEL_FRAME_TABLE_"].Update(visible=True)
                window["_MODEL_FRAME_TABLE_"].Update(
                    values=[
                        [g, p]
                        for g, p in zip(
                            pytoxo_context.model._calculate_genotypes(),
                            pytoxo_context.model.penetrances,
                        )
                    ]
                )
            except pytoxo.errors.BadFormedModelError:
                sg.popup_ok(
                    "The file contains a bad formed model. PyToxo cannot "
                    "interpret it. Revise PyToxo's file format requirements.",
                    title="File parsing error",
                    font=window_general_font,
                )
            except IOError:
                sg.popup_ok(
                    f"Error trying to open '{filename}'.",
                    title="File opening error",
                    font=window_general_font,
                )
        elif event in ("Exit", None):
            break

    window.close()


if __name__ == "__main__":
    main()
