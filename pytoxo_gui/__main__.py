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

MAX_ORDER_SUPPORTED = 12


class PyToxoContext:
    """Support container class to maintain the PyToxo stuff state during GUI
    main loop."""

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

# Model frame
model_frame = sg.Frame(
    key="_MODEL_FRAME_",
    title="Epistatic model",
    layout=[
        [
            sg.Text(
                "None model loaded",
                key="_MODEL_DISABLED_TEXT_",
                text_color="grey",
                tooltip="Use the menu to load a model from a file or enter one manually",
                visible=True,  # Pending to be disabled when a model was loaded
            ),
            sg.Table(
                key="_MODEL_TABLE_",
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

# MAFs frame text entries
mafs_entries = []
for i in range(1, MAX_ORDER_SUPPORTED + 1):
    mafs_entries.append(
        sg.InputText(
            key=f"_MAFS_INPUT_{i}_",
            size=(3, 1),
            visible=False,  # Pending to be enabled when a model was loaded
        )
    )


# Configuration frame
configuration_frame = sg.Frame(
    key="_CONFIG_FRAME_",
    title="Configuration",
    layout=[
        [
            sg.Combo(
                ("Heritability", "Prevalence"),
                key="_PREV_OR_HER_CB_",
                size=(10, 1),
                readonly=True,
            ),
            sg.InputText(key="_PREV_OR_HER_INPUT_", size=(5, 1)),
            sg.Frame(
                key="_MAFS_FRAME_",
                title="MAFs",
                layout=[
                    [
                        sg.Text(
                            "None model loaded",
                            key="_MAFS_DISABLED_TEXT_",
                            text_color="grey",
                            tooltip="You need to have set the model before setting MAFs",
                            visible=True,  # Pending to be disabled when a model was loaded
                        ),
                    ]
                    + mafs_entries
                ],
                element_justification="center",
            ),
        ],
    ],
    element_justification="left",
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
    size=(650, 600),
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
                window["_MODEL_DISABLED_TEXT_"].Update(visible=False)
                window["_MODEL_TABLE_"].Update(visible=True)
                window["_MODEL_TABLE_"].Update(
                    values=[
                        [g, p]
                        for g, p in zip(
                            pytoxo_context.model._calculate_genotypes(),
                            pytoxo_context.model.penetrances,
                        )
                    ]
                )
                # Enable MAFs
                window["_MAFS_DISABLED_TEXT_"].Update(visible=False)
                for i in range(1, pytoxo_context.model.order + 1):
                    window[f"_MAFS_INPUT_{i}_"].Update(visible=True)
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
