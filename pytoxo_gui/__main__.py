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
    key="-MODEL_FRAME-",
    title="Epistatic model",
    layout=[
        [
            sg.Text(
                "None model loaded",
                key="-MODEL_DISABLED_TEXT-",
                text_color="grey",
                tooltip="Use the menu to load a model from a file or enter one manually",
                visible=True,  # Pending to be disabled when a model was loaded
            ),
            sg.Table(
                key="-MODEL_TABLE-",
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
            key=f"-MAFS_INPUT_{i}-",
            size=(3, 1),
            visible=False,  # Pending to be enabled when a model was loaded
            enable_events=True,  # To refresh the loop and can check filled fields
        )
    )


# Configuration frame
configuration_frame = sg.Frame(
    key="-CONFIG_FRAME-",
    title="Configuration",
    layout=[
        [
            sg.Combo(
                ("Heritability", "Prevalence"),
                key="-PREV_OR_HER_CB-",
                size=(10, 1),
                readonly=True,
                enable_events=True,  # To refresh the loop and can check filled fields
            ),
            sg.InputText(
                key="-PREV_OR_HER_INPUT-",
                size=(5, 1),
                enable_events=True,  # To refresh the loop and can check filled fields
            ),
            sg.Frame(
                key="-MAFS_FRAME-",
                title="MAFs",
                layout=[
                    [
                        sg.Text(
                            "None model loaded",
                            key="-MAFS_DISABLED_TEXT-",
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
            tooltip="You need to have set all fields before calculating the table",
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
window["-MODEL_FRAME-"].expand(expand_x=True, expand_y=False)
# window["-CONFIG_FRAME-"].expand(expand_x=True, expand_y=False)
# #########################################################


# Create PyToxo context object
pytoxo_context = PyToxoContext()


# GUI event loop
def main():
    while True:
        event, values = window.read()

        # Check events
        if event == "Open model":
            filename = sg.popup_get_file(
                "Open model",
                no_window=True,  # To use a native approach
                file_types=(("Comma separated values", "*.csv"),),
            )
            if not filename:
                continue  # The operation has been canceled
            try:
                pytoxo_context.model = pytoxo.Model(filename)
                # Print model in the GUI's table
                window["-MODEL_DISABLED_TEXT-"].Update(visible=False)
                window["-MODEL_TABLE-"].Update(visible=True)
                window["-MODEL_TABLE-"].Update(
                    values=[
                        [g, p]
                        for g, p in zip(
                            pytoxo_context.model.calculate_genotypes(),
                            pytoxo_context.model.penetrances,
                        )
                    ]
                )
                # Enable MAFs
                window["-MAFS_DISABLED_TEXT-"].Update(visible=False)
                for i in range(1, pytoxo_context.model.order + 1):
                    window[f"-MAFS_INPUT_{i}-"].Update(visible=True)
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
        elif event == "-PREV_OR_HER_INPUT-":
            """Check input is a float. At the final of this loop is checked
            what fields are filled and update the GUI in consonance"""
            try:
                float(values[event])
            except ValueError:
                sg.popup_ok(
                    f"This field must be filled with a float. Try again",
                    title="Input configuration validation error",
                    font=window_general_font,
                )
                window[event].Update(value="")
        elif event == "Calculate table":
            input_mafs = []
            for i in range(1, pytoxo_context.model.order + 1):
                input_mafs.append(float(values[f"-MAFS_INPUT_{i}-"]))

            try:
                if values["-PREV_OR_HER_CB-"] == "Heritability":
                    ptable = pytoxo_context.model.find_max_prevalence_table(
                        mafs=input_mafs, h=float(values["-PREV_OR_HER_INPUT-"])
                    )
                elif values["-PREV_OR_HER_CB-"] == "Prevalence":
                    ptable = pytoxo_context.model.find_max_heritability_table(
                        mafs=input_mafs, p=float(values["-PREV_OR_HER_INPUT-"])
                    )
                # Print generated penetrance table in the GUI's table
                # noinspection PyUnboundLocalVariable
                window["-MODEL_TABLE-"].Update(
                    values=[
                        [g, p, pen]
                        for g, p, pen in zip(
                            pytoxo_context.model.calculate_genotypes(),
                            pytoxo_context.model.penetrances,
                            ptable.penetrance_values,
                        )
                    ]
                )
            except pytoxo.errors.BadFormedModelError as e:
                sg.popup_ok(
                    e.message,
                    title="Resolution error",
                    font=window_general_font,
                )
            except pytoxo.errors.UnsolvableModelError as e:
                sg.popup_ok(
                    e.message,
                    title="Unsolvable model error",
                    font=window_general_font,
                )
            except ValueError as e:
                sg.popup_ok(
                    f"{e} Check input parameters.",
                    title="Input configuration validation error",
                    font=window_general_font,
                )
        elif event in ("Exit", sg.WIN_CLOSED, None):
            break

        # Check current values state: if all configuration is filled, enable
        # calculate button, else disable it
        if window["-MODEL_TABLE-"].visible and "" not in [
            values["-PREV_OR_HER_CB-"]
        ] + [values["-PREV_OR_HER_INPUT-"]] + [
            values[f"-MAFS_INPUT_" f"{i}-"]
            for i in range(1, pytoxo_context.model.order + 1)
        ]:
            window["Calculate table"].Update(disabled=False)
        else:
            window["Calculate table"].Update(disabled=True)

    window.close()


if __name__ == "__main__":
    main()
