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
            sg.Text(
                "None model loaded",
                key="-PREV_OR_HER_DISABLED_TEXT-",
                text_color="grey",
                tooltip="You need to have set the model before setting this",
                visible=True,  # Pending to be disabled when a model was loaded
            ),
            sg.Combo(
                ("Heritability", "Prevalence"),
                key="-PREV_OR_HER_CB-",
                size=(10, 1),
                readonly=True,
                enable_events=True,  # To refresh the loop and can check filled fields
                visible=False,  # Pending to be disabled when a model was loaded
            ),
            sg.InputText(
                key="-PREV_OR_HER_INPUT-",
                size=(5, 1),
                enable_events=True,  # To refresh the loop and can check filled fields
                visible=False,  # Pending to be disabled when a model was loaded
            ),
        ],
        [
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


# ##################### GUI EVENT LOOP ####################
def main():
    while True:
        event, values = window.read()

        """First fix incomplete fields like `.2` instead of `0.2`. It must be 
        done here because then working with the current event some validations
        are done and they would fail"""
        if pytoxo_context.model and values:
            widgets_to_check = ["-PREV_OR_HER_INPUT-"] + [
                f"-MAFS_INPUT_" f"{i}-"
                for i in range(1, pytoxo_context.model.order + 1)
            ]
            for widget in widgets_to_check:
                if values[widget] == "." and event != widget:
                    window[widget].Update(value="0.0")
                    values[widget] = "0.0"
                elif values[widget].startswith(".") and event != widget:
                    window[widget].Update(value=f"0{values[widget]}")
                    values[widget] = window[widget]
                elif values[widget].endswith(".") and event != widget:
                    window[widget].Update(value=f"{values[widget]}0")
                    values[widget] = window[widget]

        # Check events
        if event in ("Exit", sg.WIN_CLOSED, None):
            break
        elif event == "Open model":
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
                for i in range(1, pytoxo_context.model.order + 1):
                    window[f"-MAFS_INPUT_{i}-"].Update(visible=True)
                window["-MAFS_DISABLED_TEXT-"].Update(visible=False)
                # Enable prevalence or heritability entry
                """Note: actually this parameter is independent of the load 
                of the model. However, thus limiting the order of filling 
                fields, it is possible to validate this field with the validator
                of the model, which otherwise might not be available. It also
                makes the user better understand the order in which the 
                interface should be used."""
                window[f"-PREV_OR_HER_CB-"].Update(visible=True)
                window[f"-PREV_OR_HER_INPUT-"].Update(visible=True)
                window["-PREV_OR_HER_DISABLED_TEXT-"].Update(visible=False)
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
        elif (
            event == "-PREV_OR_HER_INPUT-"
            and values[event] != ""
            and values[event] != "."
        ):
            """Check input is valid using 'Model' check function. At the
            final of this loop is checked what fields are filled to update
            the GUI in consonance"""
            try:
                pytoxo_context.model.check_find_table_parameters(
                    # Simply a valid value to validate th other ignoring this
                    mafs=[0.0] * pytoxo_context.model.order,
                    h_or_p=values[event],
                )
            except ValueError as e:
                sg.popup_ok(
                    f"{e} Revise this field.",
                    title="Input configuration validation error",
                    font=window_general_font,
                )
        elif (
            event.startswith("-MAFS_INPUT_")
            and values[event] != ""
            and values[event] != "."
        ):
            """Check input is valid using 'Model' check function. At the
            final of this loop is checked what fields are filled to update
            the GUI in consonance"""
            try:
                pytoxo_context.model.check_find_table_parameters(
                    mafs=values[event],
                    h_or_p=0.0,  # Simply a valid value to validate th other ignoring this
                )
            except ValueError as e:
                sg.popup_ok(
                    f"{e} Revise this field.",
                    title="Input configuration validation error",
                    font=window_general_font,
                )
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


# #########################################################


if __name__ == "__main__":
    main()
