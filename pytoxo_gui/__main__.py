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
MAX_NUMERICAL_INPUT_LEN = 20


class PyToxoContext:
    """Support container class to maintain the PyToxo stuff state during GUI
    main loop."""

    model = None  # Loaded PyToxo model
    order = 1  # Must be modified when a model is loaded. Before serves to control some workarounds with the GUI
    ptable = None  # When a penetrance table is calculated


# ####################### GUI DESIGN ######################
# General style settings
sg.theme("SystemDefaultForReal")
window_general_font = "Verdana, 13"

# Main menu
menu = sg.Menu(
    [
        ["File", ["Open model", "Close model", "Save calculated table", "Exit"]],
        # ["Edit", ["Paste", "Undo"]],  # TODO
        # ["Help", "About PyToxo GUI"],  # TODO
    ]
)

# Epistatic model table
headings = [
    " Genotype definition ",
    "Penetrance expression",
    "Calculated penetrance",
]
empty_rows = [["" for col in range(len(headings))]]

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
                values=empty_rows,
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


# ################## AUXILIARY FUNCTIONS ##################
def refresh_mafs_entries_to_check_keys(pytoxo_context: PyToxoContext) -> list[str]:
    """Auxiliary function to refresh contents used in several times within
    the GUI main loop."""
    return [f"-MAFS_INPUT_" f"{i}-" for i in range(1, pytoxo_context.order + 1)]


def refresh_text_entries_to_check_keys(
    mafs_entries_to_check_keys: list[str],
) -> list[str]:
    """Auxiliary function to refresh contents used in several times within
    the GUI main loop."""
    return ["-PREV_OR_HER_INPUT-"] + mafs_entries_to_check_keys


def refresh_text_entries_to_check_values(
    text_entries_to_check_keys: list[str], values: dict
) -> list[str]:
    """Auxiliary function to refresh contents used in several times within
    the GUI main loop."""
    return [values[k] for k in text_entries_to_check_keys]


def check_all_filled(
    window: sg.Window, values: dict, text_entries_to_check_values: list[str]
) -> bool:
    """Check current values state: if all configuration is filled, enable
    calculate button, else disable it. Also returns a bool about the all
    filled condition to add the possibility to reuse this function to a
    security check before try to calculate."""
    if (
        window["-MODEL_TABLE-"].visible
        and "" not in [values["-PREV_OR_HER_CB-"]] + text_entries_to_check_values
    ):
        window["Calculate table"].Update(disabled=False)
        return True
    else:
        window["Calculate table"].Update(disabled=True)
        return False


# #########################################################


# ##################### GUI EVENT LOOP ####################
def main():
    # Create PyToxo context object
    pytoxo_context = PyToxoContext()

    while True:
        event, values = window.read()

        # Check if exit event
        if event in ("Exit", sg.WIN_CLOSED, None):
            break

        # Load text items to check
        mafs_entries_to_check_keys = refresh_mafs_entries_to_check_keys(pytoxo_context)
        text_entries_to_check_keys = refresh_text_entries_to_check_keys(
            mafs_entries_to_check_keys
        )
        text_entries_to_check_values = refresh_text_entries_to_check_values(
            text_entries_to_check_keys, values
        )

        """Fix the input of illegal chars. In entry widgets, 
        only numerical values are allowed"""
        if pytoxo_context.model and values:
            for k in text_entries_to_check_keys:
                if len(values[k]) > 1 and (
                    values[k][-1] not in (".0123456789")  # Illegal chars
                    or (
                        "." in values[k][:-1] and values[k][-1] == "."
                    )  # Only one `.` in the entry
                    or len(values[k]) > MAX_NUMERICAL_INPUT_LEN  # Already too long
                ):
                    """Delete last char from input. The user perceives that
                    the keystroke is ignored"""
                    values[k] = values[k][:-1]
                    window[k].update(value=values[k])
                elif len(values[k]) == 1 and values[k] not in (
                    ".0123456789"
                ):  # Illegal chars
                    values[k] = ""
                    window[k].update(value=values[k])

            # Refresh text items to check
            mafs_entries_to_check_keys = refresh_mafs_entries_to_check_keys(
                pytoxo_context
            )
            text_entries_to_check_keys = refresh_text_entries_to_check_keys(
                mafs_entries_to_check_keys
            )
            text_entries_to_check_values = refresh_text_entries_to_check_values(
                text_entries_to_check_keys, values
            )

        # Beautify incomplete fields like `.2` instead of `0.2`
        if pytoxo_context.model and values:
            for k in text_entries_to_check_keys:
                if values[k] == "." and event != k:
                    # Fix `.`
                    values[k] = "0.0"
                    window[k].Update(value=values[k])
                elif values[k] != "" and event != k:
                    # Fix e.g. `00.1`, `0.4600000`, `1.` or `.23'
                    values[k] = str(float(values[k]))
                    window[k].update(value=values[k])

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
                # Load to the PyToxo context
                pytoxo_context.model = pytoxo.Model(filename)
                pytoxo_context.order = (
                    pytoxo_context.model.order
                )  # It is important to update also this

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

                # Refresh text items to check
                mafs_entries_to_check_keys = refresh_mafs_entries_to_check_keys(
                    pytoxo_context
                )
                text_entries_to_check_keys = refresh_text_entries_to_check_keys(
                    mafs_entries_to_check_keys
                )
                text_entries_to_check_values = refresh_text_entries_to_check_values(
                    text_entries_to_check_keys, values
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

        elif event == "Close model":
            # Load to the PyToxo context
            pytoxo_context = PyToxoContext()  # Refresh with a new instance

            # Remove the previous model from the GUI
            window["-MODEL_DISABLED_TEXT-"].Update(visible=True)
            window["-MODEL_TABLE-"].Update(visible=False)
            window["-MODEL_TABLE-"].Update(values=empty_rows)

            # Disable MAFs and clean values
            for k in mafs_entries_to_check_keys:
                window[k].Update(visible=False)
                values[k] = ""
                window[k].Update(value=values[k])
            window["-MAFS_DISABLED_TEXT-"].Update(visible=True)

            # Disable prevalence or heritability entry and clean values
            """Note: actually this parameter is independent of the load 
            of the model. However, thus limiting the order of filling 
            fields, it is possible to validate this field with the validator
            of the model, which otherwise might not be available. It also
            makes the user better understand the order in which the 
            interface should be used."""
            window[f"-PREV_OR_HER_CB-"].Update(visible=False)
            values[f"-PREV_OR_HER_CB-"] = ""
            window[f"-PREV_OR_HER_CB-"].Update(value=values[f"-PREV_OR_HER_CB-"])
            window[f"-PREV_OR_HER_INPUT-"].Update(visible=False)
            values[f"-PREV_OR_HER_INPUT-"] = ""
            window[f"-PREV_OR_HER_INPUT-"].Update(value=values[f"-PREV_OR_HER_INPUT-"])
            window["-PREV_OR_HER_DISABLED_TEXT-"].Update(visible=True)

            # Refresh text items to check
            mafs_entries_to_check_keys = refresh_mafs_entries_to_check_keys(
                pytoxo_context
            )
            text_entries_to_check_keys = refresh_text_entries_to_check_keys(
                mafs_entries_to_check_keys
            )
            text_entries_to_check_values = refresh_text_entries_to_check_values(
                text_entries_to_check_keys, values
            )

        elif event == "Save calculated table":
            if not pytoxo_context.ptable:
                sg.popup_ok(
                    f"There is not a calculated penetrance table. Calculate "
                    f"the table before trying to save it.",
                    title="No penetrance table",
                    font=window_general_font,
                )
            else:
                filename = sg.popup_get_file(
                    "Save calculated table",
                    save_as=True,
                    default_extension="csv",
                    no_window=True,  # To use a native approach
                    file_types=(("Comma separated values", "*.csv"),),
                )
                if not filename:
                    continue  # The operation has been canceled
                else:
                    print(filename)
        elif (
            event == "-PREV_OR_HER_INPUT-"
            and values[event] != ""
            and values[event]
            != "."  # The user is already writing, or it will solved in the next interaction
        ):
            """Check input is valid using 'Model' check function. At the
            final of this loop is checked what fields are filled to update
            the GUI in consonance"""
            try:
                pytoxo_context.model.check_find_table_parameters(
                    # Simply a valid value to validate th other ignoring this
                    mafs=[0.0] * pytoxo_context.model.order,
                    h_or_p=float(
                        values[event]
                    ),  # Try cast because actually is a string
                )
            except ValueError as e:
                # Format error message
                msg = e.__str__()
                if not msg.endswith("."):
                    msg = f"{e}."
                msg = msg.capitalize()

                sg.popup_ok(
                    f"{msg} Revise this field.",
                    title="Input configuration validation error",
                    font=window_general_font,
                )
                # Remove value
                window[event].Update(value="")

        elif (
            event.startswith("-MAFS_INPUT_")
            and values[event] != ""
            and values[event]
            != "."  # The user is already writing, or it will solved in the next interaction
        ):
            """Check input is valid using 'Model' check function. At the
            final of this loop is checked what fields are filled to update
            the GUI in consonance"""
            try:
                pytoxo_context.model.check_find_table_parameters(
                    mafs=[float(values[event])]
                    * pytoxo_context.model.order,  # Try cast because actually is a string
                    h_or_p=0.0,  # Simply a valid value to validate th other ignoring this
                )
            except ValueError as e:
                # Format error message
                msg = e.__str__()
                if not msg.endswith("."):
                    msg = f"{e}."
                msg = msg.capitalize()

                sg.popup_ok(
                    f"{msg} Revise this field.",
                    title="Input configuration validation error",
                    font=window_general_font,
                )
                # Remove value
                window[event].Update(value="")

        elif event == "Calculate table":
            """Security check to avoid some stranger cases playing with the
            input fields"""
            if check_all_filled(window, values, text_entries_to_check_values):
                input_mafs = []
                for k in mafs_entries_to_check_keys:
                    input_mafs.append(float(values[k]))

                try:
                    if values["-PREV_OR_HER_CB-"] == "Heritability":
                        ptable = pytoxo_context.model.find_max_prevalence_table(
                            mafs=input_mafs, h=float(values["-PREV_OR_HER_INPUT-"])
                        )
                    elif values["-PREV_OR_HER_CB-"] == "Prevalence":
                        ptable = pytoxo_context.model.find_max_heritability_table(
                            mafs=input_mafs, p=float(values["-PREV_OR_HER_INPUT-"])
                        )

                    # Load the table to the PyToxo context
                    # noinspection PyUnboundLocalVariable
                    pytoxo_context.ptable = ptable

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
                except pytoxo.errors.ResolutionError as e:
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

        # Finally, check if all is filled before go to the next interaction
        check_all_filled(window, values, text_entries_to_check_values)

    window.close()


# #########################################################


if __name__ == "__main__":
    main()
