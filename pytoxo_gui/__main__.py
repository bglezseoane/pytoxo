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

MAX_ORDER_SUPPORTED = 12  # The interface would need some fixes to support bigger orders
MAX_NUMERICAL_INPUT_LEN = 20
INFO_BANNER_MAX_MODEL_NAME_LEN = 20


class PyToxoContext:
    """Support container class to maintain the PyToxo stuff state during GUI
    main loop."""

    model = None  # Loaded PyToxo model
    order = 1  # Must be modified when a model is loaded. Before serves to control some workarounds with the GUI
    ptable = None  # When a penetrance table is calculated


# ####################### GUI DESIGN ######################
# General style settings
sg.theme("DarkGreen4")
window_general_font = "13"
table_font = ("Courier", "15")
main_background_color = "#044343"  # Default of `DarkGreen4` theme
table_font_color = "black"
table_headers_font = window_general_font
table_headers_font_color = "white"
table_background_color = "white"
table_alternating_background_color = "lightgrey"
text_inputs_text_color = "black"
text_inputs_background_color = "white"
selection_colors = ("black", "#c6dffc")
disabled_text_color = "grey"

# Tooltipis
tt_model_disabled_text = (
    "Use the menu to load a model from a file or enter one manually"
)
tt_prev_or_her_input_dis = "You need to have set the model before setting this"
tt_prev_or_her_input_en = "A probability value between 0 and 1"
tt_mafs_disabled_text = "You need to have set the model before setting MAFs"
tt_mafs_input = "The minor allele frequency, a value between 0 and 0.5"
tt_calculate_button_dis = "You need to have set all fields before calculating the table"
tt_calculate_button_en = "Calculate the table with the set configuration"

# Main menu
menu = sg.Menu(
    [
        [
            "File",
            ["Open model", "Close model and clean", "Save calculated table", "Exit"],
        ],
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
size_y_model_frame_component = 28  # Used by the disabled text or the table
model_frame = sg.Frame(
    key="-MODEL_FRAME-",
    title="Epistatic model",
    layout=[
        [
            sg.Table(
                key="-MODEL_TABLE-",
                headings=headings,
                values=empty_rows,
                enable_events=True,
                vertical_scroll_only=True,
                num_rows=size_y_model_frame_component - 1,  # Due to header row
                justification="center",
                display_row_numbers=True,
                hide_vertical_scroll=True,
                auto_size_columns=True,
                header_font=table_headers_font,
                header_text_color=table_headers_font_color,
                font=table_font,
                text_color=table_font_color,
                background_color=table_background_color,
                alternating_row_color=table_alternating_background_color,
                selected_row_colors=selection_colors,
            ),
        ]
    ],
    element_justification="center",
)

# Informative banner
info_banner_model_head_text = "Model: "
info_banner_order_head_text = "Order: "
info_banner_fixing_head_text = "Fixing: "
info_banner_maximizing_head_text = "Maximizing: "
info_banner_fixing_maximizing_options = ["prevalence", "heritability"]
info_banner_model_none_text = f"{info_banner_model_head_text}none"
info_banner_order_none_text = f"{info_banner_order_head_text}none"
info_banner_fixing_none_text = f"{info_banner_fixing_head_text}none"
info_banner_maximizing_none_text = f"{info_banner_maximizing_head_text}none"
info_banner = [
    [
        sg.Text(
            info_banner_model_none_text,
            key="-INFO_MODEL-",
            size=(INFO_BANNER_MAX_MODEL_NAME_LEN + len(info_banner_model_head_text), 1),
            text_color=disabled_text_color,
            justification="center",
        ),
        sg.Text(
            info_banner_order_none_text,
            key="-INFO_ORDER-",
            text_color=disabled_text_color,
            justification="center",
        ),
        sg.Text(
            info_banner_fixing_none_text,
            key="-INFO_FIXING-",
            size=(
                len(info_banner_fixing_head_text)
                + len(max(info_banner_fixing_maximizing_options)),
                1,
            ),
            text_color=disabled_text_color,
            justification="center",
        ),
        sg.Text(
            info_banner_maximizing_none_text,
            key="-INFO_MAXIMIZING-",
            size=(
                len(info_banner_maximizing_head_text)
                + len(max(info_banner_fixing_maximizing_options)),
                1,
            ),
            text_color=disabled_text_color,
            justification="center",
        ),
    ],
    sg.Text(
        "State: Ready",
        key="-INFO_STATE_READY-",
        visible=True,
        text_color="green",
        justification="center",
    ),
    sg.Text(
        "State: Calculating...",
        key="-INFO_STATE_CALCULATING-",
        visible=False,
        text_color="#f29114",
        justification="center",
    ),
]

# Prevalence or heritability frame
prev_or_her_texts = ["Heritability", "Prevalence"]
prev_or_her_frame = sg.Frame(
    key="-PREV_OR_HER_FRAME-",
    title="Fix prevalence or heritability",
    layout=[
        [
            sg.Combo(
                prev_or_her_texts,
                key="-PREV_OR_HER_CB-",
                readonly=True,
                enable_events=True,  # To refresh the loop and can check filled fields
                size=(9, 1),
                text_color=text_inputs_text_color,
                background_color=text_inputs_background_color,
            ),
            sg.InputText(
                key="-PREV_OR_HER_INPUT-",
                disabled=True,  # Pending to be disabled when a model was loaded
                enable_events=True,  # To refresh the loop and can check filled fields
                tooltip=tt_prev_or_her_input_dis,
                size=(4, 1),
                text_color=text_inputs_text_color,
                background_color=text_inputs_background_color,
            ),
        ],
    ],
    element_justification="left",
)

# MAFs frame text entries
mafs_entries = []
for i in range(1, MAX_ORDER_SUPPORTED + 1):
    mafs_entries.append(
        sg.pin(  # To fix the entries in the layout, in horizontal
            sg.InputText(
                key=f"-MAFS_INPUT_{i}-",
                visible=False,  # Pending to be enabled when a model was loaded
                enable_events=True,  # To refresh the loop and can check filled fields
                tooltip=tt_mafs_input,
                size=(4, 1),
                pad=(1, 3),  # 3 seems to be the default
                text_color=text_inputs_text_color,
                background_color=text_inputs_background_color,
            )
        )
    )

# MAFs frame
mafs_frame = sg.Frame(
    key="-MAFS_FRAME-",
    title="MAFs",
    layout=[
        [
            sg.Text(
                "None model loaded",
                key="-MAFS_DISABLED_TEXT-",
                visible=True,  # Pending to be disabled when a model was loaded
                tooltip=tt_mafs_disabled_text,
                text_color=disabled_text_color,
            ),
        ]
        + mafs_entries
    ],
    element_justification="center",
)

# Layout composition
layout = [
    [menu],
    [model_frame],
    [info_banner],
    [prev_or_her_frame, mafs_frame],
    [
        sg.Button(
            "Calculate table",
            disabled=True,
            tooltip=tt_calculate_button_dis,
        )
    ],
]

# Window composition
window = sg.Window(
    "PyToxo GUI",
    layout,
    font=window_general_font,
    finalize=True,
    size=(650, 650),
    element_justification="center",
)

# Window style patches
window["-MODEL_FRAME-"].expand(expand_x=True)


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
        window["Calculate table"].set_tooltip(tt_calculate_button_en)
        return True
    else:
        window["Calculate table"].Update(disabled=True)
        window["Calculate table"].set_tooltip(tt_calculate_button_dis)
        return False


def update_info_banner(window: sg.Window, key: str, *args: str) -> None:
    """Auxiliary function to update the informative banner. The different
    keys have different treatments. Supported keys are: `-INFO_MODEL-`,
    `-INFO_ORDER-`, `-INFO_FIXING-`, `-INFO_MAXIMIZING-`, `-INFO_STATE_READY-`,
    `-INFO_STATE_CALCULATING-` and `CLEAN`. `CLEAN` key serves to clean the
    banner as default."""
    if key == "-INFO_MODEL-":
        if len(args[0]) > INFO_BANNER_MAX_MODEL_NAME_LEN:
            name = args[0][:INFO_BANNER_MAX_MODEL_NAME_LEN]
        else:
            name = args[0]
        window[key].Update(value=f"{info_banner_model_head_text}{name}")
    elif key == "-INFO_ORDER-":
        window[key].Update(value=f"{info_banner_order_head_text}{args[0]}")
    elif key == "-INFO_FIXING-":
        window[key].Update(value=f"{info_banner_fixing_head_text}{args[0]}")
    elif key == "-INFO_MAXIMIZING-":
        window[key].Update(value=f"{info_banner_maximizing_head_text}{args[0]}")
    elif key == "-INFO_STATE_READY-":
        window[key].Update(visible=True)
        window["-INFO_STATE_CALCULATING-"].Update(visible=False)
    elif key == "-INFO_STATE_CALCULATING-":
        window[key].Update(visible=True)
        window["-INFO_STATE_READY-"].Update(visible=False)
        window.refresh()  # Needed here because it does not go through the loop
    elif key == "CLEAN":
        window["-INFO_MODEL-"].Update(value=f"{info_banner_model_none_text}")
        window["-INFO_ORDER-"].Update(value=f"{info_banner_model_none_text}")
        window.refresh()  # Needed here because it does not go through the loop
    else:
        raise ValueError(key)


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
                window[f"-PREV_OR_HER_INPUT-"].Update(disabled=False)
                window["-PREV_OR_HER_INPUT-"].set_tooltip(tt_prev_or_her_input_en)

                # Update informative banner
                update_info_banner(window, "-INFO_MODEL-", pytoxo_context.model.name)
                update_info_banner(window, "-INFO_ORDER-", str(pytoxo_context.order))

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
                window.Hide()  # Patch because popup modal function does not work in all platforms
                sg.popup_ok(
                    "The file contains a bad formed model. PyToxo cannot "
                    "interpret it. Revise PyToxo's file format requirements.",
                    title="File parsing error",
                    font=window_general_font,
                )
                window.UnHide()
            except IOError:
                window.Hide()  # Patch because popup modal function does not work in all platforms
                sg.popup_ok(
                    f"Error trying to open '{filename}'.",
                    title="File opening error",
                    font=window_general_font,
                )
                window.UnHide()

        elif event == "Close model and clean":
            # Load to the PyToxo context
            pytoxo_context = PyToxoContext()  # Refresh with a new instance

            # Remove the previous model from the GUI
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
            values[f"-PREV_OR_HER_CB-"] = ""
            window[f"-PREV_OR_HER_CB-"].Update(value=values[f"-PREV_OR_HER_CB-"])
            window[f"-PREV_OR_HER_INPUT-"].Update(disabled=True)
            window["-PREV_OR_HER_INPUT-"].set_tooltip(tt_prev_or_her_input_dis)
            values[f"-PREV_OR_HER_INPUT-"] = ""
            window[f"-PREV_OR_HER_INPUT-"].Update(value=values[f"-PREV_OR_HER_INPUT-"])

            # Update informative banner
            update_info_banner(window, "CLEAN")  # Clean all model related

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
                window.Hide()  # Patch because popup modal function does not work in all platforms
                sg.popup_ok(
                    f"There is not a calculated penetrance table. Calculate "
                    f"the table before trying to save it.",
                    title="No penetrance table",
                    font=window_general_font,
                )
                window.UnHide()
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
                    pytoxo_context.ptable.write_to_file(
                        filename=filename, overwrite=True, format="csv"
                    )

        elif event == "-MODEL_TABLE-":
            # Intercept the event to refresh window when click on the table
            pass

        elif event == "-PREV_OR_HER_CB-" and values[event] != "":
            # Update informative banner
            fixed_prev_or_her = values[event].lower()
            if fixed_prev_or_her == info_banner_fixing_maximizing_options[0]:
                maximized_prev_or_her = info_banner_fixing_maximizing_options[1]
            else:
                maximized_prev_or_her = info_banner_fixing_maximizing_options[0]
            update_info_banner(window, "-INFO_FIXING-", fixed_prev_or_her)
            update_info_banner(window, "-INFO_MAXIMIZING-", maximized_prev_or_her)

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

                window.Hide()  # Patch because popup modal function does not work in all platforms
                sg.popup_ok(
                    f"{msg} Revise this field.",
                    title="Input configuration validation error",
                    font=window_general_font,
                )
                window.UnHide()

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

                window.Hide()  # Patch because popup modal function does not work in all platforms
                sg.popup_ok(
                    f"{msg} Revise this field.",
                    title="Input configuration validation error",
                    font=window_general_font,
                )
                window.UnHide()

                # Remove value
                window[event].Update(value="")

        elif event == "Calculate table":
            """Security check to avoid some stranger cases playing with the
            input fields"""
            if check_all_filled(window, values, text_entries_to_check_values):
                input_mafs = []
                for k in mafs_entries_to_check_keys:
                    input_mafs.append(float(values[k]))

                # Update informative banner
                update_info_banner(window, "-INFO_STATE_CALCULATING-")
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
                    window.Hide()  # Patch because popup modal function does not work in all platforms
                    sg.popup_ok(
                        e.message,
                        title="Resolution error",
                        font=window_general_font,
                    )
                    window.UnHide()
                except pytoxo.errors.UnsolvableModelError as e:
                    window.Hide()  # Patch because popup modal function does not work in all platforms
                    sg.popup_ok(
                        e.message,
                        title="Unsolvable model error",
                        font=window_general_font,
                    )
                    window.UnHide()
                except ValueError as e:
                    window.Hide()  # Patch because popup modal function does not work in all platforms
                    sg.popup_ok(
                        f"{e} Check input parameters.",
                        title="Input configuration validation error",
                        font=window_general_font,
                    )
                    window.UnHide()
                finally:
                    # Update informative banner
                    update_info_banner(window, "-INFO_STATE_READY-")

        # Finally, check if all is filled before go to the next interaction
        check_all_filled(window, values, text_entries_to_check_values)

    window.close()


# #########################################################


if __name__ == "__main__":
    main()
