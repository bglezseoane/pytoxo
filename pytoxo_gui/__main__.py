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

import base64
import io
import os
import platform
import subprocess
from typing import List

import PIL.Image

import pytoxo
import pytoxo.errors

# Detect the platform where the GUI is going to be used
detected_platform = platform.system()

"""PySimpleGUI imports Tk, and this fail if Tk is not correctly installed in 
the platform"""
try:
    import PySimpleGUI as sg
except ImportError:
    print(pytoxo.errors.GUIUnsupportedPlatformError(detected_platform).message)
    exit(1)

MAX_ORDER_SUPPORTED = 12  # The interface would need some fixes to support bigger orders
MAX_NUMERICAL_INPUT_LEN = 20
INFO_BANNER_MAX_MODEL_NAME_LEN = 18


class PyToxoContext:
    """Support container class to maintain the PyToxo stuff state during GUI
    main loop."""

    model = None  # Loaded PyToxo model
    order = 1  # Must be modified when a model is loaded. Before serves to control some workarounds with the GUI
    ptable = None  # When a penetrance table is calculated


# ####################### GUI DESIGN ######################
# Adapt some patches depending of the detected platform
if detected_platform == "Darwin":  # Mac OS
    window_general_font_dep_of_platform = ("", "13")
    state_ready_font_size_dep_of_platform = 15
    state_calculating_font_size_dep_of_platform = state_ready_font_size_dep_of_platform
    table_font_size_dep_of_platform = state_ready_font_size_dep_of_platform
    window_size_dep_of_platform = (650, 800)
    model_frame_size_y_dep_of_platform = 27
    """Modal windows are not supported on Mac OS, so this patch a used to 
    emulate them hiding and un-hiding the main window"""
    hide_windows_to_emulate_modal_dep_of_platform = True
    """With Mac OS the menu is located in the upper menu bar, and it is 
    sensitive to the theme that is being used in the system, so it is consulted
    below"""
    p = subprocess.Popen(
        "defaults read -g AppleInterfaceStyle", shell=True, stdout=subprocess.PIPE
    )
    output, _ = p.communicate()
    macos_current_theme = output.decode("UTF-8").strip()
    if macos_current_theme == "Dark":
        menu_text_color_dep_of_platform = "#ffffff"
    else:
        menu_text_color_dep_of_platform = "#000000"
elif detected_platform == "Linux":
    window_general_font_dep_of_platform = ("", "10")
    state_ready_font_size_dep_of_platform = 11
    state_calculating_font_size_dep_of_platform = state_ready_font_size_dep_of_platform
    table_font_size_dep_of_platform = state_ready_font_size_dep_of_platform
    window_size_dep_of_platform = (650, 780)
    model_frame_size_y_dep_of_platform = 27
    hide_windows_to_emulate_modal_dep_of_platform = False
    menu_text_color_dep_of_platform = "#000000"
elif detected_platform == "Windows":
    window_general_font_dep_of_platform = ("", "10")
    state_ready_font_size_dep_of_platform = 11
    state_calculating_font_size_dep_of_platform = state_ready_font_size_dep_of_platform
    table_font_size_dep_of_platform = state_ready_font_size_dep_of_platform
    window_size_dep_of_platform = (650, 780)
    model_frame_size_y_dep_of_platform = 24
    hide_windows_to_emulate_modal_dep_of_platform = False
    menu_text_color_dep_of_platform = "#000000"
else:
    raise pytoxo.errors.GUIUnsupportedPlatformError(detected_platform)

# Main style settings
pytoxo_main_color = "#044343"
sg.LOOK_AND_FEEL_TABLE["PyToxoTheme"] = {
    "BACKGROUND": "#4f7b7b",
    "TEXT": "#e4e4e4",
    "INPUT": "#ffffff",
    "TEXT_INPUT": "#e4e4e4",
    "SCROLL": "#045757",
    "BUTTON": ("#e4e4e4", pytoxo_main_color),
    "PROGRESS": ("#000000", "#000000"),
    "BORDER": 1,
    "SLIDER_DEPTH": 0,
    "PROGRESS_DEPTH": 0,
    "COLOR_LIST": ["#222222", pytoxo_main_color, "#045757", "#e4e4e4"],
    "DESCRIPTION": ["#000000", "Turquoise", "Grey", "Dark"],
}  # Custom window style for PyToxo
sg.LOOK_AND_FEEL_TABLE["PyToxoLightTheme"] = {
    "BACKGROUND": "#a7bdbd",
    "TEXT": pytoxo_main_color,
    "INPUT": "#ffffff",
    "TEXT_INPUT": "#e4e4e4",
    "SCROLL": "#045757",
    "BUTTON": ("#e4e4e4", pytoxo_main_color),
    "PROGRESS": ("#000000", "#000000"),
    "BORDER": 1,
    "SLIDER_DEPTH": 0,
    "PROGRESS_DEPTH": 0,
    "COLOR_LIST": ["#222222", pytoxo_main_color, "#045757", "#e4e4e4"],
    "DESCRIPTION": ["#000000", "Turquoise", "Grey", "Dark"],
}  # Custom window style for PyToxo, light alternative

# ########## THEME SELECTION ##########
# Uncomment desired style between:
# Dark style
sg.theme("PyToxoTheme")
disabled_text_color_dep_of_style = "#dbdbdb"
state_ready_color_dep_of_style = "#1bff1B"
state_calculating_color_dep_of_style = "#e59400"
logo_pytoxo_dep_of_style = os.path.abspath("img/logo_white.gif")
logo_udc_dep_of_style = os.path.abspath("img/logo_udc.gif")
# # Light style
# sg.theme("PyToxoLightTheme")
# disabled_text_color_dep_of_style = "#656b6b"
# state_ready_color_dep_of_style = "#8cff19"
# state_calculating_color_dep_of_style = "#ee683b"
# logo_pytoxo_dep_of_style = os.path.abspath("img/logo.gif")
# logo_udc_dep_of_style =os.path.abspath("img/logo_udc_green.gif")
# #####################################

# Other style settings: fonts
window_general_font = window_general_font_dep_of_platform
table_font = ("Courier", table_font_size_dep_of_platform)
table_headers_font = window_general_font
state_ready_font = ("", state_ready_font_size_dep_of_platform, "bold")
state_calculating_font = ("", state_calculating_font_size_dep_of_platform, "bold")

# Other style settings: font colors
table_font_color = "#000000"
table_headers_font_color = "#ffffff"
text_inputs_text_color = "#000000"
state_ready_color = state_ready_color_dep_of_style
state_calculating_color = state_calculating_color_dep_of_style
disabled_text_color = disabled_text_color_dep_of_style

# Other style settings: background colors
table_background_color = "#ffffff"
table_alternating_background_color = "#d3d3d3"

# Other style settings: highlight colors
selection_colors = ("#000000", "#c6dffc")

# Tooltipis messages
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
        ["Help", "About PyToxo GUI"],
    ],
    text_color=menu_text_color_dep_of_platform,
)

# Logo preparation for the main window
logo = PIL.Image.open(logo_pytoxo_dep_of_style)
logo_size_x, logo_size_y = logo.size
logo_new_size_x = 450
logo_new_size = (logo_new_size_x, (logo_size_y * logo_new_size_x) // logo_size_x)
logo = logo.resize(logo_new_size, PIL.Image.ANTIALIAS)
buffered = io.BytesIO()
logo.save(buffered, format="GIF")  # GIF is the best format for Tkinter
logo_b64 = base64.b64encode(buffered.getvalue())

# Logo preparation for the about pop-up
logo_popup_size_x = 300
logo_popup_size = (logo_popup_size_x, (logo_size_y * logo_popup_size_x) // logo_size_x)
logo = logo.resize(logo_popup_size, PIL.Image.ANTIALIAS)
buffered = io.BytesIO()
logo.save(buffered, format="GIF")  # GIF is the best format for Tkinter
logo_b64_popup = base64.b64encode(buffered.getvalue())

# UDC's logo preparation for the about pop-up
logo_udc = PIL.Image.open(logo_udc_dep_of_style)
logo_udc_size_x, logo_udc_size_y = logo_udc.size
logo_udc_new_size_x = 300
logo_udc_new_size = (
    logo_udc_new_size_x,
    (logo_udc_size_y * logo_udc_new_size_x) // logo_udc_size_x,
)
logo_udc = logo_udc.resize(logo_udc_new_size, PIL.Image.ANTIALIAS)
buffered = io.BytesIO()
logo_udc.save(buffered, format="GIF")  # GIF is the best format for Tkinter
logo_udc_b64 = base64.b64encode(buffered.getvalue())

# Epistatic model table
headings = [
    "Genotype definition",
    "Penetrance expression",
    "Calculated penetrance",
]
empty_rows = [["" for col in range(len(headings))]]

# Model frame
model_frame_size_y = model_frame_size_y_dep_of_platform
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
                num_rows=model_frame_size_y - 1,  # Due to header row
                justification="center",
                display_row_numbers=True,
                hide_vertical_scroll=True,
                auto_size_columns=False,
                col_widths=[21] * len(headings),
                header_font=table_headers_font,
                header_text_color=table_headers_font_color,
                header_background_color=pytoxo_main_color,
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
info_banner_model_skeleton_text = "Model: {} (order {})"
info_banner_fixing_head_text = "Fixing: "
info_banner_maximizing_head_text = "Maximizing: "
info_banner_fixing_maximizing_options = ["prevalence", "heritability"]
info_banner_model_none_text = f"{info_banner_model_skeleton_text.format('none','none')}"
info_banner_fixing_none_text = f"{info_banner_fixing_head_text}none"
info_banner_maximizing_none_text = f"{info_banner_maximizing_head_text}none"
info_banner = [
    [
        sg.Text(
            info_banner_model_none_text,
            key="-INFO_MODEL-",
            size=(25, 1),
            text_color=disabled_text_color,
            justification="center",
        ),
        sg.Text(
            info_banner_fixing_none_text,
            key="-INFO_FIXING-",
            size=(20, 1),
            text_color=disabled_text_color,
            justification="center",
        ),
        sg.Text(
            info_banner_maximizing_none_text,
            key="-INFO_MAXIMIZING-",
            size=(25, 1),
            text_color=disabled_text_color,
            justification="center",
        ),
    ],
    sg.Text(
        "Ready",
        key="-INFO_STATE_READY-",
        visible=True,
        font=state_ready_font,
        text_color=state_ready_color,
        justification="center",
    ),
    sg.Text(
        "Calculating...",
        key="-INFO_STATE_CALCULATING-",
        visible=False,
        font=state_calculating_font,
        text_color=state_calculating_color,
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
            ),
            sg.InputText(
                key="-PREV_OR_HER_INPUT-",
                disabled=True,  # Pending to be disabled when a model was loaded
                enable_events=True,  # To refresh the loop and can check filled fields
                tooltip=tt_prev_or_her_input_dis,
                size=(4, 1),
                text_color=text_inputs_text_color,
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

# Available output formats
output_formats = ["GAMETES", "CSV"]

# Layout composition
layout = [
    [menu],
    [sg.Image(data=logo_b64, key="-LOGO-")],
    [model_frame],
    [info_banner],
    [prev_or_her_frame, mafs_frame],
    [
        sg.Button(
            "Calculate table",
            disabled=True,
            tooltip=tt_calculate_button_dis,
        ),
    ],
    [
        sg.Text("Save calculated table as"),
        sg.Combo(
            output_formats,
            key="-FORMATS_CB-",
            readonly=True,
            default_value=output_formats[0],
            size=(9, 1),
            text_color=text_inputs_text_color,
        ),
    ],
]

# Window composition
window = sg.Window(
    "PyToxo GUI",
    layout,
    font=window_general_font,
    finalize=True,
    size=window_size_dep_of_platform,
    element_justification="center",
)


# #########################################################


# ################## AUXILIARY FUNCTIONS ##################
def refresh_mafs_entries_to_check_keys(pytoxo_context: PyToxoContext) -> List[str]:
    """Auxiliary function to refresh contents used in several times within
    the GUI main loop."""
    return [f"-MAFS_INPUT_" f"{i}-" for i in range(1, pytoxo_context.order + 1)]


def refresh_text_entries_to_check_keys(
    mafs_entries_to_check_keys: List[str],
) -> List[str]:
    """Auxiliary function to refresh contents used in several times within
    the GUI main loop."""
    return ["-PREV_OR_HER_INPUT-"] + mafs_entries_to_check_keys


def refresh_text_entries_to_check_values(
    text_entries_to_check_keys: List[str], values: dict
) -> List[str]:
    """Auxiliary function to refresh contents used in several times within
    the GUI main loop."""
    return [values[k] for k in text_entries_to_check_keys]


def check_all_filled(
    window: sg.Window, values: dict, text_entries_to_check_values: List[str]
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
    `-INFO_FIXING-`, `-INFO_MAXIMIZING-`, `-INFO_STATE_READY-`,
    `-INFO_STATE_CALCULATING-` and `CLEAN`. `CLEAN` key serves to clean the
    banner as default."""
    if key == "-INFO_MODEL-":
        if len(args[0]) > INFO_BANNER_MAX_MODEL_NAME_LEN:
            name = f"{args[0][:INFO_BANNER_MAX_MODEL_NAME_LEN-3]}..."
        else:
            name = args[0]
        window[key].Update(
            value=f"{info_banner_model_skeleton_text.format(name,args[1])}"
        )
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
                modal=True,  # Work in all platforms (native approach)
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
                """Enable MAFs assuring to prevent from being enabled too 
                many when changing from a larger model to a smaller one. 
                Assure also that MAFs are all empty"""
                for i in range(1, pytoxo_context.model.order + 1):
                    values[f"-MAFS_INPUT_{i}-"] = ""
                    window[f"-MAFS_INPUT_{i}-"].Update(
                        value=values[f"-MAFS_INPUT_{i}-"]
                    )
                    window[f"-MAFS_INPUT_{i}-"].Update(visible=True)
                for i in range(pytoxo_context.model.order + 1, MAX_ORDER_SUPPORTED + 1):
                    window[f"-MAFS_INPUT_{i}-"].Update(visible=False)
                    values[f"-MAFS_INPUT_{i}-"] = ""
                    window[f"-MAFS_INPUT_{i}-"].Update(
                        value=values[f"-MAFS_INPUT_{i}-"]
                    )
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
                update_info_banner(
                    window,
                    "-INFO_MODEL-",
                    pytoxo_context.model.name,
                    str(pytoxo_context.order),
                )

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
                if hide_windows_to_emulate_modal_dep_of_platform:
                    window.Hide()
                sg.popup_ok(
                    "The file contains a bad formed model. PyToxo cannot "
                    "interpret it. Revise PyToxo's file format requirements.",
                    title="File parsing error",
                    modal=True,
                    font=window_general_font,
                )
                if hide_windows_to_emulate_modal_dep_of_platform:
                    window.UnHide()
            except IOError:
                if hide_windows_to_emulate_modal_dep_of_platform:
                    window.Hide()
                sg.popup_ok(
                    f"Error trying to open '{filename}'.",
                    title="File opening error",
                    modal=True,
                    font=window_general_font,
                )
                if hide_windows_to_emulate_modal_dep_of_platform:
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
                if hide_windows_to_emulate_modal_dep_of_platform:
                    window.Hide()
                sg.popup_ok(
                    f"There is not a calculated penetrance table. Calculate "
                    f"the table before trying to save it.",
                    title="No penetrance table",
                    modal=True,
                    font=window_general_font,
                )
                if hide_windows_to_emulate_modal_dep_of_platform:
                    window.UnHide()
            else:
                # Get configured output format
                output_format = values["-FORMATS_CB-"].lower()
                # Calculate default output file extension attending to output format
                if output_format == "gametes":
                    output_default_extension = ".txt"
                else:
                    output_default_extension = ".csv"

                filename = sg.popup_get_file(
                    "Save calculated table",
                    save_as=True,
                    default_extension=output_default_extension,
                    no_window=True,  # To use a native approach
                    modal=True,  # Work in all platforms (native approach)
                )
                if not filename:
                    continue  # The operation has been canceled
                else:
                    try:
                        pytoxo_context.ptable.write_to_file(
                            filename=filename, overwrite=True, format=output_format
                        )
                    except pytoxo.errors.GenericCalculationError as e:  # Improvable exception
                        # Format error message
                        msg = e.__str__()
                        if not msg.endswith("."):
                            msg = f"{e}."
                        msg = msg.capitalize()

                        if hide_windows_to_emulate_modal_dep_of_platform:
                            window.Hide()
                        sg.popup_ok(
                            f"{msg} Revise this field.",
                            title="Saving error",
                            modal=True,
                            font=window_general_font,
                        )
                        if hide_windows_to_emulate_modal_dep_of_platform:
                            window.UnHide()

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

                if hide_windows_to_emulate_modal_dep_of_platform:
                    window.Hide()
                sg.popup_ok(
                    f"{msg} Revise this field.",
                    title="Input configuration validation error",
                    modal=True,
                    font=window_general_font,
                )
                if hide_windows_to_emulate_modal_dep_of_platform:
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

                if hide_windows_to_emulate_modal_dep_of_platform:
                    window.Hide()
                sg.popup_ok(
                    f"{msg} Revise this field.",
                    title="Input configuration validation error",
                    modal=True,
                    font=window_general_font,
                )
                if hide_windows_to_emulate_modal_dep_of_platform:
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
                    if hide_windows_to_emulate_modal_dep_of_platform:
                        window.Hide()
                    sg.popup_ok(
                        e.message,
                        title="Resolution error",
                        modal=True,
                        font=window_general_font,
                    )
                    if hide_windows_to_emulate_modal_dep_of_platform:
                        window.UnHide()
                except pytoxo.errors.UnsolvableModelError as e:
                    if hide_windows_to_emulate_modal_dep_of_platform:
                        window.Hide()
                    sg.popup_ok(
                        e.message,
                        title="Unsolvable model error",
                        modal=True,
                        font=window_general_font,
                    )
                    if hide_windows_to_emulate_modal_dep_of_platform:
                        window.UnHide()
                except ValueError as e:
                    if hide_windows_to_emulate_modal_dep_of_platform:
                        window.Hide()
                    sg.popup_ok(
                        f"{e} Check input parameters.",
                        title="Input configuration validation error",
                        modal=True,
                        font=window_general_font,
                    )
                    if hide_windows_to_emulate_modal_dep_of_platform:
                        window.UnHide()
                finally:
                    # Update informative banner
                    update_info_banner(window, "-INFO_STATE_READY-")

        elif event == "About PyToxo GUI":
            """This patch is to emulate window modal behaviour on Mac OS. In
            Linux and Windows, this is innocuous, because the "about" window is
            modal and the main window will be unresponsive still close the other
            one. In Mac OS, modal windows are not supported, and previously
            approach used with popups does not work here. It is also not
            possible to use popups here, because they do not admit a layout.
            So, with this patch, at least is limited the opnening of several
            "about" windows in all platforms. In Mac OS still will be
            possible to use the main window with the "about" window opened, but
            it is not very important"""
            try:
                # In Mac OS, hide and unhide works, focus doesn't work
                # noinspection PyUnboundLocalVariable
                about_popup.Hide()
                about_popup.UnHide()
            except:
                about_popup = sg.Window(
                    "About PyToxo GUI",
                    layout=[
                        [
                            sg.Text(
                                "PyToxo GUI\nA graphical user interface for "
                                "PyToxo\n",
                                justification="center",
                            )
                        ],
                        [sg.Image(data=logo_b64_popup)],
                        [
                            sg.Text(
                                "PyToxo\nA Python library for "
                                "calculating penetrance tables of any "
                                "bivariate epistasis model\n\nCopyright 2021 Borja "
                                "González Seoane\nUniversity of A Coruña\nContact: "
                                "borja.gseoane@udc.es",
                                justification="center",
                            )
                        ],
                        [sg.Image(data=logo_udc_b64)],
                    ],
                    font=window_general_font,
                    finalize=True,
                    modal=True,
                    element_justification="center",
                )

        # Finally, check if all is filled before go to the next interaction
        check_all_filled(window, values, text_entries_to_check_values)

    window.close()


# #########################################################


if __name__ == "__main__":
    main()
